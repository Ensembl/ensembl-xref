=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Xref::Mapper::Loader - An adaptor for loading core database

=cut

=head1 DESCRIPTION

This is the base adaptor for loading xrefs into the species core database from
the xref database in production. The aim of the module is to reduce redundancy
of code and keep the SQL to a constrained set of functions

=cut

package Bio::EnsEMBL::Xref::Mapper::Loader;
use strict;
use warnings;
use Carp;
use Bio::EnsEMBL::Utils::Exception;

use parent qw( Bio::EnsEMBL::Xref::Mapper );

use DBI;

sub new {
  my( $class, $mapper ) = @_;

  my $self ={};
  bless $self,$class;
  $self->core( $mapper->core );
  $self->xref( $mapper->xref );
  $self->mapper( $mapper );
  return $self;
}


=head2 mapper


=cut

sub mapper{
  my ($self, $arg) = @_;

  if ( defined $arg ) {
    $self->{_mapper} = $arg;
  }
  return $self->{_mapper};
}


=head2 update


=cut

sub update {
  my ($self, $arg) = @_;
  # remove xref, object_xref, identity_xref, depenedent_xref, go_xref, unmapped_object, (interpro???), external_synonym, projections.

  my $verbose = $self->mapper->verbose;
  my $core_dbi = $self->core->dbc;
  my $xref_dbi = $self->xref->dbc;

  #####################################
  # first remove all the projections. #
  #####################################

  <core_tmp_functions>->delete_projected_xrefs();


  #########################################
  # Get source_id to external_db_id       #
  #########################################

  $self->name_to_external_db_id = %{ <core_tmp_functions>->get_xref_external_dbs() };

  my %source_id_to_external_db_id;
  my %source_list = $self->xref->get_valid_source_id_to_external_db_id();
  while ( my ( $name, $id ) = each %source_list ) {
    if ( defined $self->name_to_external_db_id{ $name } ) {
      $source_id_to_external_db_id{ $id } = $self->name_to_external_db_id{ $name };
    }
    elsif ( $name =~ m/notransfer$/xms) {}
    else {
      confess "ERROR: Could not find $name in external_db table please add this too continue";
    }
  }

  $self->xref->mark_mapped_xrefs_already_run(); # just incase this is being ran again


  ######################################
  # For each external_db to be updated #
  # Delete the existing ones           #
  ######################################
  while( my $base_sources = $self->xref->get_source_ids_with_xrefs() ) {
    if ( !defined $self->name_to_external_db_id{ $base_sources->name } ) {
      next;  #must end in notransfer
    }

    <core_tmp_functions>->delete_by_external_db_id(
      $self->name_to_external_db_id{$base_sources->name}
    );
  }


  ##########################################
  # Get the offsets for object_xref, xref  #
  ##########################################

  my %offsets = <core_tmp_functions>->parsing_stored_data();
  my $xref_offset = $offset{ 'xref' };
  my $object_xref_offset = $offset{ 'object_xref' };


  ####################
  # Get analysis id's
  ####################

  my %analysis_ids = $self->get_analysis();
  my $checksum_analysis_id; #do not populate until we know we need this

  print "xref offset is $xref_offset, object_xref offset is $object_xref_offset\n";


  ##########################################
  # Now add the new ones from xref to core #
  ##########################################

  $self->map_xrefs_from_xrefdb_to_coredb();


  #######################################
  # Remember to do unmapped entries
  # 1) make sure the reason exist/create them and get the ids for these.
  # 2) Process where dumped is null and type = DIRECT, DEPENDENT, SEQUENCE_MATCH, MISC seperately
  ########################################
  my %reason_id;

  # Get the cutoff values
  my $failed_sources = $self->xref->get_unmapped_regions();

  my %summary_failed = $failed_sources->summary;
  my %desc_failed    = $failed_sources->desc;

  foreach my $key (keys %desc_failed){
    my $failed_id = <core_tmp_functions>->get_unmapped_reasons( $desc_failed{$key} );

    if(!defined $failed_id ) {
      $failed_id = <core_tmp_functions>->add_unmapped_reason(
        $summary_failed{$key}, $desc_failed{$key} );
    }
    $reason_id{$key} = $failed_id;
  }

  $transaction_start_sth->execute();

  # DIRECT #
  $self->load_unmapped_direct_xref( $xref_offset );

  # MISC #
  $self->load_unmapped_misc_xref( $xref_offset );

  # DEPENDENT #
  $self->load_unmapped_dependent_xref( $xref_offset );

  # SEQUENCE_MATCH #
  $self->load_unmapped_sequence_xrefs( $xref_offset, $analysis_id );

  # WEL (What ever is left) #
  # These are those defined as dependent but the master never existed and the xref and their descriptions etc are loaded first
  # with the dependencys added later so did not know they had no masters at time of loading.
  # (e.g. EntrezGene, WikiGene, MIN_GENE, MIM_MORBID)
  $self->load_unmapped_other_xref( $xref_offset, $analysis_id );


  $transaction_end_sth->execute();

  $self->xref->insert_process_status( 'core_loaded' );

  return;
} ## end sub update


sub map_xrefs_from_xrefdb_to_coredb {
  my ( $self, $xref_offset ) = @_;

  $transaction_start_sth->execute();

  while ( my $xref_handle = $self->xref->get_dump_out_xrefs() ) {

    next if( !defined $name_to_external_db_id{ $xref_handle->name } );

    if( defined $xref_handle->where_from and $xref_handle->where_from ne q{} ){
      $xref_handle->where_from = "Generated via $where_from";
    }
    my $ex_id = $self->name_to_external_db_id{ $xref_handle->name };

    if ( $verbose ) {
      print "updating ($xref_handle->source_id) $xref_handle->name in core (for $type xrefs)\n";
    }

    my @xref_list=();  # process at end. Add synonyms and set dumped = 1;

    # dump SEQUENCE_MATCH, DEPENDENT, DIRECT, COORDINATE_OVERLAP, INFERRED_PAIR, (MISC?? same as direct come from official naming)

    ### If DIRECT ,         xref, object_xref,                  (order by xref_id)  # maybe linked to more than one?
    ### if INFERRED_PAIR    xref, object_xref
    ### if MISC             xref, object_xref

    if ( $xref_handle->type eq 'DIRECT' or $xref_handle->type eq 'INFERRED_PAIR' or
         $xref_handle->type eq 'MISC'   or $xref_handle->type eq 'SEQUENCE_MATCH' ) {
      push @xref_list, $self->load_identity_xref(
        $xref_handle->source_id, $xref_handle->type, $xref_offset, $ex_id );
    }
    elsif ($xref_handle->type eq 'CHECKSUM') {
      if(! defined $checksum_analysis_id) {
        $checksum_analysis_id = $self->get_single_analysis( 'xrefchecksum' );
      }

      push @xref_list, $self->load_checksum_xref(
        $xref_handle->source_id, $xref_handle->type, $xref_offset, $ex_id );
    }
    elsif ( $xref_handle->type eq 'DEPENDENT' ) {
      push @xref_list, $self->load_dependent_xref(
        $xref_handle->source_id, $xref_handle->type, $xref_offset, $ex_id );
    }
    else {
      print "PROBLEM:: what type is $type\n";
    }


    # Transfer data for synonym and set xref database xrefs to dumped.
    if ( @xref_list ) {
      $self->load_synonyms( @xref_list );

      $self->xref->mark_mapped_xrefs( @xref_list, 'MAPPED' );
    }

    # Update the core databases release in for source form the xref database
    if ( defined $xref_handle->release_info and $xref_handle->release_info ne q{1} ){
       my $add_release_info_sth   = $core_dbi->prepare('UPDATE external_db SET db_release = ? WHERE external_db_id = ?');
       $add_release_info_sth->execute($xref_handle->release_info, $ex_id) ||
         confess "Failed to add release info **$xref_handle->release_info** for external source $ex_id\n";
    }
  } ## end while for getting dump out xrefs
  $transaction_end_sth->execute();

  return;
}


sub load_unmapped_direct_xref {
  my ( $self, $xref_offset, $analysis_id ) = @_;

  my @xref_list = ();
  # my $analysis_id = $analysis_ids{'Transcript'};   # No real analysis here but in table it is set to not NULL
  while ( my $direct_handle = $self->xref->get_insert_direct_xref_low_priority() ) {
    if( defined $self->name_to_external_db_id{ $direct_handle->dbname } ) {
      $xref_id = $self->add_xref(
        $xref_offset,
        $direct_handle->xref_id,
        $self->name_to_external_db_id{ $direct_handle->dbname },
        $direct_handle->acc,
        $direct_handle->label,
        $direct_handle->version,
        $direct_handle->desc,
        'UNMAPPED',
        $direct_handle->info,
        $core_dbi
      );

      <core_tmp_functions>->add_unmapped_object( {
        analysis_id        => $analysis_id,
        external_db_id      => $self->name_to_external_db_id{ $direct_handle->dbname },
        identifier         => $acc,
        unmapped_reason_id => $reason_id{'NO_STABLE_ID'} } );

      push @xref_list, $xref_id;
    }
  }

  if ( @xref_list ) {
    $self->xref->mark_mapped_xrefs( @xref_list, 'UNMAPPED_NO_STABLE_ID' );
  }

  return;
}



sub load_unmapped_dependent_xref {
  my ( $self, $xref_offset, $analysis_id ) = @_;

  my @xref_list = ();

  @xref_list = ();
  my $last_acc= 0;
  # $analysis_id = $analysis_ids{'Transcript'};
  while ( my $dependent_handle = $self->xref->get_insert_dependent_xref_low_priority() ) {
    if( !defined $self->name_to_external_db_id{ $dependent_handle->dbname } ){
      next;
    }
    if($last_acc ne $acc){
      $xref_id = $self->add_xref(
        $xref_offset,
        $dependent_handle->xref_id,
        $self->name_to_external_db_id{ $dependent_handle->dbname },
        $dependent_handle->acc,
        $dependent_handle->label || $dependent_handle->acc,
        $dependent_handle->version,
        $dependent_handle->desc,
        'UNMAPPED',
        $dependent_handle->info,
        $core_dbi);
    }
    $last_acc = $acc;
    <core_tmp_functions>->add_unmapped_object( {
      analysis_id        => $analysis_id,
      external_db_id     => $self->name_to_external_db_id{ $dependent_handle->dbname },
      identifier         => $dependent_handle->acc,
      unmapped_reason_id => $reason_id{ 'MASTER_FAILED' },
      parent             => $dependent_handle->parent } );
    push @xref_list, $xref_id;
  }

  if ( @xref_list ) {
    $self->xref->mark_mapped_xrefs( @xref_list, 'UNMAPPED_MASTER_FAILED' );
  }

  return;
}



sub load_unmapped_sequence_xrefs {
  my ( $self, $xref_offset, $analysis_id ) = @_;

  @xref_list = ();
  my $last_xref = 0;
  while( my $xref_handle = $self->xref->get_insert_sequence_xref_remaining() ){
    if (
      !defined $self->name_to_external_db_id{ $dbname } ||
      ( defined $xref_handle->status and $xref_handle->status eq 'FAILED_PRIORITY')
    ) {
      next;
    }
    if($last_xref != $xref_id){
      $xref_id = $self->add_xref(
        $xref_offset,
        $xref_handle->xref_id,
        $self->name_to_external_db_id{ $dbname },
        $xref_handle->acc,
        $xref_handle->label,
        $xref_handle->version,
        $xref_handle->desc,
        'UNMAPPED',
        $xref_handle->info,
        $core_dbi);
    }
    $last_xref = $xref_id;

    if( defined $ensembl_id ) {
      $analysis_id= $analysis_ids{$xref_handle->ensembl_object_type};

      <core_tmp_functions>->add_unmapped_object( {
        analysis_id         => $analysis_id,
        external_db_id      => $self->name_to_external_db_id{ $dbname },
        identifier          => $xref_handle->acc,
        unmapped_reason_id  => $reason_id{$dbname},
        query_score         => $xref_handle->q_id,
        target_score        => $xref_handle->t_id,
        ensembl_id          => $xref_handle->ensembl_id,
        ensembl_object_type => $xref_handle->ensembl_object_type
      } );
    }
    else{
      if($xref_handle->seq_type eq 'dna'){
        $ensembl_object_type = 'Transcript';
      }
      else{
        $ensembl_object_type = 'Translation';
      }
      $analysis_id = $analysis_ids{$xref_handle->ensembl_object_type};
      <core_tmp_functions>->add_unmapped_object( {
        analysis_id         => $analysis_id,
        external_db_id      => $self->name_to_external_db_id{ $dbname },
        identifier          => $xref_handle->acc,
        unmapped_reason_id  => $reason_id{ 'MASTER_FAILED' },
        ensembl_object_type => $xref_handle->ensembl_object_type } );
    }
    push @xref_list, $xref_id;
  }


  if ( @xref_list ) {
    $self->xref->mark_mapped_xrefs( @xref_list, 'UNMAPPED_NO_MAPPING' );
  }

  return;
}


sub load_unmapped_misc_xref {
  my ( $self, $xref_offset, $analysis_id ) = @_;

  @xref_list = ();
  # $analysis_id = $analysis_ids{'Transcript'};   # No real analysis here but in table it is set to not NULL
  while( my $misc_handle = $self->xref->get_insert_misc_xref() ){
    if ( defined $name_to_external_db_id{ $dbname } ) {
      $xref_id = $self->add_xref(
        $xref_offset,
        $misc_handle->xref_id,
        $self->name_to_external_db_id{ $misc_handle->dbname },
        $misc_handle->acc,
        $misc_handle->label,
        $misc_handle->version,
        $misc_handle->desc,
        'UNMAPPED',
        $misc_handle->info,
        $core_dbi);

      <core_tmp_functions>->add_unmapped_object( {
        analysis_id        => $analysis_id,
        external_db_id     => $self->name_to_external_db_id{ $misc_handle->dbname },
        identifier         => $misc_handle->acc,
        unmapped_reason_id => $reason_id{'NO_MAPPING'} } );
      push @xref_list, $xref_id;
    }
  }

  if ( @xref_list ) {
    $self->xref->mark_mapped_xrefs( @xref_list, 'UNMAPPED_NO_MAPPING' );
  }

  return;
}



sub load_unmapped_other_xref {
  ###########################
  # WEL (What ever is left).#
  ###########################

  # These are those defined as dependent but the master never existed and the xref and their descriptions etc are loaded first
  # with the dependencys added later so did not know they had no masters at time of loading.
  # (e.g. EntrezGene, WikiGene, MIN_GENE, MIM_MORBID)

  my ( $self, $xref_offset, $analysis_id ) = @_;

  $set_unmapped_sth  =  $core_dbi->prepare("insert into unmapped_object (type, analysis_id, external_db_id, identifier, unmapped_reason_id) values ('xref', ?, ?, ?, '".$reason_id{"NO_MASTER"}."')");

  # $analysis_id = $analysis_ids{'Transcript'};   # No real analysis here but in table it is set to not NULL
  @xref_list = ();
  while( my $other_handle = $self->xref->get_insert_other_xref() ){
    if( !defined $self->name_to_external_db_id{ $other_handle->dbname } ){
      next;
    }

    $xref_id = $self->add_xref(
      $xref_offset,
      $other_handle->xref_id,
      $self->name_to_external_db_id{ $other_handle->dbname },
      $other_handle->acc,
      $other_handle->label,
      $other_handle->version,
      $other_handle->desc,
      'UNMAPPED',
      $other_handle->info,
      $core_dbi
    );

    <core_tmp_functions>->add_unmapped_object( {
      analysis_id        => $analysis_id,
      external_db_id     => $self->name_to_external_db_id{ $other_handle->dbname },
      identifier         => $other_handle->acc,
      unmapped_reason_id => $reason_id{ 'NO_MASTER' }
    } );
    push @xref_list, $xref_id;
  }

  if ( @xref_list ) {
    $self->xref->mark_mapped_xrefs( @xref_list, 'UNMAPPED_NO_MASTER' );
  }

  return;
}



sub load_identity_xref {
   my ( $self, $source_id, $type, $xref_offset, $ex_id ) = @_;

  my $last_xref = 0;
  my @xref_list = ();
  while ( my $identity_xref_handle = $self->xref->get_identity_xref( $source_id, $type ) ) {
    if( $last_xref != $identity_xref_handle->xref_id ) {
      push @xref_list, $identity_xref_handle->xref_id;
      $count++;
      $xref_id = <core_tmp_functions>->add_xref(
        $xref_offset,
        $identity_xref_handle->xref_id,
        $self->name_to_external_db_id{ $identity_xref_handle->dbname },
        $identity_xref_handle->acc,
        $identity_xref_handle->label,
        $identity_xref_handle->version,
        $identity_xref_handle->desc,
        $identity_xref_handle->type,
        $identity_xref_handle->info || $where_from);
      $last_xref = $xref_id;
    }

    $object_xref_id = <core_tmp_functions>->add_object_xref(
      $identity_xref_handle->object_xref_offset,
      $identity_xref_handle->object_xref_id,
      $identity_xref_handle->ensembl_id,
      $identity_xref_handle->ensembl_type,
      ( $xref_id + $xref_offset),
      $identity_xref_handle->analysis_ids{ $identity_xref_handle->ensembl_type },
      $core_dbi);

    if $translation_start {
      <core_tmp_functions>->$add_identity_xref( {
        object_xref_id   => ( $identity_xref_handle->object_xref_id + $object_xref_offset ),
        query_identity   => $identity_xref_handle->query_identity,
        ensembl_identity => $identity_xref_handle->target_identity,
        xref_start       => $identity_xref_handle->hit_start,
        xref_end         => $identity_xref_handle->hit_end,
        ensembl_start    => $identity_xref_handle->translation_start,
        ensembl_end      => $identity_xref_handle->translation_end,
        cigar_line       => $identity_xref_handle->cigar_line,
        score            => $identity_xref_handle->score,
        evalue           => $identity_xref_handle->evalue
      } );
    }
  }

  return @xref_list;
}



sub load_checksum_xref {
  my ( $self, $source_id, $type, $xref_offset, $ex_id ) = @_;
  my $count = 0;
  my $last_xref = 0;
  my @xref_list = ();
  while( my $checksum_xref_handle = $self->xref->get_insert_checksum_xref( $source_id, $type ) ) {
    if($last_xref != $checksum_xref_handle->xref_id) {
      push @xref_list, $checksum_xref_handle->xref_id;
      $count++;
      $xref_id = $self->add_xref(
        $xref_offset, $xref_id, $ex_id,
        $checksum_xref_handle->acc,
        $checksum_xref_handle->label,
        $checksum_xref_handle->version,
        $checksum_xref_handle->desc,
        $checksum_xref_handle->type,
        $checksum_xref_handle->info || $where_from, $core_dbi);
      $last_xref = $xref_id;
    }
    my $object_xref_id = $self->add_object_xref(
      $object_xref_offset,
      $checksum_xref_handle->object_xref_id,
      $checksum_xref_handle->ensembl_id,
      $checksum_xref_handle->ensembl_type,
      ( $checksum_xref_handle->xref_id + $xref_offset ),
      $checksum_analysis_id,
      $core_dbi);
  }
  if ( $verbose ) {
    print "CHECKSUM $count\n";
  }

  return @xref_list;
}



sub load_dependent_xref {
  my ( $self, $source_id, $type, $xref_offset, $ex_id ) = @_;

  my $count = 0;
  my $ox_count = 0;
  my @master_problems;
  my $err_master_count=0;
  # $dependent_sth->execute($source_id, $type);
  # my ($xref_id, $acc, $label, $version, $desc, $info, $object_xref_id, $ensembl_id, $ensembl_type, $master_xref_id);
  # $dependent_sth->bind_columns(\$xref_id, \$acc, \$label, \$version, \$desc, \$info, \$object_xref_id, \$ensembl_id, \$ensembl_type, \$master_xref_id);
  my $last_xref = 0;
  my $last_ensembl = 0;
  my @xref_list = ();
  while ( my $dependent_xref_handle = $self->xref->get_insert_dependent_xref( $source_id, $type ) ) {
    if($last_xref != $xref_id) {
      push @xref_list, $xref_id;
      $count++;
      $xref_id = $self->add_xref(
        $xref_offset, $xref_id, $ex_id,
        $dependent_xref_handle->acc,
        $dependent_xref_handle->label || $dependent_xref_handle->acc,
        $dependent_xref_handle->version,
        $dependent_xref_handle->desc,
        $dependent_xref_handle->type,
        $dependent_xref_handle->info || $where_from,
        $core_dbi);
      $last_xref = $xref_id;
    }
    if ( $last_xref != $xref_id or $last_ensembl != $ensembl_id ) {
      my $object_xref_id = $self->add_object_xref(
        $object_xref_offset,
        $dependent_xref_handle->object_xref_id,
        $dependent_xref_handle->ensembl_id,
        $dependent_xref_handle->ensembl_type,
        ( $dependent_xref_handle->xref_id + $xref_offset ),
        $analysis_ids{ $dependent_xref_handle->ensembl_type },
        $core_dbi
      );

      if ( defined $master_xref_id ) { # need to sort this out for FlyBase since there are EMBL direct entries from the GFF and dependent xrefs from Uniprot
        # $add_dependent_xref_sth->execute(($object_xref_id+$object_xref_offset), ($master_xref_id+$xref_offset), ($xref_id+$xref_offset) );
        <core_tmp_functions>->add_dependent_xref(
          ( $object_xref_id + $object_xref_offset ),
          ( $master_xref_id + $xref_offset ),
          ( $xref_id+$xref_offset )
        );
      }
      else {
        if ( $err_master_count < 10 ) {
          push @master_problems, $acc;
        }
        $err_master_count++;
      }
      $ox_count++;
    }
    $last_xref = $xref_id;
    $last_ensembl = $ensembl_id;
  }
  if ( @master_problems ){
    print "WARNING:: for $name $err_master_count problem master xrefs\nExamples are :-\t";
    print join ', ', @master_problems;
    print "\n";
  }

  if ( $verbose ) {
    print "DEP $count xrefs, $ox_count object_xrefs\n";
  }

  return @xref_list;
}



sub load_synonyms {
  my ( $self, $xref_list) = @_;

  my $syn_count = 0;

  my ($xref_id, $syn);

  while( my $syn_handle = $self->xref->get_synonyms_for_xref( @{ $xref_list } ) ) {
    <core_tmp_functions>->add_xref_synonym(
      $syn_handle->xref_id, $syn_handle->syn );
    $syn_count++;
  }

  if ( $syn_count ) {
    print "\tadded $syn_count synonyms\n";
  }

  return;
}



=head2 get_analysis


=cut

sub get_analysis{
  my $self = shift;
  my %type_to_logic_name = ( 'Gene'        => 'xrefexoneratedna',
                             'Transcript'  => 'xrefexoneratedna',
                             'Translation' => 'xrefexonerateprotein', );
  my %analysis_id;
  foreach my $key (qw(Gene Transcript Translation)){
    my $logic_name = $type_to_logic_name{$key};
    $analysis_id{$key} = $self->get_single_analysis($logic_name);
  }
  return %analysis_id;
} ## end sub get_analysis


=head2 get_single_analysis


=cut

sub get_single_analysis {
  my ($self, $logic_name) = @_;
  my $h = $self->core->dbc()->sql_helper();
  my $analysis_ids = $h->execute_simple(
    -SQL => 'SELECT analysis_id FROM analysis WHERE logic_name=?',
    -PARAMS => [$logic_name]
  );
  my $analysis_id;

  if(@{$analysis_ids}) {
    $analysis_id = $analysis_ids->[0];
  }
  else {
    if ( $self->verbose ) {
      print "No analysis with logic_name $logic_name found, creating ...\n";
    }
    # TODO - other fields in analysis table
    $self->core()->dbc()->sql_helper()->execute_update(
      -SQL => 'INSERT INTO analysis (logic_name, created) VALUES (?,NOW())',
      -PARAMS => [$logic_name],
      -CALLBACK => sub {
        my ($sth) = @_;
        $analysis_id = $sth->{'mysql_insertid'};
        return;
      }
    );
  }

  return $analysis_id;
} ## end sub get_single_analysis


=head2 add_xref


=cut

sub add_xref {
  my ($self, $offset, $xref_id, $external_db_id, $dbprimary_acc, $display_label, $version, $description, $info_type, $info_text, $dbc)  = @_;
  my $select_sth = $dbc->prepare("select xref_id from xref where dbprimary_acc = ? and external_db_id = ? and info_type = ? and info_text = ? and version = ?");
  my $insert_sth = $dbc->prepare("insert into xref (xref_id, external_db_id, dbprimary_acc, display_label, version, description, info_type, info_text) values (?, ?, ?, ?, ?, ?, ?, ?)");
  my $new_xref_id;
  $select_sth->execute($dbprimary_acc, $external_db_id, $info_type, $info_text, $version);
  $select_sth->bind_columns(\$new_xref_id);
  $select_sth->fetch();
  if (!$new_xref_id) {
    $insert_sth->execute(($xref_id+$offset), $external_db_id, $dbprimary_acc, $display_label, $version, $description, $info_type, $info_text);
    return $xref_id;
  }

  return $new_xref_id - $offset;
} ## end sub add_xref


=head2 add_object_xref


=cut

sub add_object_xref {
  my ($self, $offset, $object_xref_id, $ensembl_id, $ensembl_object_type, $xref_id, $analysis_id, $dbc) = @_;
  my $select_sth = $dbc->prepare("select object_xref_id from object_xref where xref_id = ? and ensembl_object_type = ? and ensembl_id = ? and analysis_id = ?");
  my $insert_sth = $dbc->prepare("insert ignore into object_xref (object_xref_id, ensembl_id, ensembl_object_type, xref_id, analysis_id) values (?, ?, ?, ?, ?)");
  my $new_object_xref_id;
  $select_sth->execute($xref_id, $ensembl_object_type, $ensembl_id, $analysis_id);
  $select_sth->bind_columns(\$new_object_xref_id);
  $select_sth->fetch();
  if (!$new_object_xref_id) {
    $insert_sth->execute(($object_xref_id+$offset), $ensembl_id, $ensembl_object_type, $xref_id, $analysis_id);
    return $object_xref_id;
  }

  return $new_object_xref_id - $offset;
} ## end sub add_object_xref

1;
