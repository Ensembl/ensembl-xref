
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

Bio::EnsEMBL::Xref::Mapper - Base mapper

=head1 SYNOPSIS

=head1 DESCRIPTION

my $mapper = Bio::EnsEMBL::Xref::Mapper->new( $xref_dba, $core_dba );
$mapper->process_file( 'file', 1 );

=cut

package Bio::EnsEMBL::Xref::Mapper;

use strict;
use warnings;

use Carp;

use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );
use Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor;

use Bio::EnsEMBL::Xref::Mapper::QC;
use Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;

=head2 new

  Description: Constructor for base Mapper.
  Returntype : Bio::EnsEMBL::Xref::Mapper
  Exceptions : none
  Caller     : general

=cut

sub new {
  my ( $caller, %args ) = @_;

  my $class = ref($caller) || $caller;
  my $self =
    bless { _xref => $args{xref_dba}, _core => $args{core_dba} },
    $class;

  confess "Required arguments missing (xref/core DBA)"
    unless defined $self->{_xref} and
    defined $self->{_core};

  assert_ref( $self->{_xref},
              'Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor' );
  assert_ref( $self->{_core}, 'Bio::EnsEMBL::DBSQL::DBAdaptor' );

  return $self;
}

=head2 xref

  Arg [1]    : Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor (optional)
  Example    : $mapper->xref($new_xref_adaptor);
  Description: Getter / Setter for the xref database.
  Returntype : Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor
  Exceptions : none

=cut

sub xref {
  my $self = shift;
  $self->{_xref} = shift if @_;

  return $self->{_xref};
}

=head2 core

  Arg [1]    : (optional)
  Example    : $mapper->core($new_core);
  Description: Getter / Setter for the core.
               info for the ensembl core database.
  Returntype : Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor
  Exceptions : none

=cut

sub core {
  my $self = shift;
  $self->{_core} = shift if @_;

  return $self->{_core};
}

=head2 farm_queue

  Arg [1]    : (optional)
  Example    : $mapper->farm_queue("long");
  Description: Getter / Setter for the farm queue.
  Returntype : string
  Exceptions : none

=cut

sub farm_queue {
  my $self = shift;
  $self->{_queue} = shift if @_;

  return $self->{_queue};
}

=head2 exonerate

  Arg [1]    : (optional)
  Example    : $mapper->exonerate("/usr/local/exonerate1.1.1");
  Description: Getter / Setter for the exonerate executable with full path.
  Returntype : string
  Exceptions : none

=cut

sub exonerate {
  my $self = shift;
  $self->{_exonerate} = shift if @_;

  return $self->{_exonerate};
}

=head2 previous_core

  Arg [1]    : (optional)
  Example    : $mapper->previous_core($old_core);
  Description: Getter / Setter for the previous release of the core db.
  Returntype : Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor
  Exceptions : none

=cut

sub previous_core {
  my $self = shift;
  $self->{_previous_core} = shift if @_;

  return $self->{_previous_core};
}

=head2 dumpcheck

  Arg [1]    : (optional)
  Example    : $mapper->dumpcheck("yes");
  Description: Getter / Setter for dumpcheck.
               If set the mapper will not dump fasta files
               if they exist already.
  Returntype : scalar
  Exceptions : none

=cut

sub dumpcheck {
  my $self = shift;
  $self->{_dumpcheck} = shift if @_;

  return $self->{_dumpcheck};
}

=head2 nofarm

=cut

sub nofarm {
  my $self = shift;
  $self->{_nofarm} = shift if @_;

  return $self->{_nofarm};
}

=head2 verbose

=cut

sub verbose {
  my $self = shift;
  $self->{_verbose} = shift if @_;

  return $self->{_verbose};
}

=head2 species_id

=cut

sub species_id {
  my $self = shift;
  $self->{_species_id} = shift if @_;

  return $self->{_species_id};
}

=head2 get_alt_alleles

=cut

sub get_alt_alleles {
  my $self = shift;

  my $dba  = $self->core->dba;
  my $aaga = Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor->new($dba);

  my $aa_list = $aaga->fetch_all();

  my $count      = scalar(@$aa_list);
  my $max_alt_id = 0;
  my %is_reference;
  my $sth;

  if ($count) {
    $self->xref->delete_alt_alleles();

    my $alt_added    = 0;
    my $num_of_genes = 0;

# Iterate through all alt-allele groups, pushing unique alleles into the xref alt allele table.
# Track the reference gene IDs.
    foreach my $aag (@$aa_list) {
      my $ref_gene = $aag->rep_Gene_id();

# Representative gene not guaranteed, try to find an alternative best fit
      unless ($ref_gene) {
        my $genes = $aag->get_all_Genes;
        foreach my $gene (@$genes) {
          if ( $gene->slice->is_reference ) {
            $ref_gene = $gene->dbID;
            last;
          }
        }
      }
      unless ($ref_gene) {
        carp
"Tried very hard but failed to select a representative gene for alt-allele-group "
          . $aag->dbID;
        next;
      }

      $is_reference{$ref_gene} = 1;
      my $others = $aag->get_all_Gene_ids('no rep');

      # Extra step in place to handle non-ref situations
      my @cleaned_others = grep { !/$ref_gene/x } @$others;

      $self->xref->add_alt_allele( $aag->dbID, $ref_gene, 1 );
      $num_of_genes++;
      $alt_added++;

      foreach my $aa (@cleaned_others) {
        $self->xref->add_alt_allele( $aag->dbID, $aa, 0 );
        $num_of_genes++;
      }

      if ( $aag->dbID > $max_alt_id ) { $max_alt_id = $aag->dbID }
    } ## end foreach my $aag (@$aa_list)

    print "$alt_added alleles found containing $num_of_genes genes\n";
  } ## end if ($count)
  else {
    print "No alt_alleles found for this species.\n";
  }

  ### LRGs added as alt_alleles in the XREF system but never added to core.

  #
  # Use $max_alt_id for new ones.
  #

  my $sql = (<<'LRG');
SELECT  ox.ensembl_id, g.gene_id
  FROM xref x, object_xref ox, external_db e, gene g
    WHERE x.xref_id = ox.xref_id AND
          e.external_db_id = x.external_db_id AND
          e.db_name like "Ens_Hs_gene" AND
          ox.ensembl_object_type = "Gene" AND
           x.display_label = g.stable_id
LRG

  $sth = $self->core->dbc->prepare($sql);
  my ( $core_gene_id, $lrg_gene_id );
  $sth->execute();
  $sth->bind_columns( \$lrg_gene_id, \$core_gene_id );

  $count = 0;

  my ( $old_count, $new_count, $lrg_count ) = ( 0, 0, 0 );

#
# If the core gene is already in an alt_allele set then use that alt_id for the LRG gene only.
# Else use a new one and add both core and LRG.
#
  while ( $sth->fetch() ) {
    my $aag = $aaga->fetch_by_gene_id($core_gene_id);
    if ($aag) {
      $self->xref->add_alt_allele( $aag->dbID, $lrg_gene_id, 0 );
      $old_count++;
    }
    else {
      $aag = $aaga->fetch_by_gene_id($lrg_gene_id);
      if ($aag) {
        $self->xref->add_alt_allele( $aag->dbID, $lrg_gene_id, 1 );
        print "LRG perculiarity\t$core_gene_id\t$lrg_gene_id\n";
        $lrg_count++;
      }
      else {
        $max_alt_id++;
        $self->xref->add_alt_allele( $max_alt_id, $lrg_gene_id,  0 );
        $self->xref->add_alt_allele( $max_alt_id, $core_gene_id, 1 );
        $new_count++;
      }
    }

    $count++;
  }

  if ($count) {
    print
"Added $count alt_allels for the lrgs. $old_count added to previous alt_alleles and $new_count new ones\n";
    print "LRG problem count = $lrg_count\n";
  }

  $self->xref->update_process_status("alt_alleles_added");

  return;
} ## end sub get_alt_alleles

#
# Default behaviour is not to do the offical naming
# Overload this method in the species file returning the
# official database name to do so.
# (ie, human-> HGNC, mouse ->MGI, zebrafisf -> ZFIN_ID)
#

=head2 get_official_name

=cut

sub get_official_name {
  return;
}

#
# Biomart insists that a source is linked to only one ensembl
# object type (Gene, Transcript, Translation). So biomart_fix
# will move $dbname entry for type1 to type 2
# i.e. move all HGNC from transcripts to Genes.
#

=head2 biomart_fix

=cut

sub biomart_fix {
  my ( $self, $db_name, $type1, $type2, $verbose ) = @_;

  print
    "$db_name is associated with both $type1 and $type2 object types\n"
    if $verbose;
  print "$db_name moved to Gene level.\n" unless $verbose;

  my ( $to, $from, $to_id, $from_id );
  if ( $type1 eq "Gene" or $type2 eq "Gene" ) {
    $to    = "Gene";
    $to_id = "gene_id";

    if ( $type1 eq "Translation" or $type2 eq "Translation" ) {
      $from    = "Translation";
      $from_id = "translation_id";
    }
    else {
      $from    = "Transcript";
      $from_id = "transcript_id";
    }
  }
  else {
    $to      = "Transcript";
    $to_id   = "transcript_id";
    $from    = "Translation";
    $from_id = "translation_id";
  }

  if ( $db_name eq 'GO' || $db_name eq 'goslim_goa' ) {
    $to      = 'Translation';
    $from    = 'Transcript';
    $to_id   = 'translation_id';
    $from_id = 'transcript_id';
  }

  print "Therefore moving all associations from $from to " . $to . "\n"
    if $verbose;

  my $sql = (<<"EOF");
  UPDATE IGNORE object_xref, gene_transcript_translation, xref, source
    SET object_xref.ensembl_object_type = "$to",
      object_xref.ensembl_id = gene_transcript_translation.$to_id 
	WHERE object_xref.ensembl_object_type = "$from" AND
	  object_xref.ensembl_id = gene_transcript_translation.$from_id AND
	    xref.xref_id = object_xref.xref_id AND
	      xref.source_id = source.source_id AND
                object_xref.ox_status = "DUMP_OUT"  AND
		  source.name = "$db_name";
EOF
  my $result = $self->xref->dbc->do($sql);

  if ( $db_name eq "GO" || $db_name eq 'goslim_goa' ) {
    $sql = (<<"EOF2");
  DELETE object_xref, identity_xref, go_xref
    FROM object_xref, xref, source, identity_xref, go_xref
      WHERE object_xref.ensembl_object_type = "$from" AND
        identity_xref.object_xref_id = object_xref.object_xref_id AND
	xref.xref_id = object_xref.xref_id AND
          go_xref.object_xref_id = object_xref.object_xref_id AND
	  xref.source_id = source.source_id AND
            object_xref.ox_status = "DUMP_OUT"  AND
	      source.name = "$db_name";
EOF2

    $result = $self->xref->dbc->do($sql);

    # Special tidying up for transcripts without translation
    # The resulting object_xref does not have an ensembl_id to map to

    $sql = (<<"EOF4");
  DELETE object_xref, identity_xref, go_xref
    FROM object_xref, xref, source, identity_xref, go_xref
      WHERE object_xref.ensembl_object_type = "$to" AND
        identity_xref.object_xref_id = object_xref.object_xref_id AND
        xref.xref_id = object_xref.xref_id AND
          go_xref.object_xref_id = object_xref.object_xref_id AND
          xref.source_id = source.source_id AND
            object_xref.ensembl_id = 0 AND
              object_xref.ox_status = "DUMP_OUT"  AND
                source.name = "$db_name";
EOF4
  } ## end if ( $db_name eq "GO" ...)
  else {
    $sql = (<<"EOF3");
  DELETE object_xref, identity_xref
    FROM xref, source, object_xref
      LEFT JOIN identity_xref
        ON identity_xref.object_xref_id = object_xref.object_xref_id
      WHERE object_xref.ensembl_object_type = "$from" AND
	xref.xref_id = object_xref.xref_id AND
	  xref.source_id = source.source_id AND
            object_xref.ox_status = "DUMP_OUT"  AND
	      source.name = "$db_name";
EOF3

    $result = $self->xref->dbc->do($sql);
  }

  # delete dependent_xref
  $sql = (<<'EOF4');
  DELETE FROM dependent_xref WHERE object_xref_id NOT IN 
   (SELECT object_xref_id FROM object_xref);
EOF4

  return;
} ## end sub biomart_fix

#
# This sub finds which source lie on multiple ensembl object types and calls biomart_fix to fix this.
#

=head2 biomart_testing

=cut

sub biomart_testing {
  my ($self) = @_;

  my $sql =
'SELECT ox.ensembl_object_type, COUNT(*), s.name  FROM xref x, object_xref ox, source s  WHERE x.xref_id = ox.xref_id AND s.source_id = x.source_id  and ox.ox_status = "DUMP_OUT" GROUP BY s.name, ox.ensembl_object_type';

  my $again = 1;
  while ($again) {
    $again = 0;

    my $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();

    my ( $type,      $count,      $name );
    my ( $last_type, $last_count, $last_name );
    $sth->bind_columns( \$type, \$count, \$name );

    $last_name = "DEFAULT";
    while ( not $again and $sth->fetch ) {
      if ( $last_name eq $name ) {
        $again = 1;
        $self->biomart_fix( $name, $last_type, $type, 1 );
      }

      $last_name  = $name;
      $last_type  = $type;
      $last_count = $count;
    }
    $sth->finish;
  }

  my $tester = Bio::EnsEMBL::Xref::Mapper::QC->new($self);
  confess 'Problems found before source_defined_move'
    if $tester->unlinked_entries;

  $self->xref->update_process_status('biomart_test_finished');

  return;
} ## end sub biomart_testing

#
# Similar to above but just reports the problems. It does not fix them
#

=head2 biomart_test

=cut

sub biomart_test {
  my $self = shift;

  my $sql =
'SELECT ox.ensembl_object_type, COUNT(*), s.name  FROM xref x, object_xref ox, source s  WHERE x.xref_id = ox.xref_id AND s.source_id = x.source_id  and ox.ox_status = "DUMP_OUT" GROUP BY s.name, ox.ensembl_object_type';

  my $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();

  my ( $type,      $count,      $name );
  my ( $last_type, $last_count, $last_name );
  $sth->bind_columns( \$type, \$count, \$name );

  $last_name = "NOTKNOWN";
  my $first = 1;
  while ( $sth->fetch ) {
    if ( $last_name eq $name ) {
      if ($first) {
        print STDERR "\nProblem Biomart test fails\n";
        $first = 0;
      }

      print STDERR "$last_name\t$last_count\t$last_type\n";
      print STDERR "$name\t$count\t$type\n";
    }

    $last_name  = $name;
    $last_type  = $type;
    $last_count = $count;
  }
  $sth->finish;

  return;
} ## end sub biomart_test

# remove a list of patterns from a string

=head2 filter_by_regexp

=cut

sub filter_by_regexp {

  my ( $self, $str, $regexps ) = @_;

  foreach my $regexp (@$regexps) {
    $str =~ s/$regexp//xig;
  }

  return $str;
}

=head2 get_species_id_from_species_name {

=cut

sub get_species_id_from_species_name {
  my ( $self, $name ) = @_;

  my $species_id;
  eval {
    $species_id = $self->xref->get_id_from_species_name($name);
    1;
    } or
    do {
    my $error = $@ || "Couldn't get ID for species name $name\n";
    $error .= "It must be one of :-\n" .
      join( "\n", @{ $self->xref->get_species_names } );
    confess $error;
    };

  return $species_id;
}

=head2 revert_to_parsing_finished

=cut

sub revert_to_parsing_finished {
  my $self = shift;

  $self->xref->clean_up();
  $self->xref->remove_mapping_data();
  $self->xref->update_process_status('parsing_finished');

  return;
}

=head2 revert_to_mapping_finished

=cut

sub revert_to_mapping_finished {
  my $self = shift;

  $self->xref->clean_up( undef, 1 );
  $self->xref->update_mapping_jobs_status('SUBMITTED');
  $self->update_process_status('mapping_finished');

  return;
}

=head2 process_alt_alleles

=cut

sub process_alt_alleles {
  my $self = shift;

  my ( $added_count, $ignored ) = $self->xref->process_alt_alleles();
  print "Added $added_count new mapping but ignored $ignored\n";

  my $tester = Bio::EnsEMBL::Xref::Mapper::QC->new($self);
  confess 'Problems found after process_alt_alleles'
    if $tester->unlinked_entries;

  $self->xref->update_process_status('alt_alleles_processed');

  return;
}

#
# Here we do the moving.
#

=head2 source_defined_move

=cut

sub source_defined_move {
  my $self = shift;

  foreach my $source ( @{ $self->xref->get_gene_specific_list() } ) {
    $self->biomart_fix( $source, "Translation", "Gene", undef, undef );
    $self->biomart_fix( $source, "Transcript",  "Gene", undef, undef );
  }

  my $tester = Bio::EnsEMBL::Xref::Mapper::QC->new($self);
  confess 'Problems found after source_defined_move'
    if $tester->unlinked_entries;

  $self->update_process_status('source_level_move_finished');

  return;
}

=head2 log_progress
  Arg [1]    : String to format and print
  Arg [2]    : params to fmt
  Description: utility method to log the mapping progress
  Return type: None
  Caller     : internal

=cut

sub log_progress {
  my ( $self, $fmt, @params ) = @_;

  return if ( ! $self->verbose );
  printf( STDERR "COORD==> %s", sprintf( $fmt, @params ) );
  return;
}


1;
