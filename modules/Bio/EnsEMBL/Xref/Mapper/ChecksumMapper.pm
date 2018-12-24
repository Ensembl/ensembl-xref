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

package Bio::EnsEMBL::Xref::Mapper::ChecksumMapper;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::SeqIO;

use parent qw(Bio::EnsEMBL::Xref::Mapper);

=head2 new

Arg [1]     : HashRef - Argument hash, regular options from Bio::EnsEMBL::Xref::Mapper
                        supplemented with those listed below
Description : Overridden constructor adds new options for checksum matching
              These are:
                external_db_name
                logic_name - the analysis name for new xrefs
                method - A coderef that can compare checksums
                object_type - Ensembl feature type, e.g. Transcript, Translation, Gene
                batch_size - how many checksums to compare with the Xref DB at once.
                             Uniparc has an awful lot of accessions
Returntype  : Bio::EnsEMBL::Xref::Mapper::ChecksumMapper

=cut

sub new {
  my ($caller, %args) = @_;

  my $self = $caller->SUPER::new(%args);

  $self->{_external_db_name} = $args{external_db_name};

  $self->{_logic_name} = $args{logic_name} // 'xrefchecksum';
  $self->{_method} = $args{method} // \&perform_mapping;
  $self->{_object_type} //= 'Translation';
  $self->{_batch_size} = $args{batch_size} // 1000;

  # Cache this so we don't have to keep getting it from the DB
  $self->{_source_id} = $self->xref->get_source_id_for_source_name($self->external_db_name);

  return $self;
}

=head2 logic_name

Description : Returns the logic name for this analysis
Returntype  : String

=cut

sub logic_name {
  my $self = shift;
  return $self->{_logic_name};
}

=head2 method

Arg [1]     : CodeRef - a function to compare checksums
Description : Getter/setter for the checksum method to use
Returntype  : CodeRef

=cut

sub method {
  my ($self, $method) = @_;
  $self->{_method} = $method if defined $method;
  return $self->{_method};
}

=head2 external_db_name

Description : Supplies the external_db_name to use as the source of checksums
              Then we are able to create the Xrefs with the correct source
Returntype  : String - matches the db_name in external_db or the source table

=cut

sub external_db_name {
  my $self = shift;
  return $self->{_external_db_name};
}

=head2 object_type

Description : The object type used for the current variety of Xref checksum matching
Returntype  : String matching the Ensembl feature types e.g. Gene, Transcript,
              Exon, Translation

=cut

sub object_type {
  my $self = shift;
  return $self->{_object_type};
}


=head2 batch_size

Arg [1]     : Int - Set batch size from default
Description : Controlled number of checksums to consider at once
Returntype  : Int - Batch size

=cut

sub batch_size {
  my ($self, $batch_size) = @_;
  $self->{_batch_size} = $batch_size if defined $batch_size;
  return $self->{_batch_size};
}

=head2 source_id

Description : Getter for the cached source ID used in the current run

=cut

sub source_id {
  my $self = shift;
  return $self->{_source_id};
}

=head2 process

Description : Trigger checksum mapping for a species
caller      : Pipeline

=cut

sub process {
  my ($self, $species_id) = @_;

  # Check has been removed here that originally tested for checksum_xrefs to decide
  # whether to run the comparison. It is removed to match the latest thinking in 
  # Bio::EnsEMBL::Production::Pipeline::Xrefs::UniParcMapping
  my $method = $self->method;
  my $results = $self->compare_checksums_versus_file();
  $self->upload($results, $species_id);

  return;
}

=head2 compare_checksums_versus_file

Description : Uses the "Core Info" protein_file to load batches of FASTA. The batches
              are then handed over to $self->method to decide how to checksum and compare
              them with those checksums stored in the Xref DB

Returntype  : ArrayRef of dbID/Uniparc ID pairings that can be uploaded as xrefs

=cut

sub compare_checksums_versus_file {
  my ($self) = @_;
  
  my $reader = Bio::SeqIO->new(-FILE => $self->core()->protein_file(), -FORMAT => 'fasta');
  my @results;
  my @tmp_list;
  my $batch_size = $self->batch_size();
  my $count = 0;

  my $checksum_method = $self->method;

  while ( my $sequence = $reader->next_seq() ) {
    push @tmp_list, $sequence;
    $count++;
    if( ($count % $batch_size) == 0) {
      my $res = $self->$checksum_method(\@tmp_list);
      push(@results, @{$res});
      $self->log_progress("Finished batch mapping of %d sequences\n", $batch_size);
      $count = 0;
      @tmp_list = ();
    }
  }
  
  #Final mapping if there were some left over
  if(@tmp_list) {
    $self->log_progress("Finishing progess\n");
    my $res = $self->$checksum_method(\@tmp_list);
    push(@results, @{$res});
    @tmp_list = ();
  }
  
  $reader->close();
  return \@results;
}

=head2 perform_mapping

Arg [1]     : Ref of sequences to checksum

Description : Overridable method that knows how to match checksums. The default matches
              Uniparc checksums against those generated by the Xref Pipeline
Returntype  : Listref of hashrefs, containing an Ensembl dbID, a Uniparc accession and
              an Ensembl object type matching Arg [3]
Caller      : compare_checksums_versus_file()

=cut

sub perform_mapping {
  my ($self, $sequences) = @_;
  
  my @final_results;
  my $source_id = $self->source_id;
  
  $self->_xref_helper->batch(
    -SQL => 'SELECT accession FROM checksum_xref WHERE checksum = ? AND source_id = ?', 
    -CALLBACK => sub {
      my ($sth) = @_;
      foreach my $sequence (@{$sequences}) {
        my $checksum = $self->md5_checksum($sequence);
        $sth->execute($checksum, $source_id);
        my $upi;
        while(my $row = $sth->fetchrow_arrayref()) {
          my ($local_upi) = @{$row};
          if(defined $upi) {
            throw sprintf('The sequence %s had a checksum of %s but this resulted in more than one UPI: [%s, %s]',
              $sequence->id(), $checksum, $upi, $local_upi
            );
          }
          $upi = $local_upi;
        }
        if (defined $upi){
          push @final_results, { id => $sequence->id(), upi => $upi };
        }
      }
      return; # exit callback
    }
  );
  
  return \@final_results;
}


=head2 upload

Arg [1]     : ArrayRef of checksum match records
              [ { id => 1, upi => 'UPI00000A' } ]
Arg [2]     : Int - the dbID for the species in case an override is required
Description : Take a list of checksum-matched Xrefs and upload them into
              the Xref database

=cut

sub upload {
  my ($self, $results, $species_id) = @_;
    
  $species_id = $self->species_id() unless defined $species_id;
   
  $self->log_progress('Deleting records from previous possible upload runs');
  $self->_delete_entries('object_xref');
  $self->_delete_entries('xref');
  
  $self->log_progress('Starting xref insertion');

  foreach my $entry (@{$results}) {
    my $upi = $entry->{upi};
    
    # No need to check if it exists, BaseAdaptor takes care of that
    my $dbID = $self->xref->add_xref({
      source_id => $self->source_id,
      acc       => $entry->{upi},
      label     => $entry->{upi},
      version   => 1,
      species_id => $species_id,
      info_type  => 'CHECKSUM'
    });
    $entry->{xref_id} = $dbID;    
  }
  
  $self->log_progress('Starting object_xref insertion');
  foreach my $entry (@{$results}) {
    $self->xref->add_object_xref({
      xref_id => $entry->{xref_id},
      ensembl_id => $entry->{id},
      object_type => $self->object_type,
      linkage_type => 'CHECKSUM',
      ox_status => 'DUMP_OUT'
    });
  }
  
  $self->log_progress('Finished insertions');
  
  return;
}

=head2 _delete_entries

Arg [1]     : String - table name that needs cleaning of Xref entries
Description : Deletes all entries of either Xrefs or ObjectXrefs for the current source

=cut

sub _delete_entries {
  my ($self, $table) = @_;
  $self->log_progress('Deleting entries from %s', $table);
  my $lookup = {
    xref => <<'SQL',
DELETE  x
FROM    xref x
WHERE   x.source_id    = ?
SQL
    object_xref => <<'SQL',
DELETE  ox
FROM    xref x,
        object_xref ox
WHERE   x.source_id    = ?
AND     ox.xref_id     = x.xref_id
SQL
  };
  
  my $sql = $lookup->{$table};
  throw "Cannot find delete SQL for the table $table" unless $sql;

  my $count = $self->_xref_helper()->execute_update(-SQL => $sql, -PARAMS => [$self->source_id]);
  my $type = ($count == 1) ? 'entry' : 'entries';
  $self->log_progress('Deleted %s %s from %s', $count, $type, $table);
  return;
}


=head2 _xref_helper

Description : Gets a core utils SqlHelper instance for the Xref DB
Returntype  : Bio::EnsEMBL::Utils::SqlHelper

=cut

sub _xref_helper {
  my $self = shift;
  return $self->xref()->dbc()->sql_helper();
}

=head2 _map_checksums

Description : Tests if there are any checksum mappings stored in the Xref DB.
              Checksum mappings indicate that we need to create Xrefs for them
Returntype  : Int - The number of rows in checksum_xref

=cut
sub _map_checksums {
  my ($self) = @_;
  my $source_id = $self->source_id;
  my $count = $self->_xref_helper->execute_single_result(
    -SQL => 'select count(*) from checksum_xref where source_id = ' . $source_id
  );
  return $count;
}

=head2 md5_checksum

Arg [1]     : String - ideally sequence
Description : Calculate the md5 for a given string
Returntype  : String - A hexademical checksum

=cut

sub md5_checksum {
  my ($self, $sequence) = @_;
  my $digest = Digest::MD5->new();
  $digest->add($sequence->seq());
  return uc $digest->hexdigest();
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
  return if (!$self->verbose);
  printf( STDERR "CHKSM==> %s\n", sprintf( $fmt, @params ) );
  return;
} ## end sub log_process

1;
