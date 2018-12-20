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



package Bio::EnsEMBL::Xref::Parser::miRBaseParser;

use strict;
use warnings;

use Carp;
use Module::Load;

use parent qw( Bio::EnsEMBL::Xref::Parser );


my $DEFAULT_LOADER_BATCH_SIZE         = 1000;
my $DEFAULT_LOADER_CHECKPOINT_SECONDS = 300;

# To save memory and processing time, when we process a record in the
# extractor we only load into memory the fields we need. Moreover, the
# same list can later on be used to confirm that we have indeed
# encountered all the mandatory fields.
# Note that care must be taken when adding new prefixes to this list
# because some of them - for instance the Rx family of fields,
# describing publications - are not compatible with the current way of
# processing.
my $mandatory_prefixes_of_interest
  = [ 'ID', 'AC', 'DE', 'SQ', q{  }, ];


=head2 run

  Arg [1]    : HashRef arguments for the parser:
                - standard list of arguments from ParseSource
                - extractor
                   - name of class used to instantiate the extractor
                - transformer
                   - name of class used to instantiate the transformer
                - loader
                   - name of class used to instantiate the loader
                - loader_batch_size
                   - how many input records to process before
                     submitting xrefs to the database
                - loader_checkpoint_seconds
                   - the maximum amount of seconds (modulo the time
                     needed to process a single entry) xrefs are
                     allowed to stay in the buffer before being
                     submitted to the database regardless of the batch
                     size
  Example    : $mirbase_parser->run({ ... });
  Description: Extract miRBase entries from text files downloaded from
               the miRBase Web site, then insert corresponding xrefs
               into the xref database.

               For each miRBase stem loop we produce a single xref,
               which will later on be sequence-matched to
               Ensembl. Creation of xrefs for mature miRNAs is
               not supported yet.

               miRBase dat files have the same structure as UniProt
               Knowledge Base dumps, which is why we use the same
               parsing engine.

  Return type: none
  Exceptions : throws on all processing errors
  Caller     : ParseSource in the xref pipeline
  Status     : Stable

=cut

sub run {
  my ( $self ) = @_;

  my $general_source_id = $self->{source_id};
  my $species_id        = $self->{species_id};
  my $species_name      = $self->{species};
  my $files             = $self->{files};
  my $release_file      = $self->{rel_file};
  my $verbose           = $self->{verbose} // 0;
  my $xref_dba          = $self->{xref_dba};

  my $loader_batch_size
    = $self->{loader_batch_size}         // $DEFAULT_LOADER_BATCH_SIZE;
  my $loader_checkpoint_seconds
    = $self->{loader_checkpoint_seconds} // $DEFAULT_LOADER_CHECKPOINT_SECONDS;

  # Try to control where ETL modules can come from, just in case
  # someone does something really weird with configuration options.
  my $extractor_class   = 'Bio::EnsEMBL::Xref::Parser::UniProtParser::';
  my $transformer_class = 'Bio::EnsEMBL::Xref::Parser::UniProtParser::';
  my $loader_class      = 'Bio::EnsEMBL::Xref::Parser::UniProtParser::';

  $extractor_class   .= $self->{extractor_class}   // 'Extractor::miRBase';
  $transformer_class .= $self->{transformer_class} // 'Transformer::miRBase';
  $loader_class      .= $self->{loader_class}      // 'Loader';

  load $extractor_class;
  load $transformer_class;
  load $loader_class;

  my $extractor = $extractor_class->new({
    'file_names'         => $files,
    'mandatory_prefixes' => $mandatory_prefixes_of_interest,
    'species_name'       => $species_name,
    'xref_dba'           => $xref_dba,
  });
  my $transformer = $transformer_class->new({
    'general_source_id' => $general_source_id,
    'species_id'        => $species_id,
    'xref_dba'          => $xref_dba,
  });
  my $loader = $loader_class->new({
    'batch_size'         => $loader_batch_size,
    'checkpoint_seconds' => $loader_checkpoint_seconds,
    'xref_dba'           => $xref_dba,
  });

 RECORD:
  while ( $extractor->get_uniprot_record() ) {

    my $extracted_record = $extractor->extract();
    if ( $extracted_record eq 'SKIP' ) {
      next RECORD;
    }

    my $transformed_data
      = $transformer->transform( $extracted_record );

    $loader->load( $transformed_data );

  }

  $extractor->finish();
  $transformer->finish();
  $loader->finish();

  return;
}


1;
