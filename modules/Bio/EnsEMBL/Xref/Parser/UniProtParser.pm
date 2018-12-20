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



package Bio::EnsEMBL::Xref::Parser::UniProtParser;

use strict;
use warnings;

use Carp;
use Module::Load;

use parent qw( Bio::EnsEMBL::Xref::Parser );


my $DEFAULT_LOADER_BATCH_SIZE         = 1000;
my $DEFAULT_LOADER_CHECKPOINT_SECONDS = 300;

my %source_name_for_section = (
  'Swiss-Prot' => 'Uniprot/SWISSPROT',
  'TrEMBL'     => 'Uniprot/SPTREMBL',
);

# List of crossreference sources for which we want to create direct or
# dependent xrefs along with matching links.
# Note that this includes sources such as 'Uniprot_gn' and 'protein_id'
# which do not actually appear in UniProt-KB files but which are
# synthesised from other data; omitting them from this list simply
# means they will not be generated.
my $crossreference_sources_of_interest = [
  'ChEMBL',
  'EMBL',
  'Ensembl',
  'MEROPS',
  'PDB',
  'Uniprot_gn',
  'protein_id',
];

# To save memory and processing time, when we process a record in the
# extractor we only load into memory the fields we need. Moreover, the
# same list can later on be used to confirm that we have indeed
# encountered all the mandatory fields.
# Note that care must be taken when adding new prefixes to this list
# because some of them - for instance the Rx family of fields,
# describing publications - are not compatible with the current way of
# processing.
my $mandatory_prefixes_of_interest
  = [ 'ID', 'AC', 'DE', 'OX', 'PE', 'SQ', q{  }, ];
my $optional_prefixes_of_interest
  = [ 'GN', 'DR', 'RG', ];


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
                   - how many UniProt entries to process before
                     submitting xrefs to the database
                - loader_checkpoint_seconds
                   - the maximum amount of seconds (modulo the time
                     needed to process a single entry) xrefs are
                     allowed to stay in the buffer before being
                     submitted to the database regardless of the batch
                     size
  Example    : $uniprot_parser->run({ ... });
  Description: Extract UniProt Knowledgebase entries from text files
               downloaded from the UniProt Web site, then insert
               corresponding xrefs and links into the xref
               database. Both Swiss-Prot and TrEMBL files are
               supported.

               There will typically be multiple xrefs and links per
               one UniProt-KB entry:
                - basic xrefs created for each entry are
                  sequence-matched to Ensembl;
                - for entries with cross-references to Ensembl
                  we create additional direct xrefs as well
                  as corresponding translation_direct_xref links;
                - for entries with cross-references to specified
                  databases (see $crossreference_sources_of_interest
                  above), we create dependent xrefs as well as
                  corresponding dependent_xref links; - if requested,
                  special dependent xrefs and corresponding mappings
                  can be created for protein IDs (extracted from
                  ChEMBL and EMBL cross-references) and UniProt gene
                  names (from dedicated fields).

               UniProt-KB dat files are record-based, with individual
               fields identified by prefixes. For details, please see
               the UniProt Knowledgebase User Manual at
               https://web.expasy.org/docs/userman.html .

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
    = $self->{loader_batch_size} // $DEFAULT_LOADER_BATCH_SIZE;
  my $loader_checkpoint_seconds
    = $self->{loader_checkpoint_seconds} // $DEFAULT_LOADER_CHECKPOINT_SECONDS;

  # Try to control where ETL modules can come from, just in case
  # someone does something really weird with configuration options.
  my $extractor_class   = 'Bio::EnsEMBL::Xref::Parser::UniProtParser::';
  my $transformer_class = 'Bio::EnsEMBL::Xref::Parser::UniProtParser::';
  my $loader_class      = 'Bio::EnsEMBL::Xref::Parser::UniProtParser::';

  $extractor_class   .= $self->{extractor_class}   // 'Extractor';
  $transformer_class .= $self->{transformer_class} // 'Transformer';
  $loader_class      .= $self->{loader_class}      // 'Loader';

  load $extractor_class;
  load $transformer_class;
  load $loader_class;

  my $extractor = $extractor_class->new({
    'file_names' => $files,
    'mandatory_prefixes' => $mandatory_prefixes_of_interest,
    'optional_prefixes'  => $optional_prefixes_of_interest,
    'species_id'         => $species_id,
    'xref_dba'           => $xref_dba,
  });
  my $transformer = $transformer_class->new({
    'accepted_crossreference_sources' => $crossreference_sources_of_interest,
    'default_direct_xref_type'        => 'Translation',  # just in case
    'species_id'                      => $species_id,
    'xref_dba'                        => $xref_dba,
  });
  my $loader = $loader_class->new({
    'batch_size'         => $loader_batch_size,
    'checkpoint_seconds' => $loader_checkpoint_seconds,
    'xref_dba'           => $xref_dba,
  });

  # Generate a map of existing dependent-xref links, which will be
  # used to prevent insertion of duplicates
  my $source_id_map = $transformer->get_source_id_map();
  while ( my ( $section, $pri_ref ) = each %{ $source_id_map } ) {
    while ( my ( $priority, $source_id ) = each %{ $pri_ref } ) {

      if ( $priority ne 'direct' ) {
        $loader->prepare_source_for_dependent_xrefs( $source_id );
      }

      # Seeing as we are already looping over all sources.
      # FIXME: do we really care about the file name here?
      if ( $verbose ) {
        print "$section $priority source id for $files->[0]: $source_id\n";
      }

    }
  }

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

  # Extract release numbers from the release file, if provided
  if ( defined $release_file ) {
    my $release_strings = $self->_get_release_strings_from_file( $release_file,
                                                                 $verbose );
    $self->_set_release_strings_on_uniprot_sources( $source_id_map,
                                                    $release_strings );
  }

  return;
}


# Extract Swiss-Prot and TrEMBL release info from the release file
sub _get_release_strings_from_file {
  my ( $self, $release_file_name, $verbose ) = @_;

  my $xref_dba = $self->{xref_dba};

  my $release_io = $xref_dba->get_filehandle( $release_file_name );
  my $release_strings = {};

  while ( my $line = $release_io->getline() ) {
    my ( $release, $section )
      = ( $line =~ m{
                      \A
                      ( UniProtKB/
                      (
                        Swiss-Prot
                      |
                        TrEMBL
                      )
                      \s+
                      Release
                      \s+
                      [^\n]+ )
                  }msx );
    if ( defined $release ) {
      $release_strings->{ $source_name_for_section{$section} } = $release;
      if ( $verbose ) {
        print "$section release is '$release'\n";
      }
    }
  }

  $release_io->close();

  return $release_strings;
}

sub _set_release_strings_on_uniprot_sources {
  my ( $self, $source_id_map, $release_strings ) = @_;

  my $xref_dba = $self->{xref_dba};

  foreach my $source ( keys %{ $source_id_map } ) {
    # Priority names are not important here, we only need source IDs
    foreach my $source_id ( values %{ $source_id_map->{$source} } ) {
      $xref_dba->set_release( $source_id, $release_strings->{$source} );
    }
  }

  return;
}

1;
