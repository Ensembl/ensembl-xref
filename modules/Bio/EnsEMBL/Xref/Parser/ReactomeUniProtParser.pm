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

Bio::EnsEMBL::Xref::Parser::ReactomeUniProtParser

=head1 DESCRIPTION

A parser class to parse the Reactome file for UniProt mappings.

-data_uri = http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt
-file_format = TSV
-columns = [uniprot_acc reactome_acc reactome_url description status species]


=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::ReactomeUniProtParser->new(
    source_id  => 85,
    species_id => 9696,
    files      => ['UniProt2Reactome_All_Levels.txt'],
    xref_dba   => $xref_dba
  );

  $parser->run();
=cut


package Bio::EnsEMBL::Xref::Parser::ReactomeUniProtParser;

use strict;
use warnings;
use Carp;
use Text::CSV;

use parent qw( Bio::EnsEMBL::Xref::Parser );

=head2 run
  Description: Runs the ReactomeUniProtParser
  Return type: N/A
  Caller     : internal
=cut

sub run {
  my ( $self ) = @_;

  my $source_id    = $self->{source_id};
  my $species_id   = $self->{species_id};
  my $species_name = $self->{species};
  my $files        = $self->{files};
  my $xref_dba     = $self->{xref_dba};

  if ( (!defined $source_id) || (!defined $species_id) || (!defined $files) ) {
    confess "Need to pass source_id, species_id and files as pairs";
  }

  my $file = shift @{$files};

  my $reactome_source_id =
    $xref_dba->get_source_id_for_source_name( "reactome_translation", "uniprot" );

#e.g.
#A0A075B6P5  R-HSA-109582  https://reactome.org/PathwayBrowser/#/R-HSA-109582  Hemostasis  TAS  Homo sapiens
#A0A075B6P5  R-HSA-1280218  https://reactome.org/PathwayBrowser/#/R-HSA-1280218  Adaptive Immune System  TAS  Homo sapiens
#A0A075B6P5  R-HSA-1280218  https://reactome.org/PathwayBrowser/#/R-HSA-1280218  Adaptive Immune System  IEA  Homo sapiens


  my $file_io = $xref_dba->get_filehandle($file);

  if ( !defined $file_io ) {
    confess "Can't open UniProt_Reactome file '$file'\n";
  }

  my $input_file = Text::CSV->new({
    sep_char       => "\t",
    binary         => 1,
    empty_is_undef => 1
  }) or confess "Cannot use file '$file': " . Text::CSV->error_diag();

  $input_file->column_names( [ 'uniprot_accession', 'accession', 'url', 'description', 'status', 'species'] );
  my (%uniprot) = 
    %{ $xref_dba->get_valid_codes( "uniprot/", $species_id ) };

  if (!%uniprot) {
    confess "No UniProt xrefs found, cannot parse dependent Reactome mappings";
  }

  my $species;
  while ( my $data = $input_file->getline_hr( $file_io ) ) {
    # Skipping non-standard characters which cannot be stored in MySQL
    if ($data->{'description'} !~ /^[A-Za-z0-9_,\(\)\/\-\.:\+'&;"\/\?%>\s\[\]]+$/x) { next; }
    $species = $data->{'species'};
    $species =~ s/\s/_/x;
    $species = lc($species);
    if ( $species ne $species_name ) {
      next;
    }
    if ( defined $uniprot{$data->{'uniprot_accession'}} ) {
      foreach my $master_xref_id (@{$uniprot{$data->{'uniprot_accession'}}} ) {
        $xref_dba->add_dependent_xref({
          acc            => $data->{'accession'},
          label          => $data->{'accession'},
          desc           => $data->{'description'},
          master_xref_id => $master_xref_id, 
          source_id      => $source_id,
          species_id     => $species_id,
          info_type      => "DEPENDENT"
        });
      }
    }
  }

  $input_file->eof or confess "Error parsing file $file: " . $input_file->error_diag();
  $file_io->close();

  return 0;
}

1;
