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

Bio::EnsEMBL::Xref::Parser::ReactomeDirectParser

=head1 DESCRIPTION

A parser class to parse the Reactome file for direct mappings to Ensembl..

-data_uri = http://www.reactome.org/download/current/Ensembl2Reactome_All_Levels.txt
-file_format = TSV
-columns = [ensembl_stable_id reactome_acc reactome_url description status species]


=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::ReactomeDirectParser->new(
    source_id  => 85,
    species_id => 9606,
    species_name => 'homo_sapiens',
    files      => ['Ensembl2Reactome_All_Levels.txt'],
    xref_dba   => $xref_dba
  );

  $parser->run();
=cut


package Bio::EnsEMBL::Xref::Parser::ReactomeDirectParser;

use strict;
use warnings;
use Carp;
use Text::CSV;

use parent qw( Bio::EnsEMBL::Xref::Parser );

=head2 run
  Description: Runs the ReactomeDirectParser
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

  my $gene_reactome_source_id =
    $xref_dba->get_source_id_for_source_name( "reactome_gene" );
  my $transcript_reactome_source_id =
    $xref_dba->get_source_id_for_source_name( "reactome_transcript" );

#e.g.
#ENSG00000000419 R-HSA-162699    https://reactome.org/PathwayBrowser/#/R-HSA-162699      Synthesis of dolichyl-phosphate mannose TAS     Homo sapiens
#ENSG00000000419 R-HSA-163125    https://reactome.org/PathwayBrowser/#/R-HSA-163125      Post-translational modification: synthesis of GPI-anchored proteins     TAS     Homo sapiens
#ENSG00000000419 R-HSA-1643685   https://reactome.org/PathwayBrowser/#/R-HSA-1643685     Disease TAS     Homo sapiens


  my $file_io = $xref_dba->get_filehandle($file);

  my $input_file = Text::CSV->new({
    sep_char       => "\t",
    binary         => 1,
    empty_is_undef => 1
  }) or confess "Cannot use file '$file': " . Text::CSV->error_diag();

  $input_file->column_names( [ 'ensembl_stable_id', 'accession', 'url', 'description', 'status', 'species'] );

  my ($reactome_source_id, $type, $species);
  while ( my $data = $input_file->getline_hr( $file_io ) ) {
    $species = $data->{'species'};
    $species =~ s/\s/_/x;
    $species = lc($species);
    # Skipping non-standard characters which cannot be stored in MySQL
    if ($data->{'description'} !~ /^[A-Za-z0-9_,\(\)\/\-\.:\+'&;"\/\?%>\s\[\]]+$/x) { next; }
    if ( $species ne $species_name ) {
      next;
    }
    if ( $data->{'ensembl_stable_id'} =~ /G[0-9]*$/x) {
      $type = 'gene';
      $reactome_source_id = $gene_reactome_source_id;
    } elsif ( $data->{'ensembl_stable_id'} =~ /T[0-9]*$/x) {
      $type = 'transcript';
      $reactome_source_id = $transcript_reactome_source_id;
    } else {
    # Stable ID cannot be identified, skip the entry
      next;
    }
    $xref_dba->add_to_direct_xrefs({
        acc        => $data->{'accession'},
        label      => $data->{'accession'},
        desc       => $data->{'description'},
        stable_id  => $data->{'ensembl_stable_id'},
        type       => $type,
        source_id  => $reactome_source_id,
        species_id => $species_id,
        info_type  => "DIRECT"
    });
  }

  $input_file->eof or confess "Error parsing file $file: " . $input_file->error_diag();
  $file_io->close();

  return 0;
}

1;
