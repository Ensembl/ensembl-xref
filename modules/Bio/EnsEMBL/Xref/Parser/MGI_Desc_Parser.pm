
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

Bio::EnsEMBL::Xref::Parser::MGI_Desc_Parser

=head1 DESCRIPTION

A parser class to parse the MGI (descriptions) source. Creates 'MISC' xref using MGI accession with description and
also creates the synonyms extracted from the pipe seperated synonym_field

-species = mus_musculus
-species_id = 10090
-data_uri = http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt
-file_format = TSV
-columns = [accession chromosome position start end strand label status marker marker_type feature_type synonym_field]

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::MGI_Desc_Parser->new(
    source_id  => 58,
    species_id => 10090,
    files      => ["MRK_List2.rpt"],
    xref_dba   => $xref_dba  # xref db adaptor
  );

  $parser->run();
=cut

package Bio::EnsEMBL::Xref::Parser::MGI_Desc_Parser;

use strict;
use warnings;
use Carp;
use Text::CSV;

use parent qw( Bio::EnsEMBL::Xref::Parser );

=head2
The run method does the actual parsing and creation of xrefs and synonyms.
Parser gets initialized as noted above and run is called from
Bio::EnsEMBL::Production::Pipeline::Xrefs::ParseSource

my $parser = Bio::EnsEMBL::Xref::Parser::MGI_Desc_Parser->new(..)
$parser->run();
=cut

sub run {

  my ( $self, $ref_arg ) = @_;
  my $source_id  = $self->{source_id};
  my $species_id = $self->{species_id};
  my $files      = $self->{files};
  my $verbose    = $self->{verbose} // 0;
  my $xref_dba   = $self->{xref_dba};

  my $file = shift @{$files};

  my $mgi_io = $xref_dba->get_filehandle($file);

  if ( !defined $mgi_io ) {
    confess "Could not open $file\n";
  }

  my $xref_count = 0;
  my $syn_count  = 0;
  my %acc_to_xref;

  my $input_file = Text::CSV->new(
     { sep_char => "\t", empty_is_undef => 1, allow_loose_quotes => 1, }
    ) or
    confess "Cannot use file $file: " . Text::CSV->error_diag();

  # expected columns
  my @expected_columns =
    qw(accession chromosome position start end strand label status marker marker_type feature_type synonym_field);

  # read header
  my $header = $input_file->getline($mgi_io);
  if ( scalar @{$header} != scalar @expected_columns ) {
    confess "input file $file has an incorrect number of columns";
  }
  $input_file->column_names( \@expected_columns );
  while ( my $data = $input_file->getline_hr($mgi_io) ) {
    my $accession = $data->{'accession'};
    my $marker =
      defined( $data->{'marker'} ) ? $data->{'marker'} : undef;
    $acc_to_xref{$accession} =
      $xref_dba->add_xref(
                           { acc        => $accession,
                             label      => $data->{'label'},
                             desc       => $marker,
                             source_id  => $source_id,
                             species_id => $species_id,
                             info_type  => "MISC" } );
    if ( $verbose && !$marker ) {
      print "$accession has no description\n";
    }
    $xref_count++;
    my @synonyms;
    if ( defined( $acc_to_xref{$accession} ) ) {
      @synonyms = split( /\|/x, $data->{'synonym_field'} )
        if ( $data->{'synonym_field'} );
      foreach my $syn (@synonyms) {
        $xref_dba->add_synonym( $acc_to_xref{$accession}, $syn );
        $syn_count++;
      }
    }
  } ## end while ( my $data = $input_file...)
  $mgi_io->eof or
    confess "Error parsing file $file: " . $input_file->error_diag();
  $mgi_io->close();

  if ($verbose) {
    print "$xref_count MGI Description Xrefs added\n";
    print "$syn_count synonyms added\n";
  }

  return 0;    #successful
} ## end sub run

1;

