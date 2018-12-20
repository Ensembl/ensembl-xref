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

Bio::EnsEMBL::Xref::Parser::VGNCParser

=head1 DESCRIPTION

A parser class to parse the VGNC source.
VGNC is the official naming source for some vertebrates species
(currently chimp, cow, dog and horse).

-data_uri = ftp://ftp.ebi.ac.uk/pub/databases/genenames/vgnc/tsv/vgnc_gene_set_All.txt.gz
-file_format = TSV
-columns = [
    taxon_id
    vgnc_id
    symbol
    name
    locus_group
    locus_type
    status
    location
    location_sortable:
    alias_symbol
    alias_name
    prev_symbol
    prev_name
    gene_family
    gene_family_id
    date_approved_reserved
    date_symbol_changed
    date_name_changed
    date_modified
    entrez_id
    ensembl_gene_id
    uniprot_ids
  ]

Only columns listed in @required_columns are mandatory.

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::VGNCParser->new(
    source_id  => 144,
    species_id => 9598,
    files      => ['VGNC/vgnc_gene_set_All.txt.gz'],
    xref_dba   => $xref_dba
  );

  $parser->run();
=cut

package Bio::EnsEMBL::Xref::Parser::VGNCParser;

use strict;
use warnings;
use Carp;
use Text::CSV;

use parent qw( Bio::EnsEMBL::Xref::Parser );

=head2 run
  Description: Runs the VGNCParser
  Return type: N/A
  Caller     : internal
=cut

sub run {
  my ( $self ) = @_;

  my $source_id    = $self->{source_id};
  my $species_id   = $self->{species_id};
  my $files        = $self->{files};
  my $xref_dba     = $self->{xref_dba};
  my $verbose      = $self->{verbose} // 0;


  if ( (!defined $source_id) || (!defined $species_id) || (!defined $files) ) {
    confess "Need to pass source_id, species_id and files as pairs";
  }

  my $file = shift @{$files};

  my $count = 0;

  my $file_io = $xref_dba->get_filehandle($file);

  if ( !defined $file_io ) {
    confess "Can't open VGNC file $file\n";
  }

  my $input_file = Text::CSV->new({
    sep_char       => "\t",
    empty_is_undef => 1
  }) or confess "Cannot use file $file: " . Text::CSV->error_diag();

  # header must contain these columns
  my @required_columns = qw(
    taxon_id
    ensembl_gene_id
    vgnc_id
    symbol
    name
    alias_symbol
    prev_symbol
  );

  # get header columns
  my @columns = @{ $input_file->getline( $file_io ) };

  # die if some required_column is not in columns
  foreach my $colname (@required_columns) {
    if ( !grep { /$colname/xms } @columns ) {
      confess "Can't find required column $colname in VGNC file $file\n";
    }
  }

  $input_file->column_names( @columns );

  while ( my $data = $input_file->getline_hr( $file_io ) ) {

    # skip data for other species
    next if ( $species_id != $data->{'taxon_id'} );

    if ( $data->{'ensembl_gene_id'} ) {              # Ensembl direct xref
      $xref_dba->add_to_direct_xrefs({
        stable_id  => $data->{'ensembl_gene_id'},
        type       => 'gene',
        acc        => $data->{'vgnc_id'},
        label      => $data->{'symbol'},
        desc       => $data->{'name'},
        source_id  => $source_id,
        species_id => $species_id
      });

      $xref_dba->add_synonyms_for_hgnc_vgnc({
        source_id  => $source_id,
        name       => $data->{'vgnc_id'},
        species_id => $species_id,
        dead       => $data->{'alias_symbol'},
        alias      => $data->{'prev_symbol'}
      });

      $count++;
    }

  }

  $input_file->eof or confess "Error parsing file $file: " . $input_file->error_diag();
  $file_io->close();

  if ($verbose) {
    print "Loaded a total of $count VGNC xrefs\n";
  }

  return 0; # successful
}


1;
