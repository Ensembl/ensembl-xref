
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

Bio::EnsEMBL::Xref::Parser::MGIParser

=head1 DESCRIPTION

A parser class to parse the MGI (official) source, 
creating a DIRECT xref between MGI accession and ensembl mouse gene stable id ENSMUSG*

-species = mus_musculus
-species_id = 10090
-data_uri = http://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt
-file_format = TSV
-columns = [accession symbol name position chrom ens_gene_stableid] ##ignore other columns

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::MGIParser->new(
    source_id  => 55,
    species_id => 10090,
    files      => ["MRK_ENSEMBL.rpt"],
    xref_dba   => $xref_dba  # xref db adaptor
  );

  $parser->run();
=cut

package Bio::EnsEMBL::Xref::Parser::MGIParser;

use strict;
use warnings;
use Carp;
use DBI;
use Text::CSV;

use parent qw( Bio::EnsEMBL::Xref::Parser );

=head2 run
  Description: Runs the MGIParser
  Return type: N/A
  Caller     : internal
=cut

sub run {

  my ( $self, $ref_arg ) = @_;
  my $source_id  = $self->{source_id};
  my $species_id = $self->{species_id};
  my $files      = $self->{files};
  my $verbose    = $self->{verbose} // 0;
  my $xref_dba   = $self->{xref_dba};

  my $file = shift @{$files};

  my $file_io = $xref_dba->get_filehandle($file);

  if ( !defined $file_io ) {
    croak "Could not open $file\n";
  }

  #synonyms, move this to SynonymAdaptor?!
  my $syn_hash = $xref_dba->get_ext_synonyms( "MGI", $xref_dba->dbi );

  #Init input file
  my $input_file = Text::CSV->new(
    {
      sep_char           => "\t",
      empty_is_undef     => 1,
      allow_loose_quotes => 1,
    }
  ) or croak "Cannot use file $file: " . Text::CSV->error_diag();

  $input_file->column_names(
    [qw(accession symbol name position chrom ens_gene_stableid)] );
  my $count     = 0;
  my $syn_count = 0;

  while ( my $data = $input_file->getline_hr($file_io) ) {
    my $acc     = $data->{'accession'};
    my $ensid   = $data->{'ens_gene_stableid'};
    my $xref_id = $xref_dba->add_xref(
      {
        acc        => $acc,
        version    => 0,
        label      => $data->{'symbol'},
        desc       => $data->{'name'},
        source_id  => $source_id,
        species_id => $species_id,
        info_type  => "DIRECT"
      }
    );

    $xref_dba->add_direct_xref( $xref_id, $ensid, "Gene", undef );
    if ( exists $syn_hash->{$acc} ) {
      foreach my $syn ( @{ $syn_hash->{$acc} } ) {
        $xref_dba->add_to_syn( $acc, $source_id, $syn, $species_id );
        $syn_count++;
      }
    }
    $count++;

  }
  $input_file->eof
    or croak "Error parsing file $file: " . $input_file->error_diag();
  $file_io->close();

  if ($verbose) {
    print "$count direct MGI xrefs added\n";
    print $syn_count. " synonyms added\n";
  }
  return 0;

}

1;
