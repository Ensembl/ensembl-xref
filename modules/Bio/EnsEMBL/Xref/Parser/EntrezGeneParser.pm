
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

=head1 NAME

Bio::EnsEMBL::Xref::Parser::EntrezGeneParser

=head1 DESCRIPTION

This parser will read and create dependent xrefs from a simple 
comma-delimited file downloaded from the EntrezGene database.

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::EntrezGeneParser->new(
    source_id  => 11,
    species_id => 9606,
    files      => [ "gene_info.gz" ],
    xref_dba   => $xref_dba
  );

  $parser->run();

=cut

package Bio::EnsEMBL::Xref::Parser::EntrezGeneParser;

use strict;
use warnings;

use Carp;
use Text::CSV;

use parent qw( Bio::EnsEMBL::Xref::Parser );

=head2 run

  Arg []     : None
  Description: Add dependent xrefs from EntrezGene to the xref database
  Return type: Int; 0 upon success 
  Caller     : An analysis step in the Xref pipeline

=cut

sub run {

  my ($self)       = @_;
  my $source_id    = $self->{source_id};
  my $species_id   = $self->{species_id};
  my $species_name = $self->{species};
  my $files        = $self->{files};
  my $verbose      = $self->{verbose} // 0;
  my $xref_dba     = $self->{xref_dba};

  unless ( defined $source_id and
           defined $species_id and
           defined $files )
  {
    confess "Need to pass source_id, species_id and files";
  }

  my $file = shift @{$files};

  my $wiki_source_id =
    $xref_dba->get_source_id_for_source_name( "WikiGene", undef );

  my $eg_io = $xref_dba->get_filehandle($file);
  confess "Could not open $file" unless defined $eg_io;

  my $input_file = Text::CSV->new(
    { sep_char => "\t", empty_is_undef => 1, allow_loose_quotes => 1 } )
    or
    confess "Cannot use file $file: " . Text::CSV->error_diag();

  # process header
  $input_file->column_names( @{ $input_file->getline($eg_io) } );

  # read data and load xrefs
  my $xref_count = 0;
  my $syn_count  = 0;
  my %seen;    # record already processed xrefs

  while ( my $data = $input_file->getline_hr($eg_io) ) {
    # species_id corresponds to the species taxonomy id, see:
    # https://github.com/Ensembl/ensembl-xref/pull/31#issuecomment-445838474
    next unless $data->{'#tax_id'} eq $species_id;

    my $acc = $data->{'GeneID'};
    next if exists $seen{$acc};

    my $symbol = $data->{'Symbol'};
    my $desc   = $data->{'description'};

    $xref_dba->add_xref(
                         { acc        => $acc,
                           label      => $symbol,
                           desc       => $desc,
                           source_id  => $source_id,
                           species_id => $species_id,
                           info_type  => "DEPENDENT" } );

    $xref_dba->add_xref(
                         { acc        => $acc,
                           label      => $symbol,
                           desc       => $desc,
                           source_id  => $wiki_source_id,
                           species_id => $species_id,
                           info_type  => "DEPENDENT" } );
    $xref_count++;

    my (@syn) = split( /\|/, $data->{'Synonyms'} );
    foreach my $synonym (@syn) {
      if ( $synonym ne "-" ) {
        $xref_dba->add_to_syn( $acc, $source_id, $synonym,
                               $species_id );
        $syn_count++;
      }
    }

    $seen{$acc} = 1;
  } ## end while ( my $data = $input_file...)

  $input_file->eof or
    confess "Error parsing file $file, should be EOF: " . $input_file->error_diag();
  $eg_io->close();

  print $xref_count.
    " EntrezGene Xrefs added with $syn_count synonyms\n"
    if $verbose;

  return 0;    # success
} ## end sub run

1;
