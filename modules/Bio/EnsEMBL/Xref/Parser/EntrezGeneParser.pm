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

Bio::EnsEMBL::Xref::Parser::EntrezGeneParser

=head1 DESCRIPTION

A parser class to parse the EntrezGene source file.

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::XenopusJamboreeParser->new(
    source_id  => ,
    species_id => ,
    files      => ["entrezgene.txt"],
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

sub run {
  my ( $self, $ref_arg ) = @_;
  my $source_id  = $self->{source_id};
  my $species_id = $self->{species_id};
  my $species_name = $self->{species_name};
  my $files      = $self->{files};
  my $verbose    = $self->{verbose} // 0;
  my $xref_dba   = $self->{xref_dba};

  my $file = shift @{$files};

  my $wiki_source_id =
    $self->get_source_id_for_source_name( "WikiGene" );

  my %species_tax_id =
    %{ $self->get_taxonomy_from_species_id( $species_id ) };
  $species_tax_id{$species_id} = $species_id if defined $species_name;
  
  my $eg_io = $xref_dba->get_filehandle($file);
  croak "Could not open $file\n" unless defined $eg_io;
  
  my $input_file = Text::CSV->new(
    {
      sep_char       => "\t",
      empty_is_undef => 1,
      allow_loose_quotes => 1
    }
  ) or croak "Cannot use file $file: " . Text::CSV->error_diag();

  $input_file->column_names( @{ $input_file->getline($eg_io) } );

  # read data and load xrefs
  my $xref_count = 0;
  my $syn_count  = 0;
  my %seen; # record already processed xrefs

  while ( my $data = $input_file->getline_hr($eg_io) ) {
    next unless exists $species_tax_id{ $data->{'#tax_id'} };

    my $acc = $data->{'GeneID'};
    next if $seen{$acc};

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
                       info_type  => "DEPENDENT" }
    );
    $xref_count++;

    my (@syn) = split( /\|/, $data->{'Synonyms'} );
    foreach my $synonym (@syn) {
      if ( $synonym ne "-" ) {
        $self->add_to_syn( $acc, $source_id, $synonym, $species_id );
        $syn_count++;
      }
    }

    $seen{$acc} = 1;
  }

  $input_file->eof or
    croak "Error parsing file $file: " . $input_file->error_diag();
  $eg_io->close();

  print $xref_count.
    " EntrezGene Xrefs added with $syn_count synonyms\n"
    if $verbose;

  return 0; # success
}

1;
