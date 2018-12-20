
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

Bio::EnsEMBL::Xref::Parser::XenopusJamboreeParser

=head1 DESCRIPTION

A parser class to parse the Xenbase source file.

-species = xenopus_tropicalis
-species_id = 8364
-data_uri = ftp://ftp.xenbase.org/pub/GenePageReports/GenePageEnsemblModelMapping.txt
-file_format = TSV
-columns = [acc label desc stable_id]


=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::XenopusJamboreeParser->new(
    source_id  => 150,
    species_id => 8364,
    files      => ["xenopusjamboree.txt"],
    xref_dba   => $xref_dba
  );

  $parser->run();
=cut

package Bio::EnsEMBL::Xref::Parser::XenopusJamboreeParser;

use strict;
use warnings;
use Carp;
use Text::CSV;

use parent qw( Bio::EnsEMBL::Xref::Parser );

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

  my $input_file =
    Text::CSV->new( { sep_char => "\t", empty_is_undef => 1, } ) or
    croak "Cannot use file $file: " . Text::CSV->error_diag();

  $input_file->column_names( [qw(acc label desc stable_id)] );
  my $count = 0;
  while ( my $data = $input_file->getline_hr($file_io) ) {

    my $desc;
    $desc = $self->parse_description( $data->{'desc'} )
      if ( defined( $data->{'desc'} ) );

    $data->{'label'} = $data->{'acc'}
      if ( $data->{'label'} eq "unnamed" );

    $xref_dba->add_to_direct_xrefs(
                                    { stable_id => $data->{'stable_id'},
                                      type      => 'gene',
                                      acc       => $data->{'acc'},
                                      label     => $data->{'label'},
                                      desc      => $desc,
                                      source_id => $source_id,
                                      species_id => $species_id } );
    $count++;
  }

  $input_file->eof or
    croak "Error parsing file $file: " . $input_file->error_diag();
  $file_io->close();

  print $count . " XenopusJamboreeParser xrefs succesfully parsed\n"
    if ($verbose);

  return 0;
} ## end sub run

sub parse_description {
  my ( $self, $desc ) = @_;

  # Remove some provenance information encoded in the description
  $desc =~ s/\[.*\]//xms;

  # Remove labels of type 5 of 14 from the description
  $desc =~ s/,\s+[0-9]+\s+of\s+[0-9]+//xms;
  return $desc;
}

1;
