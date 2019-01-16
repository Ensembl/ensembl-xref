
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

Bio::EnsEMBL::Xref::Parser::ZFIN::UniprotParser

=head1 DESCRIPTION

A specialised parser class to parse the ZFIN source Uniprot related file

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::ZFIN::UniprotParser->new(
    source_id  => 149,
    species_id => 7955,
    files      => [ 'uniprot.txt' ],
    xref_dba   => $xref_dba
  );

  $parser->run();

=cut

package Bio::EnsEMBL::Xref::Parser::ZFIN::UniprotParser;

use strict;
use warnings;

use Carp;

use Text::CSV;

use Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;

use parent qw( Bio::EnsEMBL::Xref::Parser::ZFIN );

=head2 run

  Arg []     : None
  Description: Add dependent xrefs from ZFIN to the xref database
  Return type: Int; 0 upon success 
  Caller     : An analysis step in the Xref pipeline

=cut

sub run {

  my ($self) = @_;

  my $source_id  = $self->{source_id};
  my $species_id = $self->{species_id};
  my $files      = $self->{files};
  my $xref_dba   = $self->{xref_dba};
  my $verbose    = $self->{verbose} // 0;

  unless ( defined $source_id and
           defined $species_id and
           defined $files )
  {
    confess "Need to pass source_id, species_id and files";
  }

  my $file         = shift @{$files};
  my $swissprot_io = $xref_dba->get_filehandle($file);
  confess "Could not open ZFIN uniprot/swissprot $file"
    unless defined $swissprot_io;

  my $swissprot_csv = Text::CSV->new(
            { sep_char => "\t", empty_is_undef => 1, strict => 1, } ) or
    confess "Could not use swissprot file $file: " .
    Text::CSV->error_diag();

  $swissprot_csv->column_names( [ 'zfin', 'so', 'label', 'acc', ] );

  my ( $spcount, $mismatch ) = ( 0, 0 );

  my $acc2desc = $self->description();
  my (%swiss) =
    %{ $xref_dba->get_valid_codes( "uniprot/swissprot", $species_id ) };

  while ( my $swissprot_line =
          $swissprot_csv->getline_hr($swissprot_io) )
  {
    $mismatch++ and next
      unless defined $swiss{ $swissprot_line->{'acc'} };

    foreach my $xref_id ( @{ $swiss{ $swissprot_line->{'acc'} } } ) {
      $xref_dba->add_dependent_xref(
                     { master_xref_id => $xref_id,
                       acc            => $swissprot_line->{'zfin'},
                       label          => $swissprot_line->{'label'},
                       desc => $acc2desc->{ $swissprot_line->{'zfin'} },
                       source_id  => $source_id,
                       species_id => $species_id } );
      $spcount++;
    }
  }

  $swissprot_io->close();

  print "\t$spcount xrefs from UniProt ($mismatch mismatches)\n"
    if $verbose;

  return 0;    # success

} ## end sub run

1;
