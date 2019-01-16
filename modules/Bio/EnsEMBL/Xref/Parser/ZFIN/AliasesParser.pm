
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

Bio::EnsEMBL::Xref::Parser::ZFIN::AliasesParser

=head1 DESCRIPTION

A specialised parser class to parse the ZFIN source aliases related file

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::ZFIN::AliasesParser->new(
    source_id  => 149,
    species_id => 7955,
    files      => [ 'aliases.txt' ],
    xref_dba   => $xref_dba
  );

  $parser->run();

=cut

package Bio::EnsEMBL::Xref::Parser::ZFIN::AliasesParser;

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

  my $file = shift @{$files};

  my $zfin_io = $xref_dba->get_filehandle($file);
  confess "Could not open zfin aliases file $file"
    unless defined $zfin_io;

  my $zfin_csv = Text::CSV->new(
            { sep_char => "\t", empty_is_undef => 1, strict => 1, } ) or
    confess "could not use zfin file $file: " . Text::CSV->error_diag();

  $zfin_csv->column_names(qw( acc cur_name cur_symbol syn so ));

  my $syncount = 0;
  my (%zfin) = %{ $xref_dba->get_valid_codes( "zfin", $species_id ) };

  while ( my $zfin_line = $zfin_csv->getline_hr($zfin_io) ) {
    my ( $acc, $syn ) = @{$zfin_line}{qw( acc syn )};

    next unless defined $zfin{$acc} and $syn;
    $xref_dba->add_to_syn_for_mult_sources( $acc, $self->{source_ids},
                                            $syn, $species_id );
    $syncount++;
  }

  $zfin_io->close();

  print "\t$syncount synonyms loaded\n" if $verbose;

  return 0;    # success

} ## end sub run

1;
