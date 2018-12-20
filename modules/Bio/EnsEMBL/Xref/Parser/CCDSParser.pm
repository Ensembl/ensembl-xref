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

Bio::EnsEMBL::Xref::Parser::CCDSParser

=head1 DESCRIPTION

A parser class to parse the CCDS source. Fetches the CCDS to Ensembl mappings
from the core database. The CCDS entries are added as xrefs in the xref
database as well as in the direct xref table.

-species = only used by human and mouse

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::CCDSParser->new(
    source_id  => 1,
    species_id => 9606,
    xref_dba   => $xref_dba # xref db adaptor
    dba        => $dba   # core db adaptor
  );

  $parser->run();
=cut


package Bio::EnsEMBL::Xref::Parser::CCDSParser;

use strict;
use warnings;
use Carp;

use parent qw( Bio::EnsEMBL::Xref::Parser );

=head2
The run method does the actual parsing and creation of direct xrefs.
Parser gets initialized as noted above and run is called from
Bio::EnsEMBL::Production::Pipeline::Xrefs::ParseSource

my $parser = Bio::EnsEMBL::Xref::Parser::CCDSParser-E<gt>new(..)
$parser-E<gt>run();

=cut

sub run {

  my ( $self, $ref_arg ) = @_;
  my $source_id    = $self->{source_id};
  my $species_id   = $self->{species_id};
  my $verbose      = $self->{verbose} // 0;
  my $xref_dba     = $self->{xref_dba};
  my $dba          = $self->{dba};
  my $files        = $self->{files};

    if (   ( !defined $source_id )
        or ( !defined $species_id )
        or ( !defined $dba) )
    {
        confess 'Need to pass source_id, species_id, and database';
    }

  my $sql =(<<'CCDS');
    SELECT ta.value, t.stable_id
    FROM transcript t
    INNER JOIN transcript_attrib ta ON t.transcript_id = ta.transcript_id
    INNER JOIN attrib_type a ON ta.attrib_type_id = a.attrib_type_id
    WHERE a.code = 'ccds_transcript';
CCDS

  my $xref_count = 0;
  my $sth = $dba->dbc->prepare($sql);
  $sth->execute() or confess ( $dba->errstr() );
  while ( my ($ccds_id, $ens_id) = $sth->fetchrow_array() ) {
    my ( $acc, $version ) = split ( qr{ \. }msx, $ccds_id ) ;
    $xref_dba->add_to_direct_xrefs({
      acc        => $acc,
      version    => $version,
      label      => $acc,
      stable_id  => $ens_id,
      type       => 'transcript',
      source_id  => $source_id,
      species_id => $species_id,
      info_type  => 'DIRECT'
    } );
    $xref_count++;
  }
  $sth->finish;

  print "Added $xref_count DIRECT xrefs\n" if ($verbose);
  if ( !$xref_count ) {
   confess "No CCDS xref added\n";
  }

  return 0;      # successful

}

1;
