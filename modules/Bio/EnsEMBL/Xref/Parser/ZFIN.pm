
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

Bio::EnsEMBL::Xref::Parser::ZFIN - Base parser class for ZFIN

=head1 DESCRIPTION

This is a base class providing common data and methods for specialisations
tailored to parsing UniProt, RefSeq and aliases ZFIN files.

It's organised this way since the eHive pipeline is assumed to process
one file at a time.

=cut

package Bio::EnsEMBL::Xref::Parser::ZFIN;

use strict;
use warnings;

use Carp;

use parent qw( Bio::EnsEMBL::Xref::Parser );

=head1 METHODS

=head2 new

Base constructor

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  # get the source ids for ZFIN
  $self->{source_ids} = [
    $self->{xref_dba}->get_source_ids_for_source_name_pattern("ZFIN_ID")
  ];

# map xref descriptions to accessions, store count of how many (with description) have been loaded
  $self->{loaded} = $self->_build_description_map();

  return $self;
}

=head2 loaded

Return number of loaded ZFIN xrefs

=cut

sub loaded {
  return shift->{loaded};
}

=head2 description

=cut

sub description {
  return shift->{description};

}

=head2 _build_description_map

Build hash mappinng ZFIN xref descriptions to their corresponding accession
=cut

sub _build_description_map {
  my $self = shift;

  my $sql =
"SELECT accession, label, version, description FROM xref WHERE source_id in ("
    . join( ", ", @{ $self->{source_ids} } ) . ")";
  my $sth = $self->{xref_dba}->dbi->prepare_cached($sql);
  $sth->execute();

  my ( $acc, $lab, $ver, $desc );
  $sth->bind_columns( \$acc, \$lab, \$ver, \$desc );

  my $loaded = 0;
  while ( my @row = $sth->fetchrow_array() ) {
    $self->{description}{$acc} = $desc if defined $desc;
    $loaded++;
  }

  return $loaded;

}

1;
