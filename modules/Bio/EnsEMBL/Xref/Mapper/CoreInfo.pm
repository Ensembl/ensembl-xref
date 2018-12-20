
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

Bio::EnsEMBL::Xref::Mapper::CoreInfo

=cut

=head1 DESCRIPTION

This is the base adaptor for loading xrefs into the species xref database during
the xref parsing stage in production. The aim of the module is to reduce
redundancy of code and keep the SQL to a constrained set of functions

=cut

package Bio::EnsEMBL::Xref::Mapper::CoreInfo;

use strict;
use warnings;

use Carp;

use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );

# Get info from the core database.

# Need to load tables:-
#
# gene_transcript_translation
# gene_stable_id
# transcript_stable_id
# translation_stable_id

=head2 new

=cut

sub new {
  my ( $caller, $mapper ) = @_;

  assert_ref( $mapper, 'Bio::EnsEMBL::Xref::Mapper' );

  my $class = ref($caller) || $caller;
  my $self = bless {}, $class;

  $self->mapper($mapper);

  return $self;
}

=head2 mapper

=cut

sub mapper {
  my $self = shift;
  $self->{_mapper} = shift if @_;

  return $self->{_mapper};
}

=head2 get_core_data

=cut

sub get_core_data {
  my $self = shift;

  # gene_transcript_translation
  # gene_stable_id
  # transcript_stable_id
  # translation_stable_id

  $self->load_gene_transcript_translation();
  $self->load_stable_ids();

  $self->mapper->xref->update_process_status('core_data_loaded');

  return;
}

=head2 load_gene_transcript_translation

=cut

sub load_gene_transcript_translation {
  my $self = shift;

  my $ins_sth = $self->mapper->xref->dbi->prepare(
"INSERT IGNORE INTO gene_transcript_translation (gene_id, transcript_id, translation_id) VALUES (?, ?, ?)"
  );
  my $sql =
"SELECT tn.gene_id, tn.transcript_id, tl.translation_id FROM transcript tn LEFT JOIN translation tl ON tl.transcript_id = tn.transcript_id";
  my $sth = $self->mapper->core->dbc->prepare($sql);
  $sth->execute();

  my ( $gene_id, $transcript_id, $translation_id );
  $sth->bind_columns( \$gene_id, \$transcript_id, \$translation_id );

  while ( $sth->fetch() ) {
    $ins_sth->execute( $gene_id, $transcript_id, $translation_id );
  }
  $ins_sth->finish;
  $sth->finish;

  return;
}

=head2 load_stable_ids

=cut

sub load_stable_ids {
  my $self = shift;

  my ( $id, $stable_id, $biotype );
  foreach my $table (qw(gene translation)) {
    my $sth = $self->mapper->core->dbc->prepare(
                 "SELECT " . $table . "_id, stable_id FROM " . $table );
    my $ins_sth =
      $self->mapper->xref->dbi->prepare( "INSERT IGNORE INTO " .
         $table . "_stable_id (internal_id, stable_id) VALUES (?, ?)" );
    $sth->execute();
    $sth->bind_columns( \$id, \$stable_id );

    while ( $sth->fetch ) {
      $ins_sth->execute( $id, $stable_id );
    }
    $ins_sth->finish;
    $sth->finish;
  }

  # populate transcript_stable_id table incuding the biotype column
  my $table = "transcript";
  my $sth = $self->mapper->core->dbc->prepare(
        "SELECT " . $table . "_id, stable_id, biotype FROM " . $table );
  my $ins_sth =
    $self->mapper->xref->dbi->prepare( "INSERT IGNORE INTO " . $table .
      "_stable_id (internal_id, stable_id, biotype) VALUES (?, ?, ?)" );
  $sth->execute();
  $sth->bind_columns( \$id, \$stable_id, \$biotype );

  while ( $sth->fetch ) {
    $ins_sth->execute( $id, $stable_id, $biotype );
  }
  $ins_sth->finish;
  $sth->finish;

  return;
} ## end sub load_stable_ids

1;
