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

use strict;
use warnings;

use Readonly;

=head2 get_valid_source_id_to_external_db_id
  Description: Create a hash of all the external db names and ids that have
               associated xrefs
  Return type: Hashref
  Caller     : internal

=cut

=head2 get_xref_external_dbs
  Description: Create a hash of all the external db names and ids
  Return type: Hashref
  Caller     : internal

=cut

sub get_xref_external_dbs {

  my $self = shift;

  my %externalname_to_externalid;

  my $sth = $self->dbi->prepare_cached('SELECT db_name, external_db_id FROM external_db');
  $sth->execute() or confess( $self->dbi->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $external_db_name = $row[0];
    my $external_db_id   = $row[1];
    $externalname_to_externalid{$external_db_name} = $external_db_id;
  }

  return %externalname_to_externalid;
} ## end sub get_xref_external_dbs


sub get_valid_source_id_to_external_db_id {

  my $self = shift;

  my %source_id_to_external_db_id;

  my $sql = (<<'SQL')
    SELECT s.source_id, s.name
    FROM source s, xref x
    WHERE x.source_id = s.source_id
    GROUP BY s.source_id
SQL

  my $sth = $self->dbi->prepare_cached( $sql );
  $sth->execute() or confess( $self->dbi->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $source_name = $row[0];
    my $source_id   = $row[1];
    $source_id_to_external_db_id{ $source_name } = $source_id;
  }


  return %source_id_to_external_db_id;
} ## end sub get_valid_source_id_to_external_db_id


sub delete_projected_xrefs {
  my $sql = (<<'SQL');
    DELETE es
    FROM xref x, external_synonym es
    WHERE x.xref_id = es.xref_id and x.info_type = 'PROJECTION'
SQL
  my $sth = $core_dbi->prepare($sql);
  my $affected_rows = $sth->execute();
  if ( $verbose ) {
    print "\tDeleted $affected_rows PROJECTED external_synonym row(s)\n";
  }

  # Delete all ontologies, as they are done by a separate pipeline
  $sql = (<<'SQL');
    DELETE ontology_xref, object_xref, xref, dependent_xref
    FROM ontology_xref, object_xref, xref
    LEFT JOIN dependent_xref ON xref_id = dependent_xref_id
    WHERE ontology_xref.object_xref_id = object_xref.object_xref_id AND object_xref.xref_id = xref.xref_id
SQL
  $sth = $core_dbi->prepare($sql);
  $affected_rows = $sth->execute();
  if ( $verbose ) {
  print "\tDeleted $affected_rows PROJECTED ontology_xref row(s)\n";
  }

  $sql = (<<'SQL');
    DELETE object_xref
    FROM object_xref, xref
    WHERE object_xref.xref_id = xref.xref_id AND xref.info_type = 'PROJECTION'
SQL
  $sth = $core_dbi->prepare($sql);
  $affected_rows = $sth->execute();
  if ( $verbose ) {
    print "\tDeleted $affected_rows PROJECTED object_xref row(s)\n";
  }

  $sql = (<<'SQL');
    DELETE xref
    FROM xref
    WHERE xref.info_type = 'PROJECTION'
SQL
  $sth = $core_dbi->prepare($sql);
  $affected_rows = $sth->execute();
  if ( $verbose ) {
    print "\tDeleted $affected_rows PROJECTED xref row(s)\n";
  }
}


=head2 delete_by_external_db_id
  Arg [1]    : external_db_id
  Description: Delete xrefs and associated links for a given external db ID
  Return type:
  Exceptions : confess on a failed UPDATE
  Caller     : internal

=cut

sub delete_by_external_db_id {
  my ( $self, $external_db_id ) = @_;

  my $external_synonym = (<<'SQL');
    DELETE external_synonym
    FROM external_synonym, xref
    WHERE external_synonym.xref_id = xref.xref_id AND
          xref.external_db_id = ?
SQL

  my $ontology_xref = (<<'SQL');
    DELETE ontology_xref.*
    FROM ontology_xref, object_xref, xref
    WHERE ontology_xref.object_xref_id = object_xref.object_xref_id AND
          object_xref.xref_id = xref.xref_id AND
          xref.external_db_id = ?
SQL

  my $identity = (<<'SQL');
    DELETE identity_xref
    FROM identity_xref, object_xref, xref
    WHERE identity_xref.object_xref_id = object_xref.object_xref_id AND
          object_xref.xref_id = xref.xref_id AND
          xref.external_db_id = ?
SQL

  my $object_xref = (<<'SQL');
    DELETE object_xref
    FROM object_xref, xref
    WHERE object_xref.xref_id = xref.xref_id AND
          xref.external_db_id = ?
SQL

  my $master = (<<'SQL');
    DELETE ox, d
    FROM xref mx, xref x, dependent_xref d
         LEFT JOIN object_xref ox ON ox.object_xref_id = d.object_xref_id
    WHERE mx.xref_id = d.master_xref_id AND
          dependent_xref_id = x.xref_id AND
          mx.external_db_id = ?
SQL

  my $dependent = (<<'SQL');
    DELETE d
    FROM dependent_xref d, xref x
    WHERE d.dependent_xref_id = x.xref_id and x.external_db_id = ?
SQL

  my $xref = (<<'SQL');
    DELETE FROM xref WHERE xref.external_db_id = ?
SQL

  my $unmapped = (<<'SQL');
    DELETE FROM unmapped_object WHERE type="xref" and external_db_id = ?
SQL

  %sql_hash = (
    external_synonym => $external_synonym,
    ontology_xref    => $ontology_xref,
    identity_xref    => $identity_xref,
    object_xref      => $object_xref,
    master           => $master,
    dependent        => $dependent,
    xref             => $xref,
    unmapped         => $unmapped,
  );

  my @sth;
  my $i = 0;
  foreach my $table ( qw(external_synonym ontology_xref identity_xref
                         object_xref master dependent xref unmapped) ) {
    $sth[ $i++ ] = $self->dbi->prepare_cached( $sql_hash{$table} );
  }

  my $transaction_start_sth  =  $core_dbi->prepare('start transaction');
  my $transaction_end_sth    =  $core_dbi->prepare('commit');
  $transaction_start_sth->execute();

  for Readonly::Array my $ii ( 0..7 ) {
    $sth[$ii]->execute($external_db_id) or
      confess $self->dbi->errstr() . "\n $external_db_id\n\n";
  }

  $transaction_end_sth->execute();

  return;
} ## end sub delete_by_external_db_id


=head2 parsing_stored_data
  Description: Store data needed to be able to revert to same stage as after parsing
  Return type:
  Caller     : internal

  Notes      : Store max id for

    xref                                        xref_id
    object_xref                                 object_xref_id

=cut

sub parsing_stored_data {
  my $self = shift;

  my %table_and_key = (
    'xref'        => 'SELECT MAX(xref_id) FROM xref',
    'object_xref' => 'SELECT MAX(object_xref_id) FROM object_xref' );

  my %results = (
    xref => 0
    object_xref => 0
  );

  foreach my $table ( keys %table_and_key ) {
    my $sth = $self->dbi->prepare_cached( $table_and_key{$table} );
    $sth->execute;
    my $max_val;
    $sth->bind_columns( \$max_val );
    $sth->fetch;
    $self->add_meta_pair( $table . '_offset',
                          $max_val || 0);
    $sth->finish();

    $results{ $table } = $max_val || 0;
  }
  return %results;
} ## end sub parsing_stored_data




sub add_identity_xref {
  my ( $self, $xref ) = @_;

  my $object_xref_id   = $xref->{object_xref_id};
  my $xref_identity    = $xref->{xref_identity};
  my $ensembl_identity = $xref->{ensembl_identity};
  my $xref_start       = $xref->{xref_start};
  my $xref_end         = $xref->{xref_end};
  my $ensembl_start    = $xref->{ensembl_start};
  my $ensembl_end      = $xref->{ensembl_end};
  my $cigar_line       = $xref->{cigar_line};
  my $score            = $xref->{score};
  my $evalue           = $xref->{evalue};

  my $identity_sql = (<<'SQL');
    INSERT IGNORE INTO
      identity_xref ( object_xref_id, xref_identity, ensembl_identity,
                      xref_start, xref_end, ensembl_start, ensembl_end,
                      cigar_line, score, evalue )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
SQL

  my $add_identity_xref_sth  = $core_dbi->prepare_cached( $identity_sql );

  $add_identity_xref_sth->execute(
    ($object_xref_id+$object_xref_offset), $query_identity, $target_identity,
    $hit_start, $hit_end, $translation_start, $translation_end, $cigar_line,
    $score, $evalue);

  return;
}



sub add_dependent_xref {
  my ( $self, $object_xref_id, $master_xref_id, $dependent_xref_id ) = @_;

  my $sql = (<<'SQL');
    INSERT IGNORE INTO
      dependent_xref ( object_xref_id, master_xref_id, dependent_xref_id)
    VALUES (?, ?, ?)
SQL

my $sth  = $core_dbi->prepare_cached( $sql );

  $sth->execute(
    $object_xref_id, $master_xref_id, $dependent_xref_id);

  return;
}



sub add_xref_synonym {
  my ( @self, $xref_id, $syn)
  my $add_syn_sth = $core_dbi->prepare_cached(
    'INSERT IGNORE INTO external_synonym (xref_id, synonym) VALUES (?, ?)');

  $add_syn_sth->execute( ($xref_id+$xref_offset), $syn);

  return;
}



sub get_unmapped_reason_id {
  my ( $self, $desc_failed ) = @_;

  my $sql = (<<'SQL');
    SELECT unmapped_reason_id
    FROM unmapped_reason
    WHERE full_description LIKE ?
SQL

  my $sth = $core_dbi->prepare_cached( $sql );
  $sth->execute( $desc_failed );
  my $failed_id=undef;
  $sth->bind_columns(\$failed_id);
  $sth->fetch;

  return $failed_id;
}



sub add_unmapped_reason {
  my ( $self, $summary_failed, $desc_failed ) = @_;

  my $sql = (<<'SQL');
    INSERT INTO
      unmapped_reason (summary_description, full_description)
    VALUES (?, ?)
SQL

  my $sth = $core_dbi->prepare_cached( $sql);
  $sth->execute( $summary_failed, $desc_failed );
  return $sth->{'mysql_insertid'};
}



sub add_unmapped_object {
  my ( $self, $param ) = @_;

  my %params = %{ $param };

  my $sql = (<<"SQL");
    INSERT INTO
      unmapped_object (type, analysis_id, external_db_id, identifier, unmapped_reason_id )
    VALUES ('xref', @{ [ join', ', ('?') x keys %params ] } )
SQL

  my $sth =  $core_dbi->prepare_cached( $sql );
  $sth->execute( $analysis_id, external_db_id, $accession, $reason_id );
  return;
}

1;
