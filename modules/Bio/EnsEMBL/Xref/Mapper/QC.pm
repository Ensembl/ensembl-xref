
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

Bio::EnsEMBL::Xref::Mapper::QC - Basic health checks related to xref mapping

=head1 SYNOPSIS

  my $tester = Bio::EnsEMBL::Xref::Mapper::QC->new( $mapper );
  confess 'Problems found' if $tester->unlinked_entries;

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Xref::Mapper::QC;

use strict;
use warnings;

use Bio::EnsEMBL::Xref::Mapper;

##########################################
# Testing  (may be moved to healthchecks)
##########################################

##### Unlinked entries ##############

# ERRORS
#    dependent_xref            and xref
#    primary_xref              and xref
#    transcript_direct_xref    and xref
#    translation_direct_xref   and xref
#    gene_direct_xref          and xref
#    synonym                   and xref

#    identity_xref             and object_xref
#    go_xref                   and object_xref

#    gene_transcript_translation   and gene_stable_id
#    gene_transcript_translation   and transcript_stable_id
#    gene_transcript_translation   and translation_stable_id

# WARNNGS
#    gene_direct_xref              and gene_stable_id
#    transcript
#    translation

# All object_xref of type go have a go_xref entry

##### Numbers between xref and core (xref and object_xref) are similar

##### if human or mouse check the number of gene name changes.

=head2 new

=cut

sub new {
  my ( $caller, $mapper ) = @_;

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

=head2 unlinked_entries

=cut

sub unlinked_entries {
  my $self = shift;

  $self->mapper->xref->update_process_status('tests_started');

  my $dbi = $self->mapper->xref->dbi;

  # dependent_xref and xref
  my ( $xref_id, $count );
  my $count_sql = (<<"COUNT1");
SELECT COUNT(1) 
  FROM dependent_xref d LEFT JOIN xref x 
    ON d.master_xref_id = x.xref_id
    WHERE x.xref_id IS NULL;
COUNT1

  my $sth = $dbi->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns( \$count );
  $sth->fetch();
  $sth->finish;

  my $failed = 0;
  my $sql;
  if ($count) {
    $failed = 1;

    $sql = (<<"DISTINCT1");
SELECT DISTINCT(d.master_xref_id) 
  FROM dependent_xref d LEFT JOIN xref x 
    ON d.master_xref_id = x.xref_id 
    WHERE x.xref_id IS NULL 
  LIMIT 10;
DISTINCT1

    $sth = $dbi->prepare($sql);
    $sth->execute();
    $sth->bind_columns( \$xref_id );

    print STDERR "SQL QUERY: $sql\n";
    while ( $sth->fetch ) {
      print STDERR "Problem with master xref $xref_id\n";
    }
    $sth->finish;
  }

  $count_sql = (<<"COUNT2");
SELECT COUNT(1) 
  FROM dependent_xref d LEFT JOIN xref x 
    ON d.dependent_xref_id = x.xref_id 
    WHERE x.xref_id IS NULL;
COUNT2

  $sql = (<<"DISTINCT2");
SELECT DISTINCT(d.dependent_xref_id)
  FROM dependent_xref d LEFT JOIN xref x
    ON d.dependent_xref_id = x.xref_id 
    WHERE x.xref_id IS NULL 
  LIMIT 10;
DISTINCT2

  $sth = $dbi->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns( \$count );
  $sth->fetch();
  $sth->finish;

  if ($count) {
    $failed = 1;
    $sth    = $dbi->prepare($sql);
    $sth->execute();
    $sth->bind_columns( \$xref_id );

    print STDERR "SQL QUERY: $sql\n";
    while ( $sth->fetch ) {
      print STDERR "Problem with dependent xref $xref_id\n";
    }
    $sth->finish;
  }

  $count_sql = (<<"COUNT3");
SELECT COUNT(1) 
  FROM primary_xref d LEFT JOIN xref x 
    ON d.xref_id = x.xref_id 
    WHERE x.xref_id IS NULL;
COUNT3

  $sql = (<<"DISTINCT3");
SELECT DISTINCT(d.xref_id) 
  FROM primary_xref d LEFT JOIN xref x 
    ON d.xref_id = x.xref_id 
    WHERE x.xref_id IS NULL 
  LIMIT 10;
DISTINCT3

  $sth = $dbi->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns( \$count );
  $sth->fetch();
  $sth->finish;

  if ($count) {
    $failed = 1;
    $sth    = $dbi->prepare($sql);
    $sth->execute();
    $sth->bind_columns( \$xref_id );

    print STDERR "SQL QUERY: $sql\n";
    while ( $sth->fetch ) {
      print STDERR "Problem with primary xref $xref_id\n";
    }
    $sth->finish;
  }

  foreach my $type (qw(transcript translation gene)) {
    $count_sql = (<<"COUNT4");
SELECT COUNT(1) 
  FROM "${type}_direct_xref" d LEFT JOIN xref x 
    ON d.general_xref_id = x.xref_id 
    WHERE x.xref_id IS NULL;
COUNT4

    $sql = (<<"DISTINCT4");
SELECT DISTINCT(d.general_xref_id) 
  FROM ${type}_direct_xref d LEFT JOIN xref x 
    ON d.general_xref_id = x.xref_id 
    WHERE x.xref_id IS NULL 
  LIMIT 10;
DISTINCT4

    $sth = $dbi->prepare($count_sql);
    $sth->execute();
    $sth->bind_columns( \$count );
    $sth->fetch();
    $sth->finish;

    if ($count) {
      $failed = 1;
      $sth    = $dbi->prepare($sql);
      $sth->execute();
      $sth->bind_columns( \$xref_id );

      print STDERR "SQL QUERY: $sql\n";
      while ( $sth->fetch ) {
        print STDERR "Problem with " .
          $type . "_direct_xref $xref_id\n";
      }
      $sth->finish;
    }
  } ## end foreach my $type (qw(transcript translation gene))

  $count_sql = (<<"COUNT5");
SELECT COUNT(1) 
  FROM synonym d LEFT JOIN xref x 
    ON d.xref_id = x.xref_id 
    WHERE x.xref_id IS NULL;
COUNT5

  $sql = (<<"DISTINCT5");
SELECT DISTINCT(d.xref_id) 
  FROM synonym d LEFT JOIN xref x 
    ON d.xref_id = x.xref_id 
    WHERE x.xref_id IS NULL 
  LIMIT 10;
DISTINCT5

  $sth = $dbi->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns( \$count );
  $sth->fetch();
  $sth->finish;

  if ($count) {
    $failed = 1;
    $sth    = $dbi->prepare($sql);
    $sth->execute();
    $sth->bind_columns( \$xref_id );

    print STDERR "SQL QUERY: $sql\n";
    while ( $sth->fetch ) {
      print STDERR "Problem with synonym $xref_id\n";
    }
    $sth->finish;
  }

  $count_sql = (<<"COUNT6");
SELECT COUNT(1) 
  FROM identity_xref d LEFT JOIN object_xref o 
    ON d.object_xref_id = o.object_xref_id 
    WHERE o.object_xref_id IS NULL;
COUNT6

  $sql = (<<"DISTINCT6");
SELECT DISTINCT(d.object_xref_id) 
  FROM identity_xref d LEFT JOIN object_xref o 
    ON d.object_xref_id = o.object_xref_id 
    WHERE o.object_xref_id IS NULL 
  LIMIT 10;
DISTINCT6

  $sth = $dbi->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns( \$count );
  $sth->fetch();
  $sth->finish;

  if ($count) {
    $failed = 1;
    $sth    = $dbi->prepare($sql);
    $sth->execute();
    $sth->bind_columns( \$xref_id );

    print STDERR "SQL QUERY: $sql\n";
    while ( $sth->fetch ) {
      print STDERR "Problem with object_xref $xref_id\n";
    }
    $sth->finish;
  }

  $count_sql = (<<"COUNT7");
SELECT COUNT(1) 
  FROM go_xref d LEFT JOIN object_xref o 
    ON d.object_xref_id = o.object_xref_id 
    WHERE o.object_xref_id IS NULL;
COUNT7

  $sql = (<<"DISTINCT7");
SELECT DISTINCT(d.object_xref_id) 
  FROM go_xref d LEFT JOIN object_xref o 
    ON d.object_xref_id = o.object_xref_id 
    WHERE o.object_xref_id IS NULL 
  LIMIT 10;
DISTINCT7

  $sth = $dbi->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns( \$count );
  $sth->fetch();
  $sth->finish;

  if ($count) {
    $failed = 1;
    $sth    = $dbi->prepare($sql);
    $sth->execute();
    $sth->bind_columns( \$xref_id );

    print STDERR "SQL QUERY: $sql\n";
    while ( $sth->fetch ) {
      print STDERR "Problem with object_xref $xref_id\n";
    }
    $sth->finish;
  }

  foreach my $type (qw(transcript translation gene)) {
    $count_sql = (<<"COUNT8");
SELECT COUNT(1) 
  FROM gene_transcript_translation d LEFT JOIN ${type}_stable_id x 
    ON d.${type}_id = x.internal_id 
    WHERE x.internal_id IS NULL AND 
    d.${type}_id IS NOT NULL;
COUNT8

    $sql = (<<"DISTINCT8");
SELECT DISTINCT(d.${type}_id) 
  FROM gene_transcript_translation d LEFT JOIN ${type}_stable_id x 
    ON d.${type}_id = x.internal_id 
    WHERE x.internal_id IS NULL AND 
      d.${type}_id IS NOT NULL 
    LIMIT 10;
DISTINCT8

    $sth = $dbi->prepare($count_sql);
    $sth->execute();
    $sth->bind_columns( \$count );
    $sth->fetch();
    $sth->finish;

    if ($count) {
      $failed = 1;
      $sth    = $dbi->prepare($sql);
      $sth->execute();
      $sth->bind_columns( \$xref_id );

      print STDERR "SQL QUERY: $sql\n";
      while ( $sth->fetch ) {
        print STDERR "Problem with " . $type . "_id $xref_id\n";
      }
      $sth->finish;
    }
  } ## end foreach my $type (qw(transcript translation gene))

  $count_sql = (<<"COUNT9");
SELECT COUNT(1) 
  FROM xref x, source s, object_xref o LEFT JOIN go_xref g 
    ON o.object_xref_id = g.object_xref_id 
    WHERE x.xref_id = o.xref_id AND s.source_id = x.source_id AND 
      s.name LIKE 'GO' AND 
      ox_status IN ('DUMP_OUT') AND 
      g.object_xref_id IS NULL;
COUNT9

  $sql = (<<"DISTINCT9");
SELECT DISTINCT(o.object_xref_id) 
  FROM xref x, source s, object_xref o LEFT JOIN go_xref g 
    ON o.object_xref_id = g.object_xref_id 
    WHERE x.xref_id = o.xref_id AND 
      s.source_id = x.source_id AND 
      s.name LIKE 'GO'          AND 
      ox_status IN ('DUMP_OUT') AND 
      g.object_xref_id IS NULL 
  LIMIT 10;
DISTINCT9

  $sth = $dbi->prepare($count_sql);
  $sth->execute();
  $sth->bind_columns( \$count );
  $sth->fetch();
  $sth->finish;

  if ($count) {
    $failed = 1;
    $sth    = $dbi->prepare($sql);
    $sth->execute();
    $sth->bind_columns( \$xref_id );

    print STDERR "SQL QUERY: $sql\n";
    while ( $sth->fetch ) {
      print STDERR
"Problem with object_xref $xref_id which is linked to a GO source but has no go_xref reference\n";
    }
    $sth->finish;
  }

  if ( !$failed ) {
    $self->mapper->xref->update_process_status('tests_finished');
  }
  else {
    $self->mapper->xref->update_process_status('tests_failed');
  }

  return $failed;
} ## end sub unlinked_entries

=head2 entry_number_check

=cut

sub entry_number_check {
  my $self = shift;
  my $dbi  = $self->mapper->xref->dbi;

# No point doing xrefs object_xrefs are more important and gives a better indication of wether things went okay.
  my %old_object_xref_count;
  my %new_object_xref_count;

  my $sql = (<<"SQL1");
SELECT s.name, COUNT(distinct x.xref_id, ensembl_id)
  FROM xref x, object_xref ox, source s
    WHERE ox.xref_id = x.xref_id AND 
      x.source_id = s.source_id  AND 
      ox_status = "DUMP_OUT"     AND 
      s.name NOT LIKE "AFFY%" 
  GROUP BY s.name
SQL1

  my $sth = $dbi->prepare($sql);
  $sth->execute();

  my ( $name, $count );
  $sth->bind_columns( \$name, \$count );
  while ( $sth->fetch() ) {
    $new_object_xref_count{$name} = $count;
  }
  $sth->finish;

  $sql = (<<"SQL2");
SELECT e.db_name, COUNT(*) 
  FROM xref x, object_xref ox, external_db e 
    WHERE ox.xref_id = x.xref_id          AND 
      x.external_db_id = e.external_db_id AND 
      e.db_name NOT LIKE "AFFY%"          AND 
      (x.info_type IS NULL OR x.info_type != "PROJECTION") 
  GROUP BY e.db_name
SQL2

  $sth = $self->mapper->core->dbc->prepare($sql);
  $sth->execute();
  $sth->bind_columns( \$name, \$count );

  while ( $sth->fetch() ) {
    my $change = 0;
    $old_object_xref_count{$name} = $count;
    if ( defined $new_object_xref_count{$name} ) {
      $change =
        ( ( $new_object_xref_count{$name} - $count )/$count )*100;

      if ( $change > 5 ) {    # increase of 5%
        print "WARNING: $name has increased by " .
          int($change) .
          "\% was $count now " . $new_object_xref_count{$name} . "\n"
          if $self->mapper->verbose;
      }
      elsif ( $change < -5 ) {    # decrease by 5%
        print "WARNING: $name has decreased by " .
          int($change) .
          " \% was $count now " . $new_object_xref_count{$name} . "\n"
          if $self->mapper->verbose;
      }
    }
    else {
      print
"WARNING: xrefs $name are not in the new database but $count are in the old???\n"
        if $self->mapper->verbose;
    }
  } ## end while ( $sth->fetch() )
  $sth->finish;

  foreach my $key ( keys %new_object_xref_count ) {
    print "WARNING: $key has " . $new_object_xref_count{$key} .
      " xrefs in the new database but NONE in the old\n"
      if not defined $old_object_xref_count{$key} and
      $self->mapper->verbose;
  }

  return;
} ## end sub entry_number_check

=head2 name_change_check

=cut

sub name_change_check {
  my $self = shift;
  my $dbi  = $self->mapper->xref->dbi;

  my %new_name;    # $old_name{$gene_id} = HGNC_%name
  my %id_to_stable_id;

  my $official_name = $self->mapper->get_official_name;
  return unless defined $official_name;

  my $sql = (<<"SQL1");
SELECT x.label, gsi.internal_id, gsi.stable_id 
  FROM object_xref ox, xref x, gene_stable_id gsi, source s 
    WHERE x.xref_id = ox.xref_id AND ox.ensembl_object_type = ? AND 
      gsi.internal_id = ox.ensembl_id AND 
      x.source_id = s.source_id       AND 
      s.name LIKE ?;
SQL1

  my $sth = $dbi->prepare($sql);
  $sth->execute(qw( Gene $official_name_% ));

  my ( $name, $gene_id, $stable_id, $syn );
  $sth->bind_columns( \$name, \$gene_id, \$stable_id );
  my $count = 0;

  while ( $sth->fetch() ) {
    $new_name{$gene_id}        = $name;
    $id_to_stable_id{$gene_id} = $stable_id;
    $count++;
  }
  $sth->finish;

  # Use synonyms as well.
  my %alias;
  $sql = (<<"SQL2");
SELECT x.label, sy.synonym 
  FROM xref x, synonym sy, source so 
    WHERE x.xref_id = sy.xref_id AND 
      so.source_id = x.source_id AND 
      so.name LIKE "${official_name}_%";
SQL2

  $sth = $dbi->prepare($sql);
  $sth->execute();

  $sth->bind_columns( \$name, \$syn );
  $count = 0;

  while ( $sth->fetch() ) {
    $alias{$syn} = $name;
  }
  $sth->finish;

  $sql = (<<"SQL3");
SELECT x.label, sy.synonym 
  FROM xref x, synonym sy, source so 
    WHERE x.xref_id = sy.xref_id AND 
    so.source_id = x.source_id   AND 
    so.name LIKE "EntrezGene"';
SQL3

  $sth = $dbi->prepare($sql);
  $sth->execute();
  $sth->bind_columns( \$name, \$syn );

  while ( $sth->fetch() ) {
    $alias{$syn} = $name;
  }
  $sth->finish;

  # NOTE ncRNA has higher priority
  $sql = (<<"SQL4");
SELECT x.display_label, g.gene_id 
  FROM gene g, xref x 
    WHERE g.display_xref_id = x.xref_id AND 
      biotype = 'protein_coding';
SQL4

  $sth = $self->mapper->core->dbc->prepare($sql);
  $sth->execute();
  $sth->bind_columns( \$name, \$gene_id );

  $count = 0;
  my $total_count = 0;
  while ( $sth->fetch() ) {
    $total_count++ if defined $new_name{$gene_id};

    if ( defined $new_name{$gene_id} and $new_name{$gene_id} ne $name )
    {
      if ( not defined $alias{$name} or
           $alias{$name} ne $new_name{$gene_id} )
      {
        print STDERR "WARN: gene_id ($gene_id) " .
          $id_to_stable_id{$gene_id} .
          " new = " . $new_name{$gene_id} . " old = $name\n";
        $count++;
      }
    }
  }

  print STDERR
"$count entries with different names out of $total_count protein coding gene comparisons?\n"
    if $total_count;

  return;
} ## end sub name_change_check

=head2 direct_stable_id_check

=cut

sub direct_stable_id_check {
  my $self = shift;
  my $dbi  = $self->mapper->xref->dbi;

  foreach my $type (qw(gene transcript translation)) {

    my $sql = (<<"SQL");
SELECT s.name, COUNT(*) 
  FROM source s, xref x, ${type}_direct_xref gdx LEFT JOIN ${type}_stable_id gsi 
    ON gdx.ensembl_stable_id = gsi.stable_id 
    WHERE s.source_id = x.source_id AND 
    x.xref_id = gdx.general_xref_id AND 
    gsi.stable_id IS NULL 
    GROUP BY s.name;
SQL

    my $sth = $dbi->prepare($sql);
    $sth->execute();

    my ( $name, $count );
    $sth->bind_columns( \$name, \$count );

    my $total_count = 0;
    while ( $sth->fetch() ) {
      print STDERR "WARNING $name has $count invalid stable ids in " .
        $type . "_direct_xrefs\n";
      $total_count += $count;
    }
    $sth->finish;

    print STDERR "USEFUL SQL: $sql\n" if $total_count;
  } ## end foreach my $type (qw(gene transcript translation))

  return;
} ## end sub direct_stable_id_check

1;
