
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
use Test::More;
use Test::Exception;
use Test::Warnings;
use FindBin '$Bin';
use lib "$Bin/";

use Bio::EnsEMBL::Xref::Test::TestDB;
use Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Test::MultiTestDB;
use File::Temp qw/ tempdir /;

use_ok 'Bio::EnsEMBL::Xref::Mapper';

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $core_dba = $multi_db->get_DBAdaptor('core');

my $db     = Bio::EnsEMBL::Xref::Test::TestDB->new();
my %config = %{ $db->config };

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config{host},
  dbname => $config{db},
  user   => $config{user},
  pass   => $config{pass},
  port   => $config{port}
);

use_ok 'Bio::EnsEMBL::Xref::Mapper::CoordinateMapper';

# populate the CoordinateXref with test coordinates
# NM_001008407 => Do not overlap
# NM_001008408 => Do not meet the threshold
# NM_001008409 => Meets the threshold

$db->schema->populate(
  'CoordinateXref',
  [
    [
      qw/source_id species_id accession chromosome strand txstart txend cdsstart cdsend exonstarts exonends/
    ],
    [
      50, 9606, 'NM_001008407', '19', 1, 31166507, 31196939, 31167507,
      31196839, "31166507,31193157,31183173", "31166567,31193309,31183380"
    ],
    [
      50, 9606, 'NM_001008408', '20', 1, 31166507, 31196939, 31167507,
      31196839, "31166507,31193157,31183173", "31166567,31193309,31183380"
    ],
    [
      50, 9606, 'NM_001008409', '20', 1, 31166507, 31196939, 31167507,
      31196839,
      "31166507, 31193157, 31183173, 31186274, 31180256, 31172464,31195211",
      "31166567, 31193309, 31183380, 31186395, 31180401, 31172587, 31196939"
    ],
  ]
);

# initialize cooordinatemapper with xref and core dba
my $coord_mapper =
  Bio::EnsEMBL::Xref::Mapper::CoordinateMapper->new( xref_dba => $xref_dba, core_dba => $core_dba );

# pass species_id and base_path as arguments
my $species_id = "9606";
my $base_path = tempdir( CLEANUP => 1 );
$coord_mapper->run_coordinatemapping( 1, $species_id, $base_path );

# check if base_path dir exists
ok( -d $base_path, "$base_path exists" );

# check if mapping text files are created
ok( -e $base_path . '/object_xref_coord.txt',
  'object_xref_coord file created' );
ok(
  -e $base_path . '/unmapped_object_coord.txt',
  'unmapped_object_coord file created'
);
ok(
  -e $base_path . '/unmapped_reason_coord.txt',
  'unmapped_reason_coord file created'
);
ok( -e $base_path . '/xref_coord.txt', 'xref_coord file created' );

# check if the files are non-zero size
ok(
  -s $base_path . '/object_xref_coord.txt',
  'object_xref_coord file is non-zero size'
);
ok(
  -s $base_path . '/unmapped_object_coord.txt',
  'unmapped_object_coord file is non-zero size'
);
ok(
  -s $base_path . '/unmapped_reason_coord.txt',
  'unmapped_reason_coord file is non-zero size'
);
ok( -s $base_path . '/xref_coord.txt', 'xref_coord file is non-zero size' );

my $ta         = $core_dba->get_TranscriptAdaptor();
my $tr         = $ta->fetch_by_stable_id("ENST00000217347");
my @ucsc_xrefs = @{ $tr->get_all_xrefs("UCSC") };

my $ucsc_xref = shift @ucsc_xrefs;
ok( defined $ucsc_xref, "ucsc_xref is defined" );

# check if the xref object if of the right type with the right ids set
is( $ucsc_xref->info_type, "COORDINATE_OVERLAP",
  "Got the right info_type as " . $ucsc_xref->info_type );
is( $ucsc_xref->primary_id, "NM_001008409",
  "Got the right primary_id as " . $ucsc_xref->primary_id );
is( $ucsc_xref->display_id, "NM_001008409",
  "Got the right display_id as " . $ucsc_xref->display_id );
is(
  $ucsc_xref->db_display_name,
  "UCSC Stable ID",
  "Got the right db_display_name as " . $ucsc_xref->db_display_name
);

# now check for unmapped objects
my $uoa              = $core_dba->get_UnmappedObjectAdaptor();
my @unmapped_objects = @{ $uoa->fetch_by_identifier('NM_001008408') };
my $unmapped_xref    = shift @unmapped_objects;
ok(
  $unmapped_xref->summary eq "Did not meet threshold",
  "Got the identifier that did not meet the threshold"
);

ok(
  $unmapped_xref->query_score < 0.75,
  "Got score ("
    . $unmapped_xref->query_score
    . ") less than the threshold (0.75)"
);

@unmapped_objects = @{ $uoa->fetch_by_identifier('NM_001008407') };
$unmapped_xref    = shift @unmapped_objects;
ok(
  $unmapped_xref->summary eq "No overlap",
  "Got the identifier that has no overlap"
);

ok( !defined $unmapped_xref->query_score,
  "No score. Score is undefined as there is no overlap" );

done_testing();

