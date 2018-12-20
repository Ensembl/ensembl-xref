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

use strict;
use warnings;

use Test::More;
use Test::Exception;
use Test::Warnings;
use Bio::EnsEMBL::Xref::Test::TestDB;
use Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;

# Check that the module loaded correctly
use_ok 'Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor';

my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();
my %config = %{ $db->config };

my %search_conditions;

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config{host},
  dbname => $config{db},
  user   => $config{user},
  pass   => $config{pass},
  port   => $config{port}
);

ok($db, 'TestDB ready to go');

# Ensure that the BaseAdaptor handle is returned
ok( defined $xref_dba, 'BaseAdaptor handle returned');


# upload_xref_object_graphs
my $fake_xref = {
  accession   => 'NM01234',
  version     => 1,
  label       => 'NM01234.1',
  description => 'Fake RefSeq transcript',
  species_id  => '9606',
  info_type   => 'DIRECT',
  info_text   => 'These are normally aligned',
  dumped      => 'NO_DUMP_ANOTHER_PRIORITY'
};

my @xref_array_00 = ( $fake_xref );

throws_ok { $xref_dba->upload_xref_object_graphs() } qr/Please give me some xrefs to load/, 'Throws with no arguments';
throws_ok
  { $xref_dba->upload_xref_object_graphs( \@xref_array_00 ) }
  qr/Your xref does not have an accession-number/,
  'Throws with bad accession due to lowercase characters';


# Add an example source to the db
my $source = $db->schema->resultset('Source')->create({
  name                 => 'RefSeq',
  status               => 'KNOWN',
  source_release       => '38',
  download             => 'Y',
  priority             => 1,
  priority_description => 'Like a boss',
});

ok(defined $source->source_id, 'Was the source created in the DB?');

my $xref = $source->create_related('xrefs', $fake_xref);

my $rs = $db->schema->resultset('Xref')->search(
  { accession => 'NM01234'}
);

my $matching_xref = $rs->next;
ok(defined $matching_xref,'A result was pulled from the DB');


# Test get_filehandle
throws_ok { $xref_dba->get_filehandle() } qr/No file name/, 'Throws with no arguments';
throws_ok { $xref_dba->get_filehandle('fake_file.tsv') } qr/Can not open file/, 'Throws fake file name';


# get_source_id_for_source_name
throws_ok
  { $xref_dba->get_source_id_for_source_name() }
  qr/source_name undefined/,
  'get_source_id_for_source_name - Throws with no arguments';
throws_ok
  { $xref_dba->get_source_id_for_source_name('fake_name') }
  qr/No source_id/,
  'get_source_id_for_source_name - Throws with no matching source_id';
is(
  $xref_dba->get_source_id_for_source_name('RefSeq'),
  $source->source_id, 'get_source_id_for_source_name' );

# Add a second example source to the db
my $source2 = $db->schema->resultset('Source')->create({
  name                 => 'Second_fake_source',
  status               => 'KNOWN',
  source_release       => '38',
  download             => 'Y',
  priority             => 1,
  priority_description => 'Like a boss',
});

is(
  $xref_dba->get_source_id_for_source_name('RefSeq'),
  $source->source_id, 'get_source_id_for_source_name' );


# get_source_ids_for_source_name_pattern
throws_ok
  { $xref_dba->get_source_ids_for_source_name_pattern() }
  qr/source_name undefined/,
  'get_source_ids_for_source_name_pattern - Throws with no arguments';
ok(
  defined $xref_dba->get_source_ids_for_source_name_pattern('fake_name'),
  'get_source_ids_for_source_name_pattern - Array returned' );

is(
  $xref_dba->get_source_ids_for_source_name_pattern('RefSeq'),
  $source->source_id, 'get_source_ids_for_source_name_pattern' );


# get_source_name_for_source_id
throws_ok { $xref_dba->get_source_name_for_source_id() } qr/source_id undefined/, 'Throws with no arguments';
is( $xref_dba->get_source_name_for_source_id('fake_name'), -1, 'source not present in db' );
is( $xref_dba->get_source_name_for_source_id($source->source_id), 'RefSeq', 'get_source_name_for_source_id' );


# label_to_acc
my %l2a_result = %{ $xref_dba->label_to_acc( 'RefSeq', 9606 ) };
is( $l2a_result{ 'NM01234.1' }, 1, 'label_to_acc' );


# get_valid_codes
my %gvc_result = %{ $xref_dba->get_valid_codes( 'RefSeq', 9606 ) };
is( $gvc_result{ 'NM01234' }[0], 1, 'get_valid_codes' );


# Test to ensure that the code fails if there is no SOURCE_ID defined
my $new_xref_00 = {
  ACCESSION   => 'NM01235',
  VERSION     => 1,
  LABEL       => 'NM01235.1',
  DESCRIPTION => 'Fake RefSeq transcript',
  SPECIES_ID  => '9606',
  INFO_TYPE   => 'DIRECT',
  INFO_TEXT   => 'These are normally aligned'
};

my @xref_array_01 = ($new_xref_00);

throws_ok
  { $xref_dba->upload_xref_object_graphs( \@xref_array_01 ) }
  qr/Your xref: NM01235 does not have a source-id/,
  'upload_xref_object_graphs - Throws with no source ID';


# Test that a valid xref hashref gets added to the db
my $new_xref_01 = {
  ACCESSION   => 'NM01235',
  VERSION     => 1,
  LABEL       => 'NM01235.1',
  DESCRIPTION => 'Fake RefSeq transcript',
  SPECIES_ID  => '9606',
  SOURCE_ID   => $source->source_id,
  INFO_TYPE   => 'DIRECT',
  INFO_TEXT   => 'These are normally aligned'
};

my @xref_array_02 = ( $new_xref_01 );

ok(
  !defined $xref_dba->upload_xref_object_graphs( \@xref_array_02 ),
  "upload_xref_object_graphs - NM01235 was added to the xref table" );
is(
  _check_db( $db, 'Xref', { accession => 'NM01235' } )->accession,
  'NM01235', 'upload_xref_object_graphs - Xref loaded' );

# get_xref
my $xref_id = $xref_dba->get_xref('NM01235', $source->source_id, 9606);
ok( defined $xref_id, "NM01235 (xref_id: $xref_id) was added to the xref table" );


# upload_direct_xrefs - Testing later with the other direct_xref functions


# add_meta_pair
throws_ok {
  $xref_dba->add_meta_pair()
} qr/Need to specify/, 'Throws with required arguments not provided';
throws_ok {
  $xref_dba->add_meta_pair( 'fake_key' )
} qr/Need to specify/, 'Throws with required arguments not provided';

ok( !defined $xref_dba->add_meta_pair( 'fake_key', 'fake_value' ), 'add_meta_pair' );
is( _check_db( $db, 'Meta', { meta_key => 'fake_key' } )->meta_value, 'fake_value', 'Metadata pair added' );


# get_xref_sources
my $xref_sources = $xref_dba->get_xref_sources();
ok( defined $xref_sources, 'There are sources in the db');


# get_xref_id
throws_ok
  { $xref_dba->get_xref_id() }
  qr/Need an accession for get_xref_id/,
  'get_xref_id - Throws with no arguments 1';

throws_ok {
   $xref_dba->get_xref_id( {
      acc => 'NM01235'
   } )
} qr/Need a source_id for get_xref_id/, 'get_xref_id - Throws with no arguments 2';

throws_ok {
   $xref_dba->get_xref_id( {
      acc       => 'NM01235',
      source_id => $source->source_id
   } )
} qr/Need a species_id for get_xref_id/, 'get_xref_id - Throws with no arguments 3';

is(
   $xref_dba->get_xref_id( {
      acc        => 'NM01235',
      source_id  => $source->source_id,
      species_id => 9606
   } ),
   2,
   'get_xref_id'
);

# Test for a new coordinate_xref
my $new_coordinate_xref = {
  accession    => 'ucsc.1',
  chromosome   => 5,
  strand       => 1,
  txStart      => 1,
  txEnd        => 10,
  exonStarts   => 1,
  exonEnds     => 10,
  species_id   => '9606',
  source_id    => $source->source_id
};
my $coordinate_xref_id_new = $xref_dba->add_coordinate_xref( $new_coordinate_xref );
ok( defined $coordinate_xref_id_new, "ucsc.1 (coordinate_xref_id: $coordinate_xref_id_new) was loaded to the xref table" );
is( _check_db( $db, 'CoordinateXref', { coord_xref_id => $coordinate_xref_id_new } )->accession, 'ucsc.1', 'Coordinate_xref loaded' );


# add_xref
# Test for xref that already exists
my $new_xref_02 = {
  acc          => 'NM01235',
  version      => 1,
  label        => 'NM01235.1',
  desc         => 'Fake RefSeq transcript',
  species_id   => '9606',
  source_id    => $source->source_id,
  info_type    => 'DIRECT',
  info_text    => 'These are normally aligned',
  update_label => 1,
  update_desc  => 1
};
my $xref_id_old = $xref_dba->add_xref( $new_xref_02 );
ok( defined $xref_id_old, "NM01235 (xref_id: $xref_id_old) was loaded" );
ok($xref_id_old == $xref_id, "NM01235 already existed in the db");

# Test for a new xref
my $new_xref_03 = {
  acc          => 'NM01236',
  version      => 1,
  label        => 'NM01236.1',
  desc         => 'Fake RefSeq transcript',
  species_id   => '9606',
  source_id    => $source->source_id,
  info_type    => 'DIRECT',
  info_text    => 'These are normally aligned',
  update_label => 1,
  update_desc  => 1
};
my $xref_id_new = $xref_dba->add_xref( $new_xref_03 );
ok( defined $xref_id_new, "NM01236 (xref_id: $xref_id_new) was loaded to the xref table" );
is( _check_db( $db, 'Xref', { xref_id => $xref_id_new } )->accession, 'NM01236', 'Xref loaded' );


# add_object_xref
throws_ok {
   $xref_dba->add_object_xref()
} qr/add_object_xref needs an xref_id/, 'add_object_xref - Throws with no arguments 1';

throws_ok {
   $xref_dba->add_object_xref( { xref_id => $xref_id } )
} qr/add_object_xref needs an ensembl_id/, 'add_object_xref - Throws with no arguments 2';

throws_ok {
   $xref_dba->add_object_xref( { xref_id => $xref_id, ensembl_id => 1 } )
} qr/add_object_xref needs an object_type/, 'add_object_xref - Throws with no arguments 3';

my $object_xref_id = $xref_dba->add_object_xref( { xref_id => $xref_id, ensembl_id => 1, object_type => 'Gene' } );
ok( defined $object_xref_id, "add_object_xref - Object_xref entry inserted - $object_xref_id" );


# get_object_xref
ok( !defined $xref_dba->get_object_xref( 1000000, 'ENSFAKE000000000', 'Gene' ), 'get_object_xref - Fake' );
is( $xref_dba->get_object_xref( $xref_id, 1, 'Gene' ), 1, 'get_object_xref - Real' );


# add_identity_xref
throws_ok {
   $xref_dba->add_identity_xref()
} qr/add_identity_xref needs an object_xref_id/, 'add_identity_xref - Throws with no arguments 1';

throws_ok {
   $xref_dba->add_identity_xref(
      { object_xref_id => $object_xref_id } )
} qr/add_identity_xref needs a score/, 'add_identity_xref - Throws with no arguments 2';

throws_ok {
   $xref_dba->add_identity_xref( { object_xref_id => $object_xref_id, score => 1 } )
} qr/add_identity_xref needs a target_identity/, 'add_identity_xref - Throws with no arguments 3';

throws_ok { $xref_dba->add_identity_xref(
   { object_xref_id => $object_xref_id, score => 1, target_identity => 1 }
) } qr/add_identity_xref needs a query_identity/, 'add_identity_xref - Throws with no arguments 4';

ok(
   !defined $xref_dba->add_identity_xref(
      { object_xref_id => $object_xref_id, score => 1, target_identity => 1, query_identity => 1 } ),
   "add_identity_xref - Identity xref row added" );


# add_to_direct_xrefs
my $new_xref_04 = {
  stable_id    => 'NM01236',
  type         => 'Gene',
  acc          => 'NM01236',
  version      => 1,
  label        => 'NM01236.1',
  desc         => 'Fake RefSeq transcript',
  species_id   => '9606',
  source_id    => $source->source_id,
  info_text    => 'These are normally aligned',
  update_label => 1,
  update_desc  => 1
};
throws_ok { $xref_dba->add_to_direct_xrefs() }
   qr/Need a direct_xref on which this xref linked too/, 'add_to_direct_xrefs - Throws with no arguments 1';

throws_ok { $xref_dba->add_to_direct_xrefs(
   { stable_id => 'NM01236' }
) } qr/Need a table type on which to add/, 'add_to_direct_xrefs - Throws with no arguments 2';

throws_ok { $xref_dba->add_to_direct_xrefs(
   { stable_id => 'NM01236', type => 'Gene' }
) } qr/Need an accession of this direct xref/, 'add_to_direct_xrefs - Throws with no arguments 3';

throws_ok { $xref_dba->add_to_direct_xrefs(
   { stable_id => 'NM01236', type => 'Gene', acc => 'NM01236' }
) } qr/Need a source_id for this direct xref/, 'add_to_direct_xrefs - Throws with no arguments 4';

throws_ok { $xref_dba->add_to_direct_xrefs(
   { stable_id => 'NM01236', type => 'Gene', acc => 'NM01236', source_id => $source->source_id }
) } qr/Need a species_id for this direct xref/, 'add_to_direct_xrefs - Throws with no arguments 5';

ok( !defined $xref_dba->add_to_direct_xrefs( $new_xref_04 ), 'add_to_direct_xrefs' );
is(
  _check_db( $db, 'Xref', { accession => 'NM01236' } )->accession,
  'NM01236', 'add_to_direct_xrefs - Direct xref NM01236 has been loaded' );


# add_direct_xref
# Entry has already been added, so this test should be just for failing out
ok( !defined $xref_dba->add_direct_xref( $xref_id_new, 'NM01236', 'Gene' ), 'add_direct_xref' );


# get_valid_xrefs_for_direct_xrefs
my %valid_direct_xrefs = %{ $xref_dba->get_valid_xrefs_for_direct_xrefs( 'RefSeq', ',' ) };
is( $valid_direct_xrefs{ 'NM01236' }, '3,NM01236,Gene,', 'get_valid_xrefs_for_direct_xrefs' );


# upload_direct_xrefs
my $new_xref_05 = {
  ACCESSION    => 'NM01235',
  SPECIES_ID   => '9606',
  SOURCE_ID    => $source->source_id,
  STABLE_ID    => 'NM01235',
  ENSEMBL_TYPE => 'Transcript',
  LINKAGE_XREF => 'PROBE',
  SOURCE       => 'RefSeq'
};

my @xref_array_03 = ( $new_xref_05 );
ok( !defined $xref_dba->upload_direct_xrefs( \@xref_array_03 ), 'upload_direct_xrefs' );


# add_multiple_direct_xrefs
my $new_direct_xref_00 = {
  STABLE_ID    => 'NM01236',
  ENSEMBL_TYPE => 'Gene',
  ACCESSION    => 'DX01236',
  VERSION      => 1,
  LABEL        => 'DX01236.1',
  DESCRIPTION  => 'Fake RefSeq transcript',
  SPECIES_ID   => '9606',
  SOURCE_ID    => $source->source_id,
  INFO_TEXT    => 'These are normally aligned',
};

my @xref_array_05 = ( $new_direct_xref_00 );
ok( !defined $xref_dba->add_multiple_direct_xrefs( \@xref_array_05 ), 'add_multiple_direct_xrefs' );


# add_dependent_xref
my $new_xref_06 = {
  master_xref_id => $xref_id_new,
  type           => 'Gene',
  acc            => 'XX123456',
  version        => 1,
  label          => 'DPNDT',
  desc           => 'Fake dependent xref',
  species_id     => '9606',
  source_id      => $source->source_id,
  info_text      => 'These are normally aligned',
  linkage        => $source->source_id,
  update_label   => 1,
  update_desc    => 1
};

throws_ok { $xref_dba->add_dependent_xref() }
   qr/Need a master_xref_id on which this xref linked too/, 'add_dependent_xref - Throws with no arguments 1';

throws_ok { $xref_dba->add_dependent_xref(
   { master_xref_id => $xref_id_new }
) } qr/Need an accession of this dependent xref/, 'add_dependent_xref - Throws with no arguments 2';

throws_ok { $xref_dba->add_dependent_xref(
   { master_xref_id => $xref_id_new, acc => 'XX123456' }
) } qr/Need a source_id for this dependent xref/, 'add_dependent_xref - Throws with no arguments 3';

throws_ok { $xref_dba->add_dependent_xref(
   { master_xref_id => $xref_id_new, acc => 'XX123456', source_id => $source->source_id }
) } qr/Need a species_id for this dependent xref/, 'add_dependent_xref - Throws with no arguments 4';

my $dependent_xref_id = $xref_dba->add_dependent_xref( $new_xref_06 );
ok( defined $dependent_xref_id, "add_dependent_xref - Dependent xref entry inserted - $dependent_xref_id" );


# get_valid_xrefs_for_dependencies
my %valid_dependent_xrefs = %{ $xref_dba->get_valid_xrefs_for_dependencies( 'RefSeq', ( 'RefSeq' ) ) };
is( $valid_dependent_xrefs{'XX123456'}, $xref_id_new, 'get_valid_xrefs_for_dependencies' );


# add_multiple_dependent_xrefs
my @xref_array_06 = ( $new_xref_06 );
ok( !defined $xref_dba->add_multiple_dependent_xrefs( \@xref_array_06 ), 'add_multiple_dependent_xrefs' );

# add_to_syn_for_mult_sources
my @source_array = ( $source->source_id );
ok(
   !defined $xref_dba->add_to_syn_for_mult_sources( 'NM01234', \@source_array, 'fake_multi_synonym', 9606 ),
   'add_to_syn_for_mult_sources' );


# add_to_syn
throws_ok { $xref_dba->add_to_syn(
   'XX01236', $source->source_id, 'fake_synonym', 9606
) } qr/Could not find acc XX01236/, 'Throws with no matching accession';
ok( !defined $xref_dba->add_to_syn( 'NM01236', $source->source_id, 'fake_synonym', 9606 ), 'Add to synonyms' );


# add_synonym
ok( !defined $xref_dba->add_synonym( $xref_id_new, 'fake_synonym' ), 'Add a fake synonym' );


# add_multiple_synonyms
my @multi_syn_array = ( 'fs:000', 'fs:001', 'fs:002' );
ok( !defined $xref_dba->add_multiple_synonyms( $xref_id_new, \@multi_syn_array ), 'Add multiple fake synonyms' );


# get_label_to_acc
ok( scalar $xref_dba->get_label_to_acc( 'RefSeq', 9606, 'Like a boss' ) > 0, 'get_label_to_acc' );


# get_acc_to_label
ok( scalar $xref_dba->get_acc_to_label( 'RefSeq', 9606, 'Like a boss' ) > 0, 'get_acc_to_label' );


# get_label_to_desc
ok( scalar $xref_dba->get_label_to_desc( 'RefSeq', 9606, 'Like a boss' ) > 0, 'get_label_to_desc' );


# set_release
ok( !defined $xref_dba->set_release( $source->source_id, 100 ), 'Set release without problems' );


# get_dependent_mappings
ok( !defined $xref_dba->get_dependent_mappings( $source->source_id ), 'Dependent mappings updated' );


# get_ext_synonyms
ok( scalar $xref_dba->get_ext_synonyms( 'RefSeq' ) > 0, 'get_ext_synonyms' );


# parsing_finished_store_data
ok( !defined $xref_dba->parsing_finished_store_data(), 'parsing_finished_store_data' );


# get_meta_value
my $get_meta_value = $xref_dba->get_meta_value('PARSED_xref_id');
ok( $xref_dba->get_meta_value('PARSED_xref_id') == 5, "get_meta_value ($get_meta_value)" );


# _update_xref_info_type
ok( !defined $xref_dba->_update_xref_info_type( $xref_id_new, 'PROBE' ), '_update_xref_info_type - Real xref_id' );
is( _check_db( $db, 'Xref', { xref_id => $xref_id_new } )->info_type, 'PROBE', 'info_type updated' );
ok( !defined $xref_dba->_update_xref_info_type( 1000000, 'PROBE' ), '_update_xref_info_type - Fake xref_id' );


# _add_pair
ok( !defined $xref_dba->_add_pair( $source->source_id, 'NM01236', 'NT01236' ), '_add_pair' );
ok( defined _check_db( $db, 'Pair', { source_id => $source->source_id } ), 'Data loaded' );


# _add_primary_xref
ok( defined $xref_dba->_add_primary_xref( $xref_id_new, 'GATACCA', 'dna', 'experimental' ), '_add_primary_xref');
ok( defined _check_db( $db, 'PrimaryXref', { status => 'experimental' } ), 'Data loaded' );
is( _check_db( $db, 'PrimaryXref', { xref_id => $xref_id_new } )->sequence, 'GATACCA', 'Sequence loaded' );


# _update_primary_xref_sequence
ok( !defined $xref_dba->_update_primary_xref_sequence( $xref_id_new, 'CTATGGT' ), '_update_primary_xref_sequence');
is( _check_db( $db, 'PrimaryXref', { xref_id => $xref_id_new } )->sequence, 'CTATGGT', 'Sequence updated' );


# _update_xref_label - This should have already been covered by previous tests
# Specific tests can be added later

# _update_xref_description - This should have already been covered by previous tests
# Specific tests can be added later

note 'Test methods to support base mapper';

throws_ok {
  $xref_dba->get_id_from_species_name()
} qr/Undefined/, 'Throws with no name argument';

throws_ok {
  $xref_dba->get_id_from_species_name( 'Vampire' )
} qr/No ID/, 'Throws with unavailable species';

is( $xref_dba->get_id_from_species_name('Homo sapiens'), 1, 'Species id from name' );
my $names = $xref_dba->get_species_names();
is( scalar @{ $names }, 1, 'Number of species names' );
is( $names->[0], 'Homo sapiens', 'Species name' );

throws_ok {
  $xref_dba->add_alt_allele()
} qr/Need to specify/, 'Throws with missing arguments';
  
my $gene_stable_id = $db->schema->resultset('GeneStableId')->create({
  internal_id          => 1,
  stable_id            => 'FakeStableId0001'
});

$xref_dba->add_alt_allele( 10, 1, 1 );
is( _check_db( $db, 'AltAllele', { alt_allele_id => 10 } )->gene_id, 1, 'Add alt allele' );

$xref_dba->update_process_status( 'alt_alleles_added' );
is( _check_db( $db, 'ProcessStatus', { status => 'alt_alleles_added' } )->id, 1, 'Update process status' );

is( $xref_dba->xref_latest_status(), 'alt_alleles_added', 'Latest status' );

$xref_dba->delete_alt_alleles();
ok( !_check_db( $db, 'AltAllele' ), 'No alt alleles after deletion');

throws_ok {
  $xref_dba->update_mapping_jobs_status()
} qr/not given/, 'Throws with missing arguments';

my $mapping_job = $db->schema->resultset('MappingJob')->create({
  root_dir          => '/fake_dir/',
  status            => 'SUBMITTED',
  job_id            => 1
});
is( _check_db( $db, 'MappingJob', { job_id => 1 } )->status, 'SUBMITTED', 'Added mapping job' );

$xref_dba->update_mapping_jobs_status( 'SUCCESS' );
is( _check_db( $db, 'MappingJob', { job_id => 1 } )->status, 'SUCCESS', 'Updated mapping job status' );

# get/set_species
ok( !defined($xref_dba->species), "Species not defined yet");
is($xref_dba->species("human"), "human", "Species set to ". $xref_dba->species);

# get_dba
ok( defined($xref_dba->dba), "dba defined");
isa_ok( $xref_dba->dba, 'Bio::EnsEMBL::DBSQL::DBAdaptor' );

done_testing();


sub _check_db {
   my ($db, $table, $search_conditions) = @_;

   my $rs = $db->schema->resultset( $table )->search( $search_conditions );
   return $rs->next;
}


