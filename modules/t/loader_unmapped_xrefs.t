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

use Data::Dumper;

use Test::More;
use Test::Exception;
use Test::Warnings;
use Bio::EnsEMBL::Xref::Test::TestDB;
use Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Xref::Mapper;
use Bio::EnsEMBL::Xref::Mapper::CoreInfo;
use Bio::EnsEMBL::Xref::Mapper::Loader;

use Bio::EnsEMBL::Test::MultiTestDB;

# Check that the modules loaded correctly
use_ok 'Bio::EnsEMBL::Xref::Mapper';
use_ok 'Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor';

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $dba = $multi_db->get_DBAdaptor('core');

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


my $loader_handle = Bio::EnsEMBL::Xref::Mapper::Loader->new( xref_dba => $xref_dba, core_dba => $dba );
isa_ok( $loader_handle, 'Bio::EnsEMBL::Xref::Mapper::Loader' );


## Primary Tests
#  These are tests for single operation functions that do not rely on the
#  functionality of other functions in the module


# get_analysis
my %analysis_ids = $loader_handle->get_analysis();

# get_xref_external_dbs
my %returned_external_db_ids = $loader_handle->get_xref_external_dbs();

# parsing_stored_data
my %returned_stored_data = $loader_handle->parsing_stored_data();


## Loader Tests
#  Tests for calling the loader functions

# delete_projected_xrefs
ok(
  !defined $loader_handle->delete_projected_xrefs(
    $returned_external_db_ids{'RefSeq_dna_predicted'}
  ),
  'delete_projected_xrefs'
);

## Prepare the xref db
my $source = $db->schema->resultset('Source')->create({
  name                 => 'RefSeq_dna_predicted',
  source_release       => '38',
  priority             => 1,
  priority_description => 'Like a boss',
});

is(
  $loader_handle->xref->get_source_id_for_source_name('RefSeq_dna_predicted'),
  $source->source_id, 'get_source_id_for_source_name'
);

my $source_mapping_method = $db->schema->resultset('SourceMappingMethod')->create({
  source_id => $source->source_id,
  method    => 'Alignment',
});

my $mapping = $db->schema->resultset('Mapping')->create({
  job_id => 123456,
  type    => 'dna',
  command_line => 'ls',
  percent_query_cutoff => 75,
  percent_target_cutoff => 90,
  method => 'Alignment',
  array_size => 5,
});

my $new_xref_00 = {
  ACCESSION    => 'NM04560',
  VERSION      => 1,
  LABEL        => 'NM04560.1',
  DESCRIPTION  => 'Fake RefSeq transcript',
  SPECIES_ID   => '9606',
  SOURCE_ID    => $source->source_id,
  INFO_TYPE    => 'DIRECT',
  INFO_TEXT    => 'These are normally aligned',
  update_label => 1,
  update_desc  => 1
};

my $new_xref_01 = {
  ACCESSION    => 'NM04561',
  VERSION      => 1,
  LABEL        => 'NM04561.1',
  DESCRIPTION  => 'Fake RefSeq transcript',
  SPECIES_ID   => '9606',
  SOURCE_ID    => $source->source_id,
  INFO_TYPE    => 'DIRECT',
  INFO_TEXT    => 'These are normally aligned',
  update_label => 1,
  update_desc  => 1
};

my $new_xref_02 = {
  ACCESSION    => 'NM04562',
  VERSION      => 1,
  LABEL        => 'NM04562.1',
  DESCRIPTION  => 'Fake RefSeq transcript',
  SPECIES_ID   => '9606',
  SOURCE_ID    => $source->source_id,
  INFO_TYPE    => 'SEQUENCE_MATCH',
  INFO_TEXT    => 'These are normally aligned',
  update_label => 1,
  update_desc  => 1
};

my $new_xref_03 = {
  ACCESSION    => 'NM04563',
  VERSION      => 1,
  LABEL        => 'NM04563.1',
  DESCRIPTION  => 'Fake RefSeq misc',
  SPECIES_ID   => '9606',
  SOURCE_ID    => $source->source_id,
  INFO_TYPE    => 'MISC',
  INFO_TEXT    => 'These are normally aligned',
  update_label => 1,
  update_desc  => 1
};

my @new_xref_array = ( $new_xref_00, $new_xref_01, $new_xref_02, $new_xref_03 );
$loader_handle->xref->upload_xref_object_graphs( \@new_xref_array );

my $object_xref_id_00 = $loader_handle->xref->add_object_xref(
  {
    xref_id     => $loader_handle->xref->get_xref('NM04560', $source->source_id, 9606),
    ensembl_id  => 1,
    object_type => 'Gene'
  }
);
ok( defined $object_xref_id_00, "add_object_xref - Object_xref entry inserted - $object_xref_id_00" );

ok(
   !defined $loader_handle->xref->add_identity_xref(
      { object_xref_id => $object_xref_id_00, score => 1, target_identity => 1, query_identity => 1 } ),
   "add_identity_xref - Identity xref row added" );


my $object_xref_id_01 = $loader_handle->xref->add_object_xref(
  {
    xref_id     => $loader_handle->xref->get_xref('NM04561', $source->source_id, 9606),
    ensembl_id  => 1,
    object_type => 'Gene'
  }
);
ok( defined $object_xref_id_01, "add_object_xref - Object_xref entry inserted - $object_xref_id_01" );

my $object_xref_id_02 = $loader_handle->xref->add_object_xref(
  {
    xref_id     => $loader_handle->xref->get_xref('NM04562', $source->source_id, 9606),
    ensembl_id  => 1,
    object_type => 'Gene'
  }
);
ok( defined $object_xref_id_02, "add_object_xref - Object_xref entry inserted - $object_xref_id_02" );

ok(
  defined $loader_handle->xref->_add_primary_xref(
    $loader_handle->xref->get_xref('NM04562', $source->source_id, 9606),
    'GATACCA', 'dna', 'experimental'
  ),
  '_add_primary_xref'
);

my $object_xref_id_03 = $loader_handle->xref->add_object_xref(
  {
    xref_id     => $loader_handle->xref->get_xref('NM04563', $source->source_id, 9606),
    ensembl_id  => 1,
    object_type => 'Gene'
  }
);
ok( defined $object_xref_id_03, "add_object_xref - Object_xref entry inserted - $object_xref_id_03" );

# add_dependent_xref
my $new_xref_04 = {
  master_xref_id => $object_xref_id_01,
  type           => 'Gene',
  acc            => 'NM04564',
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

my $new_xref_04_id = $loader_handle->xref->add_dependent_xref( $new_xref_04 );
ok( defined $new_xref_04_id, "add_dependent_xref - Dependent xref entry inserted - $new_xref_04_id" );

my $object_xref_id_04 = $loader_handle->xref->add_object_xref(
  {
    xref_id     => $new_xref_04_id,
    ensembl_id  => 1,
    object_type => 'Gene'
  }
);
ok( defined $object_xref_id_04, "add_object_xref - Object_xref entry inserted - $object_xref_id_04" );

my $xref_dbi = $loader_handle->xref->dbi();

my $dependent_update_sth = $xref_dbi->prepare(
  'UPDATE object_xref SET master_xref_id = ? WHERE object_xref_id = ?' );

$dependent_update_sth->execute(
  $loader_handle->xref->get_xref('NM04561', $source->source_id, 9606),
  $object_xref_id_04 );

# Set the dumping on the object_xref table
$db->schema->resultset('ObjectXref')->update({ ox_status => 'DUMP_OUT' });

## Unmapped Xrefs
# Get the cutoff values
my %failed_sources = %{ $loader_handle->xref->get_unmapped_reason() };

my %summary_failed = %{ $failed_sources{'summary'} };
my %desc_failed    = %{ $failed_sources{'desc'} };

my %reason_id;
foreach my $key (keys %desc_failed){
  my $failed_id = $loader_handle->get_unmapped_reason_id( $desc_failed{$key} );

  if(!defined $failed_id ) {
    $failed_id = $loader_handle->add_unmapped_reason(
      $summary_failed{$key}, $desc_failed{$key} );
  }
  $reason_id{$key} = $failed_id;
}

## Wrapper Tests
#  The are functions that wrap the logical calling of multiple functions
#  unmapped_xrefs_from_xrefdb_to_coredb
ok(
  !defined $loader_handle->unmapped_xrefs_from_xrefdb_to_coredb(
    $returned_stored_data{'xref'},
    $returned_stored_data{'object_xref'},
    \%analysis_ids,
    \%reason_id
  ),
  'unmapped_xrefs_from_xrefdb_to_coredb'
);


## Final wrapper function
# update



done_testing();

1;

