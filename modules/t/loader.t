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
ok( defined $analysis_ids{'Gene'}, 'get_analysis - Gene');
ok( defined $analysis_ids{'Transcript'}, 'get_analysis - Transcript');
ok( defined $analysis_ids{'Translation'}, 'get_analysis - Translation');


# get_single_analysis
is(
  $loader_handle->get_single_analysis('xrefexoneratedna'), $analysis_ids{'Gene'},
  "get_single_analysis - Gene - xrefexoneratedna ($analysis_ids{'Gene'})"
);

is(
  $loader_handle->get_single_analysis('xrefexoneratedna'), $analysis_ids{'Transcript'},
  "get_single_analysis - Transcript - xrefexoneratedna ($analysis_ids{'Transcript'})"
);

is(
  $loader_handle->get_single_analysis('xrefexonerateprotein'), $analysis_ids{'Translation'},
  "get_single_analysis - Translation - xrefexonerateprotein ($analysis_ids{'Translation'})"
);


# add_xref
my $returned_xref_id = $loader_handle->add_xref(
  10,
  123,
  1,
  'dbprimary_acc',
  'display_label',
  0,
  'description',
  'MISC',
  'info_text',
  $loader_handle->core->dbc
);
is( $returned_xref_id, 123, "add_xref - xref ($returned_xref_id)");

my $dependent_xref_id = $loader_handle->add_xref(
  10,
  124,
  1,
  'dependent_primary_acc',
  'dependent_display_label',
  0,
  'dependent_description',
  'SEQUENCE_MATCH',
  'info_text',
  $loader_handle->core->dbc
);
is( $dependent_xref_id, 124, "add_xref - dependent ($dependent_xref_id)");


# add_object_xref
my $returned_object_xref_id = $loader_handle->add_object_xref(
  10,
  5000,
  1,
  'Gene',
  $returned_xref_id + 10,
  $analysis_ids{'Gene'},
  $loader_handle->core->dbc
);
is( $returned_object_xref_id, 5000, 'add_object_xref');


# get_xref_external_dbs
my %returned_external_db_ids = $loader_handle->get_xref_external_dbs();
ok( defined $returned_external_db_ids{'GO'}, 'get_xref_external_dbs' );


# delete_projected_xrefs
ok(
  $loader_handle->delete_projected_xrefs(
    $returned_external_db_ids{'RefSeq'}
  ),
  'delete_projected_xrefs'
);


# delete_by_external_db_id
ok( !defined $loader_handle->delete_by_external_db_id(), 'delete_by_external_db_id' );

# parsing_stored_data
my %returned_stored_data = $loader_handle->parsing_stored_data();
ok(
  defined $returned_stored_data{'xref'},
  "parsing_stored_data - xref ($returned_stored_data{'xref'})"
);
ok(
  defined $returned_stored_data{'object_xref'},
  "parsing_stored_data - object_xref ($returned_stored_data{'object_xref'})"
);

# add_identity_xref
ok (
  !defined $loader_handle->add_identity_xref( {
    object_xref_id => $returned_object_xref_id,
    xref_identity => 100,
    ensembl_identity => 100,
    xref_start => 1,
    xref_end => 256,
    ensembl_start => 1,
    ensembl_end => 256,
    cigar_line => '256M',
    score => 100,
    evalue => 24_000_000
  } ),
  'add_identity_xref'
);


# add_dependent_xref
ok (
  !defined $loader_handle->add_dependent_xref(
    $returned_object_xref_id,
    $returned_xref_id,
    $dependent_xref_id
  ),
  'add_dependent_xref'
);


# add_xref_synonym
ok(
  !defined $loader_handle->add_xref_synonym( $returned_xref_id, 'TestsSynonym' ),
  'add_xref_synonym'
);


# get_unmapped_reason_id
ok(
  !defined $loader_handle->get_unmapped_reason_id('Example Reason'),
  'get_unmapped_reason_id'
);


# add_unmapped_reason
my $unmapped_reason_id = $loader_handle->add_unmapped_reason(
  'unaligned', 'Sequence did not match to known sequenecs'
);
ok( defined $unmapped_reason_id, "add_unmapped_reason - $unmapped_reason_id" );

ok(
  defined $loader_handle->get_unmapped_reason_id('Sequence did not match%'),
  'get_unmapped_reason_id'
);


# add_unmapped_object
ok(
  !defined $loader_handle->add_unmapped_object({
    analysis_id        => $analysis_ids{'Gene'},
    external_db_id     => $returned_external_db_ids{'RefSeq'},
    identifier         => 'TestIdentifier',
    unmapped_reason_id => $unmapped_reason_id
  } ),
  'add_unmapped_object - basic'
);

ok(
  !defined $loader_handle->add_unmapped_object({
    analysis_id        => $analysis_ids{'Gene'},
    external_db_id     => $returned_external_db_ids{'RefSeq'},
    identifier         => 'TestIdentifier',
    unmapped_reason_id => $unmapped_reason_id,
    query_score        => 100
  } ),
  'add_unmapped_object - basic and query_score'
);



## Loader Tests
#  Tests for calling the loader functions

# load_unmapped_direct_xref
# load_unmapped_dependent_xref
# load_unmapped_sequence_xrefs
# load_unmapped_misc_xref
# load_unmapped_other_xref
# load_identity_xref
# load_checksum_xref
# load_dependent_xref
# load_synonyms


## Wrapper Tests
#  The are functions that wrap the logical calling of multiple functions

# update
# map_xrefs_from_xrefdb_to_coredb



done_testing();


sub _check_db {
   my ($dba, $table, $search_conditions) = @_;

   my $result_set = $dba->schema->resultset( $table )->search( $search_conditions );
   return $result_set->next;
}

