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
    $returned_external_db_ids{'RefSeq_dna_predicted'}
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
    external_db_id     => $returned_external_db_ids{'RefSeq_dna_predicted'},
    identifier         => 'TestIdentifier',
    unmapped_reason_id => $unmapped_reason_id
  } ),
  'add_unmapped_object - basic'
);

ok(
  !defined $loader_handle->add_unmapped_object({
    analysis_id        => $analysis_ids{'Gene'},
    external_db_id     => $returned_external_db_ids{'RefSeq_dna_predicted'},
    identifier         => 'TestIdentifier',
    unmapped_reason_id => $unmapped_reason_id,
    query_score        => 100
  } ),
  'add_unmapped_object - basic and query_score'
);


## Loader Tests
#  Tests for calling the loader functions

## Prepare the xref db
my $source = $db->schema->resultset('Source')->create({
  name                 => 'RefSeq_dna_predicted',
  status               => 'KNOWN',
  source_release       => '38',
  download             => 'Y',
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


## Mapped Xrefs
# load_identity_xref
# get_insert_identity_xref
my $insert_identity_xrefs = $loader_handle->xref->get_insert_identity_xref(
  $source->source_id, 'DIRECT' );

while( my $insert_identity_xref_ref = $insert_identity_xrefs->() ) {
  my %insert_identity_xref = %{ $insert_identity_xref_ref };
  is(
    $insert_identity_xref{'acc'}, 'NM04560' ,
    "get_insert_identity_xref - $insert_identity_xref{'acc'}"
  );
}

my $loaded_identity_xrefs = $loader_handle->load_identity_xref(
  $source->source_id,                                # $source_id
  'DIRECT',                                          # $type
  $returned_stored_data{'xref'},                     # $xref_offset
  $returned_external_db_ids{'RefSeq_dna_predicted'}, # $ex_id
  $returned_stored_data{'object_xref'}               # $object_xref_offset
);
is( $loaded_identity_xrefs, 1, 'load_identity_xref');


# load_checksum_xref
my $insert_checksum_xrefs = $loader_handle->xref->get_insert_checksum_xref(
  $source->source_id, 'DIRECT' );

while( my $insert_checksum_xref_ref = $insert_checksum_xrefs->() ) {
  my %insert_checksum_xref = %{ $insert_checksum_xref_ref };
  ok(
    (
      $insert_checksum_xref{'acc'} eq 'NM04560' or
      $insert_checksum_xref{'acc'} eq 'NM04561'
    ),
    "get_insert_checksum_xref - $insert_checksum_xref{'acc'}"
  );
}

my $loaded_checksum_xrefs = $loader_handle->load_checksum_xref(
  $source->source_id,                                   # $source_id
  'DIRECT',                                             # $type
  $returned_stored_data{'xref'},                        # $xref_offset
  $returned_external_db_ids{'RefSeq_dna_predicted'},    # $ex_id
  $returned_stored_data{'object_xref'},                 # $object_xref_offset
  $loader_handle->get_single_analysis( 'xrefchecksum' ) # $checksum_analysis_id
);
is( $loaded_checksum_xrefs, 2, 'load_checksum_xref');


# load_dependent_xref
my $insert_dependent_xrefs = $loader_handle->xref->get_insert_dependent_xref(
  $source->source_id, 'DEPENDENT' );

while( my $insert_dependent_xref_ref = $insert_dependent_xrefs->() ) {
  my %insert_dependent_xref = %{ $insert_dependent_xref_ref };
  is(
    $insert_dependent_xref{'acc'}, 'NM04564' ,
    "get_insert_dependent_xref - $insert_dependent_xref{'acc'}"
  );
}

my $loaded_dependent_xrefs = $loader_handle->load_dependent_xref(
  $source->source_id,                                # $source_id
  'DEPENDENT',                                       # $type
  $returned_stored_data{'xref'},                     # $xref_offset
  $returned_external_db_ids{'RefSeq_dna_predicted'}, # $ex_id
  $returned_stored_data{'object_xref'},              # $object_xref_offset
);
is( $loaded_dependent_xrefs, 1, 'load_dependent_xref');


# load_synonyms
# add_multiple_synonyms
my @multi_syn_array = ( 'fs:000', 'fs:001', 'fs:002' );
ok(
  !defined $loader_handle->xref->add_multiple_synonyms(
    $loader_handle->xref->get_xref('NM04562', $source->source_id, 9606),
    \@multi_syn_array
  ),
  'Add multiple fake synonyms'
);

my @xref_with_syn_list = (
  $loader_handle->xref->get_xref('NM04562', $source->source_id, 9606),
);

ok(
  !defined $loader_handle->load_synonyms(
    \@xref_with_syn_list,
    $returned_stored_data{'xref'}
  ),
  'load_synonyms'
);

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


# load_unmapped_direct_xref
my @unmapped_direct_xrefs = $loader_handle->load_unmapped_direct_xref(
  $returned_stored_data{'xref'} + 100,
  $analysis_ids{'Transcript'},
  $unmapped_reason_id
);

foreach my $unmapped_id ( @unmapped_direct_xrefs ) {
  ok(
    ( $unmapped_id == 1 or $unmapped_id == 2 ),
    "load_unmapped_direct_xref - $unmapped_id"
  );
}


# load_unmapped_dependent_xref
my @unmapped_dependent_xrefs = $loader_handle->load_unmapped_dependent_xref(
  $returned_stored_data{'xref'} + 100,
  $analysis_ids{'Transcript'},
  $unmapped_reason_id
);

foreach my $unmapped_id ( @unmapped_dependent_xrefs ) {
  ok( $unmapped_id == 5, "load_unmapped_dependent_xref - $unmapped_id" );
}


# load_unmapped_sequence_xrefs - ensembl_id is defined
my @unmapped_sequence_xrefs = $loader_handle->load_unmapped_sequence_xrefs(
  $returned_stored_data{'xref'} + 100,
  \%analysis_ids,
  \%reason_id
);

is( $unmapped_sequence_xrefs[0], 3, 'load_unmapped_sequence_xrefs - 3' );

# load_unmapped_sequence_xrefs - ensembl_id is !defined


# load_unmapped_misc_xref
my @unmapped_misc_xrefs = $loader_handle->load_unmapped_misc_xref(
  $returned_stored_data{'xref'} + 100,
  $analysis_ids{'Transcript'},
  $reason_id{'NO_MAPPING'}
);

is( $unmapped_misc_xrefs[0], 4, 'load_unmapped_misc_xrefs - 4' );


# load_unmapped_other_xref
my @unmapped_other_xrefs = $loader_handle->load_unmapped_other_xref(
  $returned_stored_data{'xref'} + 100,
  $analysis_ids{'Transcript'},
  $reason_id{'NO_MASTER'}
);

is( $unmapped_other_xrefs[0], 5, 'load_unmapped_other_xrefs - 5' );

done_testing();

1;

