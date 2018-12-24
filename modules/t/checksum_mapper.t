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

use_ok 'Bio::EnsEMBL::Xref::Mapper::ChecksumMapper';

use Bio::EnsEMBL::Xref::Test::TestDB;
use Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Test::MultiTestDB;

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $dba = $multi_db->get_DBAdaptor('core');

my $xref_db = Bio::EnsEMBL::Xref::Test::TestDB->new();
my $source = $xref_db->schema->resultset('Source')->create({
  name => 'UniParc'
});


my %config = %{ $xref_db->config };
my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config{host},
  dbname => $config{db},
  user   => $config{user},
  pass   => $config{pass},
  port   => $config{port}
);

my $mapper = Bio::EnsEMBL::Xref::Mapper::ChecksumMapper->new(
  xref_dba => $xref_dba,
  core_dba => $dba,
  external_db_name => 'UniParc'
);

is($mapper->logic_name, 'xrefchecksum', 'Test getters');
is($mapper->object_type,'Translation', 'Defaults to translation feature type');
cmp_ok($mapper->batch_size, '==', 1000, 'Default batch size is set');
is($mapper->external_db_name, 'UniParc', 'Provided db name is set');
cmp_ok($mapper->source_id, '==', 1, 'For the purposes of this test the source_id matching');

$mapper->method('do_stuff');
is ($mapper->method(), 'do_stuff', 'Test obvious stuff');

# $db->populate();

my $translation_adaptor = $dba->get_TranslationAdaptor;
my $trans_1 = $translation_adaptor->fetch_by_stable_id('ENSP00000308980');
my $trans_2 = $translation_adaptor->fetch_by_stable_id('ENSP00000278995');

$mapper->upload([
  { id => $trans_1->dbID, upi => 'UPI00000A', object_type => 'Translation' },
  { id => $trans_2->dbID, upi => 'UPI00000B', object_type => 'Translation' }
], 1, 1);

my $xref_1 = $xref_db->schema->resultset('Xref')->find({
  accession => 'UPI00000A'
});
ok ($xref_1, 'First accession was uploaded' );

my $xref_2 = $xref_db->schema->resultset('Xref')->find({
  accession => 'UPI00000B'
});
ok ($xref_2, 'Second accession was uploaded' );

is($xref_1->object_xref->ensembl_id , $trans_1->dbID, 'Object Xref links to original translation dbID');
is($xref_2->object_xref->ensembl_id , $trans_2->dbID, 'Object Xref links to original translation dbID');

done_testing();
