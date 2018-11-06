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
use Bio::EnsEMBL::Xref::Test::TestDB;


my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();
ok($db, 'default instantiation proceeds as planned with testdb.conf');
throws_ok { Bio::EnsEMBL::Xref::Test::TestDB->new(config_file => 'not_here') }
  qr/does not exist!/,
  'TestDB grumbles differently when an invalid config file is offered';

# This auto-deploys the schema
$db = Bio::EnsEMBL::Xref::Test::TestDB->new(
  # config_file => 'testdb.conf'
  config => { 
    driver => 'SQLite',
    file => 'test.db',
    create => 1
  }
);

ok($db, 'TestDB ready to go');

my $source = $db->schema->resultset('Source')->create({
  name => 'RefSeq',
  status => 'KNOWN',
  source_release => '38',
  download => 'Y',
  priority => 1,
  priority_description => 'Like a boss',
  ordered => 10
});

ok(defined $source->source_id, 'Was the source created in the DB?');


my $xref = $source->create_related('xrefs', {
  accession => 'NM01234',
  version => 1,
  label => 'NM01234.1',
  description => 'Fake RefSeq transcript',
  species_id => '9606',
  info_type => 'DIRECT',
  info_text => 'These are normally aligned',
  dumped => 'NO_DUMP_ANOTHER_PRIORITY'
});

my $rs = $db->schema->resultset('Xref')->search(
  { accession => 'NM01234'}
);

my $matching_xref = $rs->next;
ok(defined $matching_xref,'A result was pulled from the DB');
is($matching_xref->accession, $xref->accession, 'Retrieved xref is the same as that which was stored');

is($matching_xref->source_id, $source->source_id, 'Source IDs also match');

is($matching_xref->source->name, 'RefSeq', 'Foreign "key" relation works');

done_testing;