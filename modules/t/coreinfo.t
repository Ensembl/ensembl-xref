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

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Xref::Test::TestDB;

use Bio::EnsEMBL::Xref::Mapper;
use Bio::EnsEMBL::Xref::Mapper::CoreInfo;

# Check that the module loaded correctly
use_ok 'Bio::EnsEMBL::Xref::Mapper::CoreInfo';

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $dba = $multi_db->get_DBAdaptor('core');

my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();
ok($db, 'Xref DB');
my %config = %{ $db->config };

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config{host},
  dbname => $config{db},
  user   => $config{user},
  pass   => $config{pass},
  port   => $config{port}
);

throws_ok {
  Bio::EnsEMBL::Xref::Mapper::CoreInfo->new()
} qr/undef/, 'Throws without argument';

throws_ok {
  Bio::EnsEMBL::Xref::Mapper::CoreInfo->new([])
} qr/was expected/, 'Throws with argument of incorrect class';

my $mapper = Bio::EnsEMBL::Xref::Mapper->new( xref_dba => $xref_dba, core_dba => $dba );

my $coreinfo = Bio::EnsEMBL::Xref::Mapper::CoreInfo->new( $mapper );
isa_ok( $coreinfo, 'Bio::EnsEMBL::Xref::Mapper::CoreInfo' );

done_testing();
