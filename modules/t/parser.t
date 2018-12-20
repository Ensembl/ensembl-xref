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

use Bio::EnsEMBL::Xref::Test::TestDB;
use Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Test::MultiTestDB;

use_ok 'Bio::EnsEMBL::Xref::Parser';

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $dba = $multi_db->get_DBAdaptor('core');

my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();
my %config = %{ $db->config };

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config{host},
  dbname => $config{db},
  user   => $config{user},
  pass   => $config{pass},
  port   => $config{port}
);

# check it throws without required arguments
throws_ok { Bio::EnsEMBL::Xref::Parser->new() } qr/Need to pass/, 'Throws with no arguments';
throws_ok { Bio::EnsEMBL::Xref::Parser->new(source_id => 1) } qr/Need to pass/, 'Throws with not enough arguments';
throws_ok { Bio::EnsEMBL::Xref::Parser->new(source_id => 1, species_id => 2) } qr/Need to pass/, 'Throws with not enough arguments';
throws_ok { Bio::EnsEMBL::Xref::Parser->new(source_id => 1, species_id => 2, files => 'dummy') } qr/Need to pass/, 'Throws with not enough arguments';

# check it throws with wrong argument types
throws_ok { Bio::EnsEMBL::Xref::Parser->new(source_id => 1, species_id => 2, files => 'dummy', xref_dba => 'dummy') }
  qr/check it is a reference/, 'Throws with scalar files arg';
throws_ok { Bio::EnsEMBL::Xref::Parser->new(source_id => 1, species_id => 2, files => ['dummy'], xref_dba => {}) }
  qr/was expected/, 'Throws with non DBI dbi arg';
throws_ok { Bio::EnsEMBL::Xref::Parser->new(source_id => 1, species_id => 2, files => ['dummy'], xref_dba => $xref_dba, dba => {}) }
  qr/was expected/, 'Throws with non DBAdaptor dba arg';

my $parser = Bio::EnsEMBL::Xref::Parser->new(source_id => 1, species_id => 2, files => ['dummy'], xref_dba => $xref_dba, dba => $dba);
isa_ok($parser, 'Bio::EnsEMBL::Xref::Parser');

throws_ok { $parser->run() } qr/abstract method/, 'Throws calling abstract method';

done_testing();

