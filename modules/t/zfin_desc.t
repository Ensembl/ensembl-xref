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

use_ok 'Bio::EnsEMBL::Xref::Parser';

my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();

my $config = $db->config;

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config->{host},
  dbname => $config->{db},
  user   => $config->{user},
  pass   => $config->{pass},
  port   => $config->{port}
);

use_ok 'Bio::EnsEMBL::Xref::Parser::ZFINDescParser';

my $parser = Bio::EnsEMBL::Xref::Parser::ZFINDescParser->new(
 source_id  => 149,
 species_id => 7955,
 files      => ["$Bin/test-data/zfin_desc.txt"],
 xref_dba   => $xref_dba,
 verbose    => 1
);

isa_ok( $parser, 'Bio::EnsEMBL::Xref::Parser::ZFINDescParser' );

$parser->run();

ok(
  $db->schema->resultset('Xref')->check_direct_xref({
    accession   => 'ZDB-GENE-030131-3003',
    label       => 'hnf1bb',
    description => 'HNF1 homeobox Bb',
    source_id   => 149,
    species_id  => 7955,
    info_type   => 'MISC'
  }),
 'Sample zebrafish direct Xref has been inserted'
);

# Test if all the rows were inserted
is($db->schema->resultset('Xref')->count, 6, "All 6 rows were inserted");

my $parser_no_file = Bio::EnsEMBL::Xref::Parser::ZFINDescParser->new(
 source_id  => 149,
 species_id => 7955,
 files      => [],
 xref_dba   => $xref_dba
);

throws_ok{ $parser_no_file->run() } qr/No file name/, 'No file provided throws error' ;

done_testing();

