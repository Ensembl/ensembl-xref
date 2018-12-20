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
use Data::Dumper;
use FindBin '$Bin';
use lib "$Bin/";

use Bio::EnsEMBL::Xref::Test::TestDB;
use Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;

use_ok 'Bio::EnsEMBL::Xref::Parser';

my $db     = Bio::EnsEMBL::Xref::Test::TestDB->new();
my %config = %{ $db->config };

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config{host},
  dbname => $config{db},
  user   => $config{user},
  pass   => $config{pass},
  port   => $config{port}
);

# populate the synonyms to test if
$db->schema->populate( 'Synonym',
  [ [qw/xref_id synonym/], [ 1, 'Test-Synonym1' ], [ 1, 'Test-Synonym2' ] ] );

use_ok 'Bio::EnsEMBL::Xref::Parser::MGIParser';

my $parser = Bio::EnsEMBL::Xref::Parser::MGIParser->new(
  source_id  => 55,
  species_id => 10090,
  files      => ["$Bin/test-data/MRK_ENSEMBL.rpt"],
  xref_dba   => $xref_dba
);
isa_ok( $parser, 'Bio::EnsEMBL::Xref::Parser::MGIParser' );

$parser->run();

ok(
  $db->schema->resultset('Xref')->check_direct_xref(
    {
      accession   => 'MGI:1915733',
      label       => '1110002O04Rik',
      description => 'RIKEN cDNA 1110002O04 gene',
      source_id   => 55,
      species_id  => 10090,
      info_type   => 'DIRECT'
    }
  ),
  'Sample mouse direct Xref has been inserted'
);

ok(
  $db->schema->resultset('GeneDirectXref')->find(
    {
      ensembl_stable_id => "ENSMUSG00000102531"
    }
  ),
  'Sample mouse gene direct Xref has been inserted'
);


is($db->schema->resultset('Xref')->count, 10, "All 10 rows were inserted in to Xref");
is($db->schema->resultset('GeneDirectXref')->count, 10, "All 10 rows were inserted in to GeneDirectXref");
is($db->schema->resultset('Synonym')->count, 2, "Synonym count remained the same");


done_testing();
