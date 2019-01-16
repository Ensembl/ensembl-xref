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
use Bio::EnsEMBL::Test::MultiTestDB;

use_ok 'Bio::EnsEMBL::Xref::Parser';

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba      = $multi_db->get_DBAdaptor('core');

my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();
my %config = %{ $db->config };

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config{host},
  dbname => $config{db},
  user   => $config{user},
  pass   => $config{pass},
  port   => $config{port}
);

# populate the species info as species_id2name method is called from the ArrayExpressParser 
$db->schema->populate(
  'Species',
  [
    [qw/species_id taxonomy_id name/],
    [ 9606, 9606, 'homo_sapiens' ],
  ]
);

use_ok 'Bio::EnsEMBL::Xref::Parser::ArrayExpressParser';

my $parser = Bio::EnsEMBL::Xref::Parser::ArrayExpressParser->new(
  source_id  => 1,
  species_id => 9606,
  species    => 'homo_sapiens',
  files      => [ "$Bin/test-data/homo_sapiens.ensgene.tsv" ],
  xref_dba   => $xref_dba,
  dba        => $dba
);
isa_ok( $parser, 'Bio::EnsEMBL::Xref::Parser::ArrayExpressParser' );

$parser->run();

ok(
  $db->schema->resultset('Xref')->check_direct_xref(
    {
      accession  => "ENSG00000131044",
      label      => 'ENSG00000131044',
      info_type  => 'DIRECT',
      source_id  => 1,
      species_id => 9606
    }
  ),
  'Sample Arrayexpress direct Xref has been inserted'
);

ok (
  $db->schema->resultset('GeneDirectXref')->find(
    {
       ensembl_stable_id => 'ENSG00000131044'
    }
  ),
   'Sample Arrayexpress gene direct Xref has been inserted'
 );

# Test if all the rows from gene table were inserted as the parser loads all gene stable_ids as direct xref
is( $db->schema->resultset('Xref')->count,
  21, "All 21 rows were inserted in to xref" );

is( $db->schema->resultset('GeneDirectXref')->count,
  21, "All 21 rows were inserted in to gene_direct_xref" );

done_testing();

