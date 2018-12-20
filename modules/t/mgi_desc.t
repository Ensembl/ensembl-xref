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

my $db     = Bio::EnsEMBL::Xref::Test::TestDB->new();
my %config = %{ $db->config };

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config{host},
  dbname => $config{db},
  user   => $config{user},
  pass   => $config{pass},
  port   => $config{port}
);

use_ok 'Bio::EnsEMBL::Xref::Parser::MGI_Desc_Parser';

my $parser = Bio::EnsEMBL::Xref::Parser::MGI_Desc_Parser->new(
  source_id  => 58,
  species_id => 10090,
  files      => ["$Bin/test-data/MRK_List2.rpt"],
  xref_dba   => $xref_dba
);
isa_ok( $parser, 'Bio::EnsEMBL::Xref::Parser::MGI_Desc_Parser' );

$parser->run();

# check xrefs
ok(
  $db->schema->resultset('Xref')->check_direct_xref(
    {
      accession   => 'MGI:1341858',
      label       => '03B03F',
      description => 'DNA segment, 03B03F (Research Genetics)',
      source_id   => 58,
      species_id  => 10090,
      info_type   => 'MISC'
    }
  ),
  'Sample mouse misc description Xref has been inserted'
);

is( $db->schema->resultset('Xref')->count, 10, "All 10 rows were inserted in to Xref" );

# check synonym
ok(
  $db->schema->resultset('Synonym')->check_synonym(
    {
      xref_id => 10,
      synonym => "Ecrg4"
    }
  ),
  'Synonym Ecrg4 has been inserted'
);

is( $db->schema->resultset('Synonym')->count, 2, "Inserted two synonyms" );

done_testing();
