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

use FindBin '$Bin';
use lib "$Bin/";

use Bio::EnsEMBL::Xref::Test::TestDB;
use Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;

use_ok 'Bio::EnsEMBL::Xref::Parser';

my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();
my %config = %{ $db->config };

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config{host},
  dbname => $config{db},
  user   => $config{user},
  pass   => $config{pass},
  port   => $config{port}
);

use_ok 'Bio::EnsEMBL::Xref::Parser::EntrezGeneParser';

# add EntrezGene/WikiGene source to the db
$db->schema->populate(
  'Source',
  [
   [ qw/name status ordered/ ],
   [ 'EntrezGene', 'KNOWN', 10 ],
   [ 'WikiGene', 'KNOWN', 10 ],
  ]
);

my ($source_id, $species_id) = ( 1, 9606 );

throws_ok {
  Bio::EnsEMBL::Xref::Parser::EntrezGeneParser->new(
    species_id => $species_id,
    files      => [ "$Bin/test-data/gene_info.gz" ],
    xref_dba   => $xref_dba)->run()
} qr/source_id/, 'Run without source';

throws_ok {
  Bio::EnsEMBL::Xref::Parser::EntrezGeneParser->new(
    source_id  => $source_id,
    files      => [ "$Bin/test-data/gene_info.gz" ],
    xref_dba   => $xref_dba)->run()
} qr/species_id/, 'Run without species_id';

throws_ok {
  Bio::EnsEMBL::Xref::Parser::EntrezGeneParser->new(
    source_id  => $source_id,
    species_id => $species_id,
    xref_dba   => $xref_dba)->run()
} qr/files/, 'Run without files';


my $parser = Bio::EnsEMBL::Xref::Parser::EntrezGeneParser->new(
 source_id  => $source_id,
 species_id => $species_id,
 files      => [ "$Bin/test-data/gene_info.gz" ],
 xref_dba   => $xref_dba
);

isa_ok( $parser, 'Bio::EnsEMBL::Xref::Parser::EntrezGeneParser' );

$parser->run();

my @xrefs = $db->schema->resultset('Xref')->search(
  {
   accession   => "10",
   label       => 'NAT2',
   description => 'N-acetyltransferase 2',
   species_id  => $species_id,
  }
);

# parser insert both EntrezGene and WikiGene xrefs
is( scalar @xrefs, 2, 'Sample EntrezGene/WikiGene xrefs inserted' );

# check synonyms have been inserted as well
foreach my $syn ( qw/ AAC2 NAT-2 PNAT / ) {
  ok( $db->schema->resultset('Synonym')->search( { syn => $syn } ), 'Synonym inserted' );
} 

@xrefs = $db->schema->resultset('Xref')->search(
  {
   accession   => "13",
   label       => 'AADAC',
   description => 'arylacetamide deacetylase',
   species_id  => $species_id,
  }
);

is( scalar @xrefs, 2, 'Sample EntrezGene/WikiGene xrefs inserted' );

foreach my $syn ( qw/ CES5A1 DAC / ) {
  ok( $db->schema->resultset('Synonym')->search( { syn => $syn } ), 'Synonym inserted' );
} 

# Test if all 10 entries were inserted (twice)
is($db->schema->resultset('Xref')->count, 20, "All entries inserted");

done_testing();

