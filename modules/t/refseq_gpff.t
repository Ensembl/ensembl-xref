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

use Bio::EnsEMBL::Xref::Test::TestDB;
use Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Test::MultiTestDB;

use_ok 'Bio::EnsEMBL::Xref::Parser';


my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();

my $config = $db->config;

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config->{host},
  dbname => $config->{db},
  user   => $config->{user},
  pass   => $config->{pass},
  port   => $config->{port}
);

# populate used sources
$db->schema->resultset('Source')->populate([
  {
    source_id            => 95,
    name                 => 'RefSeq_mRNA',
    priority_description => 'refseq'
  },{
    source_id            => 97,
    name                 => 'RefSeq_mRNA_predicted',
    priority_description => 'refseq'
  },{
    source_id            => 99,
    name                 => 'RefSeq_ncRNA',
    priority_description => 'refseq'
  },{
    source_id            => 101,
    name                 => 'RefSeq_ncRNA_predicted',
    priority_description => 'refseq'
  },{
    source_id            => 109,
    name                 => 'RefSeq_peptide',
    priority_description => 'refseq'
  },{
    source_id            => 110,
    name                 => 'RefSeq_peptide_predicted',
    priority_description => 'refseq'
  },{
    source_id            => 23,
    name                 => 'EntrezGene'
  },{
    source_id            => 146,
    name                 => 'WikiGene'
  }
]);

use_ok 'Bio::EnsEMBL::Xref::Parser::RefSeqGPFFParser';

my $parser_refseq_peptide = Bio::EnsEMBL::Xref::Parser::RefSeqGPFFParser->new(
 source_id  => 108,
 species_id => 9606,
 species    => 'homo_sapiens',
 files      => ["$Bin/test-data/refseq_gpff.protein.gpff"],
 rel_file   => "$Bin/test-data/refseq_gpff.release.txt",
 xref_dba   => $xref_dba,
 verbose    => 1
);

isa_ok( $parser_refseq_peptide, 'Bio::EnsEMBL::Xref::Parser::RefSeqGPFFParser' );

$parser_refseq_peptide->run();

# Test if all the rows were inserted
is($db->schema->resultset('Xref')->count, 3, "All 3 relevant RefSeq_peptide xrefs were inserted");

# Check the data of the first xref
ok(
  $db->schema->resultset('Xref')->check_direct_xref({
    accession   => 'NP_055050',
    label       => 'NP_055050.1',
    description => 'ubiquitin-like protein 4A [Homo sapiens].',
    source_id   => 109,
    species_id  => 9606
  }),
 'Sample human RefSeq_peptide direct Xref has been inserted'
);

# Check if the release file is correctly used to update source release
my $source_release = $db->schema->resultset('Source')->search({
  source_id => 95,
})->single->source_release;

is($source_release, 'NCBI Reference Sequence (RefSeq) Database Release 91, November 5, 2018', "Source release has been updated");

my $parser_refseq_dna = Bio::EnsEMBL::Xref::Parser::RefSeqGPFFParser->new(
 source_id  => 108,
 species_id => 9598,
 species    => 'pan_troglodytes',
 files      => ["$Bin/test-data/refseq_gpff.rna.gbff"],
 rel_file   => "$Bin/test-data/refseq_gpff.release.txt",
 xref_dba   => $xref_dba,
);

$parser_refseq_dna->run();

# Test if all the rows were inserted
is($db->schema->resultset('Xref')->count, 7, "All 4 relevant RefSeq_dna xrefs were inserted");

# Check the data of the first xref
ok(
  $db->schema->resultset('Xref')->check_direct_xref({
    accession   => 'NM_001009111',
    label       => 'NM_001009111.1',
    description => 'Pan troglodytes olfactory receptor OR93Ch (OR93CH), mRNA.',
    source_id   => 95,
    species_id  => 9598
  }),
 'Sample chimpanzee RefSeq_dna direct Xref has been inserted'
);

# Check the if missing rel_file throws error
my $parser_no_release = Bio::EnsEMBL::Xref::Parser::RefSeqGPFFParser->new(
 source_id  => 108,
 species_id => 9606,
 species    => 'homo_sapiens',
 files      => ["$Bin/test-data/refseq_gpff.protein.gpff"],
 xref_dba   => $xref_dba,
);

throws_ok{ $parser_no_release->run() } qr/Need to pass .+ and rel_file as pairs/, 'No release file provided throws error' ;

done_testing();

