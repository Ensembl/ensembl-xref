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
use Bio::EnsEMBL::Test::MultiTestDB;

use_ok 'Bio::EnsEMBL::Xref::Parser';


my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $multi_db->get_DBAdaptor('otherfeatures');
my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();

my $config = $db->config;

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config->{host},
  dbname => $config->{db},
  user   => $config->{user},
  pass   => $config->{pass},
  port   => $config->{port}
);

$db->schema->resultset('Species')->create({
  species_id  => 9606,
  taxonomy_id => 9606,
  name        => 'homo_sapiens'
});

# populate used sources
$db->schema->resultset('Source')->populate([
  {
    source_id            => 95,
    name                 => 'RefSeq_mRNA',
    priority_description => 'otherfeatures'
  },{
    source_id            => 98,
    name                 => 'RefSeq_mRNA_predicted',
    priority_description => 'otherfeatures'
  },{
    source_id            => 100,
    name                 => 'RefSeq_ncRNA',
    priority_description => 'otherfeatures'
  },{
    source_id            => 102,
    name                 => 'RefSeq_ncRNA_predicted',
    priority_description => 'otherfeatures'
  },{
    source_id            => 109,
    name                 => 'RefSeq_peptide',
    priority_description => 'otherfeatures'
  },{
    source_id            => 111,
    name                 => 'RefSeq_peptide_predicted',
    priority_description => 'otherfeatures'
  },{
    source_id            => 23,
    name                 => 'EntrezGene'
  }
]);

use_ok 'Bio::EnsEMBL::Xref::Parser::RefSeqCoordinateParser';

my $parser = Bio::EnsEMBL::Xref::Parser::RefSeqCoordinateParser->new(
 source_id  => 94,
 species_id => 9606,
 files      => [],
 xref_dba   => $xref_dba,
 dba        => $dba,
 verbose    => 1
);

isa_ok( $parser, 'Bio::EnsEMBL::Xref::Parser::RefSeqCoordinateParser' );

$parser->run();


# TODO

# # Test if all the rows were inserted
# is($db->schema->resultset('Xref')->count, 6, "All 6 rows were inserted");

# ok(
#   $db->schema->resultset('Xref')->check_direct_xref({
#     accession   => 'VGNC:14659',
#     label       => 'CYYR1',
#     description => 'cysteine and tyrosine rich 1',
#     source_id   => 144,
#     species_id  => 9598
#   }),
#  'Sample chimpanzee direct Xref has been inserted'
# );


my $parser_no_file = Bio::EnsEMBL::Xref::Parser::RefSeqCoordinateParser->new(
 source_id  => 94,
 species_id => 9606,
 files      => [],
 xref_dba   => $xref_dba
);

throws_ok{ $parser_no_file->run() } qr/Missing otherfeatures database adaptor/, 'No dba provided throws error' ;

done_testing();

