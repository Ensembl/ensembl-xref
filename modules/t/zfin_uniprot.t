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

use Bio::EnsEMBL::Xref::Parser::ZFINDescParser;

use_ok 'Bio::EnsEMBL::Xref::Parser::ZFIN::UniprotParser';

my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();

my $config = $db->config;

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config->{host},
  dbname => $config->{db},
  user   => $config->{user},
  pass   => $config->{pass},
  port   => $config->{port}
);

note 'Create ZFIN source and some corresponding xrefs';
my $zfin_source = $db->create_db_row('Source', {
  name => 'zfin_id',
  priority_description => 'Test ZFIN source',
  priority => 10
});

my ( $source_id, $species_id ) = ( $zfin_source->source_id, 7955 );

note 'Insert some xrefs using the zfin_desc parser';
my $zfin_desc_parser = Bio::EnsEMBL::Xref::Parser::ZFINDescParser->new(
 source_id  => $zfin_source->source_id,
 species_id => 7955,
 files      => ["$Bin/test-data/zfin_desc.txt"],
 xref_dba   => $xref_dba,
 verbose    => 1
);

$zfin_desc_parser->run();

note 'Create Uniprot source and a corresponding xref';
my $uniprot_source = $db->create_db_row('Source', {
  name => 'uniprot/swissprot',
  priority_description => 'No idea',
  priority => 20
});

# Create Uniprot xref to which ZFIN test data refers
my $uniprot_xref = $db->create_db_row('Xref', {
  accession => 'Q8UW00',
  version => 1,
  label => 'Something irrelevant',
  source_id => $uniprot_source->source_id,
  species_id => $species_id,
  info_type => 'COORDINATE_OVERLAP',
  info_text => '',
});

throws_ok {
  Bio::EnsEMBL::Xref::Parser::ZFIN::UniprotParser->new(
    species_id => $species_id,
    files      => [ "$Bin/test-data/zfin/uniprot.txt" ],
    xref_dba   => $xref_dba )->run()
} qr/source_id/, 'Run without source_id';

throws_ok {
  Bio::EnsEMBL::Xref::Parser::ZFIN::UniprotParser->new(
    source_id  => $source_id,
    files      => [ "$Bin/test-data/zfin/uniprot.txt" ],
    xref_dba   => $xref_dba )->run()
} qr/species_id/, 'Run without species_id';

throws_ok {
  Bio::EnsEMBL::Xref::Parser::ZFIN::UniprotParser->new(
    source_id  => $source_id,
    species_id => $species_id,
    xref_dba   => $xref_dba )->run()
} qr/files/, 'Run without files';

my $parser = Bio::EnsEMBL::Xref::Parser::ZFIN::UniprotParser->new(
  xref_dba => $xref_dba,
  source_id => $zfin_source->source_id,
  species_id => 7955,
  files => [ "$Bin/test-data/zfin/uniprot.txt" ]
);

isa_ok( $parser, 'Bio::EnsEMBL::Xref::Parser::ZFIN::UniprotParser' );

$parser->run();

# test ZFIN dependent xref has been assigned to refseq xref above
my $hit = $db->schema->resultset( 'DependentXref' )->fetch_dependent_xref( 'Q8UW00', 'ZDB-GENE-030131-3003' );

ok( $hit, 'Found a dependent Xref' );
is( $hit->dependent_xref->accession, 'ZDB-GENE-030131-3003', 'ZFIN xref is dependent on Uniprot accession' );

is( $db->schema->resultset('DependentXref')->count, 1, "All rows were inserted into dependent xref table" );

done_testing();
