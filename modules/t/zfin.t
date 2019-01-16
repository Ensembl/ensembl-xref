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
use FindBin '$Bin';

use Bio::EnsEMBL::Xref::Test::TestDB;
use Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::Xref::Parser::ZFINDescParser;

use_ok 'Bio::EnsEMBL::Xref::Parser::ZFIN';

my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();

my $config = $db->config;

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config->{host},
  dbname => $config->{db},
  user   => $config->{user},
  pass   => $config->{pass},
  port   => $config->{port}
);

# create ZFIN source and some corresponding xrefs
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

note 'Test the ZFIN parser base class methods';
my $parser = Bio::EnsEMBL::Xref::Parser::ZFIN->new(
  xref_dba => $xref_dba,
  source_id => $zfin_source->source_id,
  species_id => 7955,
  files => [ 'dummy.txt' ]
);

isa_ok( $parser, 'Bio::EnsEMBL::Xref::Parser::ZFIN' );

is( $parser->loaded(), 6, 'xrefs loaded' );
is_deeply( $parser->description(), { 'ZDB-GENE-030131-3003' => 'HNF1 homeobox Bb',
				     'ZDB-GENE-030131-1077' => 'hepatocyte nuclear factor 4, alpha',
				     'ZDB-GENE-040718-488' => 'WD repeat domain, phosphoinositide interacting 2',
				     'ZDB-GENE-070117-2473' => 'wirbel',
				     'ZDB-GENE-070117-2162' => 'lawrence welk',
				     'ZDB-GENE-040426-2161' => 'wntless Wnt ligand secretion mediator' }, 'description map' );

done_testing();
