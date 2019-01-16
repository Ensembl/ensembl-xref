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
use Bio::EnsEMBL::Xref::Mapper;
use Bio::EnsEMBL::Xref::Mapper::CoreInfo;

use Bio::EnsEMBL::Test::MultiTestDB;

use_ok 'Bio::EnsEMBL::Xref::Mapper';

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new;
my $dba = $multi_db->get_DBAdaptor('core');

my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();
ok($db, 'Xref DB');
my %config = %{ $db->config };

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config{host},
  dbname => $config{db},
  user   => $config{user},
  pass   => $config{pass},
  port   => $config{port}
);

throws_ok {
  Bio::EnsEMBL::Xref::Mapper->new()
} qr/Required/, 'Throws without required arguments';

throws_ok {
  Bio::EnsEMBL::Xref::Mapper->new( xref_dba => $xref_dba )
} qr/Required/, 'Throws without required arguments';

throws_ok {
  Bio::EnsEMBL::Xref::Mapper->new( core_dba => $dba )
} qr/Required/, 'Throws without required arguments';

throws_ok {
  Bio::EnsEMBL::Xref::Mapper->new( xref_dba => [], core_dba => $dba )
} qr/was expected/, 'Throws with argument not of the correct class';

throws_ok {
  Bio::EnsEMBL::Xref::Mapper->new( xref_dba => $xref_dba, core_dba => [] )
} qr/was expected/, 'Throws with argument not of the correct class';


my $mapper = Bio::EnsEMBL::Xref::Mapper->new( xref_dba => $xref_dba, core_dba => $dba );
isa_ok( $mapper, 'Bio::EnsEMBL::Xref::Mapper' );

is( $mapper->xref, $xref_dba, 'Correct xref dba assignment' );
is( $mapper->core, $dba, 'Correct core dba assignment' );

ok( !$mapper->get_official_name, 'Just returns' );

my $desc = 'BA392M18.1 (NOVEL PROTEIN SIMILAR TO TESTIS SPECIFIC PROTEIN TSPY). [Source:SPTREMBL;Acc:Q9H489]';
my $filtered_desc = $mapper->filter_by_regexp( $desc, gene_description_filter_regexps() );
is( $filtered_desc, $desc, 'Unfiltered description' );
$desc = 'weird stuff';
$filtered_desc = $mapper->filter_by_regexp( $desc, gene_description_filter_regexps() );
is( $filtered_desc, 'stuff', 'Filtered description' );

# Add an example species to the db
my $species = $db->schema->resultset('Species')->create({
  species_id           => 10,
  taxonomy_id          => 9606,
  name                 => 'Homo sapiens'
});
is($species->name, 'Homo sapiens', 'Species created in the DB');

throws_ok {
  $mapper->get_species_id_from_species_name('Mus musculus')
} qr/It must be one of/, 'Throws with unknown species name';
is($mapper->get_species_id_from_species_name('Homo sapiens'), 10, 'Species ID from name');

sub gene_description_filter_regexps {
  return [ '^BA\S+\s+\(NOVEL PROTEIN\)\.?', 'weird\s*' ];
}

done_testing();

