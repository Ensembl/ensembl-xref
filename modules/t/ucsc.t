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

use_ok 'Bio::EnsEMBL::Xref::Parser::UCSC_human_parser';

# add EntrezGene/WikiGene source to the db
$db->schema->populate(
  'Source',
  [
   [ qw/name/ ],
   [ 'UCSC_human' ],
  ]
);

my ($source_id, $species_id) = ( 1, 9606 );

my $parser = Bio::EnsEMBL::Xref::Parser::UCSC_human_parser->new(
 source_id  => $source_id,
 species_id => $species_id,
 files      => [ "$Bin/test-data/ucsc.txt" ],
 xref_dba   => $xref_dba
);

isa_ok( $parser, 'Bio::EnsEMBL::Xref::Parser::UCSCParser' );

$parser->run();

my @xrefs = $db->schema->resultset('CoordinateXref')->search(
  {
   accession   => "uc031tla.1",
   chromosome  => '1',
   txStart     => '17368',
   txEnd       => '17436',
   strand      => '-1',
   exonStarts  => '17368',
   exonEnds    => '17436',
   species_id  => $species_id,
  }
);

# Test if all 10 entries were inserted
is($db->schema->resultset('CoordinateXref')->count, 10, "All entries inserted");

done_testing();

