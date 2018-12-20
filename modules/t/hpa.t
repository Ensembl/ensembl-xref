
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

use_ok 'Bio::EnsEMBL::Xref::Parser::HPAParser';

# initialize the parser
my $parser = Bio::EnsEMBL::Xref::Parser::HPAParser->new(
    source_id  => 11,
    species_id => 9601,
    files      => ["$Bin/test-data/hpa.txt"],
    xref_dba   => $xref_dba
);

isa_ok( $parser, 'Bio::EnsEMBL::Xref::Parser::HPAParser' );

$parser->run();

ok(
    $db->schema->resultset('Xref')->check_direct_xref(
        {
            accession  => "1",
            label      => 'CAB000001',
            info_type  => 'DIRECT',
            source_id  => 11,
            species_id => 9601
        }
    ),
    'Sample hpa direct Xref has been inserted'
);

# Test if all the distinct rows were inserted in to xref
is( $db->schema->resultset('Xref')->count,
    2, "2 rows with distinct accessions were inserted" );

# Test if all the rows were inserted in to translation_direct_xref
is( $db->schema->resultset('TranslationDirectXref')->count,
    10, "10 rows were inserted in to translation_direct_xref" );

done_testing();

