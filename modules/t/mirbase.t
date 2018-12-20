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

## no critic 'RequireFilenameMatchesPackage'
package Bio::EnsEMBL::Xref::Test::Parser::miRBaseParser;

use strict;
use warnings;

use Test::More;
use Test::Exception;
use Test::Warnings 'warnings';

use English '-no_match_vars';
use FindBin '$Bin';
use Readonly;

use Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Xref::Test::TestDB;

use Bio::EnsEMBL::Xref::Parser::miRBaseParser;


Readonly my $SOURCE_ID_MIRBASE   => 169;
Readonly my $SPECIES_ID_HUMAN    => 9606;
Readonly my $SPECIES_NAME_HUMAN  => 'homo_sapiens';


my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();

my $config = $db->config();
my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  'host'   => $config->{host},
  'dbname' => $config->{db},
  'user'   => $config->{user},
  'pass'   => $config->{pass},
  'port'   => $config->{port},
);


# This is for the time being necessary for Transformer to initialise
Readonly my $SOURCE_ID_UNIPROT_SP_DIRECT      => 136;
Readonly my $SOURCE_ID_UNIPROT_SP_SEQUENCE    => 137;
Readonly my $SOURCE_ID_UNIPROT_TR_DIRECT      => 132;
Readonly my $SOURCE_ID_UNIPROT_TR_LOWEVIDENCE => 134;
Readonly my $SOURCE_ID_UNIPROT_TR_SEQUENCE    => 133;
$db->schema->populate( 'Source', [
  [ qw{ source_id name priority_description priority } ],
  [ $SOURCE_ID_UNIPROT_SP_DIRECT,      'Uniprot/SWISSPROT', 'direct',                0 ],
  [ $SOURCE_ID_UNIPROT_SP_SEQUENCE,    'Uniprot/SWISSPROT', 'sequence_mapped',       1 ],
  [ $SOURCE_ID_UNIPROT_TR_DIRECT,      'Uniprot/SPTREMBL',  'direct',                0 ],
  [ $SOURCE_ID_UNIPROT_TR_LOWEVIDENCE, 'Uniprot/SPTREMBL',  'protein_evidence_gt_2', 2 ],
  [ $SOURCE_ID_UNIPROT_TR_SEQUENCE,    'Uniprot/SPTREMBL',  'sequence_mapped',       1 ],
] );


my $parser;


subtest 'Problems with input file' => sub {

  $parser = Bio::EnsEMBL::Xref::Parser::miRBaseParser->new(
    source_id  => $SOURCE_ID_MIRBASE,
    species_id => $SPECIES_ID_HUMAN,
    species    => $SPECIES_NAME_HUMAN,
    files      => [],
    xref_dba   => $xref_dba,
  );
  throws_ok( sub { $parser->run(); },
             qr{ \A No[ ]file[ ]name[ ] }msx,
             'Throws on no file name' );

  $parser = Bio::EnsEMBL::Xref::Parser::miRBaseParser->new(
    source_id  => $SOURCE_ID_MIRBASE,
    species_id => $SPECIES_ID_HUMAN,
    species    => $SPECIES_NAME_HUMAN,
    files      => [ "$Bin/test-data/mirbase-NONEXISTENT.dat" ],
    xref_dba   => $xref_dba,
  );
  throws_ok( sub { $parser->run(); },
             qr{ \A Can[ ]not[ ]open[ ]file[ ] }msx,
             'Throws on no input file missing' );

};


subtest 'Successful inserts' => sub {

  $parser = Bio::EnsEMBL::Xref::Parser::miRBaseParser->new(
    source_id  => $SOURCE_ID_MIRBASE,
    species_id => $SPECIES_ID_HUMAN,
    species    => $SPECIES_NAME_HUMAN,
    files      => [ "$Bin/test-data/mirbase-mini.dat" ],
    xref_dba   => $xref_dba,
  );
  isa_ok( $parser, 'Bio::EnsEMBL::Xref::Parser::miRBaseParser', 'Instantiated miRBase parser' );
  lives_ok( sub { $parser->run(); }, 'Parsed miRBase data without errors' );

};


subtest 'Primary xrefs' => sub {
  my $rs;
  my $matching_xref;

  is( $db->schema->resultset('Xref')->count, 2,
  'Expected number of miRBase xrefs' );

  is( $db->schema->resultset('PrimaryXref')->count, 2,
  'Expected number of primary-xref sequences' );

  # FIXME: this method is misnamed, it checks xrefs and not *direct* xrefs
  # FIXME: this isn't informative enough - it doesn't show WHERE the mismatch is if there is one
  ok(
     $db->schema->resultset('Xref')->check_direct_xref({
       accession   => 'MI0000061',
       label       => 'hsa-let-7a-2',
       description => 'hsa-let-7a-2',
       source_id   => $SOURCE_ID_MIRBASE,
       species_id  => $SPECIES_ID_HUMAN,
       info_type   => 'SEQUENCE_MATCH',
     }),
     'Sample primary xref inserted as expected'
   );

  $rs = $db->schema->resultset('Xref')->search({
    accession => 'MI0000061',
  });
  $matching_xref = $rs->next;
  $rs = $db->schema->resultset('PrimaryXref')->search({
    xref_id => $matching_xref->xref_id,
  });
  like( $rs->next->sequence, qr{ \A AGGTTGAGGT }msx,
      'A primary xref points to the expected sequence' );

};

subtest 'Replay safety' => sub {

  my $total_xref_count = $db->schema->resultset('Xref')->count;

  $parser = Bio::EnsEMBL::Xref::Parser::miRBaseParser->new(
    source_id  => $SOURCE_ID_MIRBASE,
    species_id => $SPECIES_ID_HUMAN,
    species    => $SPECIES_NAME_HUMAN,
    files      => [ "$Bin/test-data/mirbase-mini.dat" ],
    xref_dba   => $xref_dba,
  );
  lives_ok( sub { $parser->run(); }, 'Re-parsed miRBase data without errors' );

  is( $db->schema->resultset('Xref')->count(), $total_xref_count,
      'No new xrefs inserted by the replay' );

  is( $db->schema->resultset('PrimaryXref')->count, $total_xref_count,
      'No new primary-xref sequences inserted by the replay' );


  # Ideally we would also make sure the replay has not modified
  # existing entries, no quick way of doing so though.

};


done_testing();


1;
