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
package Bio::EnsEMBL::Xref::Test::Parser::UniProtParser;

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

use Bio::EnsEMBL::Xref::Parser::UniProtParser;


Readonly my $SOURCE_ID_CHEMBL                 => 127;
Readonly my $SOURCE_ID_EMBL                   => 128;
Readonly my $SOURCE_ID_MEROPS                 => 129;
Readonly my $SOURCE_ID_PDB                    => 130;
Readonly my $SOURCE_ID_PROTEIN_ID             => 131;
Readonly my $SOURCE_ID_UNIPROT_GN             => 139;
Readonly my $SOURCE_ID_UNIPROT                => 138;
Readonly my $SOURCE_ID_UNIPROT_SP_DIRECT      => 136;
Readonly my $SOURCE_ID_UNIPROT_SP_SEQUENCE    => 137;
Readonly my $SOURCE_ID_UNIPROT_TR_DIRECT      => 132;
Readonly my $SOURCE_ID_UNIPROT_TR_LOWEVIDENCE => 134;
Readonly my $SOURCE_ID_UNIPROT_TR_SEQUENCE    => 133;
Readonly my $SPECIES_ID_HUMAN                 => 9606;

# Needed by the tests of individual ETL modules
Readonly my @EXTRACTOR_UNIPROT_MANDATORY_PREFIXES
  => ( 'ID', 'AC', 'DE', 'OX', 'PE', 'SQ', q{  }, );
Readonly my @EXTRACTOR_UNIPROT_OPTIONAL_PREFIXES
  => ( 'GN', 'DR', 'RG', );
Readonly my $LOADER_BATCH_SIZE
  => 2;


my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();

my $config = $db->config();
my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  'host'   => $config->{host},
  'dbname' => $config->{db},
  'user'   => $config->{user},
  'pass'   => $config->{pass},
  'port'   => $config->{port},
);


my $parser;


subtest 'Missing required source IDs' => sub {

  $parser = Bio::EnsEMBL::Xref::Parser::UniProtParser->new(
    source_id  => $SOURCE_ID_UNIPROT,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/uniprot-mini.dat" ],
    xref_dba   => $xref_dba,
  );
  throws_ok( sub { $parser->run(); },
             qr{ \A No[ ]source_id[ ]for[ ]source_name= }msx,
             'Throws on source IDs missing from DB' );

};


# We will need these later
$db->schema->populate( 'Source', [
  [ qw{ source_id name priority_description priority } ],
  [ $SOURCE_ID_CHEMBL,                 'ChEMBL',            q{},                     1 ],
  [ $SOURCE_ID_EMBL,                   'EMBL',              q{},                     1 ],
  [ $SOURCE_ID_MEROPS,                 'MEROPS',            q{},                     1 ],
  [ $SOURCE_ID_PDB,                    'PDB',               q{},                     1 ],
  [ $SOURCE_ID_PROTEIN_ID,             'protein_id',        q{},                     1 ],
  [ $SOURCE_ID_UNIPROT_GN,             'Uniprot_gn',        q{},                     1 ],
  [ $SOURCE_ID_UNIPROT_SP_DIRECT,      'Uniprot/SWISSPROT', 'direct',                0 ],
  [ $SOURCE_ID_UNIPROT_SP_SEQUENCE,    'Uniprot/SWISSPROT', 'sequence_mapped',       1 ],
  [ $SOURCE_ID_UNIPROT_TR_DIRECT,      'Uniprot/SPTREMBL',  'direct',                0 ],
  [ $SOURCE_ID_UNIPROT_TR_LOWEVIDENCE, 'Uniprot/SPTREMBL',  'protein_evidence_gt_2', 2 ],
  [ $SOURCE_ID_UNIPROT_TR_SEQUENCE,    'Uniprot/SPTREMBL',  'sequence_mapped',       1 ],
] );


subtest 'Problems with input file' => sub {

  $parser = Bio::EnsEMBL::Xref::Parser::UniProtParser->new(
    source_id  => $SOURCE_ID_UNIPROT,
    species_id => $SPECIES_ID_HUMAN,
    files      => [],
    xref_dba   => $xref_dba,
  );
  throws_ok( sub { $parser->run(); },
             qr{ \A No[ ]file[ ]name[ ] }msx,
             'Throws on no file name' );

  $parser = Bio::EnsEMBL::Xref::Parser::UniProtParser->new(
    source_id  => $SOURCE_ID_UNIPROT,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/uniprot-NONEXISTENT.dat" ],
    xref_dba   => $xref_dba,
  );
  throws_ok( sub { $parser->run(); },
             qr{ \A Can[ ]not[ ]open[ ]file[ ] }msx,
             'Throws on no input file missing' );

};


subtest 'Error states in ETL modules' => sub {

  subtest 'Extractor' => sub {

    my $extractor_trunc
      = Bio::EnsEMBL::Xref::Parser::UniProtParser::Extractor->new({
        file_names         => [ "$Bin/test-data/uniprot-truncated.dat" ],
        mandatory_prefixes => \@EXTRACTOR_UNIPROT_MANDATORY_PREFIXES,
        optional_prefixes  => \@EXTRACTOR_UNIPROT_OPTIONAL_PREFIXES,
        species_id         => $SPECIES_ID_HUMAN,
        xref_dba           => $xref_dba,
      });
    throws_ok( sub { $extractor_trunc->get_uniprot_record(); },
              qr{ \A Incomplete[ ]input[ ]record[ ] }msx,
               'Throws on truncated input file' );
    $extractor_trunc->finish();

    my $extractor
      = Bio::EnsEMBL::Xref::Parser::UniProtParser::Extractor->new({
        file_names         => [ "$Bin/test-data/uniprot-badData.dat" ],
        mandatory_prefixes => \@EXTRACTOR_UNIPROT_MANDATORY_PREFIXES,
        optional_prefixes  => \@EXTRACTOR_UNIPROT_OPTIONAL_PREFIXES,
        species_id         => $SPECIES_ID_HUMAN,
        xref_dba           => $xref_dba,
      });

    $extractor->get_uniprot_record();
    throws_ok( sub { $extractor->extract(); },
               # This is a string match, not a "complex regex".
               ## no critic (ProhibitComplexRegexes)
               qr{ \A One[ ]or[ ]more[ ]required[ ]fields[ ]
                  missing[ ]in[ ]record[ ] }msx,
               'Throws missing required field' );

    $extractor->get_uniprot_record();
    throws_ok( sub { $extractor->extract(); },
               qr{ \A Malformed[ ]final-option[ ]match[ ] }msx,
               'Throws on malformed DR line' );

    $extractor->get_uniprot_record();
    throws_ok( sub { $extractor->extract(); },
               qr{ found[ ]\'Synonyms\'[ ]but[ ]no[ ]\'Name\'[ ] }msx,
               'Throws on malformed GN line' );

    $extractor->get_uniprot_record();
    throws_ok( sub { $extractor->extract(); },
               qr{ \A Invalid[ ]entry[ ]status[ ] }msx,
               'Throws on invalid status in ID line' );

    $extractor->get_uniprot_record();
    throws_ok( sub { $extractor->extract(); },
               qr{ \A Invalid[ ]protein[ ]evidence[ ]level[ ] }msx,
               'Throws on invalid evidence level in PE line' );

    $extractor->get_uniprot_record();
    throws_ok( sub { $extractor->extract(); },
               qr{ \A Cannot[ ]use[ ]taxon-DB[ ]qualifier[ ] }msx,
               'Throws on unknown database qualifier in OX line' );

  };

  # There are in principle a few things one could test in Transformer
  # but at this point, it is much more convenient to do so when its
  # output can be accessed using DBIC.

  subtest 'Loader' => sub {
    my $loader;
    my $graph_node_1 = {
      ACCESSION     => 'P31337',
      INFO_TYPE     => 'SEQUENCE_MATCH',
      LABEL         => 'P31337',
      SEQUENCE      => 'GATTACA',
      SEQUENCE_TYPE => 'dna',
      SOURCE_ID     => $SOURCE_ID_UNIPROT,
      SPECIES_ID    => $SPECIES_ID_HUMAN,
    };
    my $graph_node_2 = {
      ACCESSION     => 'P90210',
      INFO_TYPE     => 'SEQUENCE_MATCH',
      LABEL         => 'P90210',
      SEQUENCE      => 'AAA',
      SEQUENCE_TYPE => 'dna',
      SOURCE_ID     => $SOURCE_ID_UNIPROT,
      SPECIES_ID    => $SPECIES_ID_HUMAN,
    };

    $loader
      = Bio::EnsEMBL::Xref::Parser::UniProtParser::Loader->new({
          batch_size         => $LOADER_BATCH_SIZE,
          xref_dba           => $xref_dba,
        });

    $loader->load( $graph_node_1 );
    is( $db->schema->resultset('Xref')->count, 0,
        'Xrefs are not immediately uploaded when in batch mode' );
    $loader->load( $graph_node_2 );
    is( $db->schema->resultset('Xref')->count, $LOADER_BATCH_SIZE,
        'All submitted xrefs uploaded when batch size reached' );

  };

};


subtest 'Successful inserts' => sub {

  $parser = Bio::EnsEMBL::Xref::Parser::UniProtParser->new(
    source_id  => $SOURCE_ID_UNIPROT,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/uniprot-mini.dat" ],
    rel_file   => "$Bin/test-data/uniprot-reldate.txt",
    xref_dba   => $xref_dba,
  );
  isa_ok( $parser, 'Bio::EnsEMBL::Xref::Parser::UniProtParser', 'Instantiated UniProt parser' );
  lives_ok( sub { $parser->run(); }, 'Parsed UniProt data without errors' );

};


subtest 'Primary xrefs' => sub {
  my $rs;
  my $matching_xref;

  is( $db->schema->resultset('Xref')->count({
    source_id => $SOURCE_ID_UNIPROT_SP_SEQUENCE,
  }), 1, 'Expected number of Swiss-Prot primary xrefs' );

  is( $db->schema->resultset('Xref')->count({
    source_id => $SOURCE_ID_UNIPROT_TR_SEQUENCE,
  }), 0, 'Expected number of TrEMBL primary xrefs' );

  is( $db->schema->resultset('Xref')->count({
    source_id => $SOURCE_ID_UNIPROT_TR_LOWEVIDENCE,
  }), 1, 'Expected number of TrEMBL low-evidence primary xrefs' );

  # Remnants of Loader tests conducted above have a different sequence
  # type so we can easily exclude them
  is( $db->schema->resultset('PrimaryXref')->count({
    sequence_type => 'peptide',
  }), 2, 'Expected number of primary-xref sequences' );

  # FIXME: this method is misnamed, it checks xrefs and not *direct* xrefs
  # FIXME: this isn't informative enough - it doesn't show WHERE the mismatch is if there is one
  ok(
     $db->schema->resultset('Xref')->check_direct_xref({
       accession   => 'P31946',
       label       => 'P31946',
       description => '14-3-3 protein beta/alpha 14-3-3 protein beta/alpha, N-terminally processed',
       source_id   => $SOURCE_ID_UNIPROT_SP_SEQUENCE,
       species_id  => $SPECIES_ID_HUMAN,
       info_type   => 'SEQUENCE_MATCH',
     }),
     'Sample primary xref inserted as expected'
   );

  $rs = $db->schema->resultset('Xref')->search({
    accession => 'A0A346RQC0',
  });
  $matching_xref = $rs->next;
  $rs = $db->schema->resultset('PrimaryXref')->search({
    xref_id => $matching_xref->xref_id,
  });
  like( $rs->next->sequence, qr{ \A VWLRRTT }msx,
      'A primary xref points to the expected sequence' );

};

subtest 'Direct xrefs' => sub {
  my $rs;
  my $matching_link;

  is( $db->schema->resultset('TranslationDirectXref')->count, 2,
      'Expected number of translation direct-xref links' );

  # There are two links but they both point to the same xref
  is( $db->schema->resultset('Xref')->count({
    source_id => $SOURCE_ID_UNIPROT_SP_DIRECT,
  }), 1, 'Expected number of Swiss-Prot direct xrefs' );

  is( $db->schema->resultset('Xref')->count({
    source_id => $SOURCE_ID_UNIPROT_TR_DIRECT,
  }), 0, 'Expected number of TrEMBL direct xrefs' );

  # FIXME: this isn't informative enough - it doesn't show WHERE the mismatch is if there is one
  ok(
     $db->schema->resultset('Xref')->check_direct_xref({
       accession   => 'P31946',
       label       => 'P31946',
       description => '14-3-3 protein beta/alpha 14-3-3 protein beta/alpha, N-terminally processed',
       source_id   => $SOURCE_ID_UNIPROT_SP_DIRECT,
       species_id  => $SPECIES_ID_HUMAN,
       info_type   => 'DIRECT',
     }),
     'Sample direct xref inserted as expected'
   );

  # Start on the link side because the corresponding xref has got two
  # links, if we started from there we would have to either force
  # specific order of returned links or check for both of them, in
  # either order.
  $rs = $db->schema->resultset('TranslationDirectXref')->search({
    ensembl_stable_id => 'ENSP00000300161',
  });
  $matching_link = $rs->next;
  $rs = $db->schema->resultset('Xref')->search({
    xref_id => $matching_link->general_xref_id,
  });
  is( $rs->next->accession, 'P31946',
      'Sample direct xref has a matching direct-xref link' );

};

subtest 'Dependent xrefs' => sub {
  my $rs;
  my $matching_link;
  my $matching_xref;
  my $matching_master;

  # Swiss-Prot entry: 1 ChEMBL, 5 EMBL, 5 PDB; 6 protein_id, 1 Uniprot_gn
  # TrEMBL entry: 2 EMBL; 2 protein_id, 1 Uniprot_gn
  is( $db->schema->resultset('DependentXref')->count, 23, 'Expected number of dependent xrefs' );

  # FIXME: this method is misnamed, it checks xrefs and not *direct* xrefs
  # FIXME: this isn't informative enough - it doesn't show WHERE the mismatch is if there is one
  ok(
     $db->schema->resultset('Xref')->check_direct_xref({
       accession   => 'A0A346RQC0',
       label       => 'amoA',
       description => undef,
       source_id   => $SOURCE_ID_UNIPROT_GN,
       species_id  => $SPECIES_ID_HUMAN,
       info_type   => 'DEPENDENT',
     }),
     'Sample dependent xref inserted as expected'
   );

  # FIXME: we might want to add an extra test to confirm we correctly
  # assemble protein_id xrefs from EMBL/ChEMBL ones. Leave that for
  # later though, there are multiple ones of each so it would take a
  # while to do it properly.

  $rs = $db->schema->resultset('Xref')->search({
    accession => 'CHEMBL3710403',
  });
  $matching_xref = $rs->next;
  $rs = $db->schema->resultset('DependentXref')->search({
    dependent_xref_id => $matching_xref->xref_id,
  });
  $matching_link = $rs->next;
  $rs = $db->schema->resultset('Xref')->search({
    xref_id => $matching_link->master_xref_id,
  });
  $matching_master = $rs->next;
  is( $matching_master->accession, 'P31946',
      'Sample dependent xref links correctly to the primary' );
  is( $matching_link->linkage_source_id, $matching_master->source_id,
      'linkage_source_id is the primary\'s source_id' );
  is( $matching_link->linkage_annotation, $matching_xref->source_id,
      'linkage_annotation is the dependent\'s source_id' );

};


subtest 'Synonyms' => sub {
  my $rs;
  my $matching_synonym;

  # Swiss-Prot entry has two synonyms, TrEMBL one has one
  is( $db->schema->resultset('Synonym')->count, 2, 'Expected number of synonyms' );

  $rs = $db->schema->resultset('Synonym')->search({
    synonym => 'E1P616',
  });
  $matching_synonym = $rs->next;
  $rs = $db->schema->resultset('Xref')->search({
    xref_id => $matching_synonym->xref_id,
  });
  is( $rs->next->accession, 'P31946',
      'Sample direct xref has a matching synonym' );

};


subtest 'Release file' => sub {
  Readonly my $SWISSPROT_RELEASE => 'UniProtKB/Swiss-Prot Release 2018_11 of 05-Dec-2018';
  Readonly my $TREMBL_RELEASE    => 'UniProtKB/TrEMBL Release 2018_11 of 05-Dec-2018';
  my $rs;
  my $matching_source;

  foreach my $source_id ( $SOURCE_ID_UNIPROT_SP_DIRECT, $SOURCE_ID_UNIPROT_SP_SEQUENCE ) {
    $rs = $db->schema->resultset('Source')->search({
      source_id => $source_id,
    });
    $matching_source = $rs->next;
    is( $matching_source->source_release, $SWISSPROT_RELEASE,
        "Expected release string set for Swiss-Prot source $source_id" );
  }
  foreach my $source_id ( $SOURCE_ID_UNIPROT_TR_DIRECT, $SOURCE_ID_UNIPROT_TR_LOWEVIDENCE, $SOURCE_ID_UNIPROT_TR_SEQUENCE ) {
    $rs = $db->schema->resultset('Source')->search({
      source_id => $source_id,
    });
    $matching_source = $rs->next;
    is( $matching_source->source_release, $TREMBL_RELEASE,
        "Expected release string set for TrEMBL source $source_id" );
  }

};



subtest 'Replay safety' => sub {

  my $total_xref_count = $db->schema->resultset('Xref')->count;
  my $dependent_xref_link_count = $db->schema->resultset('DependentXref')->count;

  $parser = Bio::EnsEMBL::Xref::Parser::UniProtParser->new(
    source_id  => $SOURCE_ID_UNIPROT,
    species_id => $SPECIES_ID_HUMAN,
    files      => [ "$Bin/test-data/uniprot-mini.dat" ],
    xref_dba   => $xref_dba,
  );
  lives_ok( sub { $parser->run(); }, 'Re-parsed UniProt data without errors' );

  is( $db->schema->resultset('Xref')->count(), $total_xref_count,
      'No new xrefs inserted by the replay' );

  is( $db->schema->resultset('PrimaryXref')->count({
    sequence_type => 'peptide',
  }), 2, 'No new peptide primary-xref sequences inserted by the replay' );

  is( $db->schema->resultset('TranslationDirectXref')->count(), 2,
      'No new translation direct-xref links inserted by the replay' );

  is( $db->schema->resultset('DependentXref')->count(), $dependent_xref_link_count,
      'No new dependent-xref links inserted by the replay' );

  is( $db->schema->resultset('Synonym')->count, 2,
      'No new synonyms inserted by the replay' );

  # Ideally we would also make sure the replay has not modified
  # existing entries, no quick way of doing so though.

};


done_testing();


1;
