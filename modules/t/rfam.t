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
use FindBin '$Bin';
use File::Spec;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Xref::Parser::RFAMParser;

use Bio::EnsEMBL::Xref::Test::TestDB;
use Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis;

use_ok 'Bio::EnsEMBL::Xref::Parser';

my $test_db     = Bio::EnsEMBL::Xref::Test::TestDB->new();
my %config = %{ $test_db->config };

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
    host   => $config{host},
    dbname => $config{db},
    user   => $config{user},
    pass   => $config{pass},
    port   => $config{port}
);

my $core_test_db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $core_dba = $core_test_db->get_DBAdaptor('core');

my $test_file = File::Spec->catfile($Bin,'test-data','rfam.seed');

my $parser = Bio::EnsEMBL::Xref::Parser::RFAMParser->new(
  source_id => 1,
  species_id => 1,
  files => [$test_file],
  xref_dba => $xref_dba,
  dba => $core_dba,
);

populate_test_supporting_features($core_dba);

my $rfam_mappings = $parser->get_rfam_supporting_features($core_dba);

is ($rfam_mappings->{RF00000000}[0], 'ENST00000000001', 'Basic retrieval of supporting features for Rfam entities');
# FIXME Needs more testing, but at least the basics work

my $meta = $parser->parse_stockholm_format($test_file);

is($meta->[0]{AC}, 'RF00001', 'First RFAM record parsed');
is($meta->[1]{AC}, 'RF00002', 'Second RFAM record parsed');


# Now that parsing and DB fetching have been tested individually, we can call run() and check the DB content

$parser->run();

my $result = $test_db->schema->resultset('Xref')->find({
  accession => 'RF00001'
});

is ($result->label,'5S_rRNA','Attributes about first RFAM match stored');
is ($result->description,'5S ribosomal RNA','Attributes about first RFAM match stored');

$result = $test_db->schema->resultset('Xref')->find({
  accession => 'RF00002'
});

ok(!$result, 'No stored xref for RFAM ID that did not have a supporting feature');

done_testing();



# Generates features and DNAAlignFeatures that will match a custom RFAM alignment
sub populate_test_supporting_features {
  my $dba = shift;

  my $slice_adaptor = $dba->get_SliceAdaptor;

  # Create an alignment and an analysis to make it findable
  my $slice = $slice_adaptor->fetch_by_region('chromosome',1,100000,101000,1);

  my $analysis_adaptor = $dba->get_AnalysisAdaptor;
  
  my $rfam_analysis = Bio::EnsEMBL::Analysis->new(
    -logic_name => 'rfamblast',
  );
  my $dbid = $analysis_adaptor->store($rfam_analysis);

  my $alignment = Bio::EnsEMBL::DnaDnaAlignFeature->new(
    -slice => $slice,
    -start => 1,
    -end => 800,
    -strand => 1,
    -hseqname => 'RF00000000',  # a synthetic RFAM ID
    -hstart => 51,
    -hend => 851,
    -hstrand => 1,
    -analysis => $rfam_analysis,
    -align_type => 'cigar', # Doesn't really matter here, but needs to be something
    -cigar_string => 'STOGIE'
  );

  my $alignment_adaptor = $dba->get_DnaAlignFeatureAdaptor;

  $alignment_adaptor->store($alignment);

  # Add some alignment features to match the Rfam.seed file

  my $rfam_matching_alignment = Bio::EnsEMBL::DnaDnaAlignFeature->new(
    -slice => $slice,
    -start => 1,
    -end => 800,
    -strand => 1,
    -hseqname => 'RF00001',  # a fake RFAM ID from the test file
    -hstart => 51,
    -hend => 851,
    -hstrand => 1,
    -analysis => $rfam_analysis,
    -align_type => 'cigar',
    -cigar_string => 'CUBAN'
  );

  $alignment_adaptor->store(
    $rfam_matching_alignment  
  );

  # Not storing RF00002, so that we have a negative case

  # Now create a non-coding transcript feature that matches coordinates with the alignment

  my $gene_analysis = Bio::EnsEMBL::Analysis->new(
    -logic_name => 'ncRNA',
  );
  $dbid = $analysis_adaptor->store($gene_analysis);


  my $exon = Bio::EnsEMBL::Exon->new(
    -slice => $slice,
    -start => 1,
    -end => 800,
    -strand => 1,
    -stable_id => 'ENSE00000000001',
    -phase => 0,
    -end_phase => 2
  );

  my $transcript = Bio::EnsEMBL::Transcript->new(
    -exons => [ $exon ],
    -stable_id => 'ENST00000000001',
    -external_name => 'TheWorst',
    -biotype => 'pseudogene',
    -source => 'testing',
    -analysis => $gene_analysis
  );

  my $transcript_adaptor = $dba->get_TranscriptAdaptor;
  $transcript_adaptor->store($transcript);

  my $gene = Bio::EnsEMBL::Gene->new(
    -start => 1,
    -end => 800,
    -strand => 1,
    -slice => $slice,
    -stable_id => 'ENSG000000001',
    -biotype => 'pseudogene',
    -source => 'testing',
    -analysis => $gene_analysis,
    -transcripts => [$transcript]
  );
  my $gene_adaptor = $dba->get_GeneAdaptor;

  $gene_adaptor->store($gene);


  # Now create a supporting feature that links the alignment to the exon
  my $sf_adaptor = $dba->get_SupportingFeatureAdaptor;
  $sf_adaptor->store($exon->dbID, [$alignment, $rfam_matching_alignment]);

  return;
};