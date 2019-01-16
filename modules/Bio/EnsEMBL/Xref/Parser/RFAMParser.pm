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

=head1 NAME

Bio::EnsEMBL::Xref::Parser::RFAMParser

=head1 DESCRIPTION

A parser class to parse the RFAM seed files in Stockholm format.

It mostly throws away the content, we use it solely as an authoritative list
of valid RFAM IDs and their description. We then rely on the Genebuild 
alignments of RFAM sequences to decide which Xrefs to create.

Can apply to all ensembl species but does not. Core DB must contain relevant
data, so many non-vertebrates are excluded.

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::RFAMParser->new(
    source_id  => 1,
    species_id => 9606,
    files => [$path.'Rfam.seed.gz']
    xref_dba   => $xref_dba,  # xref db adaptor
    core_dba   => $dba   # core db adaptor
  );

  $parser->run();
=cut


package Bio::EnsEMBL::Xref::Parser::RFAMParser;

use strict;
use warnings;
use Carp;

use parent 'Bio::EnsEMBL::Xref::Parser';

=head2 run

Description: Reimplements Bio::EnsEMBL::Xref::Parser::run . Requires source_id,
             species_id, a list containing the Rfam seed file, and database
             adaptors for core DB and xref DB
Caller     : Xref pipeline ParseSource

=cut

sub run {
  my ( $self ) = @_;
  my $source_id    = $self->{source_id};
  my $species_id   = $self->{species_id};
  my $files        = $self->{files};
  my $xref_dba     = $self->{xref_dba};
  my $core_dba     = $self->{dba};

  my $file = shift @$files;

  my $analysis_adaptor = $core_dba->get_AnalysisAdaptor();

  my $rfam_analysis = $analysis_adaptor->fetch_by_logic_name('rfamblast');

  # RFAM IDs are not present for all species. If Genebuild have not run the 
  # analysis, skip this species
  return 0 if ! defined $rfam_analysis;

  my @complete_meta = @{ $self->parse_stockholm_format($file) };

  # Iterate through all RFAM records
  my %supporting_features = %{ $self->get_rfam_supporting_features($core_dba) };

  # Now insert any matches as xrefs with direct links to Transcript IDs
  # referenced by supporting features
  foreach my $rfam (@complete_meta) {
    if ( exists $supporting_features{ $rfam->{AC} } ) {
      my $xref_id = $xref_dba->add_xref({
        acc => $rfam->{AC},
        version => 1,
        label => $rfam->{ID},
        desc => $rfam->{DE},
        source_id => $source_id,
        species_id => $species_id,
        info_type => 'DIRECT'
      });

      foreach my $stable_id ( @{ $supporting_features{ $rfam->{AC}} }) {
        $xref_dba->add_direct_xref($xref_id, $stable_id, 'Transcript');
      }
    }
  }
  return 0;
}

=head2 parse_stockholm_format

Arg 1      :  File path
Description:  Extracts all the metadata from an RFAM file
              We completely ignore all the alignments
Returntype :  Arrayref of hashrefs containing RFAM metadata

=cut

sub parse_stockholm_format {
  my ($self,$file_path) = @_;

  my %rfam_ids = ();
  # It's not a big file, we can do this in a single hit
  my $fh = $self->{xref_dba}->get_filehandle($file_path);
  local $/ = '//';

  my @complete_meta;
  while (my $record = <$fh>) {

    my %meta = $record =~/
      ^
      \#=GF\s+  # All meta keys start with #=GF
      (\w\w)    # Two letter code for meaning of key
      \s+
      (.+)      # Meta fields sometimes contain spaces
      $
    /mgx;

    unless ( exists $meta{AC} && exists $meta{ID} && exists $meta{DE}) {
      warn 'Skipping an incompletely annotated record in '.$file_path;
      next;
    }

    push @complete_meta, \%meta;
  }
 
  return \@complete_meta;
}

=head2 get_rfam_supporting_features

Arg 1      : $dba - a core DBAdaptor
Description: Access a core DB and retrieve DnaAlignFeatures with the correct
             ncRNA logic name, and match them up with transcript stable IDs
             This could perhaps be achieved through the core API, but it will
             be cumbersome unless we extend the TranscriptSupportingFeatureAdaptor
             for this specific purpose
Returntype : Hashref of listrefs - RFAM Ids and their Ensembl transcript ID pairings

=cut

sub get_rfam_supporting_features {
  my ($self,$dba) = @_;

  my $rfam_sql = q(
    SELECT DISTINCT t.stable_id, hit_name 
    FROM analysis a
    JOIN transcript t
      ON (
            a.analysis_id = t.analysis_id
        AND a.logic_name = 'ncRNA'
        AND t.biotype != 'miRNA'
      )
    JOIN exon_transcript et
      ON (
        t.transcript_id = et.transcript_id
      )
    JOIN supporting_feature sf
      ON (
            et.exon_id = sf.exon_id
        AND sf.feature_type = 'dna_align_feature' 
      )
    JOIN dna_align_feature df
      ON (
        sf.feature_id = df.dna_align_feature_id
      )
    ORDER BY hit_name
  );

  my $iterator = $dba->dbc->sql_helper()->execute(-SQL => $rfam_sql, -ITERATOR => 1);
  my %rfam_mapping;
  while ($iterator->has_next()) {
    my $row = $iterator->next;
    # RFAM sequences align to 1+ transcripts
    if ($row->[1] =~ /^(RF\d+)/) {
      my $accession = $1;
      push @{ $rfam_mapping{$accession} }, $row->[0];
    }
  }
  return \%rfam_mapping;
}

1;
