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

Bio::EnsEMBL::Xref::Parser::RefSeqCoordinateParser

=head1 DESCRIPTION

A parser class for the RefSeq_import source.
Requires an otherfeatures database adaptor and the species must have
RefSeq_import data or will be skipped.

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::RefSeqCoordinateParser->new(
    source_id  => 94,
    species_id => 1,
    files      => [],
    dba        => $otherfeatures_dba,
    xref_dba   => $xref_dba
  );

  $parser->run();

=cut


package Bio::EnsEMBL::Xref::Parser::RefSeqCoordinateParser;

use strict;
use warnings;
use Carp;

use parent qw( Bio::EnsEMBL::Xref::Parser );

# score is the fraction of exons matching out of all exons in the transcript
# Only scores higher than the threshold will be stored for transcripts
my $TRANSCRIPT_SCORE_THRESHOLD = 0.75;

# Only scores higher than the threshold will be stored translatable transcripts
my $TL_TRANSCRIPT_SCORE_THRESHOLD = 0.75;

# If Biotypes do not match, score will be multiplied with the penalty
my $PENALTY = 0.9;



=head2 run
  Description: Runs the RefSeqCoordinateParser
  Return type: N/A
  Caller     : internal
=cut

sub run {
  my ( $self ) = @_;

  my $source_id    = $self->{source_id};
  my $species_id   = $self->{species_id};
  my $species_name = $self->{species};
  my $of_dba       = $self->{dba};
  my $xref_dba     = $self->{xref_dba};
  my $verbose      = $self->{verbose} // 0;

  # initial param validation step
  if ( (!defined $source_id) || (!defined $species_id) ){
    confess "Need to pass source_id and species_id as pairs.";
  }

  # otherfeatures db param validation
  if ( !defined $of_dba ) {
    confess "Need to pass dba as pairs. Missing otherfeatures database adaptor.";
  }

  my $core_dba = $of_dba->dnadb();

  $self->{refseq_sources} = $xref_dba->get_refseq_sources;

  # get RefSeq source ids
  while (my ($source_prefix, $source_name) = each %{$self->{refseq_sources}}) {
    $self->{source_ids}->{$source_name} = $xref_dba->get_source_id_for_source_name( $source_name, 'otherfeatures' )
  }

  if ($verbose) {
    for my $source_name (sort values %{$self->{refseq_sources}}) {
      print "$source_name source ID = $self->{source_ids}->{$source_name}\n";
    }
  }

  # Cache EntrezGene IDs and source ID where available
  my $entrez_ids = $xref_dba->get_valid_codes('EntrezGene', $species_id);
  $self->{source_ids}->{EntrezGene} = $xref_dba->get_source_id_for_source_name('EntrezGene', undef);
  my $entrez_source_id = $xref_dba->get_source_id_for_source_name('EntrezGene', undef);

  my $sa = $core_dba->get_SliceAdaptor();
  my $sa_of = $of_dba->get_SliceAdaptor();
  my $chromosomes_of = $sa_of->fetch_all('toplevel', undef, 1);

  # Fetch analysis object for refseq
  my $aa_of = $of_dba->get_AnalysisAdaptor();

  # Not all species have refseq_import data, exit if not found
  if (!defined $aa_of->fetch_by_logic_name('refseq_import')) {
    print "No data found for RefSeq_import. Skipping\n";
    return;
  }

  # Iterate over chromosomes in otherfeatures database
  foreach my $chromosome_of (@{$chromosomes_of}) {
    my $chr_name = $chromosome_of->seq_region_name();

    my $genes_of = $chromosome_of->get_all_Genes('refseq_import', undef, 1);

    # For each gene in that chromosome in otherfeatures database
    foreach my $gene_of (@{$genes_of}) {
      my $transcripts_of = $gene_of->get_all_Transcripts();

      # Create a range registry for all the exons of the refseq transcript
      TRANSCRIPT_OF:
      foreach my $transcript_of (sort { $a->start <=> $b->start } @{$transcripts_of}) {
        my $id;
        # RefSeq accessions are now stored as xrefs rather than
        # stable ids as it used to be in the past. This means
        # priority is given to the display_id, and fall back to stable_id
        # for backwards compatibility.
        if (defined $transcript_of->display_xref ) {
          $id = $transcript_of->display_xref->display_id;
        } elsif (defined $transcript_of->stable_id) {
          $id = $transcript_of->stable_id;
        }

        # Skip non supported and missing accessions
        unless ( exists $self->{refseq_sources}->{substr($id, 0, 2)} ) {
          next TRANSCRIPT_OF;
        }

        my $transcript_result;
        my $tl_transcript_result;

        my $exons_of = $transcript_of->get_all_Exons();
        my $rr_exons_of = Bio::EnsEMBL::Mapper::RangeRegistry->new();
        my $tl_exons_of = $transcript_of->get_all_translateable_Exons();
        my $rr_tl_exons_of = Bio::EnsEMBL::Mapper::RangeRegistry->new();

        # register $exons_of on $rr_exons_of
        $self->compute_exons({
          exons              => $exons_of,
          check_and_register => $rr_exons_of
        });

        # register $tl_exons_of on $rr_tl_exons_of
        $self->compute_exons({
          exons              => $tl_exons_of,
          check_and_register => $rr_tl_exons_of
        });

        # Fetch slice in core database which overlaps refseq transcript
        my $chromosome = $sa->fetch_by_region('toplevel', $chr_name, $transcript_of->seq_region_start, $transcript_of->seq_region_end);
        my $transcripts = $chromosome->get_all_Transcripts(1);

        # Create a range registry for all the exons of the ensembl transcript
        TRANSCRIPT:
        foreach my $transcript(@{$transcripts}) {
          # make sure it's the same strand
          if ($transcript->strand != $transcript_of->strand) {
            next TRANSCRIPT;
          }
          my $exons = $transcript->get_all_Exons();
          my $rr_exons = Bio::EnsEMBL::Mapper::RangeRegistry->new();
          my $tl_exons = $transcript->get_all_translateable_Exons();
          my $rr_tl_exons = Bio::EnsEMBL::Mapper::RangeRegistry->new();

          # register $exons on $rr_exons, overlap with $rr_exons_of
          my $exon_match = $self->compute_exons({
            exons              => $exons,
            check_and_register => $rr_exons,
            overlap            => $rr_exons_of
          });

          # register $tl_exons on $rr_tl_exons, overlap with $rr_tl_exons_of
          my $tl_exon_match = $self->compute_exons({
            exons              => $tl_exons,
            check_and_register => $rr_tl_exons,
            overlap            => $rr_tl_exons_of
          });

          # $exons_of overlap with $rr_exons
          my $exon_match_of = $self->compute_exons({
            exons   => $exons_of,
            overlap => $rr_exons
          });

          # $tl_exons_of overlap with $rr_tl_exons
          my $tl_exon_match_of = $self->compute_exons({
            exons   => $tl_exons_of,
            overlap => $rr_tl_exons
          });

          # Comparing exon matching with number of exons to give a score
          # score is the fraction of exons matching out of all exons in the transcript
          my $score = ( ($exon_match_of + $exon_match)) / (scalar(@{$exons_of}) + scalar(@{$exons}) );
          my $tl_score = 0;
          if (scalar(@{$tl_exons_of}) > 0) {
            $tl_score = ( ($tl_exon_match_of + $tl_exon_match)) / (scalar(@{$tl_exons_of}) + scalar(@{$tl_exons}) );
          }
          if ($transcript->biotype eq $transcript_of->biotype) {
            $transcript_result->{$transcript->stable_id} = $score;
            $tl_transcript_result->{$transcript->stable_id} = $tl_score;
          } else {
            $transcript_result->{$transcript->stable_id} = $score * $PENALTY;
            $tl_transcript_result->{$transcript->stable_id} = $tl_score * $PENALTY;
          }
        }

        my ($best_id, $best_score, $best_tl_score) = $self->compute_best_scores($transcript_result, $tl_transcript_result);

        # If a best match was defined for the refseq transcript, store it as direct xref for ensembl transcript
        if ($best_id) {
          my ($acc, $version) = split(/\./x, $id);

          my $source = $self->source_id_from_acc($acc);

          next TRANSCRIPT unless defined $source;

          my $xref_id = $xref_dba->add_xref({
            acc        => $acc,
            version    => $version,
            label      => $id,
            source_id  => $source,
            species_id => $species_id,
            info_type  => 'DIRECT'
          });
          $xref_dba->add_direct_xref($xref_id, $best_id, 'Transcript', undef);

          my $entrez_id;
          # RefSeq accessions are now stored as xrefs rather than
          # stable ids as it used to be in the past. This means
          # priority is given to the display_id, and fall back to stable_id
          # for backwards compatibility.
          if (defined $gene_of->display_xref) {
            $entrez_id = $gene_of->display_xref->display_id;
          } elsif (defined $gene_of->stable_id) {
            $entrez_id = $gene_of->stable_id;
          }

          my $tl_of = $transcript_of->translation();
          my $ta = $core_dba->get_TranscriptAdaptor();
          my $t = $ta->fetch_by_stable_id($best_id);
          my $tl = $t->translation();

          # Add link between Ensembl gene and EntrezGene
          if (defined $entrez_ids->{$entrez_id} ) {
            foreach my $dependent_xref_id (@{$entrez_ids->{$entrez_id}}) {
              $xref_dba->add_dependent_xref({
                master_xref_id => $xref_id,
                acc            => $dependent_xref_id,
                source_id      => $self->source_id_from_name('EntrezGene'),
                species_id     => $species_id
              });
            }
          }

          # Also store refseq protein as direct xref for ensembl translation, if translation exists
          if (defined $tl && defined $tl_of) {
            if ($tl_of->seq eq $tl->seq) {
              my $tl_id = $tl_of->stable_id();
              my @xrefs = grep {$_->{dbname} eq 'GenBank'} @{$tl_of->get_all_DBEntries};
              if(scalar @xrefs == 1) {
                $tl_id = $xrefs[0]->primary_id();
              }
              my ($tl_acc, $tl_version) = split(/\./xms, $tl_id);

              my $tl_source = $self->source_id_from_acc($tl_acc);

              next TRANSCRIPT unless defined $tl_source;

              my $tl_xref_id = $xref_dba->add_xref({
                acc        => $tl_acc,
                version    => $tl_version,
                label      => $tl_id,
                source_id  => $tl_source,
                species_id => $species_id,
                info_type  => 'DIRECT'
              });
              $xref_dba->add_direct_xref($tl_xref_id, $tl->stable_id(), 'Translation', undef);
            }
          }
        }
      }
    }
  }
  return 0;
}



=head2 compute_exons
  Arg [1]    : Hash ref: exons, check_and_register and overlap as pairs
  Description: exons is the array ref of the exons to process.
               check_and_register may contain a range registry to check_and_register the exons there.
               overlap may contain a range registry to calculate the overlap of the exons there.
               Returns the exon_match, which is always 0 if no overlap requested.
  Return type: Scalar (number)
  Caller     : internal
=cut

sub compute_exons {
  my ($self, $params) = @_;

  my $exon_match = 0;

  foreach my $exon (@{$params->{exons}}) {
    if (defined $params->{check_and_register}) {
      $params->{check_and_register}->check_and_register( 'exon', $exon->seq_region_start, $exon->seq_region_end );
    }
    if (defined $params->{overlap}) {
      my $overlap = $params->{overlap}->overlap_size('exon', $exon->seq_region_start, $exon->seq_region_end);
      $exon_match += $overlap / ($exon->seq_region_end - $exon->seq_region_start + 1);
    }
  }

  return $exon_match;
}

=head2 compute_best_scores
  Arg [1]    : hash ref: transcript_result
  Arg [2]    : hash ref: tl_transcript_result
  Description: Provided the hashrefs of the transcripts and translatable transcripts,
               computes the best ID and the best transcript and translatable transcript score.
  Return type: Array (best_id, best_score, best_tl_score)
  Caller     : internal
=cut

sub compute_best_scores {
  my ($self, $transcript_result, $tl_transcript_result) = @_;

  my $best_score = 0;
  my $best_tl_score = 0;
  my $best_id;

  # Comparing the scores based on coding exon overlap
  # If there is a stale mate, chose best exon overlap score
  foreach my $tid (sort { $transcript_result->{$b} <=> $transcript_result->{$a} } keys(%{$transcript_result})) {
    my $score = $transcript_result->{$tid};
    my $tl_score = $tl_transcript_result->{$tid};
    if ($score > $TRANSCRIPT_SCORE_THRESHOLD || $tl_score > $TL_TRANSCRIPT_SCORE_THRESHOLD) {
      if ($tl_score > $best_tl_score) {
        $best_id = $tid;
        $best_score = $score;
        $best_tl_score = $tl_score;
      } elsif ($tl_score == $best_tl_score) {
        if ($score > $best_score) {
          $best_id = $tid;
          $best_score = $score;
        }
      } elsif ($score >= $best_score) {
        $best_id = $tid;
        $best_score = $score;
      }
    }
  }

  return ($best_id, $best_score, $best_tl_score);
}

=head2 source_id_from_name
  Arg [1]    : Scalar (string name)
  Description: Returns the source ID for provided source name.
               Prints warning message on failure.
  Return type: Scalar (number source_id)
  Caller     : internal
=cut

sub source_id_from_name {
  my ($self, $name) = @_;

  my $source_id;

  if ( exists $self->{source_ids}->{$name} ) {
    $source_id = $self->{source_ids}->{$name};
  } elsif ( $self->{verbose} ) {
    confess "Can't get source ID for name '$name'\n";
  }

  return $source_id;
}

=head2 source_id_from_acc
  Arg [1]    : Scalar (string acc)
  Description: Returns the source ID for provided RefSeq accession.
               Requires $self->{source_ids} to have been populated.
               Prints warning message on failure.
  Return type: Scalar (number source_id)
  Caller     : internal
=cut

sub source_id_from_acc {
  my ($self, $acc) = @_;

  my $source_id;
  my $prefix = substr($acc, 0, 2);

  if ( exists $self->{refseq_sources}->{$prefix} ) {
    $source_id = $self->source_id_from_name( $self->{refseq_sources}->{$prefix} );
  } elsif ( $self->{verbose} ) {
    confess "Can't get source ID for accession '$acc'\n";
  }

  return $source_id;
}

1;
