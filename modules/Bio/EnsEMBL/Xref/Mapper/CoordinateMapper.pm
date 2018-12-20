
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

Bio::EnsEMBL::Xref::Mapper::CoordinateMapper

=cut

=head1 DESCRIPTION

Mapper class with a set of subroutines used for creating Xrefs based on
coordinate overlaps.

=cut

package Bio::EnsEMBL::Xref::Mapper::CoordinateMapper;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;

use DBI qw( :sql_types );

use Carp;
use IO::File;
use File::Spec::Functions;

use parent qw( Bio::EnsEMBL::Xref::Mapper );

my $coding_weight = 2;
my $ens_weight    = 3;

my $transcript_score_threshold = 0.75; 
# 75% identity is required at the minimum for transcript coordinates to match

=head2 run_coordinatemapping
  Arg [1]    : boolean flag set to upload result to table or not
  Arg [2]    : species id
  Arg [3]    : output dir
  Description: Main routine that calculates the overlap score between the ensembl transcripts and 
               the other data source transcripts (eg: ucsc)

              For each Ensembl transcript:
                1. Register all Ensembl exons in a RangeRegistry.
                2. Find all transcripts in the external database that are within the range of
                   this Ensembl transcript
              
              For each of those external transcripts:
                3. Calculate the overlap of the exons of the external transcript with the Ensembl
                   exons using the overlap_size() method in the RangeRegistry.
                4. Register the external exons in their own RangeRegistry.
                5. Calculate the overlap of the Ensembl exons with the external exons as in step 3
                6. Calculate the match score.
                7. Decide whether or not to keep the match.

  Return type: None
  Caller     : Bio::EnsEMBL::Production::Pipeline::Xrefs::CoordinateMapping

=cut

sub run_coordinatemapping {
  my ( $self, $do_upload, $species_id, $output_dir ) = @_;

  unless ( -d $output_dir ) {
    confess "Cannot access output directory at $output_dir";
  }

  my $xref_db = $self->xref();
  my $core_db = $self->core();

  my $species = $core_db->species();

  # We only do coordinate mapping for mouse and human for now.
  unless ( $species eq 'mus_musculus' || $species eq 'homo_sapiens' ) {
    return;
  }

  my $xref_filename        = catfile( $output_dir, 'xref_coord.txt' );
  my $object_xref_filename = catfile( $output_dir, 'object_xref_coord.txt' );
  my $unmapped_reason_filename = catfile( $output_dir, 'unmapped_reason_coord.txt' );
  my $unmapped_object_filename = catfile( $output_dir, 'unmapped_object_coord.txt' );

  my $xref_dbh = $xref_db->dbc()->db_handle();
  my $core_dbh = $core_db->dbc()->db_handle();

  ######################################################################
  # Figure out the last used 'xref_id', 'object_xref_id',              #
  # 'unmapped_object_id', and 'unmapped_reason_id' from the Core       #
  # database so we can cheat on bulk insert                            #
  ######################################################################

  my $xref_id = $core_dbh->selectall_arrayref('SELECT MAX(xref_id) FROM xref')->[0][0];
  my $object_xref_id = $core_dbh->selectall_arrayref('SELECT MAX(object_xref_id) FROM object_xref')->[0][0];

  my $unmapped_object_id = $core_dbh->selectall_arrayref(
    'SELECT MAX(unmapped_object_id) FROM unmapped_object')->[0][0];
  my $unmapped_reason_id = $core_dbh->selectall_arrayref(
    'SELECT MAX(unmapped_reason_id) FROM unmapped_reason')->[0][0];

  $self->log_progress( "Last used xref_id            is %d\n", $xref_id );
  $self->log_progress( "Last used object_xref_id     is %d\n", $object_xref_id );
  $self->log_progress( "Last used unmapped_object_id is %d\n", $unmapped_object_id );
  $self->log_progress( "Last used unmapped_reason_id is %d\n", $unmapped_reason_id );

  ######################################################################
  # Get an 'analysis_id', or discover that we need to add our analysis #
  # to the 'analyis' table later.                                      #
  ######################################################################

  my $analysis_params =
    sprintf( "weights(coding,ensembl)="
      . "%.2f,%.2f;"
      . "transcript_score_threshold=" . "%.2f",
    $coding_weight, $ens_weight, $transcript_score_threshold );

  my $analysis_adaptor = $self->core->get_AnalysisAdaptor;
  my $analysis = $analysis_adaptor->fetch_by_logic_name('xrefcoordinatemapping');

  my $analysis_id;
  if (! defined $analysis ) {
    $analysis_id = $analysis_adaptor->store(
      Bio::EnsEMBL::Analysis->new(
        -logic_name => 'xrefcoordinatemapping',
        -program => 'Bio::EnsEMBL::Production::Pipeline::Xrefs::CoordinateMapping',
        -parameters => $analysis_params,
        -program_file => 'CoordinateMapper.pm'
      )
    );
  } elsif ( $analysis->parameters ne $analysis_params && $do_upload ) {

    $analysis->parameters($analysis_params);
    $analysis_id = $analysis_adaptor->update($analysis);

  }

  if ( defined($analysis_id) ) {
    $self->log_progress( "Analysis ID                  is %d\n", $analysis_id );
  }

  ######################################################################
  # Read and store available Xrefs from the Xref database.             #
  ######################################################################

  my %unmapped;
  my %mapped;

  my $xref_sql = qq(
    SELECT  coord_xref_id, source_id, accession
    FROM    coordinate_xref
    WHERE   species_id = ?
  );

  my $xref_sth = $xref_dbh->prepare($xref_sql);
  $xref_sth->bind_param( 1, $species_id, SQL_INTEGER );

  $xref_sth->execute();

  my $external_db_id = 11000;    # FIXME HARDCODED (11000 is 'UCSC')

  while ( my $xref = $xref_sth->fetchrow_hashref() ) {

    $unmapped{ $xref->{'coord_xref_id'} } = {
      'external_db_id' => $external_db_id,
      'accession'      => $xref->{'accession'},
      'reason'         => 'No overlap',
      'reason_full'    => 'No coordinate overlap with any Ensembl transcript'
    };
  }
  $xref_sth->finish();

  
  ######################################################################
  # Do coordinate matching.                                            #
  ######################################################################

  my $slice_adaptor = $core_db->get_SliceAdaptor();
  my @chromosomes   = @{ $slice_adaptor->fetch_all('Chromosome') };

  my $sql = qq(
    SELECT    coord_xref_id, accession,
              txStart, txEnd,
              cdsStart, cdsEnd,
              exonStarts, exonEnds
    FROM      coordinate_xref
    WHERE     species_id = ?
    AND       chromosome = ? AND strand   = ?
    AND       ((txStart BETWEEN ? AND ?)        -- txStart in region
    OR         (txEnd   BETWEEN ? AND ?)        -- txEnd in region
    OR         (txStart <= ? AND txEnd >= ?))   -- region is fully contained
    ORDER BY  accession
  );

  foreach my $chromosome (@chromosomes) {
    my $chr_name = $chromosome->seq_region_name();

    $self->log_progress( "Processing chromsome '%s'\n", $chr_name );

    my @genes = @{ $chromosome->get_all_Genes( undef, undef, 1 ) };

    $self->log_progress( "There are %4d genes on chromosome '%s'\n",
      scalar(@genes), $chr_name );

    while ( my $gene = shift(@genes) ) {

      my @transcripts = @{ $gene->get_all_Transcripts() };

      my %gene_result;

      foreach my $transcript ( sort { $a->start() <=> $b->start() } @transcripts )
      {
        my @exons = @{ $transcript->get_all_Exons() };

        my %transcript_result;

        # '$rr1' is the RangeRegistry holding Ensembl exons for one
        # transcript at a time.
        my $rr1 = Bio::EnsEMBL::Mapper::RangeRegistry->new();

        my $coding_transcript;
        if ( defined( $transcript->translation() ) ) {
          $coding_transcript = 1;
        }
        else {
          $coding_transcript = 0;
        }

        foreach my $exon (@exons) {

          #-------------------------------------------------------------
          # Register each exon in the RangeRegistry.  Register both the
          # total length of the exon and the coding range of the exon.
          #-------------------------------------------------------------

          $rr1->check_and_register( 'exon', $exon->start(), $exon->end() );

          if ( $coding_transcript
            && defined( $exon->coding_region_start($transcript) )
            && defined( $exon->coding_region_end($transcript) ) )
          {
            $rr1->check_and_register(
              'coding',
              $exon->coding_region_start($transcript),
              $exon->coding_region_end($transcript)
            );
          }
        }

        #---------------------------------------------------------------
        # Get hold of all transcripts from the external database that
        # overlaps with this Ensembl transcript.
        #---------------------------------------------------------------

        my $sth = $xref_dbh->prepare_cached($sql);
        $sth->bind_param( 1, $species_id,          SQL_INTEGER );
        $sth->bind_param( 2, $chr_name,            SQL_VARCHAR );
        $sth->bind_param( 3, $gene->strand(),      SQL_INTEGER );
        $sth->bind_param( 4, $transcript->start(), SQL_INTEGER );
        $sth->bind_param( 5, $transcript->end(),   SQL_INTEGER );
        $sth->bind_param( 6, $transcript->start(), SQL_INTEGER );
        $sth->bind_param( 7, $transcript->end(),   SQL_INTEGER );
        $sth->bind_param( 8, $transcript->start(), SQL_INTEGER );
        $sth->bind_param( 9, $transcript->end(),   SQL_INTEGER );

        $sth->execute();

        my ( $coord_xref_id, $accession, $txStart, $txEnd, $cdsStart,
          $cdsEnd, $exonStarts, $exonEnds );

        $sth->bind_columns(
          \(
            $coord_xref_id, $accession, $txStart,    $txEnd,
            $cdsStart,      $cdsEnd,    $exonStarts, $exonEnds
          )
        );

        while ( $sth->fetch() ) {
          my @exonStarts = split( /,\s*/, $exonStarts );
          my @exonEnds   = split( /,\s*/, $exonEnds );
          my $exonCount  = scalar(@exonStarts);

          # '$rr2' is the RangeRegistry holding exons from the external
          # transcript, for one transcript at a time.
          my $rr2 = Bio::EnsEMBL::Mapper::RangeRegistry->new();

          my $exon_match   = 0;
          my $coding_match = 0;

          my $coding_count = 0;

          for ( my $i = 0 ; $i < $exonCount ; ++$i ) {

            #-----------------------------------------------------------
            # Register the exons from the external database in the same
            # way as with the Ensembl exons, and calculate the overlap
            # of the external exons with the previously registered
            # Ensembl exons.
            #-----------------------------------------------------------

            my $overlap =
              $rr1->overlap_size( 'exon', $exonStarts[$i], $exonEnds[$i] );

            $exon_match += $overlap / ( $exonEnds[$i] - $exonStarts[$i] + 1 );

            $rr2->check_and_register( 'exon', $exonStarts[$i], $exonEnds[$i] );

            if ( !defined($cdsStart) || !defined($cdsEnd) ) {

              # Non-coding transcript. Do nothing!
            }
            else {
              my $codingStart =
                ( $exonStarts[$i] > $cdsStart ? $exonStarts[$i] : $cdsStart );
              my $codingEnd =
                ( $exonEnds[$i] < $cdsEnd ? $exonEnds[$i] : $cdsEnd );

              if ( $codingStart < $codingEnd ) {
                my $coding_overlap =
                  $rr1->overlap_size( 'coding', $codingStart, $codingEnd );

                $coding_match +=
                  $coding_overlap / ( $codingEnd - $codingStart + 1 );

                $rr2->check_and_register( 'coding', $codingStart, $codingEnd );

                ++$coding_count;
              }
            }
          } ## end for exonCount

          my $rexon_match   = 0;
          my $rcoding_match = 0;

          my $rcoding_count = 0;

          foreach my $exon (@exons) {

            #-----------------------------------------------------------
            # Calculate the overlap of the Ensembl exons with the
            # external exons.
            #-----------------------------------------------------------

            my $overlap =
              $rr2->overlap_size( 'exon', $exon->start(), $exon->end() );

            $rexon_match += $overlap / ( $exon->end() - $exon->start() + 1 );

            if ( $coding_transcript
              && defined( $exon->coding_region_start($transcript) )
              && defined( $exon->coding_region_end($transcript) ) )
            {
              my $coding_overlap = $rr2->overlap_size(
                'coding',
                $exon->coding_region_start($transcript),
                $exon->coding_region_end($transcript)
              );

              $rcoding_match +=
                $coding_overlap /
                ( $exon->coding_region_end($transcript) -
                  $exon->coding_region_start($transcript) +
                  1 );

              ++$rcoding_count;
            }
          } ## end foreach my $exon (@exons)

          #-------------------------------------------------------------
          # Calculate the match score.
          #-------------------------------------------------------------

          my $score =
            ( ( $exon_match + $ens_weight * $rexon_match ) +
              $coding_weight *
              ( $coding_match + $ens_weight * $rcoding_match ) ) /
            ( ( $exonCount + $ens_weight * scalar(@exons) ) +
              $coding_weight *
              ( $coding_count + $ens_weight * $rcoding_count ) );

          if ( !defined( $transcript_result{$coord_xref_id} )
            || $transcript_result{$coord_xref_id} < $score )
          {
            $transcript_result{$coord_xref_id} = $score;
          }

        } ## end while ( $sth->fetch() )
        $sth->finish();

        #---------------------------------------------------------------
        # Apply transcript threshold and pick the best match(es) for
        # this transcript.
        #---------------------------------------------------------------

        my $best_score;
        foreach my $coord_xref_id (
          sort( { $transcript_result{$b} <=> $transcript_result{$a} }
            keys(%transcript_result) )
          )
        {
          my $score = $transcript_result{$coord_xref_id};

          if ( $score > $transcript_score_threshold ) {
            $best_score ||= $score;

            if ( sprintf( "%.3f", $score ) eq sprintf( "%.3f", $best_score ) ) {
              if ( exists $unmapped{$coord_xref_id} ) {
                $mapped{$coord_xref_id} = $unmapped{$coord_xref_id};
                delete $unmapped{$coord_xref_id};
                $mapped{$coord_xref_id}{'reason'}      = undef;
                $mapped{$coord_xref_id}{'reason_full'} = undef;
              }

              push @{ $mapped{$coord_xref_id}{'mapped_to'} },
                {
                  'ensembl_id'          => $transcript->dbID(),
                  'ensembl_object_type' => 'Transcript'
                };

              # This is now a candidate Xref for the gene.
              if ( !defined $gene_result{$coord_xref_id} 
                || $gene_result{$coord_xref_id} < $score )
              {
                $gene_result{$coord_xref_id} = $score;
              }

            }
            elsif ( exists( $unmapped{$coord_xref_id} ) ) {
              $unmapped{$coord_xref_id}{'reason'} = 'Was not best match';
              $unmapped{$coord_xref_id}{'reason_full'} =
                sprintf( 'Did not top best transcript match score (%.2f)',
                  $best_score );
              if ( !defined $unmapped{$coord_xref_id}{'score'}
                || $score > $unmapped{$coord_xref_id}{'score'} )
              {
                $unmapped{$coord_xref_id}{'score'} = $score;
                $unmapped{$coord_xref_id}{'ensembl_id'} =
                  $transcript->dbID();
              }
            }

          }
          elsif ( exists $unmapped{$coord_xref_id}
            && $unmapped{$coord_xref_id}{'reason'} ne 'Was not best match' )
          {
            $unmapped{$coord_xref_id}{'reason'}      = 'Did not meet threshold';
            $unmapped{$coord_xref_id}{'reason_full'} = sprintf(
              'Match score for transcript lower than threshold (%.2f)',
              $transcript_score_threshold
            );
            if ( !defined $unmapped{$coord_xref_id}{'score'}
              || $score > $unmapped{$coord_xref_id}{'score'} )
            {
              $unmapped{$coord_xref_id}{'score'} = $score;
              $unmapped{$coord_xref_id}{'ensembl_id'} = 
                $transcript->dbID();
            }
          }
        } ## end foreach my $coord_xref_id (...

      } ## end foreach my $transcript ( sort...

    } ## end while ( my $gene = shift(...
  } ## end foreach my $chromosome (@chromosomes)

  # Make all dumps.  Order is important.
  $self->dump_xref( $xref_filename, $xref_id, \%mapped, \%unmapped );
  $self->dump_object_xref( $object_xref_filename, $object_xref_id, \%mapped );
  $self->dump_unmapped_reason( $unmapped_reason_filename, $unmapped_reason_id, \%unmapped );
  $self->dump_unmapped_object( $unmapped_object_filename, $unmapped_object_id, $analysis_id, \%unmapped );

  if ($do_upload) {

    # Order is important!
    $self->upload_data( 'unmapped_reason', $unmapped_reason_filename, $external_db_id );
    $self->upload_data( 'unmapped_object', $unmapped_object_filename, $external_db_id );
    $self->upload_data( 'object_xref', $object_xref_filename, $external_db_id );
    $self->upload_data( 'xref', $xref_filename, $external_db_id );
  }

  $self->biomart_fix( "UCSC", "Translation", "Gene" );
  $self->biomart_fix( "UCSC", "Transcript",  "Gene" );
  return;
} ## end sub run_coordinatemapping

=head2 dump_xref
  Arg [1]    : filename to write
  Arg [2]    : xref id
  Arg [3]    : container of mapped coordinate ids from coordinate xref table
  Arg [4]    : container of unmapped coordinate ids from coordinate xref table
  Description: Dump mapped xrefs to a txt file
  Return type: None
  Caller     : internal

=cut

sub dump_xref {
  my ( $self, $filename, $xref_id, $mapped, $unmapped ) = @_;

  my $fh = IO::File->new( '>' . $filename )
    or confess( sprintf( "Can not open '%s' for writing", $filename ) );

  $self->log_progress( "Dumping for 'xref' to '%s'\n", $filename );

  foreach my $xref ( values %{$unmapped}, values %{$mapped} ) {

    # Assign 'xref_id' to this Xref.
    $xref->{'xref_id'} = ++$xref_id;

    my $accession = $xref->{'accession'};

    my ($version) = ( $accession =~ /\.(\d+)$/ );
    $version ||= 0;

    $fh->printf(
      "%d\t%d\t%s\t%s\t%d\t%s\t%s\t%s\n",
      $xref->{'xref_id'},
      $xref->{'external_db_id'},
      $accession,
      $accession,
      $version,
      '\N',
      'COORDINATE_OVERLAP',
      ''    # FIXME (possibly)
    );
  }
  $fh->close();

  $self->log_progress("Dumping for 'xref' done\n");
  return;
} ## end sub dump_xref

=head2 dump_object_xref
  Arg [1]    : filename to write
  Arg [2]    : object xref id
  Arg [3]    : container of mapped object xref ids
  Description: Dump mapped object xrefs to a txt file
  Return type: None
  Caller     : internal

=cut

sub dump_object_xref {
  my ( $self, $filename, $object_xref_id, $mapped ) = @_;

  my $fh = IO::File->new( '>' . $filename )
    or confess( sprintf( "Can not open '%s' for writing", $filename ) );

  $self->log_progress( "Dumping for 'object_xref' to '%s'\n", $filename );

  foreach my $xref ( values %{$mapped} ) {
    foreach my $object_xref ( @{ $xref->{'mapped_to'} } ) {

      # Assign 'object_xref_id' to this Object Xref.
      $object_xref->{'object_xref_id'} = ++$object_xref_id;

      $fh->printf(
        "%d\t%d\t%s\t%d\t%s\t%s\n",   $object_xref->{'object_xref_id'},
        $object_xref->{'ensembl_id'}, $object_xref->{'ensembl_object_type'},
        $xref->{'xref_id'},           '\N',
        '0'
      );
    }
  }
  $fh->close();

  $self->log_progress("Dumping for 'object_xref' done\n");
  return;
} ## end sub dump_objexref

=head2 dump_unmapped_reason
  Arg [1]    : filename to write
  Arg [2]    : unmapped id from coordinate xref table
  Arg [3]    : container of unmapped ids
  Description: Dump unmapped reason to a txt file (eg: "Did not meet threshold")
  Return type: None
  Caller     : internal

=cut

sub dump_unmapped_reason {
  my ( $self, $filename, $unmapped_reason_id, $unmapped ) = @_;

  # Create a list of the unique reasons.
  my %reasons;

  foreach my $xref ( values %{$unmapped} ) {
    if ( !exists $reasons{ $xref->{'reason_full'} } ) {
      $reasons{ $xref->{'reason_full'} } = {
        'summary' => $xref->{'reason'},
        'full'    => $xref->{'reason_full'}
      };
    }
  }

  my $fh = IO::File->new( '>' . $filename )
    or confess( sprintf( "Can not open '%s' for writing", $filename ) );

  $self->log_progress( "Dumping for 'unmapped_reason' to '%s'\n", $filename );

  my $core_dbh = $self->core->dbc->db_handle;

  my $sth =
    $core_dbh->prepare( 'SELECT unmapped_reason_id '
      . 'FROM unmapped_reason '
      . 'WHERE full_description = ?' );

  foreach
    my $reason ( sort { $a->{'full'} cmp $b->{'full'} } values(%reasons) )
  {
    # Figure out 'unmapped_reason_id' from the core database.
    $sth->bind_param( 1, $reason->{'full'}, SQL_VARCHAR );

    $sth->execute();

    my $id;
    $sth->bind_col( 1, \$id );
    $sth->fetch();

    if ( defined $id ) {
      $reason->{'unmapped_reason_id'} = $id;
    }
    else {
      $reason->{'unmapped_reason_id'} = ++$unmapped_reason_id;
    }

    $sth->finish();

    $fh->printf( "%d\t%s\t%s\n", $reason->{'unmapped_reason_id'},
      $reason->{'summary'}, $reason->{'full'} );

  }
  $fh->close();

  $self->log_progress("Dumping for 'unmapped_reason' done\n");

  # Assign reasons to the unmapped Xrefs from %reasons.
  foreach my $xref ( values %{$unmapped} ) {
    $xref->{'reason'}      = $reasons{ $xref->{'reason_full'} };
    $xref->{'reason_full'} = undef;
  }
  return;
} ## end sub dump_unmapped_reason

=head2 dump_unmapped_object
  Arg [1]    : filename to write
  Arg [2]    : unmapped object id
  Arg [3]    : analysis id
  Arg [4]    : unmapped container
  Description: Dump unmapped object to a txt file
  Return type: None
  Caller     : internal

=cut

sub dump_unmapped_object {
  my ( $self, $filename, $unmapped_object_id, $analysis_id, $unmapped ) = @_;

  my $fh = IO::File->new( '>' . $filename )
    or confess( sprintf( "Can not open '%s' for writing", $filename ) );

  $self->log_progress( "Dumping for 'unmapped_object' to '%s'\n", $filename );

  foreach my $xref ( values %{$unmapped} ) {

    # Assign 'unmapped_object_id' to this Xref.
    $xref->{'unmapped_object_id'} = ++$unmapped_object_id;

    $fh->printf(
      "%d\t%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n",
      $xref->{'unmapped_object_id'},
      'xref',
      $analysis_id || '\N',    # '\N' (NULL) means no analysis exists
                               # and uploading this table will fail.
      $xref->{'external_db_id'},
      $xref->{'accession'},
      $xref->{'reason'}->{'unmapped_reason_id'},
      (
        defined( $xref->{'score'} )
        ? sprintf( "%.3f", $xref->{'score'} )
        : '\N'
      ),
      '\N',
      $xref->{'ensembl_id'} || '\N',
      ( defined( $xref->{'ensembl_id'} ) ? 'Transcript' : '\N' ),
      '\N'
    );
  }
  $fh->close();

  $self->log_progress("Dumping for 'unmapped_object' done\n");
  return;
} ## end sub dump_unmapped_object

=head2 upload_data
  Arg [1]    : table name to upload data
  Arg [2]    : filaname
  Arg [3]    : $external_db_id
  Arg [4]    : database handle
  Description: Upload data from a file to a table
  Return type: None
  Caller     : internal

=cut

sub upload_data {
  my ( $self, $table_name, $filename, $external_db_id ) = @_;

  if ( !-r $filename ) {
    confess( sprintf( "Can not open '%s' for reading", $filename ) );
  }

  my $cleanup_sql = '';
  if ( $table_name eq 'unmapped_reason' ) {
    $cleanup_sql = qq(
    DELETE  ur
    FROM    unmapped_object uo,
            unmapped_reason ur
    WHERE   uo.external_db_id       = ?
    AND     ur.unmapped_reason_id   = uo.unmapped_reason_id);
  }
  elsif ( $table_name eq 'unmapped_object' ) {
    $cleanup_sql = qq(
    DELETE  uo
    FROM    unmapped_object uo
    WHERE   uo.external_db_id = ?);
  }
  elsif ( $table_name eq 'object_xref' ) {
    $cleanup_sql = qq(
    DELETE  ox
    FROM    xref x,
            object_xref ox
    WHERE   x.external_db_id    = ?
    AND     ox.xref_id          = x.xref_id);
  }
  elsif ( $table_name eq 'xref' ) {
    $cleanup_sql = qq(
    DELETE  x
    FROM    xref x
    WHERE   x.external_db_id    = ?);
  }
  else {
    confess( sprintf( "Table '%s' is unknown\n", $table_name ) );
  }

  my $load_sql =
    sprintf( "LOAD DATA LOCAL INFILE ? REPLACE INTO TABLE %s", $table_name );

  $self->log_progress( "Removing old data (external_db_id = '%d') from table '%s'\n",
    $external_db_id, $table_name );

  my $dbh = $self->core->dbc->db_handle();
  my $rows = $dbh->do( $cleanup_sql, undef, $external_db_id )
    or confess( $dbh->strerr() );

  $self->log_progress( "Removed %d rows\n", $rows );

  $self->log_progress( "Uploading for '%s' from '%s'\n", $table_name, $filename );

  $rows = $dbh->do( $load_sql, undef, $filename )
    or confess( $dbh->errstr() );

  $dbh->do("OPTIMIZE TABLE $table_name") or confess( $dbh->errstr() );

  $self->log_progress( "Uploading for '%s' done (%d rows)\n", $table_name, $rows );
  return;
} ## end sub upload_data

1;
