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

Bio::EnsEMBL::Xref::Parser::UCSCParser

=head1 DESCRIPTION

A parser class to parse UCSC data for human and mouse.

-data_uri = ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz
-file_format = TSV
-columns = [
    ensembl_id
    chromosome
    strand
    txStart
    txEnd
    cdsStart
    cdsEnd
    nb_exons
    startExons
    endExons
    uniprot_accession
    ucsc_accession
  ]

Only columns listed in @required_columns are mandatory.

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::UCSCParser->new(
    source_id  => 1,
    species_id => 9606,
    files      => ['UCSC_human/knownGene.txt.gz'],
    xref_dba   => $xref_dba
  );

  $parser->run();
=cut

package Bio::EnsEMBL::Xref::Parser::UCSCParser;

use strict;
use warnings;
use Carp;
use Text::CSV;

use parent qw( Bio::EnsEMBL::Xref::Parser );

=head2 run
  Description: Runs the UCSCParser 
  Return type: N/A
  Caller     : internal
=cut

sub run {
  my ( $self ) = @_;

  my $source_id    = $self->{source_id};
  my $species_id   = $self->{species_id};
  my $files        = $self->{files};
  my $xref_dba     = $self->{xref_dba};
  my $verbose      = $self->{verbose} // 0;


  if ( (!defined $source_id) || (!defined $species_id) || (!defined $files) ) {
    confess "Need to pass source_id, species_id and files as pairs";
  }

  my $file = shift @{$files};

  my $count = 0;

  my $file_io = $xref_dba->get_filehandle($file);

  if ( !defined $file_io ) {
    confess "Can't open UCSC file $file\n";
  }

  my $input_file = Text::CSV->new({
    sep_char       => "\t",
    empty_is_undef => 1
  }) or confess "Cannot use file $file: " . Text::CSV->error_diag();

  # header must contain these columns
  my @columns = qw(
    name
    chromosome
    strand
    txStart
    txEnd
    cdsStart
    cdsEnd
    nbExons
    exonStarts
    exonEnds
    uniprot
    accession
  );

  $input_file->column_names( @columns );
  my ($chromosome, $strand, $cdsStart, $cdsEnd, $txStart, $txEnd, $exonStarts, $exonEnds);

  while ( my $data = $input_file->getline_hr( $file_io ) ) {

    # UCSC uses slightly different chromosome names, at least for
    # human and mouse, so chop off the 'chr' in the beginning.  We do
    # not yet translate the names of the special chromosomes, e.g.
    # "chr6_cox_hap1" (UCSC) into "c6_COX" (Ensembl).
    $chromosome = $data->{'chromosome'};
    $chromosome =~ s/^chr//;

    # They also use '+' and '-' for the strand, instead of -1, 0, or 1.
    $strand = $data->{'strand'};
    if    ( $strand eq '+' ) { $strand = 1 }
    elsif ( $strand eq '-' ) { $strand = -1 }
    else                     { $strand = 0 }

    # ... and non-coding transcripts have cdsStart == cdsEnd.  We would
    # like these to be stored as NULLs.
    $cdsStart = $data->{'cdsStart'};
    $cdsEnd = $data->{'cdsEnd'};
    if ( $cdsStart == $cdsEnd ) {
      undef($cdsStart);
      undef($cdsEnd);
    }

    # ... and they use the same kind of "inbetween" coordinates as e.g.
    # exonerate, so increment all start coordinates by one.
    $txStart = $data->{'txStart'};
    $txStart += 1;
    $exonStarts = $data->{'exonStarts'};
    $exonStarts = join( ',', map( { ++$_ } split( /,/, $exonStarts ) ) );
    if ( defined($cdsStart) ) { $cdsStart += 1 }

    # Cut off the last comma from $exonEnds, if it exists.  This is done
    # for $exonStarts already (above).
    $exonEnds = $data->{'exonEnds'};
    if ( substr( $exonEnds, -1, 1 ) eq ',' ) { chop($exonEnds) }

    $xref_dba->add_coordinate_xref(
      {
        accession  => $data->{'accession'},
        source_id  => $source_id,
        species_id => $species_id,
        chromosome => $chromosome,
        strand     => $strand,
        txStart    => $txStart,
        txEnd      => $data->{'txEnd'},
        cdsStart   => $cdsStart,
        cdsEnd     => $cdsEnd,
        exonStarts => $exonStarts,
        exonEnds   => $exonEnds
      }
    );

    $count++;

  }

  $input_file->eof or confess "Error parsing file $file: " . $input_file->error_diag();
  $file_io->close();

  if ($verbose) {
    print "Loaded a total of $count VGNC xrefs\n";
  }

  return 0; # successful
}


1;
