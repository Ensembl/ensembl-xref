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

Bio::EnsEMBL::Xref::Parser::HGNCParser

=head1 DESCRIPTION

A parser class to parse the HGNC source.
HGNC is the official naming source for Human.

-data_uri = https://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=gd_pub_refseq_ids&col=gd_ccds_ids&col=gd_lsdb_links&status=Approved&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit
-file_format = TSV
-columns = [
    HGNC ID
    Approved symbol
    Approved name
    Previous symbols
    Synonyms
    NCBI Gene ID
    Ensembl gene ID
    RefSeq IDs
    CCDS IDs
    Locus specific databases
  ]

A core database adaptor is required.

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::HGNCParser->new(
    source_id  => 46,
    species_id => 9606,
    files      => ['hgnc_data.tsv'],
    xref_dba   => $xref_dba,
    dba        => $core_dba,
  );

  $parser->run();

=cut

package Bio::EnsEMBL::Xref::Parser::HGNCParser;

use strict;
use warnings;
use Carp;
use Text::CSV;
use utf8;

use parent qw( Bio::EnsEMBL::Xref::Parser );


# HGNC sources to be processed
my @SOURCES = (
  'ccds',
  'entrezgene_manual',
  'refseq_manual',
  'ensembl_manual',
  'desc_only'
);



=head2 run
  Description: Runs the HGNCParser
  Return type: N/A
  Caller     : internal
=cut

sub run {
  my ( $self ) = @_;

  my $source_id  = $self->{source_id};
  my $species_id = $self->{species_id};
  my $files      = $self->{files};
  my $core_dba   = $self->{dba};
  my $xref_dba   = $self->{xref_dba};
  my $verbose    = $self->{verbose} // 0;

  if ((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    confess "Need to pass source_id, species_id and files as pairs";
  }

  if (!defined $core_dba) {
    confess "Need to pass dba as pairs, ensembl core database adaptor required\n";
  }

  my $file = shift @{$files};

  # Prepare lookup lists
  $self->{refseq} = $xref_dba->get_valid_codes( 'refseq', $species_id );
  $self->{entrezgene} = $xref_dba->get_valid_xrefs_for_dependencies( 'EntrezGene', ['refseq_peptide', 'refseq_mRNA'] );

  # Prepare sources
  my $self_source_name = $xref_dba->get_source_name_for_source_id( $source_id );

  # get RefSeq source ids
  foreach my $source_name (@SOURCES) {
    $self->{source_ids}->{$source_name} = $xref_dba->get_source_id_for_source_name( $self_source_name, $source_name );
  }
  $self->{source_ids}->{'lrg'} = $xref_dba->get_source_id_for_source_name( 'LRG_HGNC_notransfer', undef );

  # fetch CCDS data
  my $sql =(<<'CCDS');
  SELECT ta.value, t.stable_id
  FROM transcript t
  INNER JOIN transcript_attrib ta ON t.transcript_id = ta.transcript_id
  INNER JOIN attrib_type a ON ta.attrib_type_id = a.attrib_type_id
  WHERE a.code = 'ccds_transcript';
CCDS

  my $sth = $core_dba->dbc->prepare($sql);
  $sth->execute() or confess( $core_dba->dbc->errstr() );

  while ( my ($ccds_id, $ens_id) = $sth->fetchrow_array() ) {
    # Remove version
    $ccds_id =~ s/\.\d+//x;
    $self->{ccds_to_ens}->{$ccds_id} = $ens_id;
  }
  $sth->finish;



# prepare data for inserts
# COLUMNS OF FILE ARE AS FOLLOW:

# HGNC ID
# Approved symbol
# Approved name
# Previous symbols
# Synonyms
# NCBI Gene ID
# Ensembl gene ID
# RefSeq IDs
# CCDS IDs
# Locus specific databases


  my $input_file = Text::CSV->new({
    sep_char       => "\t",
    empty_is_undef => 1,
    binary         => 1,
    auto_diag      => 1
  }) or confess "Cannot use file $file: ".Text::CSV->error_diag ();

  my $disk_fh = $xref_dba->get_filehandle($file);
  if ( !defined $disk_fh ) {
    confess "Can't open HGNC file '$file'\n";
  }

  my $mem_file = do {
    undef local $/;
    <$disk_fh>
  };

  close $disk_fh;

  # make sure it's utf8
  utf8::encode($mem_file);
  # get rid of non-conventional " used in the Locus specific databases field
  $mem_file =~ s/"//xg;

  open my $fh, '<', \$mem_file or confess "Can't open HGNC in-memory file: $!\n";

  $input_file->column_names( @{ $input_file->getline( $fh ) } );

  # loop through each row
  while ( my $data = $input_file->getline_hr( $fh ) ) {
    $self->process_row( $data );
  }

  close $fh;

  if ( $verbose ){
    print "HGNC xrefs loaded:\n";
    foreach my $type (sort keys %{$self->{name_count}}){
      print "\t$self->{name_count}->{$type}\t$type\n";
    }
    $self->{mismatch} //= 0;
    print "$self->{mismatch} HGNC ids could not be associated in xrefs\n";
  }
  return 0; # successful
}



=head2 process_row
  Arg [1]    : hashref (data row)
  Description: Processes a row of data and inserts all relevant HGNC xrefs
  Return type: none
  Caller     : internal
=cut

sub process_row {
  my ($self, $data) = @_;

  my $acc              = $data->{'HGNC ID'};
  my $symbol           = $data->{'Approved symbol'};
  my $name             = $data->{'Approved name'};
  my $previous_symbols = $data->{'Previous symbols'};
  my $synonyms         = $data->{'Synonyms'};

  my $seen = 0;

  # Direct CCDS to ENS mappings
  my $ccds = $data->{'CCDS IDs'};
  my @ccds_list;

  if ( defined $ccds ) {
    @ccds_list = split( /,\s/x, $ccds );
  }

  CCDS:
  foreach my $ccds (@ccds_list) {
    my $enst_id = $self->{ccds_to_ens}->{$ccds};

    if (!defined $enst_id) {
      next CCDS;
    }

    $self->{xref_dba}->add_to_direct_xrefs({
        stable_id  => $enst_id,
        type       => 'gene',
        acc        => $acc,
        label      => $symbol,
        desc       => $name,
        source_id  => $self->{source_ids}->{'ccds'},
        species_id => $self->{species_id}
    });

    $self->{xref_dba}->add_synonyms_for_hgnc_vgnc({
        source_id  => $self->{source_ids}->{'ccds'},
        name       => $acc,
        species_id => $self->{species_id},
        dead       => $previous_symbols,
        alias      => $synonyms
    });
    $self->{name_count}->{ccds}++;
  }

  # Direct LRG to ENST mappings
  my $lrg_id = $data->{'Locus specific databases'};
  if ( defined $lrg_id && $lrg_id =~ m/(LRG_\d+)|/x ){
    $lrg_id = $1;
  }

  if ( defined $lrg_id ){
    $self->{xref_dba}->add_to_direct_xrefs({
        stable_id   => $lrg_id,
        type        => 'gene',
        acc         => $acc,
        label       => $symbol,
        desc        => $name,
        source_id   => $self->{source_ids}->{'lrg'},
        species_id  => $self->{species_id}
    });

    $self->{xref_dba}->add_synonyms_for_hgnc_vgnc({
        source_id  => $self->{source_ids}->{'lrg'},
        name       => $acc,
        species_id => $self->{species_id},
        dead       => $previous_symbols,
        alias      => $synonyms
    });
    $self->{name_count}->{lrg}++;
  }

  # Direct Ensembl mappings
  my $ensg_id = $data->{'Ensembl gene ID'};
  if ( defined $ensg_id ){
    $seen = 1;

    $self->{xref_dba}->add_to_direct_xrefs({
        stable_id  => $ensg_id,
        type       => 'gene',
        acc        => $acc,
        label      => $symbol,
        desc       => $name,
        source_id  => $self->{source_ids}->{'ensembl_manual'},
        species_id => $self->{species_id}
    });

    $self->{xref_dba}->add_synonyms_for_hgnc_vgnc({
        source_id  => $self->{source_ids}->{'ensembl_manual'},
        name       => $acc,
        species_id => $self->{species_id},
        dead       => $previous_symbols,
        alias      => $synonyms
    });
    $self->{name_count}->{ensembl_manual}++;
  }

  # RefSeq
  my $refseq_id = $data->{'RefSeq IDs'};
  if ($refseq_id) {
    if ( defined $self->{refseq}->{$refseq_id} ){
      $seen = 1;
      foreach my $xref_id ( @{$self->{refseq}->{$refseq_id}} ){
        $self->{xref_dba}->add_dependent_xref({
            master_xref_id => $xref_id,
            acc            => $acc,
            label          => $symbol,
            desc           => $name,
            source_id      => $self->{source_ids}->{'refseq_manual'},
            species_id     => $self->{species_id}
        });
        $self->{name_count}->{refseq_manual}++;
      }

      $self->{xref_dba}->add_synonyms_for_hgnc_vgnc({
          source_id  => $self->{source_ids}->{'refseq_manual'},
          name       => $acc,
          species_id => $self->{species_id},
          dead       => $previous_symbols,
          alias      => $synonyms
      });
    }
  }

  # EntrezGene
  my $entrez_id = $data->{'NCBI Gene ID'};
  if ( defined $entrez_id ){
    if ( defined $self->{entrezgene}->{$entrez_id} ){
      $seen = 1;
      $self->{xref_dba}->add_dependent_xref({
          master_xref_id => $self->{entrezgene}->{$entrez_id},
          acc            => $acc,
          label          => $symbol,
          desc           => $name,
          source_id      => $self->{source_ids}->{'entrezgene_manual'},
          species_id     => $self->{species_id}
      });

      $self->{xref_dba}->add_synonyms_for_hgnc_vgnc({
          source_id  => $self->{source_ids}->{'entrezgene_manual'},
          name       => $acc,
          species_id => $self->{species_id},
          dead       => $previous_symbols,
          alias      => $synonyms
      });
      $self->{name_count}->{entrezgene_manual}++;
    }
  }

  # Store to keep descriptions if stored yet
  if ( !$seen ){
    $self->{xref_dba}->add_xref({
        acc        => $acc,
        label      => $symbol,
        desc       => $name,
        source_id  => $self->{source_ids}->{'desc_only'},
        species_id => $self->{species_id},
        info_type  => "MISC"
    });

    $self->{xref_dba}->add_synonyms_for_hgnc_vgnc({
        source_id  => $self->{source_ids}->{'desc_only'},
        name       => $acc,
        species_id => $self->{species_id},
        dead       => $previous_symbols,
        alias      => $synonyms
    });
    $self->{mismatch}++;
  }

  return;

}


1;
