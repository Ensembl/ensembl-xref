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

Bio::EnsEMBL::Xref::Parser::RefSeqGPFFParser

=head1 DESCRIPTION

A parser class for the RefSeq_peptide and RefSeq_dna sources.
Data files are in genbank format and should be of types
*.protein.gpff or *.rna.gbff.
rel_file is mandatory and should be thge release file for the
provided data files.

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::RefSeqGPFFParser->new(
    source_id  => 108,
    species_id => 9606,
    species    => 'homo_sapiens',
    files      => ['data.protein.gpff'],
    rel_file   => 'release.txt',
    xref_dba   => $xref_dba
  );

  $parser->run();

=cut


package Bio::EnsEMBL::Xref::Parser::RefSeqGPFFParser;

use strict;
use warnings;
use Carp;
use List::Util 1.45 qw(uniq);

use parent qw( Bio::EnsEMBL::Xref::Parser );


=head2 run
  Description: Runs the RefSeqGPFFParser
  Return type: N/A
  Caller     : internal
=cut

sub run {
  my ( $self ) = @_;

  my $source_id    = $self->{source_id};
  my $species_id   = $self->{species_id};
  my $species_name = $self->{species};
  my $files        = $self->{files};
  my $release_file = $self->{rel_file};
  my $xref_dba     = $self->{xref_dba};
  my $verbose      = $self->{verbose} // 0;

  if ( (!defined $source_id) || (!defined $species_id) || (!defined $species_name) || (!defined $files) || (!defined $release_file) ){
    confess "Need to pass source_id, species_id, species, files and rel_file as pairs";
  }

  $self->{refseq_sources} = $xref_dba->get_refseq_sources;

  # get RefSeq source ids
  while (my ($source_prefix, $source_name) = each %{$self->{refseq_sources}}) {
    $self->{source_ids}->{$source_name} = $xref_dba->get_source_id_for_source_name( $source_name, undef )
  }

  # get extra source ids
  $self->{source_ids}->{EntrezGene} = $xref_dba->get_source_id_for_source_name( 'EntrezGene', undef );
  $self->{source_ids}->{WikiGene} = $xref_dba->get_source_id_for_source_name( 'WikiGene', undef );

  # Retrieve existing RefSeq mRNA
  $self->{refseq_ids} = { %{$xref_dba->get_valid_codes('RefSeq_mRNA', $species_id )},
      %{$xref_dba->get_valid_codes('RefSeq_mRNA_predicted', $species_id )} };
  $self->{entrez_ids} = $xref_dba->get_valid_codes('EntrezGene', $species_id );
  $self->{wiki_ids} = $xref_dba->get_valid_codes('WikiGene', $species_id );

  if ($verbose) {
    for my $source_name (sort values %{$self->{refseq_sources}}) {
      print "$source_name source ID = $self->{source_ids}->{$source_name}\n";
    }
  }

  # populate entrez gene id => label hash
  $self->{entrez} = $xref_dba->get_acc_to_label('EntrezGene', $species_id, undef );

  # process the source files
  foreach my $file (@{$files}) {

    # type from the file (peptide or dna)
    my $type = $self->type_from_file($file);

    # get the file handler
    my $refseq_fh = $xref_dba->get_filehandle($file);

    if ( !defined $refseq_fh ) {
      confess "Can't open RefSeqGPFF file '$file'\n";
    }

    # this will hold the array of xrefs to bulk insert
    my $xrefs;

    {
      local $/ = "//\n";
      while ( my $genbank_rec = $refseq_fh->getline() ) {
        my $xref = $self->xref_from_record( $genbank_rec, $type );
         if ($xref) {
          push @{$xrefs}, $xref;
        }
      }
    }

    $refseq_fh->close();

    # upload the xrefs
    $xref_dba->upload_xref_object_graphs( $xrefs ) if $xrefs;

  }

  # get the release file handle
  my $release_fh = $xref_dba->get_filehandle($release_file);

  # get file header
  my $release = do { local $/ = "\n*"; <$release_fh> };
  $release_fh->close();

  # replace all whitespace with a single space
  $release =~ s/\s+/ /xg;

  # (NCBI.*Release  > match and get NCBI blablabla Release
  # \s\d+)          > the release number
  # \s(.*)          > the release date
  # \sDistribution  > anchor
  if ( $release =~ m/(NCBI.*Release\s\d+)\s(.*)\sDistribution/x ) {
    my ($rel_number, $rel_date) = ($1, $2);
    my $release_string = "$rel_number, $rel_date";

    # set release info
    $xref_dba->set_release( $source_id, $release_string );
    for my $source_name (sort values %{$self->{refseq_sources}}) {
      $xref_dba->set_release( $self->{source_ids}->{$source_name}, $release_string );
    }

    print "RefSeq release: '$release_string'\n" if $verbose;

  } else {
    print "Could not set release info from release file '$release_file'\n";
  }

  return 0;
}



=head2 xref_from_record
  Arg [1]    : scalar (string genbank_record)
  Arg [2]    : scalar (string type)
  Description: Processes the provided genbank record of the provided type
               returns xref for bulk insert and creates related dependent_xrefs
  Return type: hashref
  Caller     : internal
=cut

sub xref_from_record {
  my ($self, $genbank_rec, $type) = @_;

  # Get the record species
  my ($record_species) = $genbank_rec =~ /ORGANISM\s*(.*)\n/x;
  $record_species = lc $record_species;
  $record_species =~ s/\s+/_/xg;

  # skip if species is not the required
  return unless ( defined $record_species && ($record_species eq $self->{species}) );


  my ($acc) = $genbank_rec =~ /ACCESSION\s+(\S+)/x;

  my $prefix = substr($acc, 0, 2);

  # skip if acc is not of known type
  return if ( !exists $self->{refseq_sources}->{$prefix} );


  my $acc_source_id = $self->source_id_from_acc($acc);

  my $xref = {
    ACCESSION     => $acc,
    SPECIES_ID    => $self->{species_id},
    SOURCE_ID     => $acc_source_id,
    SEQUENCE_TYPE => $type,
    INFO_TYPE     => 'SEQUENCE_MATCH'
  };

  my ($ver_acc, $ver_num) = $genbank_rec =~ /VERSION\s+(\w+)\.(\d+)/x;

  if ($acc eq $ver_acc) {
    $xref->{LABEL} = "${acc}.${ver_num}";
    $xref->{VERSION} = $ver_num;
  } else {
    print "WARNING: accession $acc mismatch with version ${acc}.${ver_num}\n" if $self->{verbose};
  }

  my ($description) = $genbank_rec =~ /DEFINITION\s+  # Find the field identifier
                                  (.*)                # get the description
                                  \s+ACCESSION/xms;   # until the next field

  # TODO remove when upload_xref_object_graphs() in BaseAdaptor is updated to do this automatically
  # remove any newlines and spaces, and make sure is within mysql limits
  $description =~ s/\n//xg;
  $description =~ s/\s+/ /xg;
  $description = substr($description, 0, 255) if (length($description) > 255);

  $xref->{DESCRIPTION} = $description;

  # sequence is multiline, each line starts with base number and has spaces all over. ends with //
  my ($seq) = $genbank_rec =~ /\s*ORIGIN\s+ # Find the field identifier
                          (.+)         # get all sequence lines
                          \/\//xms;    # until the end of the field (//)

  # get rid of the base number and the whitespace for a sequence string
  $seq =~ s/[\d\s]+//xg;

  $xref->{SEQUENCE} = $seq;


  my @protein_ids = $genbank_rec =~ /\/protein_id=\"(.+?)\"/xg;

  my $protein_id = pop @protein_ids;

  $xref->{PROTEIN} = $protein_id if defined $protein_id;


  my @coded_by_list = $genbank_rec =~ /\/coded_by=\"(.*?):/xg;

  my $coded_by = pop @coded_by_list;

  $xref->{PAIR} = $coded_by if defined $coded_by;

  my ($refseq_pair) = $genbank_rec =~ /DBSOURCE\s+REFSEQ: accession (\S+)/x;

  if (defined $refseq_pair && !exists $xref->{PAIR}) {
    $xref->{PAIR} = $refseq_pair;
  }


  my @gene_ids = $genbank_rec =~ /db_xref=\"GeneID:(.+?)\"/xg;
  @gene_ids = uniq @gene_ids;

  # process existing entrez_gene_ids as dependent xrefs
  GENEID:
  foreach my $gene_id (@gene_ids) {

    next GENEID unless (defined $self->{entrez}->{$gene_id});

    push @{$xref->{DEPENDENT_XREFS}}, {
        SOURCE_ID         => $self->source_id_from_name('EntrezGene'),
        LINKAGE_SOURCE_ID => $acc_source_id,
        ACCESSION         => $gene_id,
        LABEL             => $self->{entrez}->{$gene_id}
    };

    push @{$xref->{DEPENDENT_XREFS}}, {
        SOURCE_ID         => $self->source_id_from_name('WikiGene'),
        LINKAGE_SOURCE_ID => $acc_source_id,
        ACCESSION         => $gene_id,
        LABEL             => $self->{entrez}->{$gene_id}
    };

    next GENEID unless (defined $refseq_pair);

    # remove the version number
    $refseq_pair =~ s/\.\d*//x;

    # Add xrefs for RefSeq mRNA as well where available
    foreach my $refseq_acc (@{ $self->{refseq_ids}->{$refseq_pair} }) {
      foreach my $entrez_id (@{ $self->{entrez_ids}->{$gene_id} }) {
        $self->{xref_dba}->add_dependent_xref({
          master_xref_id => $refseq_acc,
          acc            => $entrez_id,
          source_id      => $self->source_id_from_name('EntrezGene'),
          species_id     => $self->{species_id},
        });
      }
      foreach my $wiki_id (@{ $self->{wiki_ids}->{$gene_id} }) {
        $self->{xref_dba}->add_dependent_xref({
          master_xref_id => $refseq_acc,
          acc            => $wiki_id,
          source_id      => $self->source_id_from_name('WikiGene'),
          species_id     => $self->{species_id},
        });
      }
    }
  }

  return $xref;

}



=head2 source_id_from_name
  Arg [1]    : Scalar (string name)
  Description: Provided the name of a source, returns the source_id.
               Requires $self->{source_ids} to have been populated.
  Return type: Scalar (integer source_id)
  Caller     : internal
=cut

sub source_id_from_name {
  my ($self, $name) = @_;

  my $source_id;

  if ( exists $self->{source_ids}->{$name} ) {
    $source_id = $self->{source_ids}->{$name};
  } else {
    confess "Can't get source ID for name '$name'\n";
  }

  return $source_id;
}


=head2 source_id_from_acc
  Arg [1]    : Scalar (string acc)
  Description: Provided a RefSeq accession, returns the source_id.
               Requires $self->{source_ids} to have been populated.
  Return type: Scalar (integer source_id)
  Caller     : internal
=cut

sub source_id_from_acc {
  my ($self, $acc) = @_;

  my $source_id;
  my $prefix = substr($acc, 0, 2);

  if ( exists $self->{refseq_sources}->{$prefix} ) {
    $source_id = $self->source_id_from_name( $self->{refseq_sources}->{$prefix} );
  } else {
    confess "Can't get source ID for accession '$acc'\n";
  }

  return $source_id;
}


=head2 type_from_file
  Arg [1]    : Scalar (string file)
  Description: Provided the data filename of the form *.(protein|rna).*
               returns the appropriate type:
               protein => peptide
               rna => dna
  Return type: Scalar (string type)
  Caller     : internal
=cut

sub type_from_file {
  my ($self, $file) = @_;

  my ($type) = $file =~ /\.(protein|rna)\./x;

  if ( $type eq 'protein' ) {
    $type = 'peptide';
  } elsif ( $type eq 'rna') {
    $type = 'dna';
  }

  if ( !defined $type ) {
    confess "Could not work out sequence type for file '$file'\n";
  }

  return $type;
}

1;
