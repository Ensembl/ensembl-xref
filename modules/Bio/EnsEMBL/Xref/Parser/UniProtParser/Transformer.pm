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



package Bio::EnsEMBL::Xref::Parser::UniProtParser::Transformer;

use strict;
use warnings;

use Carp;


# These sources are artificially populated using data from
# crossreferences so we cannot simply read their names from the
# input. Therefore, keep these names hardcoded in just one place
# (i.e. here) rather than everywhere we use them.
my $PROTEIN_ID_SOURCE_NAME = 'protein_id';
my $UNIPROT_GN_SOURCE_NAME = 'Uniprot_gn';

# List of all 1:1 mapping from a primary-xref sequence type to the
# Ensembl type of corresponding direct xrefs. One-to-many mappings
# (e.g. for the type 'dna') have to be handled the hard way.
my %direct_xref_type_for_seq_type = (
  'peptide' => 'Translation',
);

my $MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD = 2;
my %source_selection_criteria_for_status
  = (
     'Reviewed'   => [ 'Uniprot/SWISSPROT',
                       sub {
                         return 'sequence_mapped';
                       }, ],
     'Unreviewed' => [ 'Uniprot/SPTREMBL', sub {
                         my ( $level ) = @_;
                         return ( $level <= $MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD ) ?
                           'sequence_mapped' :
                           "protein_evidence_gt_$MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD";
                       }, ],
   );

my %protein_id_extraction_recipe_for_database
  = (
     'ChEMBL' => \&_get_protein_id_xref_from_embldb_xref,
     'EMBL'   => \&_get_protein_id_xref_from_embldb_xref,
   );
sub _get_protein_id_xref_from_embldb_xref {
  my ( $embldb_extra_info, $linkage_source_id, $source_id, $species_id ) = @_;

  # For both EMBL and ChEMBL entries protein ID immediately follows
  # their respective accessions i.e will be the first element of extra_info.
  my $protein_id = $embldb_extra_info->[0];

  # Strip the version number, if any, from the protein ID. At the same
  # time, filter out entries with no ID - in which case the ID is a
  # lone hyphen.
  my ( $unversioned_protein_id )
    = ( $protein_id =~ m{
                          \A
                          # Allow hyphens if they are not
                          # the first character
                          ( [^-.] [^.]+ )
                      }msx );

  if ( ! defined $unversioned_protein_id ) {
    return;
  }

  my $xref_link = {
                   'ACCESSION'          => $unversioned_protein_id,
                   'LABEL'              => $protein_id,
                   'LINKAGE_ANNOTATION' => $source_id,
                   'LINKAGE_SOURCE_ID'  => $linkage_source_id,
                   'SOURCE_ID'          => $source_id,
                   'SPECIES_ID'         => $species_id,
              };
  return $xref_link;
}


=head2 new

  Arg [1]    : HashRef arguments for the constructor:
                - ArrayRef accepted_crossreference_sources
                   - list of names of crossreference sources for which
                     we are to generate direct or dependent xrefs
                     along with respective links. By default said list
                     is empty, in which case only primary
                     (sequence-mapped) xrefs are created.
                     To emphasise the point: this option DOES control
                     creation of direct xrefs, not just dependent
                     ones;
                - default_direct_xref_type
                   - Ensembl type (e.g. Gene, Transcript, Translation)
                     to be assigned to direct xrefs in the event of
                     either being unable to unambiguously determine it
                     from sequence type or having failed to extract
                     the sequence type to begin with;
                - general_source_id
                   - Xref-source ID as provided by the xref
                     pipeline. Note that certain parsers, most notably
                     those distributing xrefs across several different
                     sources, do not use this information.
                - species_id
                   - Ensembl ID of the species under
                     consideration. Records pertaining to other
                     species will be quietly ignored;
                - xref_dba
                   - DBAdaptor object passed from the xref pipeline
  Description: Constructor.
  Return type: Transformer object
  Exceptions : throws on failure to load all required maps from the
               database
  Caller     : UniProtParser::run()
  Status     : Stable

=cut

sub new {
  my ( $proto, $arg_ref ) = @_;

  my $crossref_source_names
    = $arg_ref->{'accepted_crossreference_sources'} // [];

  my $self = {
    'crossref_source_whitelist' => {},
    'default_direct_xref_type'  => $arg_ref->{'default_direct_xref_type'},
    'general_source_id'         => $arg_ref->{'general_source_id'},
    'maps'                      => {},
    'species_id'                => $arg_ref->{'species_id'},
    'xref_dba'                  => $arg_ref->{'xref_dba'},
  };
  my $class = ref $proto || $proto;
  bless $self, $class;

  $self->_load_maps();

  # Store this as a hash to speed up lookups
  foreach my $source_name ( @{ $crossref_source_names } ) {
    $self->{'crossref_source_whitelist'}->{$source_name} = 1;
  }

  return $self;
}


# Destructor. Makes sure the clean-up code gets executed regardless of
# whether the user has explicitly called finish() or not.
sub DESTROY {
  my ( $self ) = @_;

  $self->finish();

  return;
}


=head2 finish

  Description: Wrap-up routine. Does nothing at present, could
               e.g. print statistics.
  Return type: none
  Exceptions : none
  Caller     : destructor, UniProtParser::run()
  Status     : Stable

=cut

sub finish {
  my ( $self ) = @_;

  return;
}


=head2 transform

  Arg [1]    : HashRef extracted_record Structured UniProt-KB record as
               provided by an extractor

  Description: Process extracted_record to create Ensembl xref
               entries, in a form which can be consumed by a
               loader. The only loader in existence so far passes its
               input unmodified to
               BaseAdaptor::upload_xref_object_graphs() so that's the
               format we produce here.

               For each UniProt-KB record matching the specified
               species, this method produces:
                - a primary xref and a corresponding sequence-match
                  xref;
                - one or more direct xrefs and corresponding links for
                  records with Ensembl cross-references;
                - one or more dependent xrefs and corresponding links
                  for records with other cross-references and/or
                  declared gene names;
                - synonyms for each of the above, as needed.
               with the exact list of cross-reference sources to
               process defined in
               $self->{crossref_source_whitelist}. It also determines
               the correct source ID for the record's evidence level.

  Return type: HashRef
  Exceptions : throws on processing errors
  Caller     : UniProtParser::run()
  Status     : Stable

=cut

sub transform {
  my ( $self, $extracted_record ) = @_;

  $self->{'extracted_record'} = $extracted_record;
  $self->{'cache'} = {};

  my ( $accession, @synonyms )
    = @{ $extracted_record->{'accession_numbers'} };
  my $record_source_id = $self->_get_record_source_id();

  my $xref_graph_node
    = {
       'ACCESSION'     => $self->_get_accession(),
       'DESCRIPTION'   => $extracted_record->{'description'},
       'INFO_TYPE'     => 'SEQUENCE_MATCH',
       'LABEL'         => $self->_prepare_label(),
       'SEQUENCE'      => $self->_prepare_sequence(),
       'SEQUENCE_TYPE' => $extracted_record->{'sequence'}->{'type'},
       'SOURCE_ID'     => $record_source_id,
       'SPECIES_ID'    => $self->{'species_id'},
       'STATUS'        => 'experimental',
       'SYNONYMS'      => $self->_get_synonyms(),
     };

  # UniProt Gene Names links come from the 'gene_names' fields
  my $genename_dependent_xrefs
    = $self->_make_links_from_gene_names( $xref_graph_node );
  # Do not assign an empty array to DEPENDENT_XREFS, current insertion code
  # doesn't like them.
  if ( scalar @{ $genename_dependent_xrefs } > 0 ) {
    push @{ $xref_graph_node->{'DEPENDENT_XREFS'} }, @{ $genename_dependent_xrefs };
  }

  # All other xref links come from crossreferences
  my ( $direct_xrefs, $dependent_xrefs )
    = $self->_make_links_from_crossreferences( $xref_graph_node );
  # Do not assign empty arrays to FOO_XREFS, current insertion code
  # doesn't like them.
  if ( scalar @{ $direct_xrefs } > 0 ) {
    push @{ $xref_graph_node->{'DIRECT_XREFS'} }, @{ $direct_xrefs };
  }
  if ( scalar @{ $dependent_xrefs } > 0 ) {
    push @{ $xref_graph_node->{'DEPENDENT_XREFS'} }, @{ $dependent_xrefs };
  }

  return $xref_graph_node;
}


=head2 get_source_id_map

  Description: Returns a map between source-name/priority pairs and
               source IDs. This map is used both in the transformer
               itself and e.g. in the steering code to handle the
               setting of release numbers on sources, and it is more
               resource-efficient to only retrieve it from the
               database once - hence this method.
  Return type: HashRef
  Exceptions : throws if the relevant map does not exist
  Caller     : UniProtParser::run()
  Status     : Stable

=cut

sub get_source_id_map {
  my ( $self ) = @_;

  # Just in case, even though we presently call _load_maps() in the
  # constructor so it shouldn't be possible to call this method before
  # maps have been loaded
  if ( ! exists $self->{'maps'}->{'named_source_ids'} ) {
    confess 'Source-ID map is missing';
  }

  return $self->{'maps'}->{'named_source_ids'};
}


sub _load_maps {
  my ( $self ) = @_;

  my $xref_dba = $self->{'xref_dba'};

  my $source_id_map
    = {
       'Uniprot/SWISSPROT'
       => {
           'direct'          => undef,
           'sequence_mapped' => undef,
         },
       'Uniprot/SPTREMBL'
       => {
           'direct'            => undef,
           "protein_evidence_gt_$MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD" => undef,
           'sequence_mapped'   => undef,
         },
     };
  while ( my ( $source_name, $pri_ref ) = each %{ $source_id_map } ) {
    foreach my $priority ( keys %{ $pri_ref } ) {
      $pri_ref->{$priority}
        = $xref_dba->get_source_id_for_source_name( $source_name,
                                                    $priority );
    }
  }
  $self->{'maps'}->{'named_source_ids'}
    = $source_id_map;

  my %dependent_source_map
    = $xref_dba->get_xref_sources();
  $self->{'maps'}->{'dependent_sources'}
    = \%dependent_source_map;

  return;
}


# Returns true if the current record describes an entry derived from
# Ensembl, false otherwise.
sub _entry_is_from_ensembl {
  my ( $self ) = @_;

  # The old parser's way of identifying proteins derived from Ensembl
  # was to search for a comment topic "CAUTION" stating "The sequence
  # shown here is derived from an Ensembl". It was confirmed with
  # UniProt in November 2018 that this hadn't been done since late
  # 2017, therefore we no longer check this.

  # As of November 2018, the official way of identifying proteins
  # derived from Ensembl is to look for a special citation identifying
  # the record as such. There are several ways in which such special
  # citations could be identified, it seems the simplest to look at
  # "reference group name" fields to see if any of them says
  # 'Ensembl'.
  my $citation_groups
    = $self->{'extracted_record'}->{'citation_groups'};
  my $from_group_names
    = List::Util::any { $_ eq 'Ensembl' } @{ $citation_groups };

  return $from_group_names;
}


# Get primary accession of the extracted record.
sub _get_accession {
  my ( $self ) = @_;

  if ( ! exists $self->{'cache'}->{'accession'} ) {
    my ( $accession, @synonyms )
      = @{ $self->{'extracted_record'}->{'accession_numbers'} };
    $self->{'cache'}->{'accession'} = $accession;
    $self->{'cache'}->{'synonyms'} = \@synonyms;
  }

  return $self->{'cache'}->{'accession'};
}


# Get synonyms (secondary accessions) of the extracted record.
sub _get_synonyms {
  my ( $self ) = @_;

  if ( ! exists $self->{'cache'}->{'synonyms'} ) {
    my ( $accession, @synonyms )
      = @{ $self->{'extracted_record'}->{'accession_numbers'} };
    $self->{'cache'}->{'accession'} = $accession;
    $self->{'cache'}->{'synonyms'} = \@synonyms;
  }

  return $self->{'cache'}->{'synonyms'};
}


# Translate quality of the extracted entry into the matching Ensembl
# source_id, optionally with an override of priority.
sub _get_source_id_from_quality {
  my ( $self, $priority_override ) = @_;

  my $source_id_map = $self->{'maps'}->{'named_source_ids'};

  my $entry_quality = $self->{'extracted_record'}->{'quality'};
  my $criteria = $source_selection_criteria_for_status{ $entry_quality->{'status'} };
  my $priority_mapper = $criteria->[1];

  my $source_name = $criteria->[0];
  my $priority = $priority_override
    // $priority_mapper->( $entry_quality->{'evidence_level'} );

  return $source_id_map->{$source_name}->{$priority};
}


# Make xrefs from 'crossreferences' entries in the extracted record,
# in a form suitable to attaching to the main xref's graph node as
# consumed by upload_xref_object_graphs(). Ensembl crossreferences
# become direct xrefs, everything else - dependent ones. If requested
# we additionally generate protein_id dependent xrefs from appropriate
# sources, i.e. EMBL and ChEMBL at present.
sub _make_links_from_crossreferences {
  my ( $self, $primary_xref ) = @_;

  my $crossreferences = $self->{'extracted_record'}->{'crossreferences'};
  my $dependent_sources = $self->{'maps'}->{'dependent_sources'};
  my %whitelisted_crossreference_sources
    = %{ $self->{'crossref_source_whitelist'} };

  my @direct_xrefs;
  my @dependent_xrefs;

 REF_SOURCE:
  while ( my ( $source, $entries ) = each %{ $crossreferences } ) {

    if ( ! $whitelisted_crossreference_sources{ $source } ) {
      next REF_SOURCE;
    }

    if ( $source eq 'Ensembl' ) {

      # If we cannot unambiguously deduce Ensembl type from the
      # sequence type of the primary xref (or if it hasn't been
      # extracted at all), fall back to the global default. Should be
      # safe enough because on the one hand if we do enable creation
      # of direct xrefs from specific input we really ought to have an
      # idea of what type they are, and on the other unless the
      # database schema changes a lot the loader should complain VERY
      # loudly when it tries to insert data into the very much
      # nonexistent table '_direct_xref' if the global default is
      # omitted.
      my $direct_xref_type
        = $direct_xref_type_for_seq_type{ $primary_xref->{'SEQUENCE_TYPE'} }
        // $self->{'default_direct_xref_type'};

    DIRECT_XREF:
      foreach my $direct_ref ( @{ $entries } ) {
        my $xref_link
          = {
             # We want translation ID and 'id' for these is TRANSCRIPT ID
             'STABLE_ID'    => $direct_ref->{'optional_info'}->[0],
             'ENSEMBL_TYPE' => $direct_xref_type,
             'LINKAGE_TYPE' => 'DIRECT',
             'SOURCE_ID'    => $self->_get_source_id_from_quality( 'direct' ),
             'ACCESSION'    => $primary_xref->{'ACCESSION'},
             'VERSION'      => $primary_xref->{'VERSION'},
             'LABEL'        => $primary_xref->{'LABEL'},
             'DESCRIPTION'  => $primary_xref->{'DESCRIPTION'},
             'SPECIES_ID'   => $primary_xref->{'SPECIES_ID'},
             'INFO_TEXT'    => $primary_xref->{'INFO_TEXT'},
           };
        push @direct_xrefs, $xref_link;
      }

    }
    elsif ( exists $dependent_sources->{$source} ) {
      my $dependent_source_id = $dependent_sources->{$source};

    DEPENDENT_XREF:
      foreach my $dependent_ref ( @{ $entries } ) {
        my $xref_link
          = {
             'ACCESSION'          => $dependent_ref->{'id'},
             'LINKAGE_ANNOTATION' => $dependent_source_id,
             'LINKAGE_SOURCE_ID'  => $primary_xref->{'SOURCE_ID'},
             'SOURCE_ID'          => $dependent_source_id,
             'SPECIES_ID'         => $primary_xref->{'SPECIES_ID'},
           };
        push @dependent_xrefs, $xref_link;

        my $protein_id_xref_maker
          = $protein_id_extraction_recipe_for_database{ $source };
        if ( $whitelisted_crossreference_sources{ $PROTEIN_ID_SOURCE_NAME }
             && ( defined $protein_id_xref_maker ) ) {

          # Entries for the source 'protein_id' are constructed from
          # crossreferences to other databases
          my $protein_id_xref
            = $protein_id_xref_maker->( $dependent_ref->{'optional_info'},
                                        $primary_xref->{'SOURCE_ID'},
                                        $dependent_sources->{$PROTEIN_ID_SOURCE_NAME},
                                        $primary_xref->{'SPECIES_ID'}
                                     );
          if ( defined $protein_id_xref ) {
            push @dependent_xrefs, $protein_id_xref;
          }

        }

      }

    }
  }

  return ( \@direct_xrefs, \@dependent_xrefs );
}


# Make Uniprot_gn dependent xrefs from 'gene_names' entries in the
# extracted record, in a form suitable to attaching to the main xref's
# graph node as consumed by upload_xref_object_graphs().
sub _make_links_from_gene_names {
  my ( $self, $primary_xref ) = @_;

  my @genename_xrefs;

  my %whitelisted_crossreference_sources
    = %{ $self->{'crossref_source_whitelist'} };

  # Are we supposed to process this xref source to begin with?
  if ( ! $whitelisted_crossreference_sources{ $UNIPROT_GN_SOURCE_NAME } ) {
    return [];
  }

  # UniProt Gene Name xrefs are dependent so in order to avoid
  # circular dependencies, do not generate them for proteins derived
  # from Ensembl.
  if ( $self->_entry_is_from_ensembl() ) {
    return [];
  }

  my $gene_names = $self->{'extracted_record'}->{'gene_names'};
  my $dependent_sources = $self->{'maps'}->{'dependent_sources'};
  my $dependent_source_id = $dependent_sources->{$UNIPROT_GN_SOURCE_NAME};

 GN_ENTRY:
  foreach my $gn_entry ( @{ $gene_names } ) {
    if ( ! exists $gn_entry->{'Name'} ) {
      next GN_ENTRY;
    }

    my $name = $gn_entry->{'Name'};
    my $xref = {
                'ACCESSION'          => $primary_xref->{'ACCESSION'},
                'LABEL'              => $name,
                'LINKAGE_ANNOTATION' => $dependent_source_id,
                'LINKAGE_SOURCE_ID'  => $primary_xref->{'SOURCE_ID'},
                'SOURCE_ID'          => $dependent_source_id,
                'SPECIES_ID'         => $primary_xref->{'SPECIES_ID'},
              };

    my $synonyms = $gn_entry->{'Synonyms'};
    if ( defined $synonyms ) {
      push @{ $xref->{'SYNONYMS'} }, @{ $synonyms };
    }

    push @genename_xrefs, $xref;

    # Regardless of how many gene-name entries a protein might have in
    # UniProt, we only want at most *one* such xref - there are bits
    # of Ensembl code which implicitly assume there can only be one
    # xref per unique (accerssion,source_id,species_id) triplet and
    # Uniprot_gn is essentially a last-resort source for when nothing
    # else works. At least we should be able to consistently output
    # the first encountered gene name (and synonyms), the old parser
    # happily ignored the fact there can be multiple gene names and
    # clobbered old names and synonyms every time a new one was
    # encountered - potentially resulting in xrefs with the last
    # encountered name but synonyms from an earlier one.
    last GN_ENTRY;
  }

  return \@genename_xrefs;
}


# UniProt-specific implementation of the source-ID getter: select the
# right one depending on quality of the entry.
sub _get_record_source_id {
  my ( $self ) = @_;

  return $self->_get_source_id_from_quality();
}


# UniProt-specific implementation of the label getter: should be the
# same as the primary accession.
sub _prepare_label {
  my ( $self ) = @_;

  return $self->_get_accession();
}


# UniProt-specific implementation of the sequence getter: no
# transformations needed.
sub _prepare_sequence {
  my ( $self ) = @_;

  return $self->{'extracted_record'}->{'sequence'}->{'seq'};
}


1;
