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



package Bio::EnsEMBL::Xref::Parser::UniProtParser::Extractor;

use strict;
use warnings;

# For non-destructive substitutions in regexps (/r flag)
require 5.014_000;

use Carp;
use List::Util;
use charnames ':full';


# While processing DR fields i.e. crossreferences, used to match the
# syntax attaching a crossreference to a specific isoform. Declares
# ONE capture group, the isoform identifier.
# Use named sequences for square brackets to avoid excessive
# escaping as well as for better readability.
my $QR_DR_ISOFORM_FIELD_PATTERN
  = qr{
        \s*
        \N{LEFT SQUARE BRACKET}
        \s*
        ( [^\N{RIGHT SQUARE BRACKET}]+ )
        \s*
        \N{RIGHT SQUARE BRACKET}
        \s*
    }msx;

# While processing ID fields, used to confirm that the status of an
# entry matches an expected value. Declares NO capture groups.
my $QR_ID_STATUS_FIELD
  = qr{
        Unreviewed | Reviewed
    }msx;

# While processing OX fields i.e. taxonomy cross-references, used to
# extract both the taxon code and the database qualifier. Declares TWO
# capture groups: the database qualifier, and the taxonomic code.
my $QR_OX_TAXON_DB_ENTRY
  = qr{
        # Database qualifier. Chances are the list of
        # allowed characters will change should DBs
        # other than NCBI ever become supported here.
        ( [A-Za-z_]+ )

        \s*  # just in case
        =
        \s*  # same

        # Taxon ID. This is almost certainly NCBI-specific.
        ( [0-9]+ )
    }msx;

# While processing OX or DE fields, i.e. taxonomy cross-references or
# descriptions, allows accounting for the fact some entries might be
# followed by evidence codes. As of October 2018, this syntax is not
# declared in UniProt-KB User Manual yet frequently encountered in
# data files.  Use named sequences for curly brackets to avoid
# excessive escaping as well as for better readability.
my $QR_OX_DE_EVIDENCE_CODE_LIST
  = qr{
        \N{LEFT CURLY BRACKET}
        \s*
        [^\N{RIGHT CURLY BRACKET}]+
        \s*
        \N{RIGHT CURLY BRACKET}
        \s*
    }msx;

# Various field separators
my $QR_SPLIT_COMMA_WITH_WHITESPACE
  = qr{ \s* , \s* }msx;
my $QR_SPLIT_SEMICOLON_WITH_WHITESPACE
  = qr{ \s* ; \s* }msx;

# The default value of sequence type is a longish random string in
# order to minimise the risk of a database schema change resulting in
# whatever we put in here suddenly becoming an accepted value, as it
# could conceivably be the case for e.g. 'unknown'. It is in a package
# constant so that it will be easy to reach by any validation routines
# we might add at some point.
my $SEQUENCE_TYPE_UNDEFINED = 'b0MkJApqLUNCJwp5U3JmXgpzQmtACmZGfScKP1h';

# Translates the sequence unit from the SQ line into Ensembl sequence
# type as used in e.g. primary xrefs. We do it in the extractor and
# not in a transformer because sequence-type identification depends on
# the input format, and if we touch it to begin with why not output it
# in form which will not require further transformation.
my %sequence_type_for_unit = (
  'AA' => 'peptide',
  'BP' => 'dna',
);

# Syntax: 0 for database qualifiers to be ignored, otherwise a
# reference to a function translating taxon codes from the given
# database into Ensembl taxonomy_ids.
my %taxonomy_ids_from_taxdb_codes
  = (
     # NCBI taxon codes and Ensembl taxonomy IDs are identical
     'NCBI_TaxID' => sub { return $_[0]; },
   );


=head2 new

  Arg [1]    : HashRef arguments for the constructor:
                - ArrayRef file_names
                   - list of names of files to process. Note that as
                     with the old UniProtParser, and indeed most other
                     parsers, only the first one will actually be used;
                - ArrayRef mandatory_prefixes
                   - list of prefixes of fields that are to be read
                     into memory and which must be present in every
                     record. Defaults to none. Takes precedence over
                     optional_prefixes; in the event of a prefix
                     having been declared in both;
                - ArrayRef optional_prefixes
                   - list of prefixes of fields that are to be read
                     into memory but which can in principle be absent
                     from a record.  Defaults to none. Overridden by
                     mandatory_prefixes in the event of a prefix
                     having been declared in both;
                - species_id
                - species_name
                   - Ensembl ID and/or name of the species under
                     consideration. Records pertaining to other
                     species will be quietly ignored.
                - xref_dba
                   - DBAdaptor object passed from the xref pipeline
  Description: Constructor.
  Return type: Extractor object
  Exceptions : throws on failure to acquire the file handle
  Caller     : UniProtParser::run()
  Status     : Stable

=cut

sub new {
  my ( $proto, $arg_ref ) = @_;
  my $file_names         = $arg_ref->{'file_names'};
  my $mandatory_prefixes = $arg_ref->{'mandatory_prefixes'} // [];
  my $optional_prefixes  = $arg_ref->{'optional_prefixes'}  // [];
  my $species_id         = $arg_ref->{'species_id'};
  my $species_name       = $arg_ref->{'species_name'};
  my $xref_dba           = $arg_ref->{'xref_dba'};

  my $filename = $file_names->[0];
  my $filehandle = $xref_dba->get_filehandle( $filename );

  # Keep the file name for possible debugging purposes, unless we can
  # somehow retrieve it from _io_handle
  my $self = {
    'input_name'           => $filename,
    '_io_handle'           => $filehandle,
    'prefixes_of_interest' => {},
    'species_id'           => $species_id,
    'species_name'         => $species_name,
    'xref_dba'             => $xref_dba,
  };
  my $class = ref $proto || $proto;
  bless $self, $class;

  # Store these as a hash to speed up lookups. Process optional
  # prefixes first so that in the event of a prefix appearing on both
  # list, it is treated as mandatory.
  foreach my $prefix ( @{ $optional_prefixes } ) {
    $self->{'prefixes_of_interest'}->{$prefix} = 0;
  }
  foreach my $prefix ( @{ $mandatory_prefixes } ) {
    $self->{'prefixes_of_interest'}->{$prefix} = 1;
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


=head2 extract

  Description: Extract information from the UniProt-KB record loaded
               into memory by get_uniprot_record(), apply basic sanity
               and metadata checks, and produce a key-value
               representation of the record suitable for further
               processing.

               Note that the purpose of this method is purely to
               provide an input format-independent representation of
               provided data regardless of whether it has come from a
               plain-text file (as it is the case here), an XML dump,
               or directly from the database; anything beyond that
               belongs to the transform stage. It also should
               generally avoid performing actual retrieval of input
               data.

  Return type: HashRef
  Exceptions : throws on processing errors
  Caller     : UniProtParser::run()
  Status     : Stable

=cut

sub extract {
  my ( $self ) = @_;

  $self->_record_has_all_needed_fields();

  my $accession_numbers = $self->_get_accession_numbers();

  # Only proceed if at least one taxon code in the entry maps
  # to the current species ID, and skip unreviewed entries.
  if ( ( ! $self->_record_species_matches() )
       || ( $self->_entry_is_unreviewed( $accession_numbers ) ) ) {
    return 'SKIP';
  }

  my $entry_object
    = {
       'accession_numbers' => $accession_numbers,
       'citation_groups'   => $self->_get_citation_groups(),
       'crossreferences'   => $self->_get_database_crossreferences(),
       'description'       => $self->_get_description() // undef,
       'gene_names'        => $self->_get_gene_names(),
       'quality'           => $self->_get_quality(),
       'sequence'          => $self->_get_sequence(),
     };

  return $entry_object;
}


=head2 finish

  Description: Wrap-up routine. At the moment all it does is close the
               file handle, could also e.g. print statistics.
  Return type: none
  Exceptions : none
  Caller     : destructor, UniProtParser::run()
  Status     : Stable

=cut

sub finish {
  my ( $self ) = @_;

  $self->{'_io_handle'}->close();

  return;
}


=head2 get_uniprot_record

  Description: Load the next record from input, if any, into memory,
               and indicate availability of new data through the
               return value.

               One thing that this method absolutely *must* do is
               indicate if the input contains another *complete*
               UniProt-KB record. Everything else, including whether
               or not the record is actually retrieved, scheduled for
               asynchronous retrieval or merely prepared for lazy
               loading may depend on the specific input type, as long
               there is a reasonable level of confidence that we
               *will* eventually get a complete record.

               In case of plain-text files processed by this
               particular implementation, we have to read the input in
               line by line because there is nothing in either the
               format or the data that would allow us to skip
               ahead. We do, however, avoid any processing other than
               stripping trailing newline characters and minimal
               "indexing" (see below), and even that is only performed
               on lines we know we will need (see
               $self->{'prefixes_of_interest'} above).

               The output in this case is a hashref assigned to the
               'record' property of the object, in which keys are
               names of fields and values are arrayrefs containing
               matching lines (minus the field name and the
               three-blank separator) in the order they appeared in
               the input file.

  Return type: boolean 1 if a record has been loaded, 0 if end of input
  Exceptions : throws on processing errors
  Caller     : UniProtParser::run()
  Status     : Stable

=cut

sub get_uniprot_record {
  my ( $self ) = @_;

  my $io = $self->{'_io_handle'};
  my $uniprot_record = {};
  my %prefixes_of_interest = %{ $self->{'prefixes_of_interest'} };

 INPUT_LINE:
  while ( my $file_line = $io->getline() ) {
    my $prefix = substr( $file_line, 0, 2 );

    if ( $prefix eq q{//} ) {
      # End of record, return what we have got so far
      $self->{'record'} = $uniprot_record;

      return 1;
    }

    # Do not waste time and memory on fields we do not need
    if ( ! exists $prefixes_of_interest{$prefix} ) {
      next INPUT_LINE;
    }

    chomp ( my $content = substr( $file_line, 5 ) );

    push @{ $uniprot_record->{$prefix} }, $content;

  }

  # If we began parsing fields but have never reached the //,
  # something is very wrong
  if ( scalar keys %{ $uniprot_record } > 0 ) {
    confess 'Incomplete input record';
  }

  # EOF
  return 0;
}


# Returns true if the current record describes an entry tagged as
# unreviewed, false otherwise.
sub _entry_is_unreviewed {
  my ( $self, $accession_numbers ) = @_;

  # This is the way the old UniProtParser identified unreviewed
  # entries. FIXME: is this still a thing? As of October 2018 there
  # are NO such entries in either the first ~1000 lines of the TrEMBL
  # file or anywhere in the SwissProt one.
  if ( lc( $accession_numbers->[0] ) eq 'unreviewed' ) {
    return 1;
  }

  return 0;
}


# Parse the AC fields of the current record and produce a list of
# UniProt accession numbers. The list will reflect the order in which
# accession numbers appeared in the record, which as of October 2018
# is: primary accession number first, then all the secondary ones in
# alphanumerical order.
sub _get_accession_numbers {
  my ( $self ) = @_;

  my $ac_fields = $self->{'record'}->{'AC'};
  if ( ! defined $ac_fields ) {
    return [];
  }

  my @numbers
    = split( $QR_SPLIT_SEMICOLON_WITH_WHITESPACE,
             join( q{}, @{ $ac_fields } ) );
  # FIXME: we should probably make this persist until a new record has
  # been loaded

  return \@numbers;
}


# Parse the RG fields of the current record, if any, and produce a
# list of names of groups which cite this protein.
# Warning: this is intentionally simplistic (for example it makes no
# attempt to associate group names for specific reference numbers) and
# will NOT work properly if more detailed citation information is
# required! If you need that, you will have to implement comprehensive
# processing of different Rx lines.
sub _get_citation_groups {
  my ( $self ) = @_;

  my $rg_fields = $self->{'record'}->{'RG'};
  if ( ! defined $rg_fields ) {
    return [];
  }

  # FIXME: we should probably make this persist until a new record has
  # been loaded
  my @citation_groups = split( $QR_SPLIT_SEMICOLON_WITH_WHITESPACE,
                               join( q{}, @{ $rg_fields } ) );

  return \@citation_groups;
}


# Parse the DR fields of the current record, break them into
# constituent parts and produce a list (or to be precise, a hash of
# arrayrefs) of database cross-references grouped by reference.
sub _get_database_crossreferences {
  my ( $self ) = @_;

  my $dr_fields = $self->{'record'}->{'DR'};
  if ( ! defined $dr_fields ) {
    return {};
  }

  # FIXME: we should probably make this persist until a new record has
  # been loaded
  my $crossreferences = {};

  foreach my $dr_line ( @{ $dr_fields } ) {
    my ( $res_abbrev, $res_id, @opts )
      = split( $QR_SPLIT_SEMICOLON_WITH_WHITESPACE, $dr_line);

    my ( $last_opt, $isoform )
      = ( $opts[-1] =~ m{
                          ( .+ )  # will grab all dots but the last one
                          [.]
                          (?:
                            $QR_DR_ISOFORM_FIELD_PATTERN
                          )?
                          \z
                      }msx );
    if ( ! defined $last_opt ) {
      confess "Malformed final-option match in:\n\t$dr_line";
    }

    # At the very least, strips the trailing dot
    $opts[-1] = $last_opt;

    my $crossref
      = {
         'id'            => $res_id,
         'optional_info' => \@opts,
       };
    if ( defined $isoform ) {
      $crossref->{'target_isoform'} = $isoform;
    }

    # There CAN be multiple cross-references to a database
    push @{ $crossreferences->{$res_abbrev} }, $crossref;
  }

  return $crossreferences;
}


# Parse the DE fields of the current record and produce a description
# string compatible with the output of the old UniProtParser, namely:
#  - the description begins with a semicolon-separated list of
#    top-level names (i.e. ones not belonging to a Contains or
#    Includes section), in the order they appear in the record;
#  - this list is followed by a space and a space-separated list of
#    names from Contains and Includes sections, again in the order
#    they appear in the record;
#  - we process both RecNames and SubNames, and both types are
#    considered of equal priority (i.e. ultimately it is their order
#    that matters);
#  - in either case we only consider full names;
#  - evidence codes, PubMed references etc. are discarded.
# Note that unlike most field parsers implemented so far, this one
# does NOT attempt to fully process the syntax of DE fields; this is
# in order to avoid messing with whitespace-defined context only to
# discard most of the field data anyway. Keep this in mind should you
# want to extend the parsing to e.g. extract EC numbers, in which case
# you will likely have to implement a full parser.
sub _get_description {
  my ( $self ) = @_;

  my $de_fields = $self->{'record'}->{'DE'};
  if ( ! defined $de_fields ) {
    return;
  }

  my @names;
  my @subdescs;

  foreach my $line ( @{ $de_fields } ) {
    my ( $indent, $content )
      = ( $line =~ m{
                      \A
                      ( \s* )  # To detect sub-entries such as Includes:
                      (?:
                        RecName | SubName
                      )
                      :
                      \s*
                      Full=
                      ( [^;]+ )
                  }msx );
    if ( defined $indent ) {
      # Strip evidence codes, if any. This could in principle be
      # integrated into the match but it makes the regex rather more
      # complex.
      $content =~ s{
                     \s*
                     $QR_OX_DE_EVIDENCE_CODE_LIST
                     \z
                 }{}msx;

      if ( $indent eq q{} ) {
        push @names, $content;
      }
      else {
        push @subdescs, $content;
      }
    }
  }

  my $description
    = join( q{ }, (
                   join( q{;}, @names ),
                   @subdescs
                 ) );
  # Make sure we do not return an empty string
  if ( length( $description ) > 0 ) {
    return $description;
  }
  return;
}


# Parse the GN fields of the current record and produce a structured
# list of names of genes coding the protein sequence in question.
sub _get_gene_names {
  my ( $self ) = @_;

  my $gn_fields = $self->{'record'}->{'GN'};
  if ( ! defined $gn_fields ) {
    return [];
  }

  my @gn_text_entries;

  # Gene-name entries can span multiple GN names but we cannot simply
  # concatenate them all because there can be multiple gene names per
  # record. In order to merge data correctly we must look for the
  # name-separator line.
  my $current_entry = q{};
  foreach my $line ( @{ $gn_fields } ) {
    # This is what the separator line looks like
    if ( $line eq 'and' ) {
      push @gn_text_entries, $current_entry;
      $current_entry = q{};
    }
    else {
      # No need for extra spaces here
      $current_entry .= $line;
    }
  }
  # Make sure the last entry makes it in as well
  push @gn_text_entries, $current_entry;

  my $gene_names = [];
  foreach my $entry ( @gn_text_entries ) {
    my $parsed_entry = {};

    my @entry_captures = ( $entry =~ m{
                                        \s*
                                        ( [^=]+ )
                                        \s*
                                        =
                                        \s*
                                        ( [^;]+ )
                                        \s*
                                        ;
                                    }gmsx );

    while ( my ( $key, $value ) = splice( @entry_captures, 0, 2 ) ) {
      # Strip possible evidence codes - again, that they can appear
      # here isn't mentioned in the User Manual but does appear in the
      # data. Leave the whitespace alone, splitting should take care
      # of it.
      $value =~ s{ $QR_OX_DE_EVIDENCE_CODE_LIST }{}gmsx;

      my @split_value = split( $QR_SPLIT_COMMA_WITH_WHITESPACE, $value );

      $parsed_entry->{$key}
        = ( $key eq 'Name' ) ? $value : \@split_value;
    }

    # UniProt-KB User Manual states a "Synonyms" token can only be
    # present if there is a "Name" token.
    if ( ( exists $parsed_entry->{'Synonyms'} )
         && ( ! exists $parsed_entry->{'Name'} ) ) {
      confess "Malformed input: found 'Synonyms' but no 'Name' in:\n\t$entry";
    }

    push @{ $gene_names }, $parsed_entry;
  }

  return $gene_names;
}


# Obtain quality information for the current record. This consists of
# two parts: status (i.e. whether the entry has been reviewed or not)
# from the ID line and evidence level from the PE line.
sub _get_quality {
  my ( $self ) = @_;

  my $entry_status;
  my $evidence_level;

  # There is at most one ID line
  my $id_line = $self->{'record'}->{'ID'}->[0];
  if ( defined $id_line ) {
    ( $entry_status )
      = ( $id_line =~ m{
                         \A

                         # UniProt name
                         [0-9A-Z_]+

                         \s+

                         ( $QR_ID_STATUS_FIELD )
                         \s*
                         ;
                     }msx );
    if ( ! defined $entry_status ) {
      confess "Invalid entry status in:\n\t$id_line";
    }
  }

  # Likewise, there is at most one PE line
  my $pe_line = $self->{'record'}->{'PE'}->[0];
  if ( defined $pe_line ) {
    ( $evidence_level )
      = ( $pe_line =~ m{
                         \A

                         ( [1-5] )
                         \s*
                         :
                     }msx );
    if ( ! defined $evidence_level ) {
      confess "Invalid protein evidence level in:\n\t$pe_line";
    }
  }

  return {
          'status'         => $entry_status,
          'evidence_level' => $evidence_level,
        };
}


# Parse the sequence-metadata ('SQ') and sequence (' ') fields of the
# current record, attempt to extract its type, and transform its
# sequence into a continuous string i.e. without decorative whitespace
# or mid-sequence counts.
sub _get_sequence {
  my ( $self ) = @_;

  my $seq_object = {
    'type' => $SEQUENCE_TYPE_UNDEFINED,
    'seq'  => undef,
  };

  my $metadata_line = $self->{'record'}->{'SQ'}->[0];
  if ( defined $metadata_line ) {

    # Different input types use different capitalisation here so just
    # do a case-insensitive match. At least so far they have all
    # followed the same basic syntax.
    my ( $unit ) = ( $metadata_line =~ m{
                                          \A
                                          sequence
                                          \s+
                                          \d+
                                          \s+
                                          ( \w+ )
                                          ;
                                      }imsx );
    if ( defined $unit ) {
      $seq_object->{'type'} = $sequence_type_for_unit{ uc( $unit ) };
    }
  }

  my $sequence_fields = $self->{'record'}->{ q{  } };
  if ( defined $sequence_fields ) {

    my $sequence = q{};

  SEQUENCE_LINE:
    while ( my $line = shift @{ $sequence_fields } ) {
      # Strip mid-sequence counts (and surrounding whitespace) from
      # the end of each line
      $line =~ s{
                  \s+ \d+
                  \s* \z
              }{}msx;

      # Get rid of all the remaining whitespace. Apart from the blanks
      # preceding mid-sequence counts, which have already been taken
      # care of, we do not expect to see more than whitespace at a
      # time - so instead of trying to match as long a string of them
      # as possible in order to minimise the number of independent
      # substitutions, we always match on one in order to avoid
      # unnecessary backtracking.
      $line =~ s{ \s }{}gmsx;

      $sequence .= $line;
    }

    $seq_object->{'seq'} = $sequence;
  }

  return $seq_object;
}



# Parse the OX field of the current record and produce a list of
# database_qualifier/taxon_code pairs.
sub _get_taxon_codes {
  my ( $self ) = @_;

  # There is at most one OX line
  my $ox_line = $self->{'record'}->{'OX'}->[0];
  if ( ! defined $ox_line ) {
    return [];
  }

  # On the one hand, according to UniProt-KB User Manual from October
  # 2018 there should only be a single taxon code per OX line and the
  # current SwissProt data file doesn't contain any records which
  # violate this. On the other hand we haven't checked this in the
  # TrEMBL file because of its size and the old UniProtParser did have
  # support for synonyms (albeit using different syntax -
  # "db=id1, id2, ...;" rather than "db1=id1; db2=id2; ...").
  # The code below assumes there might be multiple entries present, if
  # you want to force it to only ever look for one simply drop the /g
  # modifier from the regex match.

  my @ox_captures
    = ( $ox_line =~ m{
                       $QR_OX_TAXON_DB_ENTRY
                       \s*

                       # Optional things (e.g. evidence codes) we do
                       # not need now but must parse on the off chance
                       # someone decides e.g. that curly braces
                       # constitute quotes so it is okay to have a
                       # semicolon between them
                       (?:
                         $QR_OX_DE_EVIDENCE_CODE_LIST
                         | [^;]+
                       )?

                       # End of record
                       ;
                   }gmsx );

  my @extracted_taxon_codes;

  # Reminder: if you change the number of capture groups above
  # remember to adapt both the variable list *and* the third argument
  # to splice().
 TAXON_ENTRY:
  while ( my ( $db_qualifier, $taxon_code ) = splice( @ox_captures, 0, 2 ) ) {

    # FIXME: as it is, i.e. with /g above, it is not possible for
    # either $db_qualifier or $taxon_code to be undefined - if the
    # match fails for any OX entry Perl will simply not push its
    # captures to @ox_captures. We need something better... Or perhaps
    # we could just drop support for multiple taxon codes and stop
    # messing around with /g.

    if ( ( ! defined $db_qualifier )
         || ( ! exists $taxonomy_ids_from_taxdb_codes{$db_qualifier} ) ) {
      # Abort on malformed or new database qualifiers
      confess "Cannot use taxon-DB qualifier '${db_qualifier}'";
    }
    elsif ( ! $taxonomy_ids_from_taxdb_codes{$db_qualifier} ) {
      # Known but of no interest. Ignore it.
      next TAXON_ENTRY;
    }

    if ( ! defined $taxon_code ) {
      confess "Failed to extract taxon code from:\n\t${ox_line}";
    }

    push @extracted_taxon_codes, {
                                  'db_qualifier' => $db_qualifier,
                                  'taxon_code'     => $taxon_code,
                                };

  }

  return \@extracted_taxon_codes;
}


# Check if the species indicated in the current record matches the one
# for which we are running the parser. Specific implementation might
# vary depending on the input data.
sub _record_species_matches {
  my ( $self ) = @_;

  return $self->_taxon_codes_match_species_id();
}


# Translate extracted taxon codes into Ensembl taxonomy IDs, then
# check if any of them correspond the current species ID.
sub _taxon_codes_match_species_id {
  my ( $self ) = @_;

  my $taxon_codes = $self->_get_taxon_codes();

  my @taxonomy_ids;
  foreach my $taxon ( @{ $taxon_codes } ) {
    my $code_mapper
      = $taxonomy_ids_from_taxdb_codes{ $taxon->{'db_qualifier'} };
    push @taxonomy_ids, $code_mapper->( $taxon->{'taxon_code'} );
  }

  foreach my $taxonomy_id ( @taxonomy_ids ) {
    if ( $taxonomy_id == $self->{'species_id'} ) {
      return 1;
    }
  }

  return 0;
}


sub _record_has_all_needed_fields {
  my ( $self ) = @_;

  my %prefixes_of_interest = %{ $self->{'prefixes_of_interest'} };

  # Only check mandatory fields
  my @needed_fields
    = grep { $prefixes_of_interest{$_} } keys %prefixes_of_interest;

  my $has_all
    = List::Util::all { exists $self->{'record'}->{$_} } @needed_fields;
  if ( ! $has_all ) {
    confess 'One or more required fields missing in record';
  }

  return;
}


1;
