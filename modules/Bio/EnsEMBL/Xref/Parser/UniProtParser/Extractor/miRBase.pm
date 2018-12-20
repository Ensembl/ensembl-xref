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



package Bio::EnsEMBL::Xref::Parser::UniProtParser::Extractor::miRBase;

use strict;
use warnings;

use Carp;

use parent qw( Bio::EnsEMBL::Xref::Parser::UniProtParser::Extractor );


# In case of miRBase data, we want to use the generic miRNA identifier
# as the description. We extract it from the ID field because while it
# is also present in the DE field, we want the variant beginning with
# three-letter species abbreviation and DE uses human-readable species
# names.
sub _get_description {  ## no critic (ProhibitUnusedPrivateSubroutines)
  my ( $self ) = @_;

  # Only one ID line per record in miRBase files
  my $id_line = $self->{'record'}->{'ID'}->[0];
  if ( ! defined $id_line ) {
    # No ID -> no description. Should never happen in practice
    # because the parser should mark ID as a mandatory field.
    return;
  }

  my ( $mirna_id )
    = ( $id_line =~ m{
                       \A
                       ( \S+ )
                       \s+
                   }msx );
  if ( ! defined $mirna_id ) {
    # A problem with the pattern? Make sure it is noticed.
    confess "Failed to extract miRNA ID from:\n\t$id_line";
  }

  return $mirna_id;
}


# Even if there is anything counting as quality information in miRBase
# records, we do not use it.
sub _get_quality {  ## no critic (ProhibitUnusedPrivateSubroutines)
  my ( $self ) = @_;

  return {};
}


# miRBase-specific implementation of the species-match check: use the
# full species name
sub _record_species_matches {  ## no critic (ProhibitUnusedPrivateSubroutines)
  my ( $self ) = @_;

  return $self->_species_name_matches();
}


# Extract the full, human-readable species name from the DE line of
# the miRBase record and see if it matches the one provided by the
# xref pipeline.
sub _species_name_matches {
  my ( $self ) = @_;

  # Only one DE line per record in miRBase files
  my $de_line = $self->{'record'}->{'DE'}->[0];
  if ( ! defined $de_line ) {
    # No DE -> no species -> no match. Should never happen in practice
    # because the parser should mark DE as a mandatory field.
    return 0;
  }

  my ( $record_species )
    = ( $de_line =~ m{
                       \A
                       ( [^\n]+ )
                       \s+
                       \S+  # miRNA identifier
                       \s+
                       stem
                       [- ] # both are present in the data. Consistency!
                       loop
                   }msx );
  if ( ! defined $record_species ) {
    # A problem with the pattern? Make sure it is noticed.
    confess "Failed to extract species from:\n\t$de_line";
  }

  # Make extracted species name match Ensembl formatting
  $record_species = lc( $record_species );
  $record_species =~ s{[ ]}{_}msx;

  return ( $record_species eq $self->{'species_name'} );
}


1;
