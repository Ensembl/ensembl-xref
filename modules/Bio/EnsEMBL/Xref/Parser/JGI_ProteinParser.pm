
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

=head1 NAME

Bio::EnsEMBL::Xref::Parser::JGI_ProteinParser

=head1 DESCRIPTION

Parser for JGI protein files with gene description, FASTA format.
This inherits from the base parser to provide the sequence type

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::JGI_ProteinParser->new(
    source_id  => 70,
    species_id => 7719,
    files      => [ "ciona.prot.fasta.gz" ],
    xref_dba   => $xref_dba
  );

  $parser->run();

=cut

package Bio::EnsEMBL::Xref::Parser::JGI_ProteinParser;

use strict;
use warnings;

use parent qw( Bio::EnsEMBL::Xref::Parser::JGI_Parser );

=head2 get_sequence_type

  Arg []     : None
  Description: Return the type of sequences handled by this specialised parser
  Return type: Scalar; a string representing the sequence type
  Caller     : Parent class run method

=cut

sub get_sequence_type {

  my $self = shift;

  return 'peptide';
}

1;
