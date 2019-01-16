
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

Bio::EnsEMBL::Xref::Parser::JGI_Parser

=head1 DESCRIPTION

Parser for JGI protein files with gene description, FASTA format.
This is the base class the provides most functionality, subclasses
just implement the method to set sequence type

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::JGI_Parser->new(
    source_id  => 70,
    species_id => 7719,
    files      => [ "ciona.prot.fasta.gz" ],
    xref_dba   => $xref_dba
  );

  $parser->run();

=cut

package Bio::EnsEMBL::Xref::Parser::JGI_Parser;

use strict;
use warnings;

use Carp;

use parent qw( Bio::EnsEMBL::Xref::Parser );

=head2 run

  Arg []     : None
  Description: Parse input file, extract xrefs and add them to the xref DB
  Return type: Int; 0 upon success 
  Caller     : An analysis step in the Xref pipeline

=cut

sub run {

  my ($self) = @_;

  my $source_id  = $self->{source_id};
  my $species_id = $self->{species_id};
  my $files      = $self->{files};
  my $verbose    = $self->{verbose} // 0;
  my $xref_dba   = $self->{xref_dba};

  unless ( defined $source_id and
           defined $species_id and
           defined $files )
  {
    confess "Need to pass source_id, species_id and files";
  }

  my $file = shift @{$files};

  # the source name defines how to parse the header
  my $source_name =
    $xref_dba->get_source_name_for_source_id($source_id);

  my $file_io = $xref_dba->get_filehandle($file);
  confess "Could not open $file\n" unless defined $file_io;

  IO::Handle->input_record_separator("\n>");

  my @xrefs;

  while ( my $record = $file_io->getline() ) {

    next if ( $record =~ /^File:/ ); # skip header

    my ( $accession, $sequence ) = $record =~ /^>?(.+?)\n([^>]*)/xs or
      confess "Can't parse FASTA entry: $record";

    # split header in different ways according to source name
    if ( $source_name =~ m/cint_aniseed_.*v1|cint_jgi_v1/x ) {
      # header format is  >ci0100146277
      # JGI 1.0
      # we want accession 146277 from above
      $accession =~ s/\w{6}//x;

    } else {
      confess
        "The source-name specified ($source_name) is not matching\n" .
        "the different cases specified in JGI_Parser.pm - please\n" .
        "edit the parser \n";
    }

    # make sequence into one long string
    $sequence =~ s/\n//xg;

    # build the xref object and store it
    push @xrefs,
      { ACCESSION     => $accession,
        SEQUENCE      => $sequence,
        SOURCE_ID     => $source_id,
        SPECIES_ID    => $species_id,
        SEQUENCE_TYPE => $self->get_sequence_type() };

  } ## end while ( my $record = $file_io...)

  $file_io->close();
  IO::Handle->input_record_separator("\n");

  $xref_dba->upload_xref_object_graphs( \@xrefs );
  print scalar(@xrefs) . " JGI_ xrefs succesfully parsed\n" if $verbose;

  return 0;    # success
} ## end sub run

=head2 get_sequence_type

  Arg []     : None
  Description: Abstract method to guard against calling it in the base class

=cut

sub get_sequence_type {

  my $self = shift;

  confess
"Abstract method, please instantiate derived class providing an implementation";
} ## end sub get_sequence_type

1;
