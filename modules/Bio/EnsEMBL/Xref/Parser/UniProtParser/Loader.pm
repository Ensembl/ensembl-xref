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



package Bio::EnsEMBL::Xref::Parser::UniProtParser::Loader;

use strict;
use warnings;

use Carp;


=head2 new

  Arg [1]    : HashRef arguments for the constructor:
                - batch_size
                   - how many UniProt-KB records to accumulate in
                     the send buffer before pushing them to the
                     database. Defaults to 1 i.e. push every record
                     individually.
                - checkpoint_seconds
                   - how often, in seconds, to push the contents of
                     the send buffer to the database regardless of how
                     many entries it contains. Defaults to 0
                     i.e. disable this feature.
                - xref_dba
                   - DBAdaptor object passed from the xref pipeline
  Description: Constructor.
  Return type: Loader object
  Exceptions : none
  Caller     : UniProtParser::run()
  Status     : Stable

=cut

sub new {
  my ( $proto, $arg_ref ) = @_;

  my $self = {
    'batch_size'         => $arg_ref->{'batch_size'}         // 1,
    'checkpoint_seconds' => $arg_ref->{'checkpoint_seconds'} // 0,
    'xref_dba'           => $arg_ref->{'xref_dba'},
    'last_flush'         => time(),
  };
  my $class = ref $proto || $proto;
  bless $self, $class;
  $self->_clear_send_buffer();

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

  Description: Wrap-up routine. At present, pushes to the database any
               data still residing in the send buffer.
  Return type: none
  Exceptions : none
  Caller     : destructor, UniProtParser::run()
  Status     : Stable

=cut

sub finish {
  my ( $self ) = @_;

  $self->flush();

  return;
}


=head2 flush

  Description: Unconditionally push the contents of the send buffer to
               the database, then empty the buffer. Meant primarily
               for internal use but could in principle be called
               explicitly as well.
  Return type: none
  Exceptions : throws upon errors having been returned by
               $xref_dba->upload_xref_object_graphs()
  Caller     : finish()
  Status     : Stable

=cut

sub flush {
  my ( $self ) = @_;

  my $xref_dba = $self->{'xref_dba'};

  if ( $self->{'send_backlog'} > 0 ) {
    $xref_dba->upload_xref_object_graphs( $self->{'send_buffer'} );

    # FIXME: we might want to run this even in the event of a failure,
    # which would require capturing exceptions thrown above
    $self->_clear_send_buffer();
  }

  return;
}


=head2 load

  Arg [1]    : HashRef transformed_data UniProt-KB prepared by a
               transformer to be ready to be inserted into the database
  Description: Add transformed_data to the send buffer, then flush the
               buffer if either it has reached the specified size or
               the requested amount of time has passed since the last
               flush.
  Return type: none
  Exceptions : none
  Caller     : UniProtParser::run()
  Status     : Stable

=cut

sub load {
  my ( $self, $transformed_data ) = @_;

  if ( ! defined $transformed_data ) {
    return;
  }

  $self->_add_to_send_buffer( $transformed_data );

  # FIXME: time checkpointing isn't particularly useful as it is, the
  # final xref output by the transformer will just sit in the queue ad
  # infinitum regardless of how much more of the file there is
  # left. Will have to create a proper timer for this.
  my $current_time = time();
  if ( ( $self->{'send_backlog'} >= $self->{'batch_size'} )
       || ( $self->{'checkpoint_seconds'} > 0  )
       && ( $current_time - $self->{'last_flush'} > $self->{'checkpoint_seconds'} ) ) {
    $self->flush();
    $self->{'last_flush'} = $current_time;
  }

  return;
}


=head2 prepare_source_for_dependent_xrefs

  Arg [1]    : integer source_id Ensembl ID of the given xref source
  Description: A thin wrapper around
               BaseAdaptor::get_dependent_mappings(), which populates
               a $xref_dba-internal map with existing dependent_xref
               links in order to prevent insertion of duplicates in
               the event of a parser being re-run on the same input
               data. This method has been left to be called explicitly
               rather than being embedded in BaseAdaptor
               initialisation because generating such a map for all
               sources would likely be a costly operation in both
               memory and processing time.
  Return type: none
  Exceptions : none
  Caller     : UniProtParser::run()
  Status     : Stable

=cut

sub prepare_source_for_dependent_xrefs {
  my ( $self, $source_id ) = @_;

  my $xref_dba = $self->{'xref_dba'};

  $xref_dba->get_dependent_mappings( $source_id );

  return;
}


# Add an entry to the send buffer and increment the related counter
# (which we use in order not to repeatedly calculate array size, which
# is a relatively costly operation) by one.
sub _add_to_send_buffer {
  my ( $self, $entry ) = @_;

  push @{ $self->{'send_buffer'} }, $entry;
  $self->{'send_backlog'}++;

  return;
}


# Replace the send buffer with an empty array and zero out the size
# counter.
sub _clear_send_buffer {
  my ( $self ) = @_;

  $self->{'send_buffer'}  = [];
  $self->{'send_backlog'} = 0;

  return;
}


1;
