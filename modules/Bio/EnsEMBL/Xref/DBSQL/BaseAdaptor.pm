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

Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor - A db adaptor for the xref database

=cut

=head1 DESCRIPTION

This is the base adaptor for loading xrefs into the species xref database during
the xref parsing stage in production. The aim of the module is to reduce
redundancy of code and keep the SQL to a constrained set of functions

=cut

package Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;

use strict;
use warnings;

use Carp;

use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Xref::FetchFiles;
use Getopt::Long;
use IO::Uncompress::AnyUncompress;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $base_dir = File::Spec->curdir();

my %xref_dependent_mapped;

my $verbose;

=head2 new
  Arg [1]    : proto
  Arg [2]    : arguments : { host => string, dbname => string, user => string, pass => string, port => int}
  Description: Initialisation class for the dbi connection
  Return type: self
  Caller     : internal

=cut

sub new {
  my ( $proto, %args ) = @_;

  my $class = ref $proto || $proto;
  my $self = bless {}, $class;

  $self->dbc( Bio::EnsEMBL::DBSQL::DBConnection->new(
    -HOST   => $args{host},
    -DBNAME => $args{dbname},
    -USER   => $args{user},
    -PASS   => $args{pass} || '',
    -PORT => $args{port} || '3306'
  ) );
  $self->verbose( $args{verbose} // 0 );

  return $self;
} ## end sub new


=head2 dbc
  Description: Getter/Setter for the dbc object
  Return type: db_connection
  Caller     : internal

=cut

sub dbc {
  my ( $self, $arg ) = @_;
  ( defined $arg ) && ( $self->{_dbc} = $arg );
  return $self->{_dbc};
} ## end sub dbc


=head2 dbi
  Description: Getter/Setter for the dbi object
  Return type: DBI database handle
  Caller     : internal

=cut

sub dbi {
  my $self = shift;

  return $self->dbc->db_handle;
} ## end sub dbi


=head2 dba
  Description: Getter for the dba object
  Return type: Bio::EnsEMBL::DBSQL::DBAdaptor instance
  Caller     : internal
=cut

sub dba {
    my $self = shift;
    return Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbconn => $self->dbc, -species => $self->species);
}


=head2 species
  Arg [1]    : (optional) string $arg
               The new value of the species
  Example    : $species = $dba->species()
  Description: Getter/Setter for the current species
  Returntype : string
  Exceptions : none
  Caller     : new
=cut

sub species{
  my ($self, $arg) = @_;

  if ( defined $arg ) {
    $self->{_species} = $arg;
  }
  return $self->{_species};
}


=head2 get_filehandle
  Arg [1]    : file name
  Description: Given a file name, returns a IO::Handle object.  Supports most common
               compression formats, e.g. zip, gzip, bzip2, lzma, xz.  If the given
               file name doesn't correspond to an existing file, the routine will
               try to add '.gz' to the file name or to remove any .'Z' or '.gz' and
               try again.  Throws on failure.
  Return type: filehandle
  Exceptions : confesses if not found
  Caller     : internal

=cut

sub get_filehandle {
  my ( $self, $file_name ) = @_;

  my $io = undef;

  if ( !( defined $file_name ) || $file_name eq '' ) {
    confess "No file name";
  }
  my $alt_file_name = $file_name;
  $alt_file_name =~ s/\.(gz|Z)$//x;

  if ( $alt_file_name eq $file_name ) {
    $alt_file_name .= '.gz';
  }

  if ( !-e $file_name ) {
    print "File $file_name does not exist, " .
          "will try $alt_file_name\n";
    $file_name = $alt_file_name;
  }

  # 'Transparent' lets IO::Uncompress modules read uncompressed input.
  # It should be on by default but set it just in case.
  $io = IO::Uncompress::AnyUncompress->new($file_name,
                                           'Transparent' => 1 )
    || confess("Can not open file '$file_name'");

  if ($verbose) {
    print "Reading from '$file_name'...\n";
  }

  return $io;
} ## end sub get_filehandle


=head2 get_source_id_for_source_name
  Arg [1]    : source name
  Arg [2]    : priority description
  Description: Gets the source ID for a given source name
  Return type: integer
  Exceptions : confesses if not found
  Caller     : internal

=cut

sub get_source_id_for_source_name {
  my ( $self, $source_name, $priority_desc) = @_;

  if ( !defined $source_name ) {
    confess 'source_name undefined';
  }

  my $low_name = lc $source_name;
  my $sql = 'SELECT source_id FROM source WHERE LOWER(name)=?';

  my @sql_params = ( $low_name );
  if ( defined $priority_desc ) {
    $sql .= " AND LOWER(priority_description)=?";
    push @sql_params, lc $priority_desc;

    $source_name .= " ($priority_desc)";
  }
  my $sth = $self->dbi->prepare_cached($sql);
  $sth->execute(@sql_params);

  my @row = $sth->fetchrow_array();
  my $source_id;
  if (@row) {
    $source_id = $row[0];
  }
  else {
    my $msg = "No source_id for source_name='${source_name}'";
    if ( defined $priority_desc ) {
      $msg .= "priority_desc='${priority_desc}'";
    }
    confess $msg;
  }
  $sth->finish();

  return $source_id;
} ## end sub get_source_id_for_source_name


=head2 get_source_ids_for_source_name_pattern
  Arg [1]    : source name
  Description: Get a set of source IDs matching a source name pattern
               Adds % to each end of the source name and doe a like query to find
               all the matching source names source_ids.
  Return type: array
  Caller     : internal

=cut

sub get_source_ids_for_source_name_pattern {

  my ( $self, $source_name) = @_;

  if ( !defined $source_name ) {
    confess 'source_name undefined';
  }

  my $sth = $self->dbi->prepare_cached(
    'SELECT source_id FROM source WHERE UPPER(name) LIKE ?'
  );

  my $big_name = uc $source_name;
  $sth->execute("%${big_name}%");

  my @sources;
  while ( my @row = $sth->fetchrow_array() ) {
    push @sources, $row[0];
  }

  return @sources;
} ## end sub get_source_ids_for_source_name_pattern


=head2 get_source_name_for_source_id
  Arg [1]    : source id
  Description: Gets the source name for a given source ID
  Return type: string
  Caller     : internal

=cut

sub get_source_name_for_source_id {
  my ( $self, $source_id ) = @_;

  if ( !defined $source_id ) {
    confess 'source_id undefined';
  }

  my $source_name;

  my $sql = 'SELECT name FROM source WHERE source_id= ?';
  my $sth = $self->dbi->prepare_cached($sql);
  $sth->execute($source_id);
  my @row = $sth->fetchrow_array();
  if (@row) {
    $sth->finish();
    $source_name = $row[0];
  }
  else {
    print "There is no entity with source-id  $source_id  in the source-table of the \n";
    print "xref-database. The source-id and the name of the source-id is hard-coded in populate_metadata.sql\n";
    print "and in the parser\n";
    print "Couldn't get source name for source ID $source_id\n";
    $source_name = '-1';
  }

  return $source_name;
} ## end sub get_source_name_for_source_id


=head2 get_valid_xrefs_for_dependencies
  Arg [1]    : dependent name
  Arg [2]    : reverse ordered source name list
  Description: Get a hash to go from accession of a dependent xref to master_xref_id
               for all of source names given
  Return type: Hashref
  Caller     : internal

=cut

sub get_valid_xrefs_for_dependencies {
  my ( $self, $dependent_name, @reverse_ordered_source_list ) = @_;

  my %dependent_2_xref;

  my $sql = 'SELECT source_id FROM source WHERE LOWER(name) =?';
  my $sth = $self->dbi->prepare_cached($sql);
  my @dependent_sources;
  $sth->execute( lc $dependent_name );
  while ( my @row = $sth->fetchrow_array() ) {
    push @dependent_sources, $row[0];
  }

  my @sources;
  foreach my $name (@reverse_ordered_source_list) {
    $sth->execute( lc $name );
    while ( my @row = $sth->fetchrow_array() ) {
      push @sources, $row[0];
    }
  }

  my $dep_sql = (<<'DSS');
    SELECT d.master_xref_id, x2.accession
    FROM dependent_xref d, xref x1, xref x2
    WHERE x1.xref_id = d.master_xref_id AND
          x1.source_id = ? AND
          x2.xref_id = d.dependent_xref_id AND
          x2.source_id = ?
DSS

  $sth = $self->dbi->prepare_cached($dep_sql);
  foreach my $d (@dependent_sources) {
    foreach my $s (@sources) {
      $sth->execute( $s, $d );
      while ( my @row = $sth->fetchrow_array() ) {
        $dependent_2_xref{ $row[1] } = $row[0];
      }
    }
  }

  return \%dependent_2_xref;
} ## end sub get_valid_xrefs_for_dependencies


=head2 get_valid_xrefs_for_direct_xrefs
  Arg [1]    : direct name
  Arg [2]    : separator
  Description: Get a hash to go from accession of a direct xref to master_xref_id
               for all of source names given
  Return type: Hashref
  Caller     : internal

=cut

sub get_valid_xrefs_for_direct_xrefs {
  my ( $self, $direct_name, $separator ) = @_;

  my %direct_2_xref;

  my $sql = 'SELECT source_id FROM source WHERE name LIKE ?';
  my $sth = $self->dbi->prepare($sql);
  my @direct_sources;
  $sth->execute("${direct_name}%");
  while ( my @row = $sth->fetchrow_array() ) {
    push @direct_sources, $row[0];
  }

  my $gene_sql = <<"GDS";
    SELECT d.general_xref_id,
           d.ensembl_stable_id,
           'Gene',
           d.linkage_xref,
           x1.accession
    FROM gene_direct_xref d, xref x1
    WHERE x1.xref_id = d.general_xref_id AND
          x1.source_id=?
GDS

  my $transcript_sql = <<"TDS";
    SELECT d.general_xref_id,
           d.ensembl_stable_id,
           'Transcript',
           d.linkage_xref,
           x1.accession
    FROM transcript_direct_xref d, xref x1
    WHERE x1.xref_id = d.general_xref_id AND
          x1.source_id=?
TDS

  my $translation_sql = <<"PDS";
    SELECT d.general_xref_id,
           d.ensembl_stable_id,
           'TYPE',
           d.linkage_xref,
           x1.accession
    FROM translation_direct_xref d, xref x1
    WHERE x1.xref_id = d.general_xref_id AND
          x1.source_id=?
PDS

  my %sql_hash = (
    Gene => $gene_sql,
    Transcript => $transcript_sql,
    Translation => $translation_sql
  );

  my @sth;
  my $i = 0;
  foreach my $type (qw(Gene Transcript Translation)) {
    $sth[ $i++ ] = $self->dbi->prepare_cached( $sql_hash{$type} );
  }

  foreach my $d (@direct_sources) {
    for my $ii ( 0 .. 2 ) {
      $sth[$ii]->execute($d);
      while ( my ( $gen_xref_id, $stable_id, $type, $link, $acc ) =
              $sth[$ii]->fetchrow_array() )
      {
        if ( !defined $link ) {
          $link = '';
        }
        $direct_2_xref{$acc} =
          $gen_xref_id . $separator .
          $stable_id . $separator . $type . $separator . $link;
      }
    }
  }

  return \%direct_2_xref;
} ## end sub get_valid_xrefs_for_direct_xrefs


=head2 label_to_acc
  Arg [1]    : source name
  Arg [2]    : species id
  Description: Get a hash of label to acc for a particular source name and
               species_id
  Return type: Hashref
  Caller     : internal

=cut

sub label_to_acc {

  my ( $self, $source_name, $species_id ) = @_;

  # First cache synonyms so we can quickly add them later
  my %synonyms;
  my $syn_sth = $self->dbi->prepare_cached('SELECT xref_id, synonym FROM synonym');
  $syn_sth->execute();

  my ( $xref_id, $synonym );
  $syn_sth->bind_columns( \$xref_id, \$synonym );
  while ( $syn_sth->fetch() ) {
    push @{ $synonyms{$xref_id} }, $synonym;
  }

  my %valid_codes;
  my @sources;

  my $big_name = uc $source_name;
  my $sql =
    'SELECT source_id FROM source WHERE UPPER(name) LIKE ?';
  my $sth = $self->dbi->prepare_cached($sql);
  $sth->execute("${big_name}%");
  while ( my @row = $sth->fetchrow_array() ) {
    push @sources, $row[0];
  }

  foreach my $source (@sources) {
    $sql =
      'SELECT label, xref_id FROM xref WHERE species_id = ? AND source_id = ?';
    $sth = $self->dbi->prepare_cached($sql);
    $sth->execute( $species_id, $source );
    while ( my @row = $sth->fetchrow_array() ) {
      $valid_codes{ $row[0] } = $row[1];
      # add any synonyms for this xref as well
      foreach my $syn ( @{ $synonyms{ $row[1] } } ) {
        $valid_codes{$syn} = $row[1];
      }
    }
  }

  return \%valid_codes;
} ## end sub label_to_acc


=head2 get_valid_codes
  Arg [1]    : source name
  Arg [2]    : species id
  Description: Hash of accession to array of xrefs.
               This is an array becouse more than one entry can exist. i.e. for
               uniprot and refseq we have direct and sequence match sets and we
               need to give both.
  Return type: Hashref
  Caller     : internal

=cut

sub get_valid_codes {

  my ( $self, $source_name, $species_id ) = @_;

  my %valid_codes;
  my @sources;

  my $big_name = uc $source_name;
  my $sql =
    'SELECT source_id FROM source WHERE UPPER(name) LIKE ?';
  my $sth = $self->dbi->prepare_cached($sql);
  $sth->execute("%$big_name%");
  while ( my @row = $sth->fetchrow_array() ) {
    push @sources, $row[0];
  }

  foreach my $source (@sources) {
    $sql = 'SELECT accession, xref_id FROM xref WHERE species_id = ? AND source_id = ?';
    $sth = $self->dbi->prepare_cached($sql);
    $sth->execute( $species_id, $source );
    while ( my @row = $sth->fetchrow_array() ) {
      push @{ $valid_codes{ $row[0] } }, $row[1];
    }
  }

  return \%valid_codes;
} ## end sub get_valid_codes


=head2 upload_xref_object_graphs
  Arg [1]    : Array of xrefs
  Description: Upload xrefs to the database and associated data to the database
  Return type:
  Exceptions : confess if accession of source ID are not provided
  Caller     : internal

=cut

sub upload_xref_object_graphs {
  my ( $self, $rxrefs ) = @_;

  if ( !defined $rxrefs || !( scalar @{$rxrefs} ) ) {
    confess "Please give me some xrefs to load";
  }

  foreach my $xref ( @{$rxrefs} ) {
    if ( !( defined $xref->{ACCESSION} ) ) {
      confess "Your xref does not have an accession-number\n";
    }
    if ( !( defined $xref->{SOURCE_ID} ) ) {
      confess "Your xref: $xref->{ACCESSION} does not have a source-id\n";
    }

    # Create entry in xref table and note ID
    my $xref_id = $self->add_xref( {
      "acc"          => $xref->{ACCESSION},
      "version"      => $xref->{VERSION} // 0,
      "label"        => $xref->{LABEL}   // $xref->{ACCESSION},
      "desc"         => $xref->{DESCRIPTION},
      "source_id"    => $xref->{SOURCE_ID},
      "species_id"   => $xref->{SPECIES_ID},
      "info_type"    => $xref->{INFO_TYPE},
      "update_label" => 1, "update_desc" => 1 } );

    # If there are any direct_xrefs, add these to the relevant tables
    $self->add_multiple_direct_xrefs( @{ $xref->{DIRECT_XREFS} } );

    # create entry in primary_xref table with sequence; if this is a "cumulative"
    # entry it may already exist, and require an UPDATE rather than an INSERT
    if ( defined $xref->{SEQUENCE} ) {
      if ( $self->primary_xref_id_exists($xref_id) ) {
        $self->_update_primary_xref_sequence( $xref->{SEQUENCE}, $xref_id );

      }
      else {
        $self->_add_primary_xref(
          $xref_id, $xref->{SEQUENCE}, $xref->{SEQUENCE_TYPE}, $xref->{STATUS}
        );
      }
    }

    # if there are synonyms, add entries in the synonym table
    $self->add_multiple_synonyms( $xref_id, $xref->{SYNONYMS} );

    # if there are dependent xrefs, add xrefs and dependent xrefs for them
    $self->add_multiple_dependent_xrefs( $xref_id, $xref->{DEPENDENT_XREFS} );

    # Add the pair data. refseq dna/pep pairs usually
    if ( defined $xref->{PAIR} ) {
      $self->_add_pair( $xref->{SOURCE_ID}, $xref->{ACCESSION}, $xref->{PAIR} );
    }

  } # foreach xref

  return;
} ## end sub upload_xref_object_graphs


=head2 upload_direct_xrefs
  Arg [1]    : Array of direct xrefs
  Description: Add direct xref to the table XXX_direct_xref. (XXX -E<gt> Gene,
               Transcript or Translation. Xref has to exist already, this module
               just adds to the direct_xref table.
               $direct_xref is a reference to an array of hash objects.
  Return type:
  Caller     : internal

=cut

sub upload_direct_xrefs {
  my ( $self, $direct_xref ) = @_;
  for my $dr ( @{$direct_xref} ) {

    # Find the xref_id for this accession and source
    my $general_xref_id =
      $self->get_xref( $dr->{ACCESSION}, $dr->{SOURCE_ID},
                       $dr->{SPECIES_ID});

    # If found add the direct xref else write error message
    if ($general_xref_id) {
      $self->add_direct_xref( $general_xref_id,
                              $dr->{STABLE_ID},
                              $dr->{ENSEMBL_TYPE},
                              $dr->{LINKAGE_XREF}
                            );
    }
    else {
      print {*STDERR} 'Problem Could not find accession ' .
        $dr->{ACCESSION} . ' for source ' .
        $dr->{SOURCE} . ' so not able to add direct xref to ' .
        $dr->{STABLE_ID} . "\n";
    }
  } ## end for my $dr ( @{$direct_xref...})
  return;
} ## end sub upload_direct_xrefs


=head2 add_meta_pair
  Arg [1]    : key
  Arg [2]    : value
  Description: Insert into the meta table the key and value.
  Return type:
  Caller     : internal

=cut

sub add_meta_pair {

  my ( $self, $key, $value ) = @_;

  my $sth = $self->dbi->prepare_cached(
    'INSERT INTO meta (meta_key, meta_value, date) VALUES (?, ?, NOW())' );
  $sth->execute( $key, $value );

  return;
} ## end sub add_meta_pair


=head2 get_xref_sources
  Description: Create a hash of all the source names for xrefs
  Return type: Hashref
  Caller     : internal

=cut

sub get_xref_sources {

  my $self = shift;

  my %sourcename_to_sourceid;

  my $sth = $self->dbi->prepare_cached('SELECT name, source_id FROM source');
  $sth->execute() or croak( $self->dbi->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $source_name = $row[0];
    my $source_id   = $row[1];
    $sourcename_to_sourceid{$source_name} = $source_id;
  }


  return %sourcename_to_sourceid;
} ## end sub get_xref_sources


=head2 species_id2taxonomy
  Description: Create and return a hash that that goes from species_id to taxonomy_id
  Return type: Hashref
  Caller     : internal

=cut

sub species_id2taxonomy {

  my $self = shift;

  my %species_id2taxonomy;

  my $sth =
    $self->dbi->prepare_cached('SELECT species_id, taxonomy_id FROM species');
  $sth->execute() or croak( $self->dbi->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $species_id  = $row[0];
    my $taxonomy_id = $row[1];
    if ( defined $species_id2taxonomy{$species_id} ) {
      push @{ $species_id2taxonomy{$species_id} }, $taxonomy_id;
    }
    else {
      $species_id2taxonomy{$species_id} = [$taxonomy_id];
    }
  }

  return %species_id2taxonomy;
} ## end sub species_id2taxonomy


=head2 species_id2name
  Description: Create and return a hash that that goes from species_id to species name
  Return type:
  Caller     : internal

=cut

sub species_id2name {
  my $self = shift;

  my %species_id2name;

  my $sth = $self->dbi->prepare_cached('SELECT species_id, name FROM species');
  $sth->execute() or croak( $self->dbi->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $species_id = $row[0];
    my $name       = $row[1];
    $species_id2name{$species_id} = [$name];
  }

  ##############################################
  # Also populate the hash with all the aliases.
  ##############################################
  $sth = $self->dbi->prepare_cached('SELECT species_id, aliases FROM species');
  $sth->execute() or croak( $self->dbi->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $species_id = $row[0];
    foreach my $name ( split /,\s*/xms, $row[1] ) {
      $species_id2name{$species_id} ||= [];
      push @{ $species_id2name{$species_id} }, $name;
    }
  }

  return %species_id2name;
} ## end sub species_id2name


=head2 get_xref_id
  Arg [1]    : xref entry
  Description: Get the xref internal id
               If there was an error, an xref with the same acc & source already
               exists. If so, find its ID, otherwise get ID of xref just inserted
  Return type: integer
  Caller     : internal

=cut

sub get_xref_id {
  my ( $self, $arg_ref ) = @_;
  my $acc = $arg_ref->{acc} ||
    confess 'Need an accession for get_xref_id';
  my $source = $arg_ref->{source_id} ||
    confess 'Need a source_id for get_xref_id';
  my $species = $arg_ref->{species_id} ||
    confess 'Need a species_id for get_xref_id';

  my $id = $self->get_xref( $acc, $source, $species );

  return $id;
} ## end sub get_xref_id


=head2 primary_xref_id_exists
  Arg [1]    : xref_id
  Description: If primary xref already exists for a partiuclar xref_id return 1
               else return 0
  Return type: integer
  Caller     : internal

=cut

sub primary_xref_id_exists {

  my ( $self, $xref_id ) = @_;

  my $exists = 0;

  my $sth =
    $self->dbi->prepare_cached('SELECT xref_id FROM primary_xref WHERE xref_id=?');
  $sth->execute($xref_id) or confess $self->dbi->errstr();
  my @row    = $sth->fetchrow_array();
  my $result = $row[0];
  if ( defined $result ) { $exists = 1; }

  $sth->finish();

  return $exists;
} ## end sub primary_xref_id_exists


=head2 get_taxonomy_from_species_id
  Arg [1]    : species ID
  Description: Get the taxon id for a particular species id
  Return type: Hashref
  Caller     : internal

=cut

sub get_taxonomy_from_species_id {
  my ( $self, $species_id) = @_;
  my %hash;

  my $sth = $self->dbi->prepare_cached(
      "SELECT taxonomy_id FROM species WHERE species_id = ?");
  $sth->execute() or croak( $self->dbi->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $hash{ $row[0] } = 1;
  }

  return \%hash;
} ## end sub get_taxonomy_from_species_id


=head2 get_direct_xref
  Arg [1]    : stable ID
  Arg [2]    : ensembl object type (Gene, Transcript, Translation)
  Arg [3]    : linked xref id
  Description: xref_ids for a given stable id and linkage_xref
  Return type: array or integer
  Caller     : internal

=cut

sub get_direct_xref {
  my ( $self, $stable_id, $type, $link ) = @_;

  $type = lc $type;

  my %sql_hash = (
    gene =>
      "SELECT general_xref_id FROM gene_direct_xref d WHERE ensembl_stable_id = ? AND linkage_xref ",
    transcript =>
      "SELECT general_xref_id FROM transcript_direct_xref d WHERE ensembl_stable_id = ? AND linkage_xref ",
    translation =>
      "SELECT general_xref_id FROM translation_direct_xref d WHERE ensembl_stable_id = ? AND linkage_xref "
  );

  my $sql = $sql_hash{ $type };
  my @sql_params = ( $stable_id );
  if ( defined $link ) {
    $sql .= '= ?';
    push @sql_params, $link;
  }
  else {
    $sql .= 'IS NULL';
  }
  my $direct_sth = $self->dbi->prepare_cached($sql);

  $direct_sth->execute(@sql_params) || confess( $self->dbi->errstr() );
  if ( wantarray() ) {
    # Generic behaviour

    my @results;

    my $all_rows = $direct_sth->fetchall_arrayref();
    foreach my $row_ref ( @{$all_rows} ) {
      push @results, $row_ref->[0];
    }

    return @results;
  }
  else {
    # Backwards-compatible behaviour. FIXME: can we get rid of it?
    # There seem to be no parsers present relying on the old behaviour
    # any more
    if ( my @row = $direct_sth->fetchrow_array() ) {
      $direct_sth->finish();
      return $row[0];
    }
  }

  return;
} ## end sub get_direct_xref


=head2 get_xref
  Arg [1]    : accession
  Arg [2]    : source ID
  Arg [3]    : species ID
  Description: return the xref_id for a particular accession, source and species
               if not found return undef
  Return type: integer
  Caller     : internal

=cut

sub get_xref {
  my ( $self, $acc, $source, $species_id ) = @_;

  # If the statement handle does nt exist create it.
  my $sql =
    'SELECT xref_id FROM xref WHERE accession = ? AND source_id = ? AND species_id = ?';
  my $get_xref_sth = $self->dbi->prepare_cached($sql);

  # Find the xref_id using the sql above
  $get_xref_sth->execute( $acc, $source, $species_id ) or
    croak( $self->dbi->errstr() );
  if ( my @row = $get_xref_sth->fetchrow_array() ) {

    # Calling finish() as we only require the first row only
    $get_xref_sth->finish();

    return $row[0];
  }

  return;
} ## end sub get_xref


=head2 get_object_xref
  Arg [1]    : xref ID
  Arg [2]    : ensembl ID
  Arg [3]    : ensembl object type (Gene, Transcript, Translation)
  Description: return the object_xref_id for a particular xref_id, ensembl_id
               and ensembl_object_type. If not found return undef
  Return type: integer
  Caller     : internal

=cut

sub get_object_xref {
  my ( $self, $xref_id, $ensembl_id, $object_type ) = @_;

  my $sql =
'SELECT object_xref_id FROM object_xref WHERE xref_id = ? AND ensembl_object_type = ? AND ensembl_id = ?';
  my $get_object_xref_sth = $self->dbi->prepare_cached($sql);

  #
  # Find the object_xref_id using the sql above
  #
  $get_object_xref_sth->execute( $xref_id, $object_type, $ensembl_id )
    or
    croak( $self->dbi->errstr() );
  if ( my @row = $get_object_xref_sth->fetchrow_array() ) {
    $get_object_xref_sth->finish();
    return $row[0];
  }

  return;
} ## end sub get_object_xref


=head2 add_xref
  Arg [1]    : xref
  Description: Create an xref
               If the xref already exists, return the matching xref_id else create
               and return a new xref_id.
  Return type: integer
  Caller     : internal

=cut

sub add_xref {
  my ( $self, $arg_ref ) = @_;

  my $acc          = $arg_ref->{acc} || confess 'add_xref needs an acc';
  my $source_id    = $arg_ref->{source_id} || confess 'add_xref needs a source_id';
  my $species_id   = $arg_ref->{species_id} || confess 'add_xref needs a species_id';
  my $label        = $arg_ref->{label} // $acc;
  my $description  = $arg_ref->{desc};
  my $version      = $arg_ref->{version} // 0;
  my $info_type    = $arg_ref->{info_type} // 'MISC';
  my $info_text    = $arg_ref->{info_text} // q{};
  my $update_label = $arg_ref->{update_label};
  my $update_desc  = $arg_ref->{update_desc};

  # See if it already exists. It so return the xref_id for this one.
  my $xref_id = $self->get_xref( $acc, $source_id, $species_id );
  if ( defined $xref_id ) {
    if ( $update_label ) {
      $self->_update_xref_label( $xref_id, $label );
    }
    if ( $update_desc ) {
      $self->_update_xref_description( $xref_id, $description );
    }
    return $xref_id;
  }

    my $sql = (<<"AXS");
  INSERT INTO xref (
    accession, version, label, description,
    source_id, species_id, info_type, info_text)
  VALUES (?,?,?,?,?,?,?,?)
AXS

  my $add_xref_sth =
    $self->dbi->prepare_cached( $sql );

  # If the description is more than 255 characters, chop it off and add
  # an indication that it has been truncated to the end of it.
  if ( defined $description && ( ( length $description ) > 255 ) ) {
    my $truncmsg = ' /.../';
    substr $description, 255 - ( length $truncmsg ), length $truncmsg,
      $truncmsg;
  }

  # Add the xref and confess if it fails
  $add_xref_sth->execute( $acc, $version // 0, $label,
                          $description, $source_id, $species_id,
                          $info_type,   $info_text ) or
    confess $self->dbi->errstr() . "\n $acc\t$label\t\t$source_id\t$species_id\n";

  return $add_xref_sth->{'mysql_insertid'};
} ## end sub add_xref


=head2 add_object_xref
  Arg [1]    : xref
  Description: Create an object_xref
               If it already exists return the object_xref_id otherwise create it
               and return the new object_xref_id
  Return type: integer
  Caller     : internal

=cut

sub add_object_xref {
  my ( $self, $arg_ref ) = @_;

  my $xref_id = $arg_ref->{xref_id}         || confess 'add_object_xref needs an xref_id';
  my $ensembl_id = $arg_ref->{ensembl_id}   || confess 'add_object_xref needs an ensembl_id';
  my $object_type = $arg_ref->{object_type} || confess 'add_object_xref needs an object_type';

  # See if it already exists. It so return the xref_id for this one.
  my $object_xref_id =
    $self->get_object_xref( $xref_id, $ensembl_id, $object_type );
  if ( defined $object_xref_id ) {
    return $object_xref_id;
  }

  my $add_object_xref_sth = $self->dbi->prepare_cached(
    'INSERT INTO object_xref (ensembl_id, ensembl_object_type, xref_id) VALUES (?,?,?)'
  );

  # Add the object_xref and confess if it fails
  $add_object_xref_sth->execute( $ensembl_id, $object_type, $xref_id )
    or
    confess $self->dbi->errstr() . "\n $ensembl_id\t$object_type\t\t$xref_id\n";

  return $add_object_xref_sth->{'mysql_insertid'};
} ## end sub add_object_xref


=head2 add_identity_xref
  Arg [1]    : xref object
  Description: Create an identity_xref
  Return type:
  Exceptions : Throw if execution fails
  Caller     : internal

=cut

sub add_identity_xref {
  my ( $self, $arg_ref ) = @_;

  my $object_xref_id = $arg_ref->{object_xref_id} ||
    confess 'add_identity_xref needs an object_xref_id';
  my $score = $arg_ref->{score} ||
    confess 'add_identity_xref needs a score';
  my $target_identity = $arg_ref->{target_identity} ||
    confess 'add_identity_xref needs a target_identity';
  my $query_identity = $arg_ref->{query_identity} ||
    confess 'add_identity_xref needs a query_identity';

  my $add_identity_xref_sth =
    $self->dbi->prepare_cached( 'INSERT INTO identity_xref ' .
           '(object_xref_id, score, query_identity, target_identity) ' .
           'VALUES(?,?,?,?)' );

  # Add the object_xref and confess if it fails
  $add_identity_xref_sth->execute( $object_xref_id, $score,
                                   $query_identity, $target_identity
    ) or
    confess(
      $self->dbi->errstr() . "\n$object_xref_id\t$score\t\t$query_identity\t$target_identity\n");

  return;
} ## end sub add_identity_xref


=head2 add_to_direct_xrefs
  Arg [1]    : xref object
  Description: Create new xref if needed and add as a direct xref to a stable_id
               Note that a corresponding method for dependent xrefs is called
               add_dependent_xref()
  Return type:
  Caller     : internal

=cut

sub add_to_direct_xrefs {
  my ( $self, $arg_ref ) = @_;

  my $stable_id  = $arg_ref->{stable_id}  || confess('Need a direct_xref on which this xref linked too');
  my $type       = $arg_ref->{type}       || confess('Need a table type on which to add');
  my $acc        = $arg_ref->{acc}        || confess('Need an accession of this direct xref');
  my $source_id  = $arg_ref->{source_id}  || confess('Need a source_id for this direct xref');
  my $species_id = $arg_ref->{species_id} || confess('Need a species_id for this direct xref');
  my $version    = $arg_ref->{version}    // 0;
  my $label      = $arg_ref->{label}      // $acc;
  my $desc       = $arg_ref->{desc};
  my $linkage    = $arg_ref->{linkage};
  my $info_text  = $arg_ref->{info_text}  // q{};


  my $direct_xref_id = $self->add_xref( {
    acc        => $acc,
    source_id  => $source_id,
    species_id => $species_id,
    label      => $label,
    desc       => $desc,
    version    => $version,
    info_type  => 'DIRECT',
    info_text  => $info_text
  } );

  # Now add the direct info
  $self->add_direct_xref( $direct_xref_id, $stable_id, $type, $linkage );
  return;
} ## end sub add_to_direct_xrefs


=head2 add_direct_xref
  Arg [1]    : xref ID
  Arg [2]    : ensembl stable ID
  Arg [3]    : ensembl object type (Gene, Transcript, Translation)
  Arg [4]    : linkage type
  Arg [5]    : updated info type
  Description: Add a single record to the direct_xref table.
               Note that an xref must already have been added to the xref table
               Note that a corresponding method for dependent xrefs is called
               add_dependent_xref_maponly()
  Return type:
  Caller     : internal

=cut

sub add_direct_xref {
  my ( $self, $general_xref_id, $ensembl_stable_id, $ensembl_type,
       $linkage_type, $update_info_type )
    = @_;

  # Check if such a mapping exists yet. Make sure get_direct_xref() is
  # invoked in list context, otherwise it will fall back to legacy
  # behaviour of returning a single xref_id even when multiple ones
  # match.
  # Note: get_direct_xref() does not currently cache its output,
  # consider changing this should performance become an issue
  my @existing_xref_ids =
    $self->get_direct_xref( $ensembl_stable_id, $ensembl_type,
                            $linkage_type );
  if ( scalar grep { $_ == $general_xref_id } @existing_xref_ids ) {
    return;
  }

  my %sql_hash = (
    gene => 'INSERT INTO gene_direct_xref VALUES (?,?,?)',
    transcript => 'INSERT INTO transcript_direct_xref VALUES (?,?,?)',
    translation => 'INSERT INTO translation_direct_xref VALUES (?,?,?)'
  );

  my $add_direct_xref_sth = $self->dbi->prepare_cached( $sql_hash{ lc $ensembl_type } );

  $add_direct_xref_sth->execute( $general_xref_id, $ensembl_stable_id,
                                 $linkage_type );

  if ( defined $update_info_type ) {
    $self->_update_xref_info_type( $general_xref_id, 'DIRECT');
  }

  return;
} ## end sub add_direct_xref


=head2 add_multiple_direct_xrefs
  Arg [1]    : xref object
  Description: Add multiple records to the direct_xref table.
  Return type:
  Caller     : internal

=cut

sub add_multiple_direct_xrefs {
  my ( $self, $direct_xrefs ) = @_;

  foreach my $direct_xref ( @{ $direct_xrefs } ) {
    $self->add_to_direct_xrefs( {
      stable_id  => $direct_xref->{STABLE_ID},
      type       => $direct_xref->{ENSEMBL_TYPE},
      linkage    => $direct_xref->{LINKAGE_TYPE},
      acc        => $direct_xref->{ACCESSION},
      version    => $direct_xref->{VERSION} // 0,
      label      => $direct_xref->{LABEL}   // $direct_xref->{ACCESSION},
      desc       => $direct_xref->{DESCRIPTION},
      source_id  => $direct_xref->{SOURCE_ID},
      species_id => $direct_xref->{SPECIES_ID},
      info_text  => $direct_xref->{INFO_TEXT},
      info_type  => $direct_xref->{LINKAGE_TYPE} } );
  }

  return;
} ## sub add_multiple_direct_xrefs


=head2 add_dependent_xref
  Arg [1]    : dependent xref
  Description: Create/Add xref and add it as a dependency of the master
               Note that a corresponding method for direct xrefs is called
               add_to_direct_xrefs()
  Return type: integer
  Caller     : internal

=cut

sub add_dependent_xref {
  my ( $self, $arg_ref ) = @_;

  my $master_xref = $arg_ref->{master_xref_id}  || confess('Need a master_xref_id on which this xref linked too');
  my $acc         = $arg_ref->{acc}        || confess('Need an accession of this dependent xref');
  my $source_id   = $arg_ref->{source_id}  || confess('Need a source_id for this dependent xref');
  my $species_id  = $arg_ref->{species_id} || confess('Need a species_id for this dependent xref');
  my $version     = $arg_ref->{version}    // 0;
  my $label       = $arg_ref->{label}      // $acc;
  my $desc        = $arg_ref->{desc};
  my $linkage     = $arg_ref->{linkage};
  my $info_text   = $arg_ref->{info_text}  // q{};


  my $dependent_xref_id = $self->add_xref( {
    acc        => $acc,
    source_id  => $source_id,
    species_id => $species_id,
    label      => $label,
    desc       => $desc,
    version    => $version,
    info_type  => 'DEPENDENT',
    info_text  => $info_text
  } );

  # Now add the dependency mapping
  $self->add_dependent_xref_maponly( $dependent_xref_id, $source_id,
                                     $master_xref, $linkage );

  return $dependent_xref_id;
} ## end sub add_dependent_xref


=head2 add_dependent_xref_maponly
  Arg [1]    : dependent xref ID
  Arg [2]    : dependent source ID
  Arg [3]    : master xref ID
  Arg [4]    : master source ID
  Arg [5]    : update info type
  Description: Add a single record to the dependent_xref table.
               Note that an xref must already have been added to the xref table
               Note that a corresponding method for direct xrefs is called add_direct_xref()
  Return type:
  Caller     : internal

=cut

sub add_dependent_xref_maponly {
  my ( $self, $dependent_id, $dependent_source_id, $master_id,
       $master_source_id, $update_info_type )
    = @_;

  my $sql = (<<'ADX');
INSERT INTO dependent_xref
  (master_xref_id,dependent_xref_id,linkage_annotation,linkage_source_id)
  VALUES (?,?,?,?)
ADX
  my $add_dependent_xref_sth = $self->dbi->prepare_cached($sql);

  # If the dependency cannot be found in %xref_dependent_mapped,
  # i.e. has not been set yet, add it
  if (( !defined $xref_dependent_mapped{"$master_id|$dependent_id"} ) ||
      $xref_dependent_mapped{"$master_id|$dependent_id"} ne
      $master_source_id )
  {

    $add_dependent_xref_sth->execute( $master_id, $dependent_id,
                            $master_source_id, $dependent_source_id ) ||
      confess(
        $self->dbi->errstr() . "\n$master_id\t$dependent_id\t$master_source_id\t$dependent_source_id" );

    $xref_dependent_mapped{"$master_id|$dependent_id"} =
      $master_source_id;
  }

  if ( defined $update_info_type ) {
    $self->_update_xref_info_type( $dependent_id, 'DEPENDENT' );
  }

  return;
} ## end sub add_dependent_xref_maponly


=head2 add_multiple_dependent_xrefs
  Arg [1]    : xref ID
  Arg [2]    : xref
  Description: Add dependent xrefs to the xref table along with dependent xref mappings
  Return type:
  Caller     : internal

=cut

sub add_multiple_dependent_xrefs {
  my ( $self, $xref_id, $dependent_xrefs ) = @_;

  foreach my $depref ( @{ $dependent_xrefs } ) {
    my %dep = %{$depref};

    # Insert the xref
    my $dep_xref_id = $self->add_xref( (
      "acc"        => $dep{ACCESSION},
      "version"    => $dep{VERSION}     // 1,
      "label"      => $dep{LABEL}       // $dep{ACCESSION},
      "desc"       => $dep{DESCRIPTION},
      "source_id"  => $dep{SOURCE_ID},
      "species_id" => $dep{SPECIES_ID},
      "info_type"  => 'DEPENDENT' ) );

    # Add the linkage_annotation and source id it came from
    $self->add_dependent_xref_maponly(
      $dep_xref_id,
      $dep{LINKAGE_SOURCE_ID},
      $xref_id,
      $dep{LINKAGE_ANNOTATION} );

    # if there are synonyms, add entries in the synonym table
    foreach my $syn ( @{ $dep{SYNONYMS} } ) {
      $self->add_synonym( $dep_xref_id, $syn );
    }  # foreach syn
  }  # foreach dep

  return;
} ## end sub add_multiple_dependent_xrefs


=head2 add_to_syn_for_mult_sources
  Arg [1]    : accession
  Arg [2]    : source IDs
  Arg [3]    : synonym
  Arg [4]    : species ID
  Description: Add synonyms for a particular accession for one or more sources.
               This is for priority xrefs where we have more than one source
               but want to write synonyms for each with the same accession
  Return type:
  Caller     : internal

=cut

sub add_to_syn_for_mult_sources {
  my ( $self, $acc, $sources, $syn, $species_id ) = @_;

  foreach my $source_id ( @{$sources} ) {
    my $xref_id =
      $self->get_xref( $acc, $source_id, $species_id );
    if ( defined $xref_id ) {
      $self->add_synonym( $xref_id, $syn );
    }
  }

  return;
} ## end sub add_to_syn_for_mult_sources


=head2 add_to_syn
  Arg [1]    : accession
  Arg [2]    : source id
  Arg [3]    : synonym
  Arg [4]    : species ID
  Description: Add synomyn for an xref given by accession and source_id
  Return type:
  Caller     : internal

=cut

sub add_to_syn {
  my ( $self, $acc, $source_id, $syn, $species_id ) = @_;

  my $xref_id = $self->get_xref( $acc, $source_id, $species_id );
  if ( defined $xref_id ) {
    $self->add_synonym( $xref_id, $syn );
  }
  else {
    confess( "Could not find acc $acc in " .
             "xref table source = $source_id of species $species_id\n" );
  }

  return;
} ## end sub add_to_syn


=head2 add_synonym
  Arg [1]    : xref ID
  Arg [2]    : synonym
  Description: Add synomyn for an xref given by xref_id
  Return type:
  Caller     : internal

=cut

sub add_synonym {
  my ( $self, $xref_id, $syn ) = @_;
  my $add_synonym_sth =
    $self->dbi->prepare_cached(
      'INSERT IGNORE INTO synonym ( xref_id, synonym ) VALUES(?,?)');
  $add_synonym_sth->execute( $xref_id, $syn ) or
    confess ( $self->dbi->errstr() . "\n $xref_id\n $syn\n\n" . $add_synonym_sth->errstr() . "\n" );

  return;
} ## sub add_synonym


=head2 add_multiple_synonyms
  Arg [1]    : xref ID
  Arg [2]    : Listref : synonyms
  Description: Add multiple synomyns for an xref given by xref_id
  Return type:
  Caller     : internal

=cut

sub add_multiple_synonyms {
  my ( $self, $xref_id, $synonyms ) = @_;

  foreach my $syn ( @{ $synonyms } ) {
    $self->add_synonym( $xref_id, $syn );
  }

  return
} ## sub add_multiple_synonyms


=head2 get_label_to_acc
  Arg [1]    : description
  Arg [2]    : species ID
  Arg [3]    : source priority description
  Description: Create a hash that uses the label as a key and the acc as the
               value. Also add synonyms for these as keys.
  Return type: Hashref
  Caller     : internal

=cut

sub get_label_to_acc {
  my ( $self, $name, $species_id, $prio_desc ) = @_;
  my %hash1 = ();

  my $sql = (<<"GLA");
  SELECT  xref.accession, xref.label
  FROM xref, source
  WHERE
    source.name LIKE ? AND
    xref.source_id = source.source_id
GLA
  my @sql_params = ($name . '%');
  if ( defined $prio_desc ) {
    $sql .= " AND source.priority_description LIKE ?";
    push @sql_params, $prio_desc;
  }
  if ( defined $species_id ) {
    $sql .= " AND xref.species_id = ?";
    push @sql_params, $species_id;
  }
  my $sub_sth = $self->dbi->prepare_cached($sql);

  $sub_sth->execute(@sql_params);
  while ( my @row = $sub_sth->fetchrow_array() ) {
    $hash1{ $row[1] } = $row[0];
  }

  # Remember synonyms
  $sql = (<<"GLS");
  SELECT  xref.accession, synonym.synonym
  FROM xref, source, synonym
  WHERE
    synonym.xref_id = xref.xref_id AND
    source.name LIKE ? AND
    xref.source_id = source.source_id
GLS

  @sql_params = ($name . '%');
  if ( defined $prio_desc ) {
    $sql .= " AND source.priority_description LIKE ?";
    push @sql_params, $prio_desc;
  }
  if ( defined $species_id ) {
    $sql .= " AND xref.species_id  = ?";
    push @sql_params, $species_id;
  }
  $sub_sth = $self->dbi->prepare_cached($sql);

  $sub_sth->execute(@sql_params);
  while ( my @row = $sub_sth->fetchrow_array() ) {
    $hash1{ $row[1] } = $row[0];
  }

  return \%hash1;
} ## end sub get_label_to_acc


=head2 get_acc_to_label
  Arg [1]    : description
  Arg [2]    : species ID
  Arg [3]    : source priority description
  Description: Create a hash that uses the accession as a key and the label as
               the value.
  Return type: Hashref
  Caller     : internal

=cut

sub get_acc_to_label {
  my ( $self, $name, $species_id, $prio_desc ) = @_;
  my %hash1 = ();

  my $sql = (<<"GLA");
  SELECT  xref.accession, xref.label
  FROM xref, source
  WHERE
    source.name LIKE ? AND
    xref.source_id = source.source_id
GLA

  my @sql_params = ($name . '%');
  if ( defined $prio_desc ) {
    $sql .= " AND source.priority_description LIKE ?";
    push @sql_params, $prio_desc;
  }
  if ( defined $species_id ) {
    $sql .= " AND xref.species_id  = ?";
    push @sql_params, $species_id;
  }
  my $sub_sth = $self->dbi->prepare_cached($sql);

  $sub_sth->execute(@sql_params);
  while ( my @row = $sub_sth->fetchrow_array() ) {
    $hash1{ $row[0] } = $row[1];
  }

  return \%hash1;
} ## end sub get_acc_to_label


=head2 get_label_to_desc
  Arg [1]    : description
  Arg [2]    : species ID
  Arg [3]    : source priority description
  Description: Create a hash that uses the label as a key and the desc as the
               value. Also add synonyms for these as keys.
  Return type: Hashref
  Caller     : internal

=cut

sub get_label_to_desc {
  my ( $self, $name, $species_id, $prio_desc ) = @_;
  my %hash1 = ();

  my $sql = (<<"GDH");
  SELECT xref.description, xref.label
  FROM xref, source
  WHERE
    source.name LIKE ? AND
    xref.source_id = source.source_id
GDH

  my @sql_params = ($name . '%');
  if ( defined $prio_desc ) {
    $sql .= " AND source.priority_description LIKE ?";
    push @sql_params, $prio_desc;
  }
  if ( defined $species_id ) {
    $sql .= " and xref.species_id  = ?";
    push @sql_params, $species_id;
  }
  my $sub_sth = $self->dbi->prepare_cached($sql);

  $sub_sth->execute(@sql_params);
  while ( my @row = $sub_sth->fetchrow_array() ) {
    $hash1{ $row[1] } = $row[0];
  }

  # Also include the synonyms
  my $syn_sql = (<<"GDS");
  SELECT xref.description, synonym.synonym
  FROM xref, source, synonym
  WHERE
    synonym.xref_id = xref.xref_id AND
    source.name LIKE ? AND
    xref.source_id = source.source_id
GDS

  @sql_params = ($name . '%');
  if ( defined $prio_desc ) {
    $syn_sql .= " AND source.priority_description LIKE ?";
    push @sql_params, $prio_desc;
  }
  if ( defined $species_id ) {
    $syn_sql .= " AND xref.species_id  = ?";
    push @sql_params, $species_id;
  }
  $sub_sth = $self->dbi->prepare_cached($syn_sql);

  $sub_sth->execute(@sql_params);
  while ( my @row = $sub_sth->fetchrow_array() ) {
    $hash1{ $row[1] } = $row[0];
  }

  return \%hash1;
} ## end sub get_label_to_desc


=head2 set_release
  Arg [1]    : source ID
  Arg [2]    : source release
  Description: Set release for a particular source_id.
  Return type:
  Caller     : internal

=cut

sub set_release {
  my ( $self, $source_id, $s_release ) = @_;

  my $sth = $self->dbi->prepare_cached(
                'UPDATE source SET source_release=? WHERE source_id=?');

  if ($verbose) {
    print "Setting release to '$s_release' for source ID '$source_id'\n";
  }

  $sth->execute( $s_release, $source_id );

  return;
} ## end sub set_release


=head2 get_dependent_mappings
  Arg [1]    : source_id
  Description: create a hash of all the dependent mapping that exist for a given
               source_id of the format {master_xref_id|dependent_xref_id}
  Return type:
  Caller     : internal

=cut

sub get_dependent_mappings {
  my $self      = shift;
  my $source_id = shift;

  my $sql = (<<"GDM");
  SELECT  d.master_xref_id, d.dependent_xref_id, d.linkage_annotation
  FROM dependent_xref d, xref x
  WHERE
    x.xref_id = d.dependent_xref_id AND
    x.source_id = ?
GDM
  my $sth = $self->dbi->prepare_cached($sql);
  $sth->execute( $source_id );
  my $master_xref;
  my $dependent_xref;
  my $linkage;
  $sth->bind_columns( \$master_xref, \$dependent_xref, \$linkage );
  while ( $sth->fetch ) {
    $xref_dependent_mapped{"$master_xref|$dependent_xref"} = $linkage;
  }

  return;
} ## end sub get_dependent_mappings


=head2 get_ext_synonyms
  Arg [1]    : source name
  Description: Create a hashref that uses the accession and labels for keys and an
               array of the synonyms as the values
  Return type: Hashref
  Caller     : internal

=cut

sub get_ext_synonyms {
  my $self        = shift;
  my $source_name = shift;

  my %ext_syns;
  my %seen; # can be in more than once fro each type of external source.
  my $separator = qw{:};

  my $sql = (<<"GES");
  SELECT  x.accession, x.label, sy.synonym
  FROM xref x, source so, synonym sy
  WHERE x.xref_id = sy.xref_id AND
    so.source_id = x.source_id AND
    so.name LIKE ?
GES
  my $sth = $self->dbi->prepare_cached($sql);

  $sth->execute( '$source_name' );
  my ( $acc, $label, $syn );
  $sth->bind_columns( \$acc, \$label, \$syn );

  my $count = 0;
  while ( $sth->fetch ) {
    if ( !( defined $seen{ $acc . $separator . $syn } ) ) {
      push @{ $ext_syns{$acc} },   $syn;
      push @{ $ext_syns{$label} }, $syn;
      $count++;
    }
    $seen{ $acc . $separator . $syn } = 1;
  }

  return \%ext_syns;

} ## end sub get_ext_synonyms


=head2 parsing_finished_store_data
  Description: Store data needed to be able to revert to same stage as after parsing
  Return type:
  Caller     : internal

  Notes      : Store max id for

    gene/transcript/translation_direct_xref     general_xref_id  #Does this change??
    xref                                        xref_id
    dependent_xref                              object_xref_id is all null
    go_xref                                     object_xref_id
    object_xref                                 object_xref_id
    identity_xref                               object_xref_id

=cut

sub parsing_finished_store_data {
  my $self = shift;

  my %table_and_key = (
    'xref'        => 'SELECT MAX(xref_id) FROM xref',
    'object_xref' => 'SELECT MAX(object_xref_id) FROM object_xref' );

  foreach my $table ( keys %table_and_key ) {
    my $sth = $self->dbi->prepare_cached( $table_and_key{$table} );
    $sth->execute;
    my $max_val;
    $sth->bind_columns( \$max_val );
    $sth->fetch;
    $self->add_meta_pair( 'PARSED_' . $table . '_id',
                          $max_val || 1);
    $sth->finish();
  }
  return;
} ## end sub parsing_finished_store_data


=head2 get_meta_value
  Arg [1]    : key
  Description: Return metadata value
  Return type: integer or string
  Caller     : internal

=cut

sub get_meta_value {
  my ( $self, $key ) = @_;

  my $sth = $self->dbi->prepare_cached(
    'SELECT meta_value FROM meta WHERE meta_key LIKE ? ORDER BY meta_id' );
  $sth->execute( $key );
  my $value;
  $sth->bind_columns( \$value );
  while ( $sth->fetch ) {    # get the last one
  }

  return $value;
} ## end sub get_meta_value


=head2 _update_xref_info_type
  Arg [1]    : xref ID
  Arg [2]    : info type
  Description: Update info_type of an existing xref
  Return type:
  Exceptions : confess if UPDATE fails
  Caller     : internal

=cut

sub _update_xref_info_type {
  my ( $self, $xref_id, $info_type ) = @_;

  my $sth = $self->dbi->prepare_cached(
    'UPDATE xref SET info_type=? WHERE xref_id=?');
  if ( !$sth->execute( $info_type, $xref_id ) ) {
    confess $self->dbi->errstr() . "\n $xref_id\n $info_type\n\n";
  }

  return;
} ## end sub _update_xref_info_type


=head2 _add_pair
  Arg [1]    : source ID
  Arg [2]    : accession
  Arg [3]    : refseq dna/pep pair
  Description: Create a pair entry. refseq dna/pep pairs usually
  Return type:
  Exceptions : confess if INSERT fails
  Caller     : internal

=cut

sub _add_pair {
  my ( $self, $source_id, $accession, $pair ) = @_;

  my $pair_sth = $self->dbi->prepare_cached(
    'INSERT INTO pairs (source_id, accession1, accession2) VALUES(?,?,?)');

  # Add the pair and confess if it fails
  $pair_sth->execute( $source_id, $accession, $pair ) or
    confess $self->dbi->errstr() . "\n $source_id\t$\t$accession\t$pair\n";

  return;
} ## end sub _add_pair


=head2 _add_primary_xref
  Arg [1]    : xref ID
  Arg [2]    : sequence
  Arg [3]    : sequence type
  Arg [4]    : status
  Description: Create an primary_xref entry.
  Return type: integer
  Exceptions :confess if INSERT fails
  Caller     : internal

=cut

sub _add_primary_xref {
  my ( $self, $xref_id, $sequence, $sequence_type, $status ) = @_;

  my $add_primary_xref_sth =
    $self->dbi->prepare_cached( 'INSERT INTO primary_xref VALUES (?,?,?,?)' );

  # Add the xref and confess if it fails
  $add_primary_xref_sth->execute( $xref_id, $sequence, $sequence_type, $status ) or
    confess $self->dbi->errstr() . "\n $xref_id\t$sequence_type\t$status\n";

  return $add_primary_xref_sth->{'mysql_insertid'};
} ## end sub _add_primary_xref


=head2 _update_primary_xref_sequence
  Arg [1]    : xref ID
  Arg [2]    : sequence
  Description: Update primary_xref sequence for matching xref_id
  Return type:
  Caller     : internal

=cut

sub _update_primary_xref_sequence {
  my ( $self, $xref_id, $sequence ) = @_;

  my $sth = $self->dbi->prepare_cached(
    'UPDATE primary_xref SET sequence=? WHERE xref_id=?');

  $sth->execute( $sequence, $xref_id ) or
    confess $self->dbi->errstr() . "\n $xref_id\n $sequence\n\n";

  return;
} ## sub _update_primary_xref_sequence


=head2 _update_xref_label
  Arg [1]    : xref ID
  Arg [2]    : label
  Description: Update xref label for matching xref_id
  Return type:
  Exceptions : confess on a failed UPDATE
  Caller     : internal

=cut

sub _update_xref_label {
  my ( $self, $xref_id, $label ) = @_;

  my $sth = $self->dbi->prepare_cached(
    'UPDATE xref SET label=? WHERE xref_id=?');

  $sth->execute( $label, $xref_id ) or
    confess $self->dbi->errstr() . "\n $xref_id\n $label\n\n";

  return;
} ## sub _update_xref_label


=head2 _update_xref_description
  Arg [1]    : xref ID
  Arg [2]    : description
  Description: Update xref dfescription for matching xref_id
  Return type:
  Exceptions : confess on a failed UPDATE
  Caller     : internal

=cut

sub _update_xref_description {
  my ( $self, $xref_id, $description ) = @_;

  my $sth = $self->dbi->prepare_cached(
    'UPDATE xref SET description=? WHERE xref_id=?');

  $sth->execute( $description, $xref_id ) or
    confess $self->dbi->errstr() . "\n $xref_id\n $description\n\n";

  return;
} ## sub _update_xref_description


=head2 _update_xref_description
  Description: Refseq sources to consider. Prefixes not in this list will be ignored.
               To be used by RefSeq parsers
  Return type: hashref
  Caller     : internal

=cut

sub get_refseq_sources {
  return {
      NM => 'RefSeq_mRNA',
      NR => 'RefSeq_ncRNA',
      XM => 'RefSeq_mRNA_predicted',
      XR => 'RefSeq_ncRNA_predicted',
      NP => 'RefSeq_peptide',
      XP => 'RefSeq_peptide_predicted',
  };
}



1;

