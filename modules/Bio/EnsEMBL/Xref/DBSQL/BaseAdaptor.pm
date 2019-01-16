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

use Data::Dumper;

use Bio::EnsEMBL::Utils::Exception;
use Getopt::Long;
use IO::Uncompress::AnyUncompress '$AnyUncompressError';

use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $base_dir = File::Spec->curdir();

my %xref_dependent_mapped;

my $verbose;

=head2 new
  Arg [1]    : proto
  Arg [2]    : arguments : { host => string, dbname => string, user => string, pass => string, port => int, group => string }
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
    -PASS   => $args{pass} // q{},
    -PORT => $args{port} // '3306',
    -GROUP => $args{group} // q{}
  ) );
  $self->verbose( $args{verbose} // 0 );

  return $self;
}

=head2 dbc
  Description: Getter/Setter for the dbc object
  Return type: db_connection
  Caller     : internal

=cut

sub dbc {
  my ( $self, $arg ) = @_;
  ( defined $arg ) && ( $self->{_dbc} = $arg );
  return $self->{_dbc};
}

=head2 dbi
  Description: Getter/Setter for the dbi object
  Return type: DBI database handle
  Caller     : internal

=cut

sub dbi {
  my $self = shift;

  return $self->dbc->db_handle;
}

=head2 dba
  Description: Getter for the dba object
  Return type: Bio::EnsEMBL::DBSQL::DBAdaptor instance
  Caller     : internal
=cut

sub dba {
  my $self = shift;
  return
    Bio::EnsEMBL::DBSQL::DBAdaptor->new( -dbconn  => $self->dbc,
                                         -species => $self->species );
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

sub species {
  my ( $self, $arg ) = @_;

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

  if ( !( defined $file_name ) || $file_name eq q{} ) {
    confess 'No file name';
  }
  my $alt_file_name = $file_name;
  $alt_file_name =~ s/\.(gz|Z)$//x;

  if ( $alt_file_name eq $file_name ) {
    $alt_file_name .= '.gz';
  }

  if ( !-e $file_name ) {
    print "File $file_name does not exist, will try $alt_file_name\n";
    $file_name = $alt_file_name;
  }

  # 'Transparent' lets IO::Uncompress modules read uncompressed input.
  # It should be on by default but set it just in case.
  $io = IO::Uncompress::AnyUncompress->new($file_name,
                                           'Transparent' => 1 )
    || confess("Can not open file '$file_name' because: $AnyUncompressError");

  if ( $verbose ) {
    print "Reading from '$file_name'...\n";
  }

  return $io;
} ## end sub get_filehandle

=head2 get_id_from_species_name

  Arg [1]    : speciesname
  Description: Gets the species ID for a given species name
  Return type: integer
  Exceptions : confesses if not found
  Caller     : Bio::EnsEMBL::Xref::Mapper

=cut

sub get_id_from_species_name {
  my ( $self, $name ) = @_;

  if ( !defined $name ) {
    confess 'Undefined species name';
  }

  my $sql = 'SELECT species_id FROM species WHERE name=?';

  my @sql_params = ($name);
  my $sth        = $self->dbi->prepare_cached($sql);
  $sth->execute(@sql_params);

  my @row = $sth->fetchrow_array();
  my $species_id;
  if (@row) {
    $species_id = $row[0];
  }
  else {
    my $msg = "No ID for species name = '${name}'";
    confess $msg;
  }

  return $species_id;
} ## end sub get_id_from_species_name

=head2 get_species_names

=cut

sub get_species_names {
  my $self = shift;

  my $sth = $self->dbi->prepare_cached('SELECT name FROM species');
  $sth->execute();

  return $sth->fetchrow_arrayref;

}

=head2 get_source_id_for_source_name
  Arg [1]    : source name
  Arg [2]    : priority description
  Description: Gets the source ID for a given source name
  Return type: integer
  Exceptions : confesses if not found
  Caller     : internal

=cut

sub get_source_id_for_source_name {
  my ( $self, $source_name, $priority_desc ) = @_;

  if ( !defined $source_name ) {
    confess 'source_name undefined';
  }

  my $low_name = lc $source_name;
  my $sql      = 'SELECT source_id FROM source WHERE LOWER(name)=?';

  my @sql_params = ($low_name);
  if ( defined $priority_desc ) {
    $sql .= ' AND LOWER(priority_description)=?';
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

  my ( $self, $source_name ) = @_;

  if ( !defined $source_name ) {
    confess 'source_name undefined';
  }

  my $sth = $self->dbi->prepare_cached('SELECT source_id FROM source WHERE UPPER(name) LIKE ?');

  my $big_name = uc $source_name;
  $sth->execute("%${big_name}%");

  my @sources;
  while ( my @row = $sth->fetchrow_array() ) {
    push @sources, $row[0];
  }

  return @sources;
}

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
    my $error_msg = <<"ERROR";
There is no entity with source-id  $source_id  in the source-table of the
xref-database. The source-id and the name of the source-id is hard-coded in
populate_metadata.sql and in the parser.
Couldn't get source name for source ID $source_id
ERROR

    confess $error_msg;
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
  my $sql = 'SELECT source_id FROM source WHERE UPPER(name) LIKE ?';
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
    confess 'Please give me some xrefs to load';
  }

  foreach my $xref ( @{$rxrefs} ) {
    if ( !( defined $xref->{ACCESSION} ) ) {
      confess "Your xref does not have an accession-number\n";
    }
    if ( !( defined $xref->{SOURCE_ID} ) ) {
      confess
        "Your xref: $xref->{ACCESSION} does not have a source-id\n";
    }

    # Create entry in xref table and note ID
    my $xref_id = $self->add_xref( {
      'acc'          => $xref->{ACCESSION},
      'version'      => $xref->{VERSION} // 0,
      'label'        => $xref->{LABEL}   // $xref->{ACCESSION},
      'desc'         => $xref->{DESCRIPTION},
      'source_id'    => $xref->{SOURCE_ID},
      'species_id'   => $xref->{SPECIES_ID},
      'info_type'    => $xref->{INFO_TYPE},
      'info_text'    => $xref->{INFO_TEXT},
      'update_label' => 1,
      'update_desc' => 1
    } );

    # If there are any direct_xrefs, add these to the relevant tables
    $self->add_multiple_direct_xrefs( $xref->{DIRECT_XREFS} );

    # create entry in primary_xref table with sequence; if this is a "cumulative"
    # entry it may already exist, and require an UPDATE rather than an INSERT
    if ( defined $xref->{SEQUENCE} ) {
      if ( $self->primary_xref_id_exists($xref_id) ) {
        $self->_update_primary_xref_sequence( $xref->{SEQUENCE}, $xref_id );
      }
      else {
        $self->_add_primary_xref(
          $xref_id, $xref->{SEQUENCE},
          $xref->{SEQUENCE_TYPE},
          $xref->{STATUS} );
      }
    }

    # if there are synonyms, add entries in the synonym table
    $self->add_multiple_synonyms( $xref_id, $xref->{SYNONYMS} );

  # if there are dependent xrefs, add xrefs and dependent xrefs for them
    $self->add_multiple_dependent_xrefs( $xref_id, $xref->{DEPENDENT_XREFS} );

    # Add the pair data. refseq dna/pep pairs usually
    if ( defined $xref->{PAIR} ) {
      $self->_add_pair( $xref->{SOURCE_ID}, $xref->{ACCESSION},
                        $xref->{PAIR} );
    }

  }    # foreach xref

  return;
} ## end sub upload_xref_object_graphs


=head2 add_alt_allele

  Arg [1]    : Scalar; the allele id
  Arg [2]    : Scalar; the gene id
  Arg [3]    : Boolean; if it's reference
  Description: Insert into the alt allele table
  Return type: None
  Caller     : Bio::EnsEMBL::Xref::Mapper

=cut

sub add_alt_allele {

  my ( $self, $allele_id, $gene_id, $is_reference ) = @_;

  if ( !defined $allele_id || !defined $gene_id || !defined $is_reference ) {
    confess 'Need to specify (allele_id, gene_id, is_reference)';
  }

  my $sth = $self->dbi->prepare_cached(
    'INSERT INTO alt_allele (alt_allele_id, gene_id, is_reference) VALUES (?, ?, ?)'
  );
  $sth->execute( $allele_id, $gene_id, $is_reference );

  return;
}

=head2 delete_alt_alleles

  Arg [ ]    : none
  Description: Delete entries in alt allele table
  Return type: none
  Caller     : Bio::EnsEMBL::Xref::Mapper

=cut

sub delete_alt_alleles {
  my $self = shift;

  my $sth = $self->dbi->prepare_cached('TRUNCATE alt_allele');
  $sth->execute;

  return;
}

=head2 get_xref_sources
  Description: Create a hash of all the source names for xrefs
  Return type: Hashref
  Caller     : internal

=cut

sub get_xref_sources {

  my $self = shift;

  my %sourcename_to_sourceid;

  my $sth =
    $self->dbi->prepare_cached('SELECT name, source_id FROM source');
  $sth->execute() or croak( $self->dbi->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $source_name = $row[0];
    my $source_id   = $row[1];
    $sourcename_to_sourceid{$source_name} = $source_id;
  }

  return %sourcename_to_sourceid;
}

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

  my $sth = $self->dbi->prepare_cached(
                    'SELECT xref_id FROM primary_xref WHERE xref_id=?');
  $sth->execute($xref_id) or confess $self->dbi->errstr();
  my @row    = $sth->fetchrow_array();
  my $result = $row[0];
  if ( defined $result ) { $exists = 1; }

  $sth->finish();

  return $exists;
} ## end sub primary_xref_id_exists


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
      'SELECT general_xref_id FROM gene_direct_xref d WHERE ensembl_stable_id = ? AND linkage_xref ',
    transcript =>
      'SELECT general_xref_id FROM transcript_direct_xref d WHERE ensembl_stable_id = ? AND linkage_xref ',
    translation =>
      'SELECT general_xref_id FROM translation_direct_xref d WHERE ensembl_stable_id = ? AND linkage_xref ',
  );

  my $sql        = $sql_hash{$type};
  my @sql_params = ($stable_id);
  if ( defined $link ) {
    $sql .= '= ?';
    push @sql_params, $link;
  }
  else {
    $sql .= 'IS NULL';
  }
  my $direct_sth = $self->dbi->prepare_cached($sql);

  $direct_sth->execute(@sql_params) || confess( $self->dbi->errstr() );
  if ( wantarray ) {
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
}

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

  my $sql = (<<'SQL');
    SELECT object_xref_id
    FROM object_xref
    WHERE xref_id = ? AND
          ensembl_object_type = ? AND
          ensembl_id = ?
SQL

  my $get_object_xref_sth = $self->dbi->prepare_cached( $sql );

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

=head2 add_coordinate_xref
  Arg [1]    : coordinate_xref
  Description: Create an entry in the coordinate_xref table
               and return a new coordinate_xref_id.
  Return type: integer
  Caller     : internal

=cut

sub add_coordinate_xref {
  my ( $self, $arg_ref ) = @_;

  my $accession    = $arg_ref->{accession}  || confess 'add_coordinate_xref needs an accession';
  my $source_id    = $arg_ref->{source_id}  || confess 'add_coordinate_xref needs a source_id';
  my $species_id   = $arg_ref->{species_id} || confess 'add_coordinate_xref needs a species_id';
  my $chromosome   = $arg_ref->{chromosome} || confess 'add_coordinate_xref needs a chromosome';
  my $strand       = $arg_ref->{strand}     || confess 'add_coordinate_xref needs a strand';
  my $tx_start     = $arg_ref->{txStart}    || confess 'add_coordinate_xref needs a txStart';
  my $tx_end       = $arg_ref->{txEnd}      || confess 'add_coordinate_xref needs a txEnd';
  my $cds_start    = $arg_ref->{cdsStart};
  my $cds_end      = $arg_ref->{cdsEnd};
  my $exon_starts  = $arg_ref->{exonStarts} || confess 'add_coordinate_xref needs exonStarts';
  my $exon_ends    = $arg_ref->{exonEnds}   || confess 'add_coordinate_xref needs exonEnds';

  my $sql = (<<"SQL");
  INSERT INTO coordinate_xref
  ( source_id,  species_id,
    accession,
    chromosome, strand,
    txStart,    txEnd,
    cdsStart,   cdsEnd,
    exonStarts, exonEnds )
  VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
SQL

  my $add_coordinate_xref_sth =
    $self->dbi->prepare_cached( $sql );

  # Add the xref and confess if it fails
  $add_coordinate_xref_sth->execute( $source_id, $species_id, $accession, $chromosome,
                          $strand, $tx_start, $tx_end, $cds_start, $cds_end,
                          $exon_starts, $exon_ends ) or
    confess $self->dbi->errstr() . "\n$accession\t$source_id\t$species_id\n";

  return $add_coordinate_xref_sth->{'mysql_insertid'};

} ## end sub add_coordinate_xref


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

  my $acc          = $arg_ref->{acc}        || confess 'add_xref needs an acc';
  my $source_id    = $arg_ref->{source_id}  || confess 'add_xref needs a source_id';
  my $species_id   = $arg_ref->{species_id} || confess 'add_xref needs a species_id';
  my $label        = $arg_ref->{label}      // $acc;
  my $description  = $arg_ref->{desc};
  my $version      = $arg_ref->{version}    // 0;
  my $info_type    = $arg_ref->{info_type}  // 'MISC';
  my $info_text    = $arg_ref->{info_text}  // q{};
  my $update_label = $arg_ref->{update_label};
  my $update_desc  = $arg_ref->{update_desc};

  # See if it already exists. It so return the xref_id for this one.
  my $xref_id = $self->get_xref( $acc, $source_id, $species_id );
  if ( defined $xref_id ) {
    if ($update_label) {
      $self->_update_xref_label( $xref_id, $label );
    }
    if ($update_desc) {
      $self->_update_xref_description( $xref_id, $description );
    }
    return $xref_id;
  }

  my $sql = (<<'SQL');
  INSERT INTO xref (
    accession, version, label, description,
    source_id, species_id, info_type, info_text)
  VALUES (?,?,?,?,?,?,?,?)
SQL

  my $add_xref_sth = $self->dbi->prepare_cached($sql);

  # If the description is more than 255 characters, chop it off and add
  # an indication that it has been truncated to the end of it.
  if ( defined $description && ( ( length $description ) > 255 ) ) {
    my $truncmsg = ' /.../';
    substr $description, 255 - ( length $truncmsg ), length $truncmsg,
      $truncmsg;
  }

  # Add the xref and confess if it fails
  $add_xref_sth->execute( $acc,         $version // 0, $label,
                          $description, $source_id,    $species_id,
                          $info_type,   $info_text ) or
    confess $self->dbi->errstr() .
    "\n $acc\t$label\t\t$source_id\t$species_id\n";

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

  my $xref_id = $arg_ref->{xref_id} ||
    confess 'add_object_xref needs an xref_id';
  my $ensembl_id = $arg_ref->{ensembl_id} ||
    confess 'add_object_xref needs an ensembl_id';
  my $object_type = $arg_ref->{object_type} ||
    confess 'add_object_xref needs an object_type';

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
    confess $self->dbi->errstr() .
    "\n $ensembl_id\t$object_type\t\t$xref_id\n";

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

  my $object_xref_id = $arg_ref->{object_xref_id}   || confess 'add_identity_xref needs an object_xref_id';
  my $score = $arg_ref->{score}                     || confess 'add_identity_xref needs a score';
  my $target_identity = $arg_ref->{target_identity} || confess 'add_identity_xref needs a target_identity';
  my $query_identity = $arg_ref->{query_identity}   || confess 'add_identity_xref needs a query_identity';

  my $sql = (<<'SQL');
    INSERT INTO identity_xref
      (object_xref_id, score, query_identity, target_identity)
    VALUES (?,?,?,?)
SQL

  my $add_identity_xref_sth =
    $self->dbi->prepare_cached( $sql );

  # Add the object_xref and confess if it fails
  $add_identity_xref_sth->execute( $object_xref_id, $score,
                                   $query_identity, $target_identity
    ) or
    confess( $self->dbi->errstr() .
      "\n$object_xref_id\t$score\t\t$query_identity\t$target_identity\n"
    );

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
    gene        => 'INSERT INTO gene_direct_xref VALUES (?,?,?)',
    transcript  => 'INSERT INTO transcript_direct_xref VALUES (?,?,?)',
    translation => 'INSERT INTO translation_direct_xref VALUES (?,?,?)',
  );

  my $add_direct_xref_sth =
    $self->dbi->prepare_cached( $sql_hash{ lc $ensembl_type } );

  $add_direct_xref_sth->execute( $general_xref_id, $ensembl_stable_id,
                                 $linkage_type );

  if ( defined $update_info_type ) {
    $self->_update_xref_info_type( $general_xref_id, 'DIRECT' );
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

  foreach my $direct_xref ( @{$direct_xrefs} ) {
    $self->add_to_direct_xrefs(
          { stable_id  => $direct_xref->{STABLE_ID},
            type       => $direct_xref->{ENSEMBL_TYPE},
            linkage    => $direct_xref->{LINKAGE_TYPE},
            acc        => $direct_xref->{ACCESSION},
            version    => $direct_xref->{VERSION} // 0,
            label      => $direct_xref->{LABEL} // $direct_xref->{ACCESSION},
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

  my $master_xref = $arg_ref->{master_xref_id} ||
    confess('Need a master_xref_id on which this xref linked too');
  my $acc = $arg_ref->{acc} ||
    confess('Need an accession of this dependent xref');
  my $source_id = $arg_ref->{source_id} ||
    confess('Need a source_id for this dependent xref');
  my $species_id = $arg_ref->{species_id} ||
    confess('Need a species_id for this dependent xref');
  my $version   = $arg_ref->{version} // 0;
  my $label     = $arg_ref->{label} // $acc;
  my $desc      = $arg_ref->{desc};
  my $linkage   = $arg_ref->{linkage};
  my $info_text = $arg_ref->{info_text} // q{};

  my $dependent_xref_id =
    $self->add_xref( {
      acc        => $acc,
      source_id  => $source_id,
      species_id => $species_id,
      label      => $label,
      desc       => $desc,
      version    => $version,
      info_type  => 'DEPENDENT',
      info_text  => $info_text } );

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

  my $sql = (<<'SQL');
    INSERT INTO dependent_xref
      (master_xref_id, dependent_xref_id, linkage_annotation, linkage_source_id)
    VALUES (?,?,?,?)
SQL
  my $add_dependent_xref_sth = $self->dbi->prepare_cached($sql);

  # If the dependency cannot be found in %xref_dependent_mapped,
  # i.e. has not been set yet, add it
  if (( !defined $xref_dependent_mapped{"$master_id|$dependent_id"} ) ||
      $xref_dependent_mapped{"$master_id|$dependent_id"} ne
      $master_source_id )
  {

    $add_dependent_xref_sth->execute( $master_id, $dependent_id,
                            $master_source_id, $dependent_source_id ) ||
      confess( $self->dbi->errstr() .
"\n$master_id\t$dependent_id\t$master_source_id\t$dependent_source_id"
      );

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

  foreach my $depref ( @{$dependent_xrefs} ) {
    my %dep = %{$depref};

    # Insert the xref
    my $dep_xref_id = $self->add_xref( {
      'acc'        => $dep{ACCESSION},
      'version'    => $dep{VERSION}     // 1,
      'label'      => $dep{LABEL}       // $dep{ACCESSION},
      'desc'       => $dep{DESCRIPTION},
      'source_id'  => $dep{SOURCE_ID},
      'species_id' => $dep{SPECIES_ID},
      'info_type'  => 'DEPENDENT' } );

    # Add the linkage_annotation and source id it came from
    $self->add_dependent_xref_maponly(
                                  $dep_xref_id, $dep{LINKAGE_SOURCE_ID},
                                  $xref_id,     $dep{LINKAGE_ANNOTATION}
    );

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
    my $xref_id = $self->get_xref( $acc, $source_id, $species_id );
    if ( defined $xref_id ) {
      $self->add_synonym( $xref_id, $syn );
    }
  }

  return;
}

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
}

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
      'INSERT IGNORE INTO synonym ( xref_id, synonym ) VALUES( ?, ? )');
  $add_synonym_sth->execute( $xref_id, $syn ) or
    confess( $self->dbi->errstr() .
         "\n $xref_id\n $syn\n\n" . $add_synonym_sth->errstr() . "\n" );

  return;
}    ## sub add_synonym

=head2 add_multiple_synonyms
  Arg [1]    : xref ID
  Arg [2]    : Listref : synonyms
  Description: Add multiple synomyns for an xref given by xref_id
  Return type:
  Caller     : internal

=cut

sub add_multiple_synonyms {
  my ( $self, $xref_id, $synonyms ) = @_;

  foreach my $syn ( @{$synonyms} ) {
    $self->add_synonym( $xref_id, $syn );
  }

  return;
}    ## sub add_multiple_synonyms

=head2 add_synonyms_for_hgnc_vgnc
  Arg [1]    : hashref : source_id, name, species_id, dead, alias
  Description: Specialized class to add synonyms from HGNC and VGNC data
  Return type: N/A
  Caller     : internal
=cut

sub add_synonyms_for_hgnc_vgnc {
  my ( $self, $ref_arg ) = @_;

  my $source_id    = $ref_arg->{source_id};
  my $name         = $ref_arg->{name};
  my $species_id   = $ref_arg->{species_id};
  my $dead_string  = $ref_arg->{dead};
  my $alias_string = $ref_arg->{alias};

  # dead name, add to synonym
  if ( defined $dead_string ) {
    $dead_string =~ s/"//xg;
    my @dead_array = split( ',\s', $dead_string );
    foreach my $dead (@dead_array) {
      $self->add_to_syn( $name, $source_id, $dead, $species_id );
    }
  }

  # alias name, add to synonym
  if ( defined $alias_string ) {
    $alias_string =~ s/"//xg;
    my @alias_array = split( ',\s', $alias_string );
    foreach my $alias (@alias_array) {
      $self->add_to_syn( $name, $source_id, $alias, $species_id );
    }
  }

  return;
} ## end sub add_synonyms_for_hgnc_vgnc


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

  my $sql = (<<'SQL');
    SELECT  xref.accession, xref.label
    FROM xref, source
    WHERE
      source.name LIKE ? AND
      xref.source_id = source.source_id
SQL

  my @sql_params = ($name . q{%});

  if ( defined $prio_desc ) {
    $sql .= ' AND source.priority_description LIKE ?';
    push @sql_params, $prio_desc;
  }
  if ( defined $species_id ) {
    $sql .= ' AND xref.species_id  = ?';
    push @sql_params, $species_id;
  }
  my $sub_sth = $self->dbi->prepare_cached($sql);

  $sub_sth->execute(@sql_params);
  while ( my @row = $sub_sth->fetchrow_array() ) {
    $hash1{ $row[0] } = $row[1];
  }

  return \%hash1;
} ## end sub get_acc_to_label


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
}

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

  my $sql = (<<'SQL');
    SELECT  d.master_xref_id, d.dependent_xref_id, d.linkage_annotation
    FROM dependent_xref d, xref x
    WHERE
      x.xref_id = d.dependent_xref_id AND
      x.source_id = ?
SQL
  my $sth = $self->dbi->prepare_cached($sql);
  $sth->execute($source_id);
  my $master_xref;
  my $dependent_xref;
  my $linkage;
  $sth->bind_columns( \$master_xref, \$dependent_xref, \$linkage );
  while ( $sth->fetch ) {
    $xref_dependent_mapped{"$master_xref|$dependent_xref"} = $linkage;
  }

  return;
}

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

  my $sql = (<<'SQL');
    SELECT  x.accession, x.label, sy.synonym
    FROM xref x, source so, synonym sy
    WHERE x.xref_id = sy.xref_id AND
      so.source_id = x.source_id AND
      so.name LIKE ?
SQL
  my $sth = $self->dbi->prepare_cached($sql);

  $sth->execute( $source_name );
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


=head2 get_meta_value

  Arg [1]    : key
  Description: Return metadata value
  Return type: integer or string
  Exceptions : If argument is not provided
  Caller     : Bio::EnsEMBL::Xref::Mapper

=cut

sub get_meta_value {
  my ( $self, $key ) = @_;
  if ( !defined $key ) {
    confess 'Need to specify the key';
  }

  my $sth = $self->dbi->prepare_cached(
    'SELECT meta_value FROM meta WHERE meta_key LIKE ? ORDER BY meta_id'
  );
  $sth->execute($key);
  my $value;
  $sth->bind_columns( \$value );
  while ( $sth->fetch ) {    # get the last one
  }

  return $value;
}

=head2 update_process_status

  Arg [1]    : Scalar; the new status
  Description: Update process status
  Returntype : none
  Exceptions : If insertion fails or argument is not provided

=cut

sub update_process_status {
  my ( $self, $value ) = @_;
  if ( !defined $value ) {
    confess 'Need to specify a value';
  }

  my $sth = $self->dbi->prepare_cached(
    'INSERT INTO process_status (status, date) VALUES(?, NOW())');
  $sth->execute($value) or confess $self->dbi->errstr() . "\n $value\n";

  return;
}


=head2 clean_up

=cut

sub clean_up {
  my ( $self, $stats, $keep_core_data ) = @_;

  # remove all object_xref, identity_xref  entries
  my $sth = $self->dbi->prepare_cached('TRUNCATE table object_xref');
  $sth->execute();

  $sth = $self->dbi->prepare_cached('TRUNCATE table go_xref');
  $sth->execute();

  $sth = $self->dbi->prepare_cached('TRUNCATE table identity_xref');
  $sth->execute();

  # remove all xrefs after PARSED_xref_id
  # set dumped to NULL fro all xrefs.
  my $max_xref_id = $self->get_meta_value('PARSED_xref_id');
  if ($max_xref_id) {
    $sth = $self->dbi->prepare_cached(
                       "DELETE from xref where xref_id > $max_xref_id");
    $sth->execute();
  }

  $sth = $self->dbi->prepare_cached('UPDATE xref set dumped = null');
  $sth->execute();

  $sth = $self->dbi->prepare_cached('TRUNCATE display_xref_priority');
  $sth->execute();

  $sth = $self->dbi->prepare_cached('TRUNCATE gene_desc_priority');
  $sth->execute();

  unless ($keep_core_data) {
    # remove all from core_info tables
    #   gene_transcript_translation
    #   [gene/transcript/translation]_stable_id
    #
    $sth = $self->dbi->prepare_cached(
                                'TRUNCATE gene_transcript_translation');
    $sth->execute();

    $sth = $self->dbi->prepare_cached('TRUNCATE gene_stable_id');
    $sth->execute();

    $sth = $self->dbi->prepare_cached('TRUNCATE transcript_stable_id');
    $sth->execute();

    $sth = $self->dbi->prepare_cached('TRUNCATE translation_stable_id');
    $sth->execute();
  }

  return;
} ## end sub clean_up

=head2 remove_mapping_data

=cut

sub remove_mapping_data {
  my $self = shift;

  my $sth = $self->dbi->prepare_cached('TRUNCATE mapping_jobs');
  $sth->execute();

  $sth = $self->dbi->prepare_cached('TRUNCATE mapping');
  $sth->execute();

  $sth = $self->dbi->prepare_cached('TRUNCATE alt_allele');
  $sth->execute();

  $sth = $self->dbi->prepare_cached('TRUNCATE source_mapping_method');
  $sth->execute();

  return;
}

=head2 update_mapping_jobs_status

=cut

sub update_mapping_jobs_status {
  my ( $self, $status ) = @_;
  confess 'Status not given' unless defined $status;

  my $sth =
    $self->dbi->prepare_cached('UPDATE mapping_jobs set status = ?');
  $sth->execute($status) or
    confess $self->dbi->errstr() . "\n $status\n";

  return;
}

=head2 process_alt_alleles

=cut

sub process_alt_alleles {
  my $self = shift;

# ALL are on the Gene level now. This may change but for now it is okay.
  my ( $alt_to_ref, $ref_to_alts ) = $self->_get_alt_allele_hashes();

#
# Move the xrefs on to the reference Gene.
# NOTE: Igonore used as the xref might already be on this Gene already and we do not want it to crash
#
  my $move_sql = (<<'MOVE');
UPDATE IGNORE object_xref ox, xref x, source s
  SET ox.ensembl_id = ?
    WHERE x.source_id = s.source_id AND
          ox.xref_id = x.xref_id AND
          ox.ensembl_id = ? AND
          ox.ensembl_object_type = 'Gene' AND
          ox.ox_status = 'DUMP_OUT' AND
          s.name in (
MOVE
  $move_sql .=
    "'" . join( "', '", @{ $self->get_gene_specific_list() } ) . "')";

  print "MOVE SQL\n$move_sql\n";

#
# Now where it was already on the Gene the ignore will have stopped the move
# so we now want to just remove those ones as they already exist.
#
  my $del_ix_sql = (<<'DIX');
DELETE ix
  FROM identity_xref ix, object_xref ox, xref x, source s
    WHERE x.source_id = s.source_id AND
          ox.object_xref_id = ix.object_xref_id AND
          ox.xref_id = x.xref_id AND
          ox.ensembl_id = ? AND
          ox.ensembl_object_type = 'Gene' AND
          ox.ox_status = 'DUMP_OUT' AND
           s.name in (
DIX
  $del_ix_sql .=
    "'" . join( "', '", @{ $self->get_gene_specific_list() } ) . "')";

  my $del_sql = (<<'DEL');
DELETE ox
  FROM object_xref ox, xref x, source s
    WHERE x.source_id = s.source_id AND
          ox.xref_id = x.xref_id AND
          ox.ensembl_id = ? AND
          ox.ensembl_object_type = 'Gene' AND
          ox.ox_status = 'DUMP_OUT' AND
           s.name in (
DEL
  $del_sql .=
    "'" . join( "', '", @{ $self->get_gene_specific_list() } ) . "')";

  my $move_sth = $self->dbi->prepare_cached($move_sql) or
    confess "$move_sql cannot be prepared";
  my $del_ix_sth = $self->dbi->prepare_cached($del_ix_sql) or
    confess "$del_ix_sql cannot be prepared";
  my $del_sth = $self->dbi->prepare_cached($del_sql) or
    confess "$del_sql cannot be prepared";

  my ( $move_count, $del_ix_count, $del_ox_count ) = ( 0, 0, 0 );
  foreach my $key ( keys %$alt_to_ref ) {
    $move_sth->execute( $alt_to_ref->{$key}, $key );
    $move_count += $move_sth->rows;

    $del_ix_sth->execute($key);
    $del_ix_count += $del_ix_sth->rows;

    $del_sth->execute($key);
    $del_ox_count += $del_sth->rows;
  }
  $move_sth->finish;
  $del_sth->finish;
  $del_ix_sth->finish;

  print
"Number of rows:- moved = $move_count, identitys deleted = $del_ix_count, object_xrefs deleted = $del_ox_count\n";
#
# Now we have all the data on the reference Gene we want to copy all the data
# onto the alt alleles.
#

  my $get_data_sql = (<<'GET');
SELECT ox.object_xref_id, ox.ensembl_object_type, ox.xref_id, ox.linkage_annotation,
       ox.linkage_type, ox.ox_status, ox.unused_priority, ox.master_xref_id,
       ix.query_identity, ix.target_identity, ix.hit_start, ix.hit_end,
       ix.translation_start, ix.translation_end, ix.cigar_line, ix.score, ix.evalue
  FROM xref x, source s, object_xref ox
    LEFT JOIN identity_xref ix ON ox.object_xref_id =ix.object_xref_id
      WHERE  x.source_id = s.source_id AND
             ox.xref_id = x.xref_id AND
             ox.ensembl_id = ? AND
             ox.ox_status = 'DUMP_OUT' AND
             ox.ensembl_object_type = 'Gene' AND
              s.name in (
GET
  $get_data_sql .=
    "'" . join( "', '", @{ $self->get_gene_specific_list() } ) . "')";
  my $get_data_sth = $self->dbi->prepare_cached($get_data_sql) or
    confess "Could not prepare $get_data_sql";

  my $insert_object_xref_sql = (<<'INO');
INSERT INTO object_xref (object_xref_id, ensembl_id, ensembl_object_type, xref_id, linkage_annotation,
            linkage_type, ox_status, unused_priority, master_xref_id)
       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
INO
  my $insert_ox_sth =
    $self->dbi->prepare_cached($insert_object_xref_sql) or
    confess "Could not prepare $insert_object_xref_sql";

  my $insert_identity_xref_sql = (<<'INI');
INSERT INTO identity_xref (object_xref_id, query_identity, target_identity, hit_start, hit_end,
            translation_start, translation_end, cigar_line, score, evalue )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
INI
  my $insert_ix_sth =
    $self->dbi->prepare_cached($insert_identity_xref_sql) or
    confess "Could not prepare $insert_identity_xref_sql";

  my $max_object_xref_id;
  my $sth = $self->dbi->prepare_cached(
                         'SELECT MAX(object_xref_id) FROM object_xref');
  $sth->execute();
  $sth->bind_columns( \$max_object_xref_id );
  $sth->fetch;
  confess 'Problem getting max object_xref_id'
    if not defined $max_object_xref_id or
    not $max_object_xref_id;

  $max_object_xref_id++;

  my ( $added_count, $ignored ) = ( 0, 0 );
  foreach my $key ( keys %$ref_to_alts ) {
    $get_data_sth->execute($key);
    my ( $object_xref_id,    $ensembl_object_type,
         $xref_id,           $linkage_annotation,
         $linkage_type,      $ox_status,
         $unused_priority,   $master_xref_id,
         $query_identity,    $target_identity,
         $hit_start,         $hit_end,
         $translation_start, $translation_end,
         $cigar_line,        $score,
         $evalue );

    $get_data_sth->bind_columns( \$object_xref_id,
                                 \$ensembl_object_type,
                                 \$xref_id,
                                 \$linkage_annotation,
                                 \$linkage_type,
                                 \$ox_status,
                                 \$unused_priority,
                                 \$master_xref_id,
                                 \$query_identity,
                                 \$target_identity,
                                 \$hit_start,
                                 \$hit_end,
                                 \$translation_start,
                                 \$translation_end,
                                 \$cigar_line,
                                 \$score,
                                 \$evalue );

    while ( $get_data_sth->fetch() ) {
      foreach my $alt ( @{ $ref_to_alts->{$key} } ) {
        $max_object_xref_id++;
        $insert_ox_sth->execute( $max_object_xref_id,  $alt,
                                 $ensembl_object_type, $xref_id,
                                 $linkage_annotation,  $linkage_type,
                                 $ox_status,           $unused_priority,
                                 $master_xref_id ) or
          confess 'Could not insert object_xref data';

        # ONLY add identity xref if object_xref was added successfully.
        if ( $insert_ox_sth->rows ) {
          $added_count++;
          $insert_ix_sth->execute( $max_object_xref_id,
                                   $query_identity,
                                   $target_identity,
                                   $hit_start,
                                   $hit_end,
                                   $translation_start,
                                   $translation_end,
                                   $cigar_line,
                                   $score,
                                   $evalue ) or
            confess 'Could not insert identity_xref data';
        }
        else {
          $ignored++;
        }
      } ## end foreach my $alt ( @{ $ref_to_alts...})
    } ## end while ( $get_data_sth->fetch...)
  } ## end foreach my $key ( keys %$ref_to_alts)

  return ( $added_count, $ignored );
} ## end sub process_alt_alleles

#
# These sources should be on the gene, even if they are mapped transcript or translation.
# We define which ones are to be moved here
#

=head2 get_gene_specific_list

=cut

sub get_gene_specific_list {
  my $self = shift;

  my @list =
    qw(DBASS3 DBASS5 EntrezGene miRBase RFAM TRNASCAN_SE RNAMMER UniGene Uniprot_gn WikiGene MIM_GENE MIM_MORBID HGNC MGI ZFIN_ID FlyBaseName_gene RGD SGD_GENE VGNC wormbase_gseqname wormbase_locus Xenbase);

  # Check the sources are used in the database considered
  my ( @used_list, $sql, $sth, $count );
  foreach my $source (@list) {
    $sql =
"SELECT COUNT(*) FROM xref x, source s WHERE s.source_id = x.source_id AND s.name = ?";
    $sth = $self->dbi->prepare_cached($sql);
    $sth->execute($source);

    $sth->bind_columns( \$count );
    $sth->fetch();

    push @used_list, $source if $count > 0;
  }

  return \@used_list;
}

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
}

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
'INSERT INTO pairs (source_id, accession1, accession2) VALUES(?,?,?)' );

  # Add the pair and confess if it fails
  $pair_sth->execute( $source_id, $accession, $pair ) or
    confess $self->dbi->errstr() .
    "\n $source_id\t$\t$accession\t$pair\n";

  return;
}

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

  my $add_primary_xref_sth = $self->dbi->prepare_cached(
                           'INSERT INTO primary_xref VALUES (?,?,?,?)');

  # Add the xref and confess if it fails
  $add_primary_xref_sth->execute( $xref_id,       $sequence,
                                  $sequence_type, $status ) or
    confess $self->dbi->errstr() .
    "\n $xref_id\t$sequence_type\t$status\n";

  return $add_primary_xref_sth->{'mysql_insertid'};
}

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
}    ## sub _update_primary_xref_sequence

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
}    ## sub _update_xref_label

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
}    ## sub _update_xref_description

#
# In case we have alt alleles with xefs, these will be direct ones
# we need to move all xrefs on to the reference
#

=head2 _get_alt_allele_hashes

=cut

sub _get_alt_allele_hashes {
  my $self = shift;

  my %alt_to_ref;
  my %ref_to_alts;

  my $sql = (<<'SQL');
    SELECT alt_allele_id, gene_id, is_reference
    FROM alt_allele
    ORDER BY alt_allele_id, is_reference DESC
SQL

  my $sth = $self->dbi->prepare_cached( $sql );
  $sth->execute();

  my ( $alt_allele_id, $gene_id, $is_ref );
  $sth->bind_columns( \$alt_allele_id, \$gene_id, \$is_ref );

  my $last_alt_allele = 0;
  my $ref_gene;
  while ( $sth->fetch() ) {
    if ( $alt_allele_id != $last_alt_allele ) {
# use the first non-reference gene if there is no reference gene in an alt_allele
      $ref_gene = $gene_id;
    }
    else {
      $alt_to_ref{$gene_id} = $ref_gene;
      push @{ $ref_to_alts{$ref_gene} }, $gene_id;
    }
    $last_alt_allele = $alt_allele_id;
  }
  $sth->finish;

  return \%alt_to_ref, \%ref_to_alts;
} ## end sub _get_alt_allele_hashes


=head2 _update_xref_description
  Description: Return a map from RefSeq prefix symbols to their corresponding
               sources.
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


=head2 get_valid_source_id_to_external_db_id
  Description: Create a hash of all the external db names and ids that have
               associated xrefs
  Return type: Hashref
  Caller     : internal

=cut

sub get_valid_source_id_to_external_db_id {

  my $self = shift;

  my %source_id_to_external_db_id;

  my $sql = (<<'SQL');
    SELECT s.source_id, s.name
    FROM source s, xref x
    WHERE x.source_id = s.source_id
    GROUP BY s.source_id
SQL

  my $sth = $self->dbi->prepare_cached( $sql );
  $sth->execute() or confess( $self->dbi->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $source_name = $row[0];
    my $source_id   = $row[1];
    $source_id_to_external_db_id{ $source_name } = $source_id;
  }


  return %source_id_to_external_db_id;
} ## end sub get_valid_source_id_to_external_db_id


=head2 get_source_ids_with_xrefs
  Description: Get all sources with xrefs
  Return type: Iterator - Hashref
  Example    : my $sources = $ba_handle->get_source_ids_with_sources();
               while( my $source_ref = $sources->() ) {
                 my %source = %{ $source_ref };
                 print "$source{'name'} - $source{'count'}";
               }

=cut

sub get_source_ids_with_xrefs {
  my $self = shift;

  my $sql = (<<'SQL');
    SELECT s.name, COUNT(s.name)
    FROM xref x, object_xref ox, source s
    WHERE ox.xref_id = x.xref_id AND
          x.source_id = s.source_id
    GROUP BY s.name
SQL

  my $sth = $self->dbi->prepare_cached( $sql );
  $sth->execute() or confess( $self->dbi->errstr() );

  my ( $source_name, $source_count );
  $sth->bind_columns( \$source_name, \$source_count);

  return sub {
    if ( $sth->fetch() ) {
      return {
        name  => $source_name,
        count => $source_count
      }
    }
  }
} ## end sub get_source_ids_with_xrefs


=head2 get_dump_out_xrefs
  Description: Get all sources that have associated xrefs
  Return type: Iterator - Hashref
  Example    : my $sources = $ba_handle->get_dump_out_sources();
               while( my $source_ref = $sources->() ) {
                 my %source = %{ $source_ref };
                 print "$source{'name'} - $source{'count'}";
               }

=cut

sub get_dump_out_xrefs {
  my $self = shift;

  my $sql = (<<'SQL');
    SELECT s.name, s.source_id, COUNT(*), x.info_type, s.priority_description, s.source_release
    FROM xref x, object_xref ox, source s
    WHERE ox.xref_id = x.xref_id  AND
          x.source_id = s.source_id AND
          ox_status = 'DUMP_OUT'
    GROUP BY s.name, s.source_id, x.info_type
SQL

  my $sth = $self->dbi->prepare_cached( $sql );
  $sth->execute();

  my ( $name, $type, $source_id, $count, $where_from, $release_info );
  $sth->bind_columns(\$name,\$source_id, \$count, \$type, \$where_from, \$release_info);

  return sub {
    if ( $sth->fetch() ) {

      return {
        type         => $type,
        source_id    => $source_id,
        count        => $count,
        name         => $name,
        where_from   => $where_from,
        release_info => $release_info,
      }
    }
  }
}


=head2 get_insert_identity_xref
  Arg [1]    : integer - $source_id
  Arg [2]    : string - $type
  Description: Get all identity xrefs
  Return type: Iterator - Hashref
  Example    : my $xref_identities = $ba_handle->get_insert_identity_xref();
               while( my $xref_identities_ref = $xref_identities->() ) {
                 my %xref_identity = %{ $xref_identities_ref };
                 print "$xref_identity{'xref_id'} - $xref_identity{'acc'}";
               }

=cut

sub get_insert_identity_xref {
  my ( $self, $source_id, $type ) = @_;

  my $sql = (<<'SQL');
    SELECT x.xref_id, x.accession, x.label, x.version, x.description, x.info_text,
           ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type,
           i.query_identity, i.target_identity, i.hit_start, i.hit_end,
           i.translation_start, i.translation_end, i.cigar_line, i.score, i.evalue
    FROM xref x, object_xref ox, identity_xref i
    WHERE ox.ox_status = 'DUMP_OUT' AND
          i.object_xref_id = ox.object_xref_id AND
          ox.xref_id = x.xref_id AND
          x.source_id = ? AND
          x.info_type = ?
    ORDER BY x.xref_id
SQL

  my $sth = $self->dbi->prepare_cached($sql);

  my $count = 0;
  $sth->execute($source_id, $type);

  my ( $xref_id, $acc, $label, $version, $desc, $info, $object_xref_id,
       $ensembl_id, $ensembl_type );

  my ( $query_identity, $target_identity, $hit_start, $hit_end,
       $translation_start, $translation_end, $cigar_line, $score, $evalue);

  $sth->bind_columns(
    \$xref_id, \$acc, \$label, \$version, \$desc, \$info, \$object_xref_id,
    \$ensembl_id, \$ensembl_type, \$query_identity, \$target_identity,
    \$hit_start, \$hit_end, \$translation_start, \$translation_end,
    \$cigar_line, \$score, \$evalue);

  return sub {
    if ( $sth->fetch() ) {
      return {
        xref_id          => $xref_id,
        acc              => $acc,
        label            => $label,
        version          => $version,
        desc             => $desc,
        info             => $info,
        ensembl_id       => $ensembl_id,
        ensembl_type     => $ensembl_type,
        object_xref_id   => $object_xref_id,
        query_identity   => $query_identity,
        ensembl_identity => $target_identity,
        xref_start       => $hit_start,
        xref_end         => $hit_end,
        ensembl_start    => $translation_start,
        ensembl_end      => $translation_end,
        cigar_line       => $cigar_line,
        score            => $score,
        evalue           => $evalue
      }
    }
  }
} ## end sub get_insert_identity_xref


=head2 get_insert_checksum_xref
  Arg [1]    : integer - $source_id
  Arg [2]    : string - $type
  Description: Get all CHECKSUM xrefs
  Return type: Iterator - Hashref
  Example    : my $xref_checksums = $ba_handle->get_insert_checksum_xref();
               while( my $xref_checksums_ref = $xref_checksums->() ) {
                 my %xref_checksum = %{ $xref_checksums_ref };
                 print "$xref_checksum{'xref_id'} - $xref_checksum{'acc'}";
               }

=cut

sub get_insert_checksum_xref {
  my ( $self, $source_id, $type ) = @_;

  my $sql = (<<'SQL');
    SELECT x.xref_id, x.accession, x.label, x.version, x.description, x.info_text,
           ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type
    FROM xref x, object_xref ox
    WHERE ox.ox_status = 'DUMP_OUT' AND
          ox.xref_id = x.xref_id AND
          x.source_id = ? AND
          x.info_type = ?
    ORDER BY x.xref_id
SQL

  my $sth = $self->dbi->prepare( $sql );

  my $count = 0;
  $sth->execute($source_id, $type);
  my $last_xref = 0;
  my ($xref_id, $acc, $label, $version, $desc, $info, $object_xref_id, $ensembl_id, $ensembl_type);

  $sth->bind_columns(
    \$xref_id, \$acc, \$label, \$version, \$desc, \$info, \$object_xref_id,
    \$ensembl_id, \$ensembl_type);

  return sub {
    if ( $sth->fetch() ) {
      return {
        xref_id          => $xref_id,
        acc              => $acc,
        label            => $label,
        version          => $version,
        desc             => $desc,
        info             => $info,
        ensembl_id       => $ensembl_id,
        ensemb_type      => $ensembl_type,
        object_xref_id   => $object_xref_id
      }
    }
  }
} ## end sub get_insert_checksum_xref


=head2 get_insert_dependent_xref
  Arg [1]    : integer - $source_id
  Arg [2]    : string - $type
  Description: Get all DEPENDENT xrefs
  Return type: Iterator - Hashref
  Example    : my $xref_dependents = $ba_handle->get_insert_dependent_xref();
               while( my $xref_dependents_ref = $xref_dependents->() ) {
                 my %xref_dependent = %{ $xref_dependents_ref };
                 print "$xref_dependent{'xref_id'} - $xref_dependent{'acc'}";
               }

=cut

sub get_insert_dependent_xref {
  my ( $self, $source_id, $type ) = @_;

  my $sql = (<<'SQL');
    SELECT x.xref_id, x.accession, x.label, x.version, x.description, x.info_text,
           ox.object_xref_id, ox.ensembl_id, ox.ensembl_object_type, ox.master_xref_id
    FROM xref x, object_xref ox
    WHERE ox.ox_status = 'DUMP_OUT' AND
          ox.xref_id = x.xref_id AND
          x.source_id = ? AND
          x.info_type = ?
    ORDER BY x.xref_id, ox.ensembl_id
SQL

  my $sth = $self->dbi->prepare( $sql );
  $sth->execute($source_id, $type);

  my ( $xref_id, $acc, $label, $version, $desc, $info, $object_xref_id,
       $ensembl_id, $ensembl_type, $master_xref_id );
  $sth->bind_columns(
    \$xref_id, \$acc, \$label, \$version, \$desc, \$info, \$object_xref_id,
    \$ensembl_id, \$ensembl_type, \$master_xref_id);

  return sub {
    if ( $sth->fetch() ) {
      return {
        xref_id          => $xref_id,
        acc              => $acc,
        label            => $label,
        version          => $version,
        desc             => $desc,
        info             => $info,
        ensembl_id       => $ensembl_id,
        ensembl_type     => $ensembl_type,
        object_xref_id   => $object_xref_id,
        master_xref_id   => $master_xref_id
      }
    }
  }

} ## end sub get_insert_dependent_xref


=head2 get_synonyms_for_xref
  Arg [1]    : Arrayref - \@xref_list
  Description: For a list of xref IDs retrieve the list of matching synonyms.
  Return type: Iterator - Hashref
  Example    : my $xref_synonyms = $ba_handle->get_synonyms_for_xref();
               while( my $xref_synonyms_ref = $xref_synonyms->() ) {
                 my %xref_synonym = %{ $xref_synonyms_ref };
                 print "$xref_synonym{'xref_id'} - $xref_synonym{'syn'}";
               }

=cut

sub get_synonyms_for_xref {
  my ( $self, $xref_list ) = @_;

  my ($xref_id, $syn);
  my $syn_sql = (<<"SQL");
    SELECT xref_id, synonym
    FROM synonym
    WHERE xref_id IN ( @{ [ join',', ('?') x @{ $xref_list } ] } )
SQL
  my $sth = $self->dbi->prepare( $syn_sql );
  $sth->execute( @{ $xref_list } );
  $sth->bind_columns(\$xref_id, \$syn);

  return sub {
    if ( $sth->fetch() ) {
      return {
        xref_id => $xref_id,
        syn     => $syn
      }
    }
  }
} ## end sub get_synonyms_for_xref


# just incase this is being run again
=head2 mark_mapped_xrefs_already_run
  Description: Set dumped to NULL if the value in dumped is not set to
              NO_DUMP_ANOTHER_PRIORITY
  Return type: undef

=cut

sub mark_mapped_xrefs_already_run {
  my ( $self ) = @_;

  my $sql = (<<"SQL");
    UPDATE xref
    SET dumped = NULL
    WHERE dumped != 'NO_DUMP_ANOTHER_PRIORITY'
SQL

  my $sth = $self->dbi->prepare( $sql );
  $sth->execute() || confess 'Could not set dumped status';

  return;
} ## end sub mark_mapped_xrefs_already_run


=head2 mark_mapped_xrefs
  Arg [1]    : Arrayref - \@xref_list
  Arg [2]    : string - $status
  Description: Set the dumped column for a list of xref IDs with $status
  Return type: undef

=cut

sub mark_mapped_xrefs {
  my ( $self, $xref_list, $status ) = @_;

  my $sql = (<<"SQL");
    UPDATE xref
    SET dumped = ?
    WHERE xref_id IN ( @{[join',', ('?') x @{$xref_list}]} )
SQL

  my $xref_dumped_sth = $self->dbi->prepare( $sql );
  $xref_dumped_sth->execute( $status, @{ $xref_list } ) ||
    confess 'Could not set dumped status';

  return;
} ## end sub mark_mapped_xrefs


=head2 insert_process_status
  Arg [1]    : string - $status
  Description: Insert the current status for a process
  Return type: undef

=cut

sub insert_process_status {
  my ( $self, $status ) = @_;

  my $sql = (<<"SQL");
    INSERT INTO process_status (status, date) VALUES ( ?, now() )
SQL

  my $sth = $self->dbi->prepare( $sql );
  $sth->execute( $status ) || confess 'Could not set dumped status';

  return;
} ## end sub insert_process_status


=head2 get_unmapped_reason
  Description: Get the names and descriptions of valid reasons why an xref could
               could be unmapped. These shoudld be used to describe the xrefs
               when they are loaded into the core DB.
  Return type: Hashref

=cut

sub get_unmapped_reason {
  my ( $self ) = @_;

  my %summary_failed;
  my %desc_failed;

  # Get the cutoff values
  my $sql = (<<'SQL');
    SELECT DISTINCT
      s.name, m.percent_query_cutoff, m.percent_target_cutoff
    FROM source s,
         source_mapping_method sm,
         mapping m
    WHERE sm.source_id = s.source_id AND
          sm.method = m.method
SQL

  my $sth = $self->dbi->prepare_cached( $sql );
  $sth->execute();
  my ( $source_name, $q_cut, $t_cut );
  $sth->bind_columns( \$source_name, \$q_cut, \$t_cut );

  while ( $sth->fetch ) {
    $summary_failed{$source_name} = 'Failed to match at thresholds';
    $desc_failed{$source_name}    = "Unable to match at the thresholds of $q_cut\% for the query or $t_cut\% for the target";
  }

  $summary_failed{'NO_STABLE_ID'} = 'Failed to find Stable ID';
  $desc_failed{'NO_STABLE_ID'}    = 'Stable ID that this xref was linked to no longer exists';

  $summary_failed{'FAILED_MAP'} = 'Failed to match';
  $desc_failed{'FAILED_MAP'}    = 'Unable to match to any ensembl entity at all';

  $summary_failed{'NO_MAPPING'} = 'No mapping done';
  $desc_failed{'NO_MAPPING'}    = 'No mapping done for this type of xref';

  $summary_failed{'MASTER_FAILED'} = 'Master failed';
  $desc_failed{'MASTER_FAILED'}    = 'The dependent xref was not matched due to the master xref not being mapped';

  $summary_failed{'NO_MASTER'} = 'No Master';
  $desc_failed{'NO_MASTER'}    = 'The dependent xref was not matched due to there being no master xref';

  return {
    summary => \%summary_failed,
    desc    => \%desc_failed
  };
} ## end sub get_unmapped_reason


=head2 get_insert_direct_xref_low_priority
  Description: Iteratively retrieve all xrefs that are have info_type as DIRECT
               and ox_status as FAILED_PRIORITY and dumped as NULL
  Return type: Iterator - Hashref
  Example    : my $xrefs = $ba_handle->get_insert_direct_xref_low_priority();
               while( my $xref_ref = $xrefs->() ) {
                 my %xref = %{ $xref_ref };
                 print $xref{'acc'};
               }

=cut

sub get_insert_direct_xref_low_priority {
  my ( $self ) = @_;

  my $sql = (<<'SQL');
    SELECT x.xref_id, x.accession, x.version, x.label, x.description,
           x.info_type, x.info_text, s.name
    FROM source s,xref x
         LEFT JOIN  object_xref ox ON ox.xref_id = x.xref_id
    WHERE x.source_id = s.source_id AND
          x.dumped IS NULL AND
          ox.ox_status != 'FAILED_PRIORITY' AND
          x.info_type = 'DIRECT'
SQL

  my $sth = $self->dbi->prepare_cached( $sql );

  $sth->execute();

  my ( $xref_id, $acc, $version, $label, $desc, $type, $info, $dbname );
  $sth->bind_columns(
    \$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname);

  return sub {
    if ( $sth->fetch() ) {
      return {
        xref_id => $xref_id,
        acc     => $acc,
        label   => $label,
        version => $version,
        desc    => $desc,
        info    => $info,
        type    => $type,
        dbname  => $dbname
      }
    }
  }

} ## end sub get_insert_direct_xref_low_priority


=head2 get_insert_dependent_xref_low_priority
  Description: Iteratively retrieve all xrefs that are have info_type as DEPENDENT
               and ox_status as FAILED_PRIORITY and dumped is NULL
  Return type: Iterator - Hashref
  Example    : my $xrefs = $ba_handle->get_insert_dependent_xref_low_priority();
               while( my $xref_ref = $xrefs->() ) {
                 my %xref = %{ $xref_ref };
                 print $xref{'acc'};
               }

=cut

sub get_insert_dependent_xref_low_priority {
  my ( $self ) = @_;

  my $sql = (<<'SQL');
    SELECT DISTINCT
      x.xref_id, x.accession, x.version, x.label, x.description,
      x.info_type, x.info_text, s.name, mx.accession
    FROM xref mx, source s, xref x
         LEFT JOIN dependent_xref dx ON  dx.dependent_xref_id = x.xref_id
         LEFT JOIN object_xref ox ON ox.xref_id = x.xref_id
    WHERE x.source_id = s.source_id AND
          dx.master_xref_id = mx.xref_id AND
          x.dumped IS NULL AND
          ox.ox_status != 'FAILED_PRIORITY' AND
          x.info_type = 'DEPENDENT'
    ORDER BY s.name, x.accession
SQL

  my $sth = $self->dbi->prepare_cached( $sql );

  $sth->execute();

  my ( $xref_id, $acc, $version, $label, $desc, $type, $info, $dbname, $parent );
  $sth->bind_columns(
    \$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname, \$parent);

  return sub {
    if ( $sth->fetch() ) {
      return {
        xref_id => $xref_id,
        acc     => $acc,
        label   => $label,
        version => $version,
        desc    => $desc,
        info    => $info,
        type    => $type,
        dbname  => $dbname
      }
    }
  }

} ## end sub get_insert_dependent_xref_low_priority


=head2 get_insert_sequence_xref_remaining
  Description: Iteratively retrieve all xrefs that are have info_type as SEQUENCE_MATCH
               and dumped is NULL
  Return type: Iterator - Hashref
  Example    : my $xrefs = $ba_handle->get_insert_sequence_xref_remaining();
               while( my $xref_ref = $xrefs->() ) {
                 my %xref = %{ $xref_ref };
                 print $xref{'acc'};
               }

=cut

sub get_insert_sequence_xref_remaining {
  my ( $self ) = @_;

  my $sql = (<<'SQL');
    SELECT  x.xref_id, x.accession, x.version, x.label, x.description, x.info_type, x.info_text,
            s.name, px.sequence_type,
            ox.ensembl_object_type, ox.ensembl_id,
            ix.query_identity, ix.target_identity, ox.ox_status
    FROM source s, primary_xref px, xref x
         LEFT JOIN object_xref ox ON ox.xref_id = x.xref_id
         LEFT JOIN identity_xref ix ON ix.object_xref_id = ox.object_xref_id
    WHERE x.source_id = s.source_id AND
          px.xref_id = x.xref_id AND
          x.dumped IS NULL AND
          x.info_type = 'SEQUENCE_MATCH'
    ORDER BY x.xref_id
SQL

  my $sth = $self->dbi->prepare_cached( $sql );

  $sth->execute();

  my ( $xref_id, $acc, $version, $label, $desc, $type, $info, $dbname,
       $seq_type, $ensembl_object_type, $ensembl_id, $q_id, $t_id, $status );
  $sth->bind_columns(
    \$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname,
    \$seq_type, \$ensembl_object_type, \$ensembl_id, \$q_id, \$t_id, \$status);

  return sub {
    if ( $sth->fetch() ) {
      return {
        xref_id             => $xref_id,
        acc                 => $acc,
        label               => $label,
        version             => $version,
        desc                => $desc,
        info                => $info,
        type                => $type,
        dbname              => $dbname,
        seq_type            => $seq_type,
        ensembl_object_type => $ensembl_object_type,
        ensembl_id          => $ensembl_id,
        q_id                => $q_id,
        t_id                => $t_id,
        status              => $status
      }
    }
  }

} ## end sub get_insert_sequence_xref_remaining


=head2 get_insert_misc_xref
  Description: Iteratively retrieve all xrefs that are have info_type as MISC
               and dumped is NULL
  Return type: Iterator - Hashref
  Example    : my $xrefs = $ba_handle->get_insert_misc_xref();
               while( my $xref_ref = $xrefs->() ) {
                 my %xref = %{ $xref_ref };
                 print $xref{'acc'};
               }

=cut

sub get_insert_misc_xref {
  my ( $self ) = @_;

  my $sql =(<<'SQL');
  SELECT x.xref_id, x.accession, x.version,
         x.label, x.description, x.info_type,
         x.info_text, s.name
  FROM xref x, source s
  WHERE x.source_id = s.source_id AND
        x.dumped IS NULL AND
        x.info_type = 'MISC'
SQL

  my $sth = $self->dbi->prepare( $sql );
  $sth->execute();

  my ( $xref_id, $acc, $version, $label, $desc, $type, $info, $dbname );
  $sth->bind_columns(
    \$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname);

  return sub {
    if ( $sth->fetch() ) {
      return {
        xref_id => $xref_id,
        acc     => $acc,
        label   => $label,
        version => $version,
        desc    => $desc,
        info    => $info,
        type    => $type,
        dbname  => $dbname
      }
    }
  }
} ## end sub get_insert_misc_xref


=head2 get_insert_other_xref
  Description: WEL (Whatever is left)
               These are those defined as dependent but the master never existed
               and the xref and their descriptions etc are loaded first with the
               dependency's added later so did not know they had no masters at
               time of loading. (e.g. EntrezGene, WikiGene, MIN_GENE, MIM_MORBID)
  Return type: Iterator - HashRef
  Example    : my $xrefs = $ba_handle->get_insert_other_xref();
               while( my $xref_ref = $xrefs->() ) {
                 my %xref = %{ $xref_ref };
                 print $xref{'acc'};
               }
=cut

sub get_insert_other_xref {
  my ( $self ) = @_;

  my $sql =(<<'SQL');
    SELECT DISTINCT x.xref_id, x.accession, x.version,
                    x.label, x.description, x.info_type,
                    x.info_text, s.name
    FROM source s, xref x
    WHERE x.source_id = s.source_id AND
          x.dumped IS NULL AND
          x.info_type = 'DEPENDENT'
SQL

  my $sth = $self->dbi->prepare( $sql );
  $sth->execute();

  my ( $xref_id, $acc, $version, $label, $desc, $type, $info, $dbname );
  $sth->bind_columns(
    \$xref_id, \$acc, \$version, \$label, \$desc, \$type, \$info, \$dbname);

  return sub {
    if ( $sth->fetch() ) {
      return {
        xref_id => $xref_id,
        acc     => $acc,
        label   => $label,
        version => $version,
        desc    => $desc,
        info    => $info,
        type    => $type,
        dbname  => $dbname
      }
    }
  }
} ## end sub get_insert_other_xref

1;

