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

Bio::EnsEMBL::Xref::Parser::ArrayExpressParser

=head1 DESCRIPTION

A parser class to parse the ArrayExpress source. Fetches the stable ids from core database
and populates the xref database's accession and label columns as DIRECT xref,
provided the species is active and declared in ArrayExpress.

-species = MULTI (Used by all ensembl species)

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::ArrayExpressParser->new(
    source_id  => 1,
    species_id => 9606,
    files      => ["project=>ensembl"],
    xref_dba   => $xref_dba,  # xref db adaptor
    dba        => $dba   # core db adaptor
  );

  $parser->run();
=cut


package Bio::EnsEMBL::Xref::Parser::ArrayExpressParser;

use strict;
use warnings;
use Carp;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Xref::FetchFiles;
use URI::ftp;

use parent qw( Bio::EnsEMBL::Xref::Parser );

my $default_ftp_server = 'ftp.ebi.ac.uk';
my $default_ftp_dir =
  'pub/databases/microarray/data/atlas/bioentity_properties/ensembl';

sub run {

  my ( $self, $ref_arg ) = @_;
  my $source_id    = $self->{source_id};
  my $species_id   = $self->{species_id};
  my $species_name = $self->{species};
  my $files        = $self->{files};
  my $verbose      = $self->{verbose} // 0;
  my $xref_dba     = $self->{xref_dba};
  my $dba          = $self->{dba};

  my $file = shift @{$files};

  # project could be ensembl or ensemblgenomes
  my $project;
  if ( $file =~ /project[=][>](\S+?)[,]/ ) {
    $project = $1;
  }

  my %species_id_to_names = $xref_dba->species_id2name();

  # the species_name passed could be an alias that is different from what is stored in species table,
  # so it is better to store it as well in the hash. We later check if the species is active in ArrayExpress
  # using all these speices names

  if ( defined $species_name ) {
    push @{ $species_id_to_names{$species_id} }, $species_name;
  }

  if ( !exists $species_id_to_names{$species_id} ) { return 0; }

  my $names               = $species_id_to_names{$species_id};
  my $species_lookup      = $self->_get_species($verbose);
  my $active = $self->_is_active_species( $species_lookup, $names, $verbose );

  if ( !$active ) {
    return 0;
  }

  $species_name = $species_id_to_names{$species_id}[0];

  #get stable_ids from core and create xrefs
  my $gene_adaptor = $self->_get_gene_adaptor($project, $species_name, $dba);

  print "Finished loading the gene_adaptor\n" if $verbose;

  my @stable_ids = map { $_->stable_id } @{ $gene_adaptor->fetch_all() };

  my $xref_count = 0;
  foreach my $gene_stable_id (@stable_ids) {

    my $xref_id = $xref_dba->add_xref(
      {
        acc        => $gene_stable_id,
        label      => $gene_stable_id,
        source_id  => $source_id,
        species_id => $species_id,
        info_type  => "DIRECT"
      }
    );

    $xref_dba->add_direct_xref( $xref_id, $gene_stable_id, 'gene', '' );
    if ($xref_id) {
      $xref_count++;
    }
  }

  print "Added $xref_count DIRECT xrefs\n" if ($verbose);
  if ( !$xref_count ) {
   croak "No arrayexpress xref added\n";
  }

  return 0;      # successful

}

sub _get_gene_adaptor {
  my ( $self, $project, $species_name, $dba ) = @_;

  my $registry = 'Bio::EnsEMBL::Registry';

  my ($gene_adaptor);

  if ( defined $project && $project eq 'ensembl' ) {
    $registry->load_registry_from_multiple_dbs(
      {
        '-host' => 'mysql-ens-sta-1',
        '-port' => 4519,
        '-user' => 'ensro',
      },
    );
    $gene_adaptor = $registry->get_adaptor( $species_name, 'core', 'Gene' );
  }
  elsif ( defined $project && $project eq 'ensemblgenomes' ) {
    $registry->load_registry_from_multiple_dbs(
      {
        '-host' => 'mysql-eg-staging-1.ebi.ac.uk',
        '-port' => 4160,
        '-user' => 'ensro',
      },
      {
        '-host' => 'mysql-eg-staging-2.ebi.ac.uk',
        '-port' => 4275,
        '-user' => 'ensro',
      },
    );
    $gene_adaptor = $registry->get_adaptor( $species_name, 'core', 'Gene' );
  }
  elsif ( defined $dba ) {
    $gene_adaptor = $dba->get_GeneAdaptor();
  }
  else {
    die( "Missing or unsupported project value. Supported values: ensembl, ensemblgenomes" );
  }

  return $gene_adaptor;
}

sub _get_species {
  my ( $self, $verbose ) = @_;
  $verbose = ( defined $verbose ) ? $verbose : 0;

  my $ff = Bio::EnsEMBL::Xref::FetchFiles->new();
  my $ftp = $ff->get_ftp(URI::ftp->new("ftp://".$default_ftp_server));

  if ( !defined $ftp) {
    croak "Failed to get FTP connection to $default_ftp_server\n";
  }

  $ftp->cwd($default_ftp_dir);
  my @files = $ftp->ls() or confess "Cannot change to $default_ftp_dir: $@";
  $ftp->quit;

  my %species_lookup;
  foreach my $file (@files) {
    my ($species) = split( /\./, $file );
    $species_lookup{$species} = 1;
  }
  return \%species_lookup;
}

# checks if the species is still active in ArrayExress
sub _is_active_species {
  my ( $self, $species_lookup, $names, $verbose ) = @_;

  #Loop through the names and aliases first. If we get a hit then great
  my $active = 0;
  foreach my $name ( @{$names} ) {
    if ( $species_lookup->{$name} ) {
      printf( 'Found ArrayExpress has declared the name "%s". This was an alias' . "\n", $name ) if $verbose;
      $active = 1;
      last;
    }
  }
  return $active;
}

1;
