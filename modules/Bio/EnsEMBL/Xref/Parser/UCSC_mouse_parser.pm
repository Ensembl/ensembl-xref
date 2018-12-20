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

Bio::EnsEMBL::Xref::Parser::UCSC_mouse_parser

=head1 DESCRIPTION

A parser class to parse UCSC data for mouse.
This module replicates the generic UCSC parser for mouse specific data
This prevents cross-mapping between species by treating each species as a separate source

-data_uri = ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/knownGene.txt.gz
-file_format = TSV
-columns = [
    ensembl_id
    chromosome
    strand
    txStart
    txEnd
    cdsStart
    cdsEnd
    nb_exons
    startExons
    endExons
    uniprot_accession
    ucsc_accession
  ]

Only columns listed in @required_columns are mandatory.

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::UCSC_mouse_parser->new(
    source_id  => 1,
    species_id => 10090,
    files      => ['UCSC_mouse/knownGene.txt.gz'],
    xref_dba   => $xref_dba
  );

  $parser->run();
=cut

package Bio::EnsEMBL::Xref::Parser::UCSC_mouse_parser;

use strict;
use warnings;

use parent qw( Bio::EnsEMBL::Xref::Parser::UCSCParser );

1;
