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

=head1 POD COVERAGE

Test to ensure at minimum all public functions have descriptions

=cut

use strict;
use warnings;

use Test::More;
use Test::Warnings;

my $package_available = eval { use Test::Pod::Coverage };
plan skip_all => 'Test::Pod::Coverage required' if $package_available;

foreach my $mod ( all_modules('modules') ) {
  pod_coverage_ok($mod);
}

done_testing();
