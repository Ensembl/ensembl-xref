requires 'Config::IniFiles';
requires 'DBD::mysql';
requires 'DBI';
requires 'DBIx::Class';
requires 'HTTP::Tiny';
requires 'Moose';
requires 'SQL::Translator', '>= 0.11018';
requires 'Text::CSV';
requires 'Text::Glob';
requires 'URI::Escape';
requires 'XML::LibXML';
requires 'XML::Simple';

test_requires 'Config::General';
test_requires 'Devel::Cycle';
test_requires 'Devel::Peek';
test_requires 'Error';
test_requires 'IO::String';
test_requires 'PadWalker';
test_requires 'Perl::Critic::Moose';
test_requires 'Perl::Critic::Utils';
test_requires 'Test::Differences';
test_requires 'Test::Exception';
test_requires 'Test::Perl::Critic';
test_requires 'Test::Pod::Coverage';
test_requires 'Test::Warnings';
