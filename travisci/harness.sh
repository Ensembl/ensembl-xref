#!/bin/bash

ENSDIR="${ENSDIR:-$PWD}"

export PERL5LIB=$ENSDIR/bioperl-live:$ENSDIR/ensembl-test/modules:$PWD/modules:$ENSDIR/ensembl/modules
export TEST_AUTHOR=$USER

echo "Running test suite"
if [ "$COVERALLS" = 'true' ]; then
  PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl-test,+ignore,ensembl' perl $ENSDIR/ensembl-test/scripts/runtests.pl -verbose modules/t
else
  perl $ENSDIR/ensembl-test/scripts/runtests.pl modules/t
fi

rt=$?
if [ $rt -eq 0 ]; then
  if [ "$COVERALLS" = 'true' ]; then
    echo "Running Devel::Cover coveralls report"
    cover --nosummary -report coveralls
  fi
  exit $?
else
  exit $rt
fi
