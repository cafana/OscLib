#!/bin/bash

set +ex
env

if [[ $QUALIFIER != *:n313* && $QUALIFIER != *:n315* && $QUALIFIER != *:n316* ]]
then
    echo Unspecified nutools version in qualifier $QUALIFIER -- must be n313, n315, or n316
    exit 1
fi

if [[ $QUALIFIER != *e20* && $QUALIFIER != *e26* && $QUALIFIER != *c7* $QUALIFIER != *c14* ]]
then
    echo Unknown compiler in qualifier $QUALIFIER -- must be e20, e26, c7, or c14
    exit 1
fi

if [[ $QUALIFIER != *debug* && $QUALIFIER != *prof* ]]
then
    echo Unknown optimization level in qualifier $QUALIFIER -- must be debug or prof
    exit 1
fi

if [[ x$STAN != *stan* ]]
then
    echo Must specify stan or stanfree in STAN variable $STAN
    exit 1
fi

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh || exit 1


# Looping over lines is a total pain in bash. Easier to just send it to a file
TMPFILE=`mktemp`
# Expect to be run in the directory one above....
jenkins/dependencies.sh $QUALIFIER:$STAN | sed 's/^/setup /' > $TMPFILE
cat $TMPFILE
source $TMPFILE


make clean # don't trust my build system
time make -j || exit 2


mkdir -p OscLib/ups
jenkins/make_table.sh > OscLib/ups/osclib.table
