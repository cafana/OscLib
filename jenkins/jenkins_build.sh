#!/bin/bash

set +ex
env

source /cvmfs/nova.opensciencegrid.org/externals/setup || exit 1

if [[ $QUALIFIER == *e19* ]]
then
    # DUNE lblpwgtools versions
    setup root v6_18_04d -q ${QUALIFIER}:py2 || exit 1
    setup boost v1_70_0 -q $QUALIFIER || exit 1
else
    # NOvA versions
    setup root v6_16_00 -q $QUALIFIER || exit 1
    setup boost v1_66_0a -q $QUALIFIER || exit 1
fi

setup eigen v3.3.5 || exit 1

if [ $STAN == stan ]
then
    if [[ $QUALIFIER == *e19* ]]
    then
        # This is silly...
        setup stan_math v2_18_0 -q $QUALIFIER || exit 1
    else
        setup stan_math v2.18.0 -q $QUALIFIER || exit 1
    fi
fi

make clean # don't trust my build system
time make -j || exit 2
