#!/bin/bash

set +ex
env

if [[ $QUALIFIER == *e20* || $QUALIFIER == *e19* || $QUALIFIER == *c7* ]]
then
    # DUNE lblpwgtools versions
    source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh || exit 1
    # Add SBN products path to get additional products
    export PRODUCTS=/cvmfs/sbn.opensciencegrid.org/products/sbn/:$PRODUCTS
    setup root v6_22_06a -q ${QUALIFIER}:p383b || exit 1
    setup boost v1_73_0 -q $QUALIFIER || exit 1
    setup eigen v3_3_9a || exit 1
    if [ $STAN == stan ]; then setup stan_math v2_18_0a -q $QUALIFIER || exit 1; fi
else
    # NOvA versions
    source /cvmfs/nova.opensciencegrid.org/externals/setup || exit 1
    setup root v6_16_00 -q $QUALIFIER || exit 1
    setup boost v1_66_0a -q $QUALIFIER || exit 1
    setup eigen v3.3.5 || exit 1
    if [ $STAN == stan ]; then setup stan_math v2.18.0 -q $QUALIFIER || exit 1; fi
fi

make clean # don't trust my build system
time make -j || exit 2
