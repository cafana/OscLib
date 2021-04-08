#!/bin/bash

set +ex
env

if [[ $QUALIFIER == *e20* || $QUALIFIER == *c7* ]]
then
    # These are the versions for nutools v3_09_02 (latest larsoft)
    source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh || exit 1
    # Add SBN products path to get additional products
    export PRODUCTS=/cvmfs/sbn.opensciencegrid.org/products/sbn/:$PRODUCTS

    setup root v6_22_06a -q ${QUALIFIER}:p383b || exit 1
    setup boost v1_73_0 -q $QUALIFIER || exit 1
    setup eigen v3_3_9a || exit 1

    if [ $STAN == stan ]; then setup stan_math v4_0_1 -q $QUALIFIER || exit 1; fi
fi

if [[ $QUALIFIER == *e19* ]]
then
    # NOvA versions (nutools v3_08_00)
    source /cvmfs/nova.opensciencegrid.org/externals/setup || exit 1

    # These are for nutools v3_08_00 (current nova)
    setup root v6_18_04d -q $QUALIFIER || exit 1
    setup boost v1_70_0 -q $QUALIFIER || exit 1
    setup eigen v3_3_9a || exit 1

    if [ $STAN == stan ]; then setup stan_math v4_0_1 -q $QUALIFIER || exit 1; fi
fi

if [[ $QUALIFIER == *e17* ]]
then
    # old NOvA versions
    source /cvmfs/nova.opensciencegrid.org/externals/setup || exit 1
    setup root v6_16_00 -q $QUALIFIER || exit 1
    setup boost v1_66_0a -q $QUALIFIER || exit 1
    setup eigen v3.3.5 || exit 1
    if [ $STAN == stan ]; then
        echo e17 build no longer supports stan
        exit 1
    fi
fi

make clean # don't trust my build system
time make -j || exit 2
