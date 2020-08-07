#!/bin/bash

set +ex
env

source /cvmfs/nova.opensciencegrid.org/externals/setup || exit 1
setup root v6_16_00 -q $QUALIFIER || exit 1
setup eigen v3.3.5 || exit 1
setup boost v1_66_0a -q $QUALIFIER || exit 1

if [ $STAN == stan ]
then
    setup stan_math v2.18.0 -q $QUALIFIER || exit 1
fi

make clean # don't trust my build system
time make -j || exit 2
