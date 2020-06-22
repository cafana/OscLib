#!/bin/bash

set +ex
env

source /cvmfs/nova.opensciencegrid.org/externals/setup || exit 1
setup root v6_16_00 -q prof:e17 || exit 1
setup eigen v3.3.5 || exit 1
setup boost v1_66_0a -q prof:e17 || exit 1

make clean # don't trust my build system
time make -j || exit 2
