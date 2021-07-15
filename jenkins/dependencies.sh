#!/bin/bash

# Centralize all the logic about what versions to depend on in this one file

if [ $# != 1 ]
then
    echo Usage: dependencies QUALIFIER >&2
    exit 1
fi

QUAL=$1

if [[ $QUAL == *:n308* ]]; then NQUAL=n308; QUAL=${QUAL/:n308/}; fi
if [[ $QUAL == *:n309* ]]; then NQUAL=n309; QUAL=${QUAL/:n309/}; fi

WANTSTAN=yes
if [[ $QUAL == *:stanfree ]]; then WANTSTAN=no; QUAL=${QUAL/:stanfree/}; else QUAL=${QUAL/:stan/}; fi

if [[ $NQUAL == n308 ]]
then
    # These are the current (Apr 2021) nova versions (nutools v3_08_00)
    echo root v6_18_04d -q$QUAL
    echo boost v1_70_0 -q$QUAL
else
    # These were the sbn versions in Apr 2021 (nutools v3_09_02)
#    echo root v6_22_06a -q${QUAL}:p383b
#    echo boost v1_73_0 -q$QUAL

    # These are the current (July 2021) sbn versions (nutools v3_09_04)
    echo root v6_22_08b -q${QUAL}:p383b
    echo boost v1_73_0 -q$QUAL
fi

echo eigen v3_3_9a

if [ $WANTSTAN == yes ]; then echo stan_math v4_0_1 -q$QUAL; fi
