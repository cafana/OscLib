#!/bin/bash

# Centralize all the logic about what versions to depend on in this one file

if [ $# != 1 ]
then
    echo Usage: dependencies QUALIFIER >&2
    exit 1
fi

QUAL=$1

if [[ $QUAL == *:n313* ]]; then NQUAL=n313; QUAL=${QUAL/:n313/}; fi
if [[ $QUAL == *:n316* ]]; then NQUAL=n316; QUAL=${QUAL/:n316/}; fi
if [[ $QUAL == *:n319* ]]; then NQUAL=n319; QUAL=${QUAL/:n319/}; fi


WANTSTAN=yes
if [[ $QUAL == *:stanfree ]]
then
  WANTSTAN=no
  QUAL=${QUAL/:stanfree/}
else
  if [[ $QUAL == *:stanthread ]]
  then
    QUAL=${QUAL/:stanthread/}
  else
    QUAL=${QUAL/:stan/}
  fi
fi

if [[ $NQUAL == n313 ]]
then
    # These are the older e20-era versions (nutools v3_13_04b)
    echo root v6_22_08d -q${QUAL}:p392
    echo boost v1_75_0 -q$QUAL
    echo eigen v3_4_0
elif [[ $NQUAL == n316 ]]
then
    # nulite v3_16_6
    echo root v6_28_12 -q${QUAL}:p3915
    echo boost v1_82_0 -q$QUAL
    echo eigen v23_08_01_66e8f
else
    # nulite v3_19_01
    echo root v6_28_12 -q${QUAL}:p3915
    echo boost v1_82_0 -q$QUAL
    echo eigen v23_08_01_66e8f
fi

if [ $WANTSTAN == yes ]
then
    echo stan_math v4_9_0a
    echo sundials v6_1_1
fi
