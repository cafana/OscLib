#!/bin/bash

TAG=`git describe --tags`

echo FILE=TABLE
echo PRODUCT=osclib
echo VERSION=$TAG
echo
echo

for EXPT in n313 n315 n316 n319
do
    for OPT in debug prof
    do
        for COMPILER in e20 e26 c7 c14
        do
            for STAN in stan stanfree
            do

                echo FLAVOR=ANY
                echo QUALIFIERS=\"${OPT}:${COMPILER}:${EXPT}:${STAN}\"
                echo
                echo 'ACTION=SETUP'
                echo '  setupEnv()'
                echo '  proddir()'
                echo
                echo '  # for get-directory-name'
                echo '  setupRequired(cetpkgsupport)'
                # Expect to be run in the directory one above....
                jenkins/dependencies.sh ${OPT}:${COMPILER}:${EXPT}:${STAN} | while read line
                do
                    echo '  setupRequired('$line')'
                done
                echo
                echo '  EnvSet(OSCLIB_VERSION, ${UPS_PROD_VERSION} )'
                echo '  EnvSet(OSCLIB_INC, ${UPS_PROD_DIR}/include )'
                echo '  EnvSet(OSCLIB_FQ_DIR, ${OSCLIB_DIR}`get-directory-name subdir`.`echo ${UPS_PROD_QUALIFIERS} | tr ":" "."` )'
                echo '  EnvSet(OSCLIB_LIB, ${OSCLIB_FQ_DIR}/lib )'
                echo '  EnvSet(OSCLIB_BIN, ${OSCLIB_FQ_DIR}/bin )'
                echo '  pathPrepend(LD_LIBRARY_PATH, ${OSCLIB_LIB})'
                echo
                echo
            done
        done
    done
done
