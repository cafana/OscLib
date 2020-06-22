#!/bin/bash

TAG=v00.01
rm -rf $TAG # clean out any previous build
mkdir $TAG
mkdir $TAG/include
mkdir $TAG/src

ls

for olddir in 'OS='*
do
    newdir=$TAG/$olddir
    newdir=${newdir/OS=SLF6/slf6.x86_64}
    newdir=${newdir/OS=SLF7/slf7.x86_64}
    newdir=${newdir/QUALIFIER=/}
    newdir=${newdir/STAN=/}
    newdir=${newdir//,/.}
    newdir=${newdir//:/.}
    echo mv $olddir $newdir
    mv $olddir $newdir
    mv $newdir/OscLib/bin $newdir/$bin
    mv $newdir/OscLib/lib $newdir/$lib
    # will overwrite each other but should all be identical
    for k in `find $TAG/$newdir/OscLib -name '*.h'`
    do
        fname=${k/$newdir/}
        mkdir -p $TAG/include/`dirname $fname`
        cp $k $TAG/include/$fname
    done
    for k in `find $TAG/$newdir/OscLib`
    do
        fname=${k/$newdir/}
        mkdir -p $TAG/src/`dirname $fname`
        mv $k $TAG/src/$fname
    done
done

mkdir $TAG/ups
touch $TAG/ups/osclib.table
