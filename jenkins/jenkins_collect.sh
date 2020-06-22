#!/bin/bash

TAG=v00.01
rm -rf $TAG # clean out any previous build
mkdir $TAG

ls

for oldname in 'OS='*
do
    newname=$oldname
    newname=${newname/OS=SLF6/slf6.x86_64}
    newname=${newname/OS=SLF7/slf7.x86_64}
    newname=${newname/QUALIFIER=/}
    newname=${newname/STAN=/}
    newname=${newname/,/.}
    newname=${newname/;/.}
    echo mv $oldname $TAG/$newname
    mv $oldname $TAG/$newname
    mkdir $TAG/$newname/bin
    mkdir $TAG/$newname/lib
done

cd $TAG
mkdir src
mkdir include
mkdir ups
touch ups/osclib.table
