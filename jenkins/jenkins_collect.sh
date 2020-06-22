#!/bin/bash

TAG=v00.01
mkdir $TAG

ls

for oldname in 'os='*
do
    newname=$oldname
    newname=${newname/os=SLF6/slf6.x86_64}
    newname=${newname/qualifier=/}
    newname=${newname/stan=/}
    newname=${newname/,/.}
    newname=${newname/;/.}
    mv $oldname $TAG/$newname
    mkdir $TAG/$newname/bin
    mkdir $TAG/$newname/lib
done

cd $TAG
mkdir src
mkdir include
mkdir ups
touch ups/osclib.table
