#!/bin/sh -e


module purge
module load gcc/7.2 cmake/3.9.6

mkdir -p bld
rm -rf bld/*

pushd bld
    cmake ../src/
    make -j6
popd

module purge
