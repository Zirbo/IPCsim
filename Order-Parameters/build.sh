#!/bin/bash
pushd src
    make
    make clean
    mv crystal_analysis ..
popd
