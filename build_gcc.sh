#!/bin/bash                                                     
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                  
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/x86_64-linux-gnu/c++/4.*
 
mkdir $DIR/build
mkdir $DIR/build/gcc

mkdir $DIR/build/gcc/debug
pushd $DIR/build/gcc/debug
cmake -DCMAKE_BUILD_TYPE=Debug $@ -DCMAKE_INSTALL_PREFIX=install ../../..
make
make install
popd

mkdir $DIR/build/gcc/release
pushd $DIR/build/gcc/release
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=install $@ ../../..
make
make install
popd

source data.sh gcc
julia runtests.jl
