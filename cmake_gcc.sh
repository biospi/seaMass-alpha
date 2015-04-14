#!/bin/bash                                                     
. /opt/intel/bin/iccvars.sh intel64
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                  
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/x86_64-linux-gnu/c++/4.*
 
mkdir $DIR/build
mkdir $DIR/build/gcc

mkdir $DIR/build/gcc/debug
pushd $DIR/build/gcc/debug
cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icc -DCMAKE_BUILD_TYPE=Debug $@ ../../..
popd

mkdir $DIR/build/gcc/release
pushd $DIR/build/gcc/release
cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icc -DCMAKE_BUILD_TYPE=Release $@ ../../..
popd
