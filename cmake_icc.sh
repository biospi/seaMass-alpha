#!/bin/bash                                                     
. /opt/intel/bin/iccvars.sh intel64
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                  
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/x86_64-linux-gnu/c++/4.*
 
mkdir $DIR/build
mkdir $DIR/build/icc

mkdir $DIR/build/icc/debug
pushd $DIR/build/icc/debug
cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icc -DCMAKE_BUILD_TYPE=Debug $@ ../../..
popd

mkdir $DIR/build/icc/release
pushd $DIR/build/icc/release
cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icc -DCMAKE_BUILD_TYPE=Release $@ ../../..
popd
