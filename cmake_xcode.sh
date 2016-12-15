#!/bin/bash                                                     
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                  
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/x86_64-linux-gnu/c++/4.*
 
mkdir $DIR/build
mkdir $DIR/build/xcode

mkdir $DIR/build/xcode/debug
pushd $DIR/build/xcode/debug
cmake -G Xcode -DCMAKE_BUILD_TYPE=Debug $@ ../../..
popd

mkdir $DIR/build/xcode/release
pushd $DIR/build/xcode/release
cmake -G Xcode -DCMAKE_BUILD_TYPE=Release $@ ../../..
popd
