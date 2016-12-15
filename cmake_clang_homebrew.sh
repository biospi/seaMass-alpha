#!/bin/bash                                                     
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                  
export LIBRARY_PATH=/usr/local/opt/llvm/lib
 
mkdir $DIR/build
mkdir $DIR/build/clang

mkdir $DIR/build/clang/debug
pushd $DIR/build/clang/debug
cmake -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++ -DCMAKE_BUILD_TYPE=Debug $@ ../../..
popd

mkdir $DIR/build/clang/release
pushd $DIR/build/clang/release
cmake -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++ -DCMAKE_BUILD_TYPE=Release $@ ../../..
popd
