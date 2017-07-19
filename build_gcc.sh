#!/bin/bash                                                     
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                  
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/x86_64-linux-gnu/c++/4.*
 
mkdir $DIR/build
mkdir $DIR/build/gcc

mkdir $DIR/build/gcc/debug
pushd $DIR/build/gcc/debug
cmake -DCMAKE_BUILD_TYPE=Debug $@ -DCMAKE_INSTALL_PREFIX=install ../../..
popd

mkdir $DIR/build/gcc/release
pushd $DIR/build/gcc/release
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=install $@ ../../..
make -j $1
make install
popd

source data.sh gcc HYE124_TTOF6600_64var_lgillet_I150211_008__index_59994 p-55-227433333333 6
source data.sh gcc P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500 p 6

julia runtests.jl
