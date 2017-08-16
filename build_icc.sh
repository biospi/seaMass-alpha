#!/bin/bash                                                     
. /opt/intel/bin/iccvars.sh intel64
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                  
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/x86_64-linux-gnu/c++/4.*
 
mkdir $DIR/build
mkdir $DIR/build/icc

mkdir $DIR/build/icc/debug
pushd $DIR/build/icc/debug
cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icc -DCMAKE_CXX_FLAGS="-std=c++11" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=install $@ ../../..
popd

mkdir $DIR/build/icc/release
pushd $DIR/build/icc/release
cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icc -DCMAKE_CXX_FLAGS="-std=c++11" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=install $@ ../../..
make -j $1
make install
popd

source data.sh icc HYE124_TTOF6600_64var_lgillet_I150211_008__index_59994 p-55-227433333333 6
source data.sh icc P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500 p 6

julia runtests.jl
