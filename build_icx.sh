#!/bin/bash                                                     
source /opt/intel/oneapi/setvars.sh
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                  
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/x86_64-linux-gnu/c++/4.*
 
mkdir $DIR/build
mkdir $DIR/build/icx

mkdir $DIR/build/icx/debug
pushd $DIR/build/icx/debug
cmake -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCMAKE_CXX_FLAGS="-std=c++11" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=install $@ ../../..
popd

mkdir $DIR/build/icx/release
pushd $DIR/build/icx/release
cmake -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCMAKE_CXX_FLAGS="-std=c++11" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=install $@ ../../..
make -j $1
make install
popd

### commented out as seamass broken atm (only asrl works)
#source data.sh icx HYE124_TTOF6600_64var_lgillet_I150211_008__index_59994 p-55-227433333333 6
#source data.sh icx P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500 p 6

#julia runtests.jl
