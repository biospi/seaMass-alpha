#!/bin/bash                                                     
. /opt/intel/bin/iccvars.sh intel64
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                   
mkdir $DIR/build
pushd $DIR/build
cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icc ..
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/x86_64-linux-gnu/c++/4.8
make
popd
