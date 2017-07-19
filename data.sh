#!/bin/bash                                                     
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                  

mkdir $DIR/data/out

export MKL_NUM_THREADS=4
export OMP_NUM_THREADS=4

mkdir $DIR/data/out/$2
pushd $DIR/data/out/$2


mkdir 1.seamass
pushd 1.seamass
../../../../build/$1/release/commandline/mzmlb2smb ../../../$2.mzMLb
../../../../build/$1/release/commandline/seamass $2.$3.smb -w$4 -d1
popd

mkdir 2.seamass-restore
pushd 2.seamass-restore
../../../../build/$1/release/commandline/seamass-restore ../1.seamass/$2.$3.smv -d1
../../../../build/$1/release/commandline/smb2mzmlb ../../../$2.mzMLb -i .
popd

mkdir 3.seamass-restore_--deconvolve
pushd 3.seamass-restore_--deconvolve
../../../../build/$1/release/commandline/seamass-restore ../1.seamass/$2.$3.smv -v -d1
../../../../build/$1/release/commandline/smb2mzmlb ../../../$2.mzMLb -i .
popd

mkdir 4.seamass-restore_--reconstruct
pushd 4.seamass-restore_--reconstruct
../../../../build/$1/release/commandline/seamass-restore ../1.seamass/$2.$3.smv -r -d1
../../../../build/$1/release/commandline/smb2mzmlb ../../../$2.mzMLb -i .
popd

mkdir 5.seamass-restore_--centroid
pushd 5.seamass-restore_--centroid
../../../../build/$1/release/commandline/seamass-restore ../1.seamass/$2.$3.smv -c -d1
../../../../build/$1/release/commandline/smb2mzmlb ../../../$2.mzMLb -i .
popd


popd
