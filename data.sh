#!/bin/bash                               
set -e                   
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                  

mkdir $DIR/data/out || true

export MKL_NUM_THREADS=4
export OMP_NUM_THREADS=4

mkdir $DIR/data/out/$3 || true
pushd $DIR/data/out/$3


mkdir 1.seamass || true
pushd 1.seamass
../../../../build/$1/$2/commandline/mzmlb2smb ../../../$3.mzMLb
../../../../build/$1/$2/commandline/seamass $3.$4.smb -w$5 -d11
popd

mkdir 2.seamass-restore || true
pushd 2.seamass-restore
../../../../build/$1/$2/commandline/seamass-restore ../1.seamass/$3.$4.smv -d1
../../../../build/$1/$2/commandline/smb2mzmlb ../../../$3.mzMLb -i .
popd

mkdir 3.seamass-restore_--deconvolve || true
pushd 3.seamass-restore_--deconvolve
../../../../build/$1/$2/commandline/seamass-restore ../1.seamass/$3.$4.smv -v -d1
../../../../build/$1/$2/commandline/smb2mzmlb ../../../$3.mzMLb -i .
popd

mkdir 4.seamass-restore_--reconstruct || true
pushd 4.seamass-restore_--reconstruct
../../../../build/$1/$2/commandline/seamass-restore ../1.seamass/$3.$4.smv -r -d1
../../../../build/$1/$2/commandline/smb2mzmlb ../../../$3.mzMLb -i .
popd

mkdir 5.seamass-restore_--centroid || true
pushd 5.seamass-restore_--centroid
../../../../build/$1/$2/commandline/seamass-restore ../1.seamass/$3.$4.smv -c -d1
../../../../build/$1/$2/commandline/smb2mzmlb ../../../$3.mzMLb -i .
popd

mkdir 6.seamass-restore_--deconvolve_--centroid || true
pushd 6.seamass-restore_--deconvolve_--centroid
../../../../build/$1/$2/commandline/seamass-restore ../1.seamass/$3.$4.smv -v -c -d1
../../../../build/$1/$2/commandline/smb2mzmlb ../../../$3.mzMLb -i .
popd


popd
