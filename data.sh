#!/bin/bash                                                     
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                  

rm -rf $DIR/data/out
mkdir $DIR/data/out

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

mkdir $DIR/data/out/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500__index3163
pushd $DIR/data/out/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500__index3163


mkdir 1.mzmlb2smb
pushd 1.mzmlb2smb
../../../../build/$1/debug/commandline/mzmlb2smb ../../../P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500__index3163.mzMLb -d1
popd

mkdir 2.seamass
pushd 2.seamass
../../../../build/$1/debug/commandline/seamass ../1.mzmlb2smb/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500__index3163.p-411-25-40-27.smb -d1
popd

mkdir 3.seamass-restore
pushd 3.seamass-restore
../../../../build/$1/debug/commandline/seamass-restore ../2.seamass/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500__index3163.p-411-25-40-27.smv -d1
popd

mkdir 4.smb2mzmlb
pushd 4.smb2mzmlb
../../../../build/$1/debug/commandline/smb2mzmlb ../../../P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500__index3163.mzMLb -i ../3.seamass-restore -d1
popd

mkdir 5.seamass-restore_--centroid
pushd 5.seamass-restore_--centroid
../../../../build/$1/debug/commandline/seamass-peak ../2.seamass/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500__index3163.p-411-25-40-27.smv -d1
popd

mkdir 6.smb2mzmlb
pushd 6.smb2mzmlb
../../../../build/$1/debug/commandline/smb2mzmlb ../../../P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500__index3163.mzMLb -i ../5.seamass-restore_--centroid -d1
popd


mkdir $DIR/data/out/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500
pushd $DIR/data/out/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500

mkdir 1.mzmlb2smb
pushd 1.mzmlb2smb
../../../../build/$1/debug/commandline/mzmlb2smb ../../../P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500.mzMLb -d1
popd

mkdir 2.seamass
pushd 2.seamass
../../../../build/$1/debug/commandline/seamass ../1.mzmlb2smb/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500.p.smb -d1
popd

mkdir 3.seamass-restore
pushd 3.seamass-restore
../../../../build/$1/debug/commandline/seamass-restore ../2.seamass/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500.p.smv -d1
popd

mkdir 4.smb2mzmlb
pushd 4.smb2mzmlb
../../../../build/$1/debug/commandline/smb2mzmlb ../../../P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500.mzMLb -i ../3.seamass-restore -d1
popd

mkdir 5.seamass-restore_--centroid
pushd 5.seamass-restore_--centroid
../../../../build/$1/debug/commandline/seamass-peak ../2.seamass/P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500.p.smv -d1
popd

mkdir 6.smb2mzmlb
pushd 6.smb2mzmlb
../../../../build/$1/debug/commandline/smb2mzmlb ../../../P02U_Swath_1__mzWindow_602_605__scanTime_2300_3500.mzMLb -i ../5.seamass-restore_--centroid -d1
popd


popd


popd
