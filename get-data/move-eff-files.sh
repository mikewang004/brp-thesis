#!/bin/sh

module load jpp/master
declare -i startdata=14399
declare -i enddata=14400
for (( i = $startdata; i <= $enddata; i++ ))
do
    echo "$i"
    cd /
    cd /sps/km3net/repo/data/calibration/KM3NeT_00000133/v8.0_PMTeff_new/calibration/
    cp KM3NeT_00000133_000%i.v8.0_PMTeff_new.ToT.QE.PMTeff.txt  ../../../../../../../pbs/home/m/mwang/zee-symfonie/get-data/KM3NeT_00000133_000%i.v8.0_PMTeff_new.ToT.QE.PMTeff.txt
done

