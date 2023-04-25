\declare -i startdata=14423
declare -i enddata=14457
for (( i = $startdata; i <= $enddata; i++ ))
do
    declare  filename="KM3NeT_00000133_000${i}.v8.0_PMTeff_new.ToT.QE.PMTeff.txt"
    cp /../../../../../../../sps/km3net/repo/data/calibration/KM3NeT_00000133/v8.0_PMTeff_new/calibration/$filename ./$filename
done

