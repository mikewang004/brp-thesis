module load jpp/master
declare -i startdata=14440
declare -i enddata=14442
z=$(( enddata - startdata ))
for (( i = 0; i <= $z; i++ ))
do
    JRunAnalyzer -f root://ccxroot:1999//hpss/in2p3.fr/group/km3net/data/raw/sea/KM3NeT_00000133/14/KM3NeT_00000133_000$z.root -a /pbs/throng/km3net/detectors/KM3NeT_00000133_20221025.detx -o jra_133_$z.root -@ "p0=1"
done
