module load jpp/master
declare -i startdata=14440
declare -i enddata=14442
z=$(( enddata - startdata ))
for (( i = 0; i <= $z; i++ ))
do
    echo $((i + $startdata))
done
