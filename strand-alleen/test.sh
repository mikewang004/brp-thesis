\declare -i startdata=3092
declare -i enddata=3103
for (( i = $startdata; i <= $enddata; i++ ))
do
    declare  filename="mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.${i}.root"
    cp /../../../../../../../sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/$filename ./showerdata/$filename
done

