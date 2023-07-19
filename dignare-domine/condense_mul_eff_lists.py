import numpy as np

path = "../get-data/"
run_numbers = np.arange(14413,14440, 1)

test = np.loadtxt(path + "KM3NeT_00000133_00014413.v8.0_PMTeff_new.ToT.QE.PMTeff.txt", skiprows = 148, usecols=[1,2,3])
length = test.shape[0]
del test
eff_array = np.zeros([len(run_numbers), int(length), 3])
for i in range(0, len(run_numbers)):
    str_effs = path + "KM3NeT_00000133_000%i.v8.0_PMTeff_new.ToT.QE.PMTeff.txt" %(run_numbers[i])
    eff_array[i, :, :] = np.loadtxt(str_effs, skiprows = 148, usecols=[1,2,3])
    
eff_array_mean = np.mean(eff_array, axis = 0)
np.savetxt("runs14413-14440.txt", eff_array_mean)