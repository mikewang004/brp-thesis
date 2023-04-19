import numpy as np

non_valid_runs = np.loadtxt("non-valid-runs.txt", delimiter=",")

runs = np.arange(14399, 14422)

mask = np.in1d(runs, non_valid_runs)



if mask.any() == True:
    print(runs[mask])
    raise Exception("do not use this data-set")
else:
    with open("download_data.sh" , 'r+') as file:
        filedata = file.readlines()
        print(filedata[1])
        filedata[1] = "declare -i startdata=%i\n" %(runs[0]) 
        filedata[2] = "declare -i enddata=%i\n" %(runs[-1])
        file.seek(0)
        file.writelines(filedata)


