import numpy as np

non_valid_runs = np.loadtxt("non-valid-runs.txt", delimiter=",")

runs = np.array([13670, 14661, 14662])

mask = np.in1d(runs, non_valid_runs)

print(mask)

if mask.any() == True:
    print(runs[mask])
    raise Exception("do not use this data-set")
else:
    print("all ok!")