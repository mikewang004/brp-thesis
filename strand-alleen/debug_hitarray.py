import numpy as np 
import matplotlib.pyplot as plt 

hit_array = np.load("muon_hit_data_real-001375x.npy")


hit_array_avg = np.mean(hit_array, axis =1 )

print(hit_array_avg)
print(hit_array_avg.shape)

np.savetxt("debug-pmt-id-map.txt", hit_array_avg, fmt="%0.10f")