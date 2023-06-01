import numpy as np 
import wrapper_muonsim

hit_array = np.load("muon_hit_data-debug.npy")
for i in range(0, hit_array.shape[0]):
    print(hit_array[i, 0, :])
