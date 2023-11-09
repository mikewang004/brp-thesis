import numpy as np

hit_array = np.load("muon_hit_data.npy")


trial_map = np.unique(hit_array[:, 0, 0])

trial_map_2 = np.zeros([trial_map.shape[0], 2])
rand = np.random.randint(0, 31, size = trial_map.shape[0])
trial_map_2[:, 0] = trial_map
for i in range(0, trial_map.shape[0]):
    trial_map_2[:, 1] = rand

print(trial_map_2)

np.savetxt("randmap.txt", trial_map_2)