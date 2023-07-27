from efficiency9 import *

#start_index = 6; stop_index = 17
start_index = 7; stop_index = 16 
apply_shadow_mask = False

data_real_old = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 0)
data_real_new = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 1)

eff_all = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 5)

hits_real_all = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, apply_shadow_mask= apply_shadow_mask)
hits_sim_all = speedrun_heatmap(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number, apply_shadow_mask= apply_shadow_mask)
hits_ratio_all = heatmap(hits_real_all.heatmap/hits_sim_all.heatmap)

# hits_real_all.plot_heatmap(indices, pmt_letters, "Mean hits per DOM", save_map = "plots/hits/rings")
#hits_ratio_all.plot_heatmap(indices, pmt_letters, "Ratio of data vs MC hits per DOM", save_map = "plots/data-mc/rings", 
#    include_mean=True)
heatmap = hits_ratio_all.plot_heatmap_summarised_ring(indices, pmt_letters, "Ratio of data vs MC hits per string seperated by PMT ring", 
           save_map = "plots/data-mc", string=True, include_mean = True, include_mean_of_mean = True)

print(np.mean(heatmap[-1, :-1]))
print(np.std(heatmap[-1, :-1]))
print((heatmap[-1, 2] - np.mean(heatmap[-1, :-1]))/np.std(heatmap[-1, :-1]))