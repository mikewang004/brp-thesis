from efficiency9 import *

#start_index = 6; stop_index = 17
start_index = 7; stop_index = 16 
apply_shadow_mask = False

data_real_old = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 0)
data_real_new = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 1)

eff_all = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 5).heatmap

hits_real_all = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, apply_shadow_mask= apply_shadow_mask)
hits_sim_all = speedrun_heatmap(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number, apply_shadow_mask= apply_shadow_mask)
hits_ratio_all = heatmap(hits_real_all.heatmap/hits_sim_all.heatmap)

eff_all[:, 6:16, 3] = np.nan
eff_all[:, :, 4] = np.nan

print(eff_all.shape)
# hits_real_all.plot_heatmap(indices, pmt_letters, "Mean hits per DOM", save_map = "plots/hits/rings")
#hits_ratio_all.plot_heatmap(indices, pmt_letters, "Ratio of data vs MC hits per DOM", save_map = "plots/data-mc/rings", 
#    include_mean=True)
# heatmap = hits_ratio_all.plot_heatmap_summarised_ring(indices, pmt_letters, "Ratio of data vs MC hits per string seperated by PMT ring", 
#            save_map = "plots/data-mc", string=True, include_mean = True, include_mean_of_mean = True)
# print(np.mean(heatmap[-1, :-1]))
# print(np.std(heatmap[-1, :-1]))
# print((heatmap[-1, 2] - np.mean(heatmap[-1, :-1]))/np.std(heatmap[-1, :-1]))


eff_lower = heatmap(eff_all[:2, :, :])
eff_mid = heatmap(eff_all[2:4, :, :])
eff_upper = heatmap(eff_all[4:, :, :])

zmin, zmax = get_zminmax(eff_all)

# eff_lower.plot_heatmap_summarised_ring(indices, pmt_letters[:4], "Mean of PMT efficiencies, ring A-D only", save_map = "plots/efficiencies", 
#     string = True, include_mean = True, include_mean_of_mean = True, zmin = zmin, zmax = zmax)
# eff_upper.plot_heatmap_summarised_ring(indices, pmt_letters[4:], "Mean of PMT efficiencies, ring E-F only", save_map = "plots/efficiencies", 
#     string = True, include_mean = True, include_mean_of_mean = True, zmin = zmin, zmax = zmax)


#Plot also these distributions
dist_plots_lower = dist_plots(eff_lower.heatmap)
dist_plots_mid = dist_plots(eff_mid.heatmap)
dist_plots_upper = dist_plots(eff_upper.heatmap)
#dist_plots_all_seperate = dist_plots(eff_all)


dist_plots_lower.plot_dist_barebones(label = "rings A-B")
dist_plots_mid.plot_dist_barebones(label = "rings C-D")
dist_plots_upper.plot_dist_barebones(label = "rings E-F")
#dist_plots_all_seperate.plot_dist_forloop(pmt_letters)
plt.xlabel("efficiency")
plt.ylabel("count")
plt.legend()
dist_plots_lower.plot_dist_save("Distribution of efficiencies, all strings except string 9 and 12", "plots/efficiencies")