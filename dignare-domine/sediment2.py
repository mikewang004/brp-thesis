from efficiency9 import *

# rates_or_eff = 7

# hits_real_all = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = rates_or_eff, apply_shadow_mask= False, new_versions = 1)
# hits_sim_all = speedrun_heatmap(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number,int_rates_or_eff = rates_or_eff, apply_shadow_mask= False, new_versions = 1)
# hits_ratio_all = heatmap(hits_real_all.heatmap/hits_sim_all.heatmap)

# hits_real_all.plot_heatmap(indices, pmt_letters, "Location of the new PMTs", save_map = "plots/new-vs-old", include_mean = False)


# hits_ratio_dist = heatmap(speedrun_heatmap_31(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = rates_or_eff, apply_shadow_mask= False).heatmap/
#     speedrun_heatmap_31(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = rates_or_eff, apply_shadow_mask= False).heatmap)

# hits_ratio_upper = heatmap(np.nanmean(hits_ratio_dist.heatmap[18:, :, :], axis=0)); hits_ratio_lower = heatmap(np.nanmean(hits_ratio_dist.heatmap[:18, :, :], axis = 0));
# ratio_upper_lower = (hits_ratio_upper.heatmap/hits_ratio_lower.heatmap)


# print(hits_ratio_all)

#plot_heatmap_ultra_basic(ratio_upper_lower, "Ratio of upper half DOMs divided by lower half, test map", stringlist, floorlist, save_map = "plots/sediment", include_mean = True)

num_bins = 50
rates_or_eff = 5 #4 for hits 5 for efficiencies 8 for rates 
dist_hits_heatmap = speedrun_heatmap_31(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = rates_or_eff, new_versions= 0)
dist_hits_pmts = speedrun_heatmap_31(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = rates_or_eff, new_versions = 1)
dist_hits = dist_plots(dist_hits_heatmap.heatmap); dist_hits_new = dist_plots(dist_hits_pmts.heatmap)



#print(dist_hits_pmts.heatmap[:, : 1])
# test = (0, 12000); limits = (0.3, 1.3)
# for i in range(0, 21):
#     plt.figure()
#     if dist_hits_heatmap.heatmap[:, :, i].size != 0:
#         dist_hits.plot_dist_barebones(num_bins = num_bins, heatmap = dist_hits_heatmap.heatmap[:, :, i], label = "old PMTs",range = limits)
#     if dist_hits_pmts.heatmap[:, :, i].size != 0:
#         dist_hits_new.plot_dist_barebones(num_bins = num_bins, heatmap = dist_hits_pmts.heatmap[:, :, i], label = "new PMTs",range = limits)
#     plt.xlabel("efficiencies")
#     plt.ylabel("count")
#     plt.legend()
#     dist_hits.plot_dist_save("Distributions of efficiencies, old and new PMTs, string %i" %(stringlist[i]), "plots/new-vs-old/distributions")
#     plt.close()


limits = (0.3, 1.3)
plt.figure()
dist_hits.plot_dist_barebones(num_bins = 200, heatmap = dist_hits_heatmap.heatmap[:, :, :], label = "old PMTs", range = limits)
dist_hits_new.plot_dist_barebones(num_bins = 200, heatmap = dist_hits_pmts.heatmap[:, :, :], label = "new PMTs", range = limits)
plt.xlabel("efficiencies")
plt.ylabel("count")
plt.legend()
dist_hits.plot_dist_save("Distributions of efficiencies, old and new PMTs", "plots/new-vs-old/")
plt.close()

print(stats.ks_2samp(dist_hits.counts, dist_hits_new.counts))