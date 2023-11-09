from efficiency9 import *

#Collar shadows PMTs C2, C5, E2, E5. Furthermore more dark rates for PMTs A1, B5, B6 and equator tape shadowing for PMTs D, E.
#TODO create filter for said PMTs

eff_all = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 5, apply_shadow_mask=False)
eff_all_shadowed = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 5, apply_shadow_mask=True)

eff_all.plot_heatmap(indices, pmt_letters, "Efficiencies per string, all except shadowed PMTs", save_map="plots/shadowing/rings", include_mean=True, zmax_array = 1.2, zmin_array = 0.6)
eff_all_shadowed.plot_heatmap(indices, pmt_letters, "Efficiencies per string, shadowed PMTs", save_map="plots/shadowing/rings", include_mean=True, zmax_array = 1.2, zmin_array = 0.6)

eff_ratio = heatmap(eff_all_shadowed.heatmap/eff_all.heatmap)
eff_ratio.plot_heatmap(indices, pmt_letters, "Shadowed PMTs over non-shadowed PMTs", save_map="plots/shadowing/rings", include_mean=True)

muon_hit_data_ratio_else = (speedrun_heatmap_31(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 4, apply_shadow_mask=False).heatmap[:, :, :] / 
    speedrun_heatmap_31(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 4, apply_shadow_mask=False).heatmap[:, :, :])
muon_hit_data_ratio_shadow = (speedrun_heatmap_31(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 4, apply_shadow_mask=True).heatmap[:, :, :] / 
    speedrun_heatmap_31(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 4, apply_shadow_mask=True).heatmap[:, :, :])

#eff_all.plot_heatmap_summarised_ring(indices, pmt_letters, title="Efficiencies per string seperated by PMT ring, all except shadowed PMTs", 
#    save_map = "plots/shadowing", include_mean = True, include_mean_of_mean = True)


#Look now at data/mc ratios 
apply_shadow_mask = True
hits_real_shadowed = speedrun_heatmap_31(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, apply_shadow_mask= apply_shadow_mask)
hits_sim_shadowed = speedrun_heatmap_31(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number, apply_shadow_mask= apply_shadow_mask)
hits_ratio_shadowed = heatmap(hits_real_shadowed.heatmap/hits_sim_shadowed.heatmap)

hits_ratio_else = heatmap(speedrun_heatmap_31(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, apply_shadow_mask= False).heatmap/
    speedrun_heatmap_31(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number, apply_shadow_mask= False).heatmap)

# Compare differencence between shadowed and else heatmaps 

#hits_ratio_ratio = heatmap(((hits_ratio_else.heatmap - hits_ratio_shadowed.heatmap)/hits_ratio_shadowed.heatmap)*100)

#print(np.nanmean(hits_ratio_ratio.heatmap[4:, :12, 18:]))

#hits_ratio_ratio.plot_heatmap(indices, pmt_letters, "Differences between shadowed and non-shadowed PMTs", save_map="plots/shadowing/rings/data-mc", include_mean=True, zmax_array =40, zmin_array = -40)


#hits_ratio_shadowed.plot_heatmap(indices, pmt_letters, "Ratio of data over MC hits, shadowed PMTs only", save_map="plots/shadowing/rings/data-mc", include_mean=True, zmax_array =1.0, zmin_array = 0.65)
#hits_ratio_else.plot_heatmap(indices, pmt_letters, "Ratio of data over MC hits, non-shadowed PMTs only", save_map="plots/shadowing/rings/data-mc", include_mean=True, zmax_array =1.0, zmin_array = 0.65)

# #Create distributions 
#dist_all = dist_plots(hits_ratio_shadowed.heatmap)
#dist_shadowed = dist_plots(hits_ratio_.heatmap)
dist_non_shadowed = dist_plots(muon_hit_data_ratio_else)
dist_lower_shadowed = dist_plots(muon_hit_data_ratio_shadow[:18, :, :])
dist_upper_shadowed = dist_plots(muon_hit_data_ratio_shadow[18:, :, :])
dist_lower_clean = dist_plots(muon_hit_data_ratio_else[:18, :, :])
dist_upper_clean = dist_plots(muon_hit_data_ratio_else[18:, :, :])
bins = 150; range = (0.4, 1.3)
#dist_all.plot_dist_barebones(label = "all PMTs", num_bins = bins, range = range)
#dist_shadowed.plot_dist_barebones(label = "all shadowed PMTs", num_bins = bins, range = range)
#dist_non_shadowed.plot_dist_barebones(label = "non-shadowed PMTs", num_bins = bins, range = range)
dist_lower_clean.plot_dist_barebones(label = "rings A-D, non-shadowed", num_bins = bins, range = range)
dist_upper_clean.plot_dist_barebones(label = "rings E-F, non-shadowed", num_bins = bins, range = range)
dist_lower_shadowed.plot_dist_barebones(label = "C2, C5", num_bins = bins, range = range)
dist_upper_shadowed.plot_dist_barebones(label = "E2, E5", num_bins = bins, range = range)
plt.xlabel("data over MC ratio")
plt.ylabel("count")
plt.legend(loc="upper left", fontsize="x-small")
dist_non_shadowed.plot_dist_save("Distribution of data over MC ratio, shadowed and non-shadowed PMTs", "plots/shadowing")

# print(stats.ks_2samp(dist_lower_shadowed.counts, dist_lower_clean.counts))
# print(stats.ks_2samp(dist_upper_shadowed.counts, dist_upper_clean.counts))