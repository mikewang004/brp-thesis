from efficiency9 import *

rates_or_eff = 4

hits_real_all = speedrun_heatmap_31(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = rates_or_eff, apply_shadow_mask= False, new_versions = None)

upper_lower_ratio = (np.nanmean(11/30*hits_real_all.heatmap[18:, :, :], axis = 0)/np.nanmean(18/30*hits_real_all.heatmap[:18, :, :], axis = 0))

plot_heatmap_ultra_basic(upper_lower_ratio, "Ratio of upper PMT hits vs lower PMT hits", stringlist, floorlist, save_map = "plots/gel-issue", include_mean = True)