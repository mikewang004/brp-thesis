from efficiency9 import *
hits_real_all = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, apply_shadow_mask= False)
hits_sim_all = speedrun_heatmap(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number, apply_shadow_mask= False)
hits_ratio_all = heatmap(hits_real_all.heatmap/hits_sim_all.heatmap)

hits_ratio_dist = heatmap(speedrun_heatmap_31(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, apply_shadow_mask= False).heatmap/
    speedrun_heatmap_31(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number, apply_shadow_mask= False).heatmap)

x