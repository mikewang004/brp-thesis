from efficiency9 import *




#data_real_old = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 0)
#data_real_new = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 1)

hits_real_old= speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 0)
hits_real_new = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 1)

hits_real_old.plot_heatmap_summarised_ring(indices, pmt_letters, "Hits per string seperated by PMT ring, old PMTs", 
           save_map = "plots/new-vs-old", string=True)

hits_real_new.plot_heatmap_summarised_ring(indices, pmt_letters, "Hits per string seperated by PMT ring, new PMTs", 
          save_map = "plots/new-vs-old", string=True)



summarised_heatmap_ratio(hits_real_new,  hits_real_old, "Ratio of old vs new PMTs hits per floor seperated by PMT ring", 
    indices, pmt_letters, string=False, save_map = "plots/new-vs-old", start_index = 6, stop_index = 17)

