from efficiency9 import *




#data_real_old = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 0)
#data_real_new = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 1)

hits_real_old= speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 0)
hits_real_new = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 1)
eff_old = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 0, int_rates_or_eff=5)
eff_new = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 1, int_rates_or_eff=5)

hits_real_old.plot_heatmap_summarised_ring(indices, pmt_letters, "Hits per string seperated by PMT ring, old PMTs", 
           save_map = "plots/new-vs-old", string=False)

hits_real_new.plot_heatmap_summarised_ring(indices, pmt_letters, "Hits per string seperated by PMT ring, new PMTs", 
          save_map = "plots/new-vs-old", string=False)



summarised_heatmap_ratio(hits_real_new,  hits_real_old, "Ratio of old vs new mean PMTs hits per floor seperated by PMT ring, strings 13-27 only", 
    indices, pmt_letters, plot_string=False, slice_string = True, save_map = "plots/new-vs-old", start_index = 6, stop_index = 17)

summarised_heatmap_ratio(hits_real_new,  hits_real_old, "Ratio of old vs new mean PMTs hits per floor seperated by PMT ring, all strings", 
    indices, pmt_letters, plot_string=False, slice_string = True, save_map = "plots/new-vs-old")

summarised_heatmap_ratio(eff_new,  eff_old, "Ratio of old vs new mean PMTs efficiencies per floor seperated by PMT ring, strings 13-27 only", 
    indices, pmt_letters, plot_string=False, slice_string = True, save_map = "plots/new-vs-old", start_index = 6, stop_index = 17)

summarised_heatmap_ratio(eff_new,  eff_old, "Ratio of old vs new mean PMTs efficiencies per floor seperated by PMT ring, all strings", 
    indices, pmt_letters, plot_string=False, slice_string = True, save_map = "plots/new-vs-old")
