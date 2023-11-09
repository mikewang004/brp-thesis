from efficiency9 import *

#start_index = 6; stop_index = 17
start_index = 7; stop_index = 16 


data_real_old = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 0)
data_real_new = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 1)

hits_real_old= speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 0)
hits_real_new = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 1)
eff_old = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 0, int_rates_or_eff=5)
eff_new = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 1, int_rates_or_eff=5)
eff_all = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 5)

hits_real_all = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number)
hits_sim_all = speedrun_heatmap(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number)
hits_ratio_all = heatmap(hits_real_all.heatmap/hits_sim_all.heatmap)

# hits_real_all.plot_heatmap(indices, pmt_letters, "Mean hits per DOM", save_map = "plots/hits/rings")
hits_ratio_all.plot_heatmap(indices, pmt_letters, "Ratio of data vs MC hits per DOM", save_map = "plots/data-mc/rings", 
    include_mean=True)
# hits_real_old.plot_heatmap_summarised_ring(indices, pmt_letters, "Mean hits per string seperated by PMT ring, old PMTs", 
#            save_map = "plots/hits", string=False)
# hits_real_new.plot_heatmap_summarised_ring(indices, pmt_letters, "Mean hits per string seperated by PMT ring, new PMTs", 
#            save_map = "plots/hits", string=False)


# eff_old.plot_heatmap(indices, pmt_letters, title="Efficiencies per string seperated by PMT ring, old PMTs", 
#     save_map = "plots/data-mc/rings", include_mean= True)
# eff_new.plot_heatmap(indices, pmt_letters, title="Efficiencies per string seperated by PMT ring, new PMTs", 
#     save_map = "plots/data-mc/rings", include_mean = True)

eff_all.plot_heatmap(indices, pmt_letters, title="Efficiencies per string seperated by PMT ring, all PMTs", 
    save_map = "plots/efficiencies/rings", include_mean = True)
# eff_old.plot_heatmap_summarised_ring(indices, pmt_letters, "Efficiencies per string seperated by PMT ring, old PMTs",
#     save_map = "plots/efficiencies", include_mean = True, include_mean_of_mean= True)
# eff_new.plot_heatmap_summarised_ring(indices, pmt_letters, "Efficiencies per string seperated by PMT ring, new PMTs",
#     save_map = "plots/efficiencies", include_mean = True, include_mean_of_mean=True)



# summarised_heatmap_ratio(hits_real_new,  hits_real_old, "Ratio of old vs new mean PMTs hits per floor seperated by PMT ring, strings 14-26 only", 
#     indices, pmt_letters, plot_string=False, slice_string = True, save_map = "plots/new-vs-old", 
#     start_index = start_index, stop_index = stop_index)

# summarised_heatmap_ratio(hits_real_new,  hits_real_old, "Ratio of old vs new mean PMTs hits per floor seperated by PMT ring, all strings", 
#     indices, pmt_letters, plot_string=False, slice_string = True, save_map = "plots/new-vs-old")

# summarised_heatmap_ratio(eff_new,  eff_old, "Ratio of old vs new mean PMTs efficiencies per floor seperated by PMT ring, strings 14-26 only", 
#     indices, pmt_letters, plot_string=False, slice_string = True, save_map = "plots/new-vs-old", 
#     start_index = start_index, stop_index = stop_index)

# summarised_heatmap_ratio(eff_new,  eff_old, "Ratio of old vs new mean PMTs efficiencies per floor seperated by PMT ring, all strings", 
#     indices, pmt_letters, plot_string=False, slice_string = True, save_map = "plots/new-vs-old")


# summarised_heatmap_ratio(hits_real_all,  hits_sim_all, "Ratio of real vs simulated PMTs hits per string seperated by PMT ring", 
#     indices, pmt_letters, plot_string=True, slice_string = False, save_map = "plots/data-mc")

# summarised_heatmap_ratio(hits_real_all,  hits_sim_all, "Ratio of real vs simulated PMTs hits per floor seperated by PMT ring", 
#     indices, pmt_letters, plot_string=False, slice_string = False, save_map = "plots/data-mc")