from efficiency9 import * 


eff_all = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 5)


dist_plots.plot_dist(eff_all)