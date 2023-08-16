from efficiency9 import *

#Collar shadows PMTs C2, C5, E2, E5. Furthermore more dark rates for PMTs A1, B5, B6 and equator tape shadowing for PMTs D, E.
#TODO create filter for said PMTs

eff_all = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 5, apply_shadow_mask=False)



eff_all.plot_heatmap_summarised_ring(indices, pmt_letters, title="Efficiencies per string seperated by PMT ring, all except shadowed PMTs", 
    save_map = "plots/gel-issue", include_mean = True, include_mean_of_mean = True)


