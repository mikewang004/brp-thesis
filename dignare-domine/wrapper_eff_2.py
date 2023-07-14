#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed July  12  2023

@author: mike

This file should also take into account the two different PMT versions and include a mapping per ring in the DOM. 
"""

from efficiency8 import *

data_real_old = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 0)
data_real_new = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 1)
data_sim = map_hit_data(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number)
data_real = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number)c

data_sim_old = map_hit_data(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number, new_versions = 0)
data_sim_new = map_hit_data(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number, new_versions=1)






real_eff_map, __, __, = data_real.export_heatmap(indices, int_rates_or_eff = 5)
sim_map, __, __ = data_sim.export_heatmap(indices)
real_map, __, __ =data_real.export_heatmap(indices)
sim_eff_map, __, __ = data_sim.export_heatmap(indices, int_rates_or_eff =5)

real_map_old_pmts, floorlist, stringlist = data_real_old.export_heatmap(indices)
real_map_new_pmts, __, __ = data_real_new.export_heatmap(indices)

sim_old_map, __, __ = data_sim_old.export_heatmap(indices)
sim_new_map, __, __ = data_sim_new.export_heatmap(indices)





sim_ratio_map = heatmap(calc_heatmap_ratio(real_map, sim_map), floorlist, stringlist)

sim_heatmap = heatmap(sim_map, floorlist, stringlist)
real_heatmap = heatmap(real_map, floorlist, stringlist)
eff_heatmap = heatmap(sim_eff_map, floorlist, stringlist)
real_eff_heatmap = heatmap(real_eff_map, floorlist, stringlist)

real_old_heatmap = heatmap(real_map, floorlist, stringlist)
real_new_heatmap = heatmap(real_map_new_pmts, floorlist, stringlist)

sim_old_heatmap = heatmap(sim_old_map, floorlist, stringlist)
sim_new_heatmap = heatmap(sim_new_map, floorlist, stringlist)

real_over_sim_heatmap = heatmap(real_map/sim_map, floorlist, stringlist)
eff_real_over_sim_heatmap = heatmap(real_eff_map/sim_eff_map, floorlist, stringlist)

#real_heatmap.plot_heatmap(indices, pmt_letters, "Sum of all hits per DOM for indicated PMT group, real data, old PMTs", save_map = "plots")
#real_new_heatmap.plot_heatmap(indices, pmt_letters, "Sum of all hits per DOM for indicated PMT group, real data, new PMTs", save_map = "plots")
# zmax_array, zmin_array = sim_heatmap.plot_heatmap(indices, pmt_letters, 
#         "Sum of all hits per DOM for indicated PMT group, sim. data", save_map = "plots/all-pmts", include_mean= True)
# real_heatmap.plot_heatmap(indices, pmt_letters, 
#         "Sum of all hits per DOM for indicated PMT group, real data", 
#         save_map = "plots/all-pmts", zmax_array= zmax_array, zmin_array = zmin_array, include_mean= True)


real_heatmap.plot_heatmap_summarised_ring(indices, pmt_letters, "Mean of the sum of all hits per DOM per string, real data", 
        save_map = "plots/floor-string-mean")
real_heatmap.plot_heatmap_summarised_ring(indices, pmt_letters, "Mean of the sum of all hits per DOM per floor, real data", 
         save_map = "plots/floor-string-mean", string = False)



real_eff_heatmap.plot_heatmap_summarised_ring(indices, pmt_letters, "Mean of the efficiencies per DOM per string", 
         save_map = "plots/floor-string-mean")
real_eff_heatmap.plot_heatmap_summarised_ring(indices, pmt_letters, "Mean of the efficiencies per DOM per floor", 
         save_map = "plots/floor-string-mean", string = False)



# #Now plot the efficiencies vs the mean sum in one scatter and try to draw something of conclusions from there 


# plt.figure()
# real_summarised_heatmap, __ = real_heatmap.summarise_per_ring(indices, pmt_letters)
# eff_summarised_heatmap, __ = eff_heatmap.summarise_per_ring(indices, pmt_letters)
# for i in range(0, len(indices)-1):
#         plt.scatter(real_summarised_heatmap[i, :], eff_summarised_heatmap[i, :], label=pmt_letters[i])
# plt.title("Efficiencies vs hits per ring")
# plt.xlabel("mean hits per string")
# plt.ylabel("mean efficiency per string")
# plt.legend()
# plt.show()


# Try to draw some conclusions from upper rings vs lower rings 