#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon July  3  2023

@author: mike

This file should also take into account the two different PMT versions and include a mapping per ring in the DOM. 
"""

from efficiency8 import *

data_real = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 0)
data_real_new = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 1)
data_sim = map_hit_data(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number)

#data_real.apply_mask_return_floor_str_hit()




real_eff_map, __, __, = data_real.export_heatmap(indices, int_rates_or_eff = 5)
sim_map, __, __ = data_sim.export_heatmap(indices)
sim_eff_map, __, __ = data_sim.export_heatmap(indices, int_rates_or_eff =5)

real_map, floorlist, stringlist = data_real.export_heatmap(indices)
real_map_new_pmts, __, __ = data_real_new.export_heatmap(indices)







sim_ratio_map = heatmap(calc_heatmap_ratio(real_map, sim_map), floorlist, stringlist)
real_heatmap = heatmap(real_map, floorlist, stringlist)
eff_heatmap = heatmap(sim_eff_map, floorlist, stringlist)
real_new_heatmap = heatmap(real_map_new_pmts, floorlist, stringlist)

#eff_heatmap.plot_heatmap(indices, "Map of the efficiencies")
#sim_eff_heatmap.plot_heatmap(indices, "Average efficiencies of DOMs, simulated data")

#sim_ratio_map.compare_upper_lower_pmts_heatmap(indices, "Ratio of upper/lower PMTs ratio of simulated/real rates")
#sim_ratio_map.plot_heatmap(indices, "Ratio of simulated vs real rates of PMTs per DOM")

#real_heatmap.plot_heatmap(indices, pmt_letters, "Sum of all hits per DOM for indicated PMT group, real data", save_map = "plots")
#real_new_heatmap.plot_heatmap(indices, pmt_letters, "Sum of all hits per DOM for indicated PMT group, real data, new PMTs", save_map = "plots")

#Create new dataset by laying the efficiency map over the ratio map 


print(real_heatmap.get_avg_std(indices))
print(real_new_heatmap.get_avg_std(indices))
sim_ratio_eff_map = np.zeros([2,sim_eff_map.shape[0], sim_eff_map.shape[1], sim_eff_map.shape[2]])
sim_ratio_eff_map[0, :, :, :] = sim_ratio_map.heatmap[:, :, :]
sim_ratio_eff_map[1, :, :, :] = eff_heatmap.heatmap[:, :, :]

letters_ints = np.arange(0, len(pmt_letters))


old_pmt_mean, old_pmt_std = real_heatmap.get_avg_std(indices)
new_pmt_mean, new_pmt_std = real_new_heatmap.get_avg_std(indices)
#plt.errorbar(letters_ints, old_pmt_mean, yerr = old_pmt_std, fmt="o")
#plt.errorbar(letters_ints, new_pmt_mean, yerr = new_pmt_std, fmt=".")
#plt.xticks(np.arange(0, len(letters_ints)), pmt_letters)
#plt.show()

#Now create subplots 

fig, axs = plt.subplots(2, 1, sharey=True)
plt.setp(axs, xticks = letters_ints, xticklabels = pmt_letters, ylabel= "mean hit rate")
axs[0].errorbar(letters_ints, old_pmt_mean, yerr = old_pmt_std, fmt="o")
axs[1].errorbar(letters_ints, new_pmt_mean, yerr = new_pmt_std, fmt="o")

fig.suptitle("Mean of the sum of the DOM hit rate seperated by ring") 
axs[0].set_title("Old PMTs (model R12199)")
axs[1].set_title("New PMTs (model R14374)")
plt.subplots_adjust(hspace=0.4)
plt.savefig("plots/mean-sum-doms-pmt-models.pdf")
plt.show()

