#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon July  3  2023

@author: mike

This file should also take into account the two different PMT versions and include a mapping per ring in the DOM. 
"""

from efficiency8 import *

data_real_old = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 0)
data_real_new = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, new_versions = 1)
data_sim = map_hit_data(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number)
data_real = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number)

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

sim_heatmap.plot_heatmap_summarised_ring(indices, pmt_letters, "Mean of the sum of all hits per DOM per string, simulated data", 
        save_map = "plots/floor-string-mean")
sim_heatmap.plot_heatmap_summarised_ring(indices, pmt_letters, "Mean of the sum of all hits per DOM per floor, simulated data", 
         save_map = "plots/floor-string-mean", string = False)

real_over_sim_heatmap.plot_heatmap_summarised_ring(indices, pmt_letters, "Mean of the sum of all hits per DOM per string, real over MC ratio", 
         save_map = "plots/floor-string-mean")
real_over_sim_heatmap.plot_heatmap_summarised_ring(indices, pmt_letters, "Mean of the sum of all hits per DOM per floor, real over MC ratio", 
         save_map = "plots/floor-string-mean", string = False)

real_eff_heatmap.plot_heatmap_summarised_ring(indices, pmt_letters, "Mean of the efficiencies per DOM per string", 
         save_map = "plots/floor-string-mean")
real_eff_heatmap.plot_heatmap_summarised_ring(indices, pmt_letters, "Mean of the efficiencies per DOM per floor", 
         save_map = "plots/floor-string-mean", string = False)

#eff_real_over_sim_heatmap.plot_heatmap_summarised_ring(indices, pmt_letters, "Mean of the efficiencies per DOM per string, real over MC ratio", 
#         save_map = "plots/floor-string-mean")
#eff_real_over_sim_heatmap.plot_heatmap_summarised_ring(indices, pmt_letters, "Mean of the efficiencies per DOM per floor, real over MC ratio", 
#         save_map = "plots/floor-string-mean", string = False)
        
#Create new dataset by laying the efficiency map over the ratio map 

sim_ratio_map.plot_heatmap(indices, pmt_letters, "Ratio of real vs simulated rates of PMTs per DOM", save_map = "plots/sim-vs-real")
eff_heatmap.plot_heatmap(indices, pmt_letters, "Efficiencies of PMTs per DOM, real data", save_map = "plots/sim-vs-real")


sim_ratio_eff_map = np.zeros([2,sim_eff_map.shape[0], sim_eff_map.shape[1], sim_eff_map.shape[2]])
sim_ratio_eff_map[0, :, :, :] = sim_ratio_map.heatmap[:, :, :]
sim_ratio_eff_map[1, :, :, :] = eff_heatmap.heatmap[:, :, :]

letters_ints = np.arange(0, len(pmt_letters))


old_pmt_mean, old_pmt_std = real_old_heatmap.get_avg_std(indices)
new_pmt_mean, new_pmt_std = real_new_heatmap.get_avg_std(indices)

sim_old_pmt_mean, sim_old_pmt_std = sim_old_heatmap.get_avg_std(indices)
sim_new_pmt_mean, sim_new_pmt_std = sim_new_heatmap.get_avg_std(indices)

#plot_ratio_eff_per_group(sim_ratio_eff_map, indices, "Efficiencies vs ratio of real vs sim rates of PMTs per DOM", 
#        stringlist, save_map = "plots/eff-vs-hits")

#Now create subplots 

# fig, axs = plt.subplots(2, 2, figsize = (20, 20))
# plt.setp(axs, xticks = letters_ints, xticklabels = pmt_letters, ylabel= "mean hit rate", xlabel="PMT ring", ylim = (-750, 3200))
# axs[0, 0].errorbar(letters_ints, old_pmt_mean, yerr = old_pmt_std, fmt="o")
# axs[1, 0].errorbar(letters_ints, new_pmt_mean, yerr = new_pmt_std, fmt="o")
# axs[0, 1].errorbar(letters_ints, sim_old_pmt_mean, yerr = sim_old_pmt_std, fmt="o", color="orange")
# axs[1, 1].errorbar(letters_ints, sim_new_pmt_mean, yerr = sim_new_pmt_std, fmt="o", color="orange")

# fig.suptitle("Mean of the sum of the DOM hit rate seperated by ring") 
# axs[0, 0].set_title("Real data, old PMTs (model R12199)")
# axs[0, 1].set_title("Simulated data, old PMTs (model R12199)")
# axs[1, 0].set_title("Real data, new PMTs (model R14374)")
# axs[1, 1].set_title("Simulated data, new PMTs (model R14374)")
# plt.subplots_adjust(hspace=0.4)
# plt.savefig("plots/mean-sum-doms-pmt-models.pdf")
# plt.clf()
# plt.close()

# #Now do the same to show difference between the two types PMTs 
# fig, axs = plt.subplots(2, 2, figsize = (20, 20))
# plt.setp(axs, xticks = letters_ints, xticklabels = pmt_letters, ylabel= "mean hit rate", xlabel="PMT ring", ylim = (-2000, 3200))
# axs[0, 0].errorbar(letters_ints, np.abs(old_pmt_mean - new_pmt_mean), yerr = err_prop_sum(old_pmt_std,new_pmt_std),  fmt="o")
# axs[0, 1].errorbar(letters_ints, np.abs(old_pmt_mean - sim_old_pmt_mean), yerr = err_prop_sum(old_pmt_std, sim_old_pmt_std), fmt="o")
# axs[1, 0].errorbar(letters_ints, np.abs(new_pmt_mean - sim_new_pmt_mean), yerr = err_prop_sum(sim_new_pmt_std, sim_new_pmt_std), fmt="o")
# axs[1, 1].errorbar(letters_ints, np.abs(sim_old_pmt_mean - sim_new_pmt_mean), yerr = err_prop_sum(sim_old_pmt_std,sim_new_pmt_std),  fmt="o")

# fig.suptitle("Differences between the mean of the sum of the DOM hit rate seperated by ring") 
# axs[0, 0].set_title("Real data, old and new PMTs compared")
# axs[0, 1].set_title("Real and simulated data compared, old PMTs")
# axs[1, 0].set_title("Real and simulated data compared, new PMTs")
# axs[1, 1].set_title("Simulated data, old and new PMTs compared")
# plt.subplots_adjust(hspace=0.4)
# plt.savefig("plots/mean-diff-sum-doms-pmt-models.pdf")

