from efficiency9 import *
import numpy.ma as ma

eff_all = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 5).heatmap

eff_all_31 = speedrun_heatmap_31(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 5).heatmap


#DOMs with different gel should have lower efficiencies and we suspect in any case that string 10.7-15 and 12 suffer from this issue. 
#So filter first for those. 

# heatmap, floorlist, stringlist = eff_all.plot_heatmap(indices, pmt_letters, title="Efficiencies per string seperated by PMT ring, all PMTs", 
#     save_map = "plots/efficiencies/rings", include_mean = True)

# eff_all.plot_heatmap_summarised_ring(indices, pmt_letters, "Efficiencies per floor seperated by PMT ring", 
#     save_map = "plots/efficiencies", string=False, include_mean = True, include_mean_of_mean = True)




print(eff_all.shape)
heatmap_all_pmts = np.nanmean(eff_all[:4, :, :], axis = 0) #2d heatmap of rings E-F 

print(stringlist, floorlist)
# plot_heatmap_ultra_basic(heatmap_all_pmts, "Efficiencies per DOM, upper half PMTs only", 
# stringlist, floorlist, save_map = "plots/efficiencies", include_mean=True)


print(np.mean(heatmap_all_pmts[6:15, 2])) #should be 0.8448
print(np.nanmean(heatmap_all_pmts))
avg_mean_sus_strings = (np.nanmean(heatmap_all_pmts[6:15, 2]) + np.nanmean(heatmap_all_pmts[:, 4]))/2 + 2*(np.sqrt(0.25 * np.nanstd(heatmap_all_pmts[6:15, 2])**2 + 0.25 * np.nanstd(heatmap_all_pmts[:, 4])**2))
print(avg_mean_sus_strings)
gel_mask = ma.masked_greater(heatmap_all_pmts, avg_mean_sus_strings)
gel_mask = ma.filled(gel_mask, np.nan)
inv_gel_mask = ma.masked_less(heatmap_all_pmts, avg_mean_sus_strings)
inv_gel_mask = ma.filled(inv_gel_mask, np.nan)



# plot_heatmap_ultra_basic(gel_mask, "Efficiencies per no gel-issue DOM, lower half PMTs only", 
#     stringlist, floorlist, save_map = "plots/gel-issue", include_mean=True)
lower_bound = 0.75; upper_bound = 0.91

#Now plot the distribution 
eff_all_27 = eff_all_31[:, :, 17]; cond_27 = (eff_all_27 > lower_bound) & (eff_all_27 < upper_bound)
eff_all_27 = eff_all_27[cond_27]
dist_mask_10 = dist_plots(eff_all_31[:18, 6:15, 2])
dist_mask_12= dist_plots(eff_all_31[:18, :, 4])
dist_mask_27= dist_plots(eff_all_27)
#dist_gel_mask = dist_plots(gel_mask); dist_inv_gel_mask = dist_plots(inv_gel_mask)

bins = 20; range = (lower_bound, upper_bound)
dist_mask_10.plot_dist_barebones(num_bins = bins, label = "string 10", range = range)
dist_mask_12.plot_dist_barebones(num_bins = bins, label = "string 12", range = range)
dist_mask_27.plot_dist_barebones(num_bins = bins, label = "string 27", range = range)
plt.xlabel("efficiency")
plt.ylabel("count")
plt.legend()
dist_mask_27.plot_dist_save("Distribution of efficiencies, suspected gel issue DOMs", "plots/gel-issue")
#dist_mask_27.plot_dist_save("Distribution of efficiencies, strings 10 and 27", "plots/gel-issue")


#Perform K-S test

test = stats.ks_2samp(dist_mask_27.counts, dist_mask_12.counts)
print(test)
print(stats.ks_2samp(dist_mask_27.counts, dist_mask_10.counts))

#TODO correlate these results to rate vs efficiency plots. 