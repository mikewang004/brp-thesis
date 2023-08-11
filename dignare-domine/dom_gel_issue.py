from efficiency9 import *
import numpy.ma as ma

eff_all = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 5).heatmap




#DOMs with different gel should have lower efficiencies and we suspect in any case that string 10.7-15 and 12 suffer from this issue. 
#So filter first for those. 

# heatmap, floorlist, stringlist = eff_all.plot_heatmap(indices, pmt_letters, title="Efficiencies per string seperated by PMT ring, all PMTs", 
#     save_map = "plots/efficiencies/rings", include_mean = True)

# eff_all.plot_heatmap_summarised_ring(indices, pmt_letters, "Efficiencies per floor seperated by PMT ring", 
#     save_map = "plots/efficiencies", string=False, include_mean = True, include_mean_of_mean = True)




print(eff_all.shape)
heatmap_all_pmts = np.nanmean(eff_all[4:, :, :], axis = 0) #2d heatmap of rings E-F 

print(stringlist, floorlist)
# plot_heatmap_ultra_basic(heatmap_all_pmts, "Efficiencies per DOM, upper half PMTs only", 
# stringlist, floorlist, save_map = "plots/efficiencies", include_mean=True)


print(np.mean(heatmap_all_pmts[6:15, 2]))
print(np.nanmean(heatmap_all_pmts))
avg_mean_sus_strings = (np.nanmean(heatmap_all_pmts[6:15, 2]) + np.nanmean(heatmap_all_pmts[:, 4]))/2 + np.sqrt(0.25 * np.nanstd(heatmap_all_pmts[6:15, 2])**2 + 0.25 * np.nanstd(heatmap_all_pmts[:, 4])**2)
#print(avg_mean_sus_strings)
gel_mask = ma.masked_greater(heatmap_all_pmts, avg_mean_sus_strings)
gel_mask = ma.filled(gel_mask, np.nan)
inv_gel_mask = ma.masked_less(heatmap_all_pmts, avg_mean_sus_strings)
inv_gel_mask = ma.filled(inv_gel_mask, np.nan)


# plot_heatmap_ultra_basic(gel_mask, "Efficiencies per suspected gel-issue DOM, upper half PMTs only", 
#     stringlist, floorlist, save_map = "plots/gel-issue", include_mean=True)

#Now plot the distribution 

dist_eff_all = dist_plots(eff_all)
dist_gel_mask = dist_plots(gel_mask); dist_inv_gel_mask = dist_plots(inv_gel_mask)
dist_eff_all.plot_dist(xlabel = "efficiency", title = "Distribution efficiencies upper PMTs", save_map = "plots/gel-issue")
#dist_gel_mask.plot_dist(xlabel = "efficiency", title = "Distribution efficiencies upper PMTs, suspected gel issue", 
# save_map = "plots/gel-issue")
dist_inv_gel_mask.plot_dist(xlabel = "efficiency", title = "Distribution efficiencies upper PMTs, no suspected gel issue", 
    save_map = "plots/gel-issue")


#TODO correlate these results to rate vs efficiency plots. 