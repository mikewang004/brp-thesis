from efficiency9 import *
import numpy.ma as ma

eff_all = speedrun_heatmap(muon_hit_data_real, modid_map, pmt_serial_map, magic_number, int_rates_or_eff = 5)


#DOMs with different gel should have lower efficiencies and we suspect in any case that string 10.7-15 and 12 suffer from this issue. 
#So filter first for those. 

heatmap, floorlist, stringlist = eff_all.plot_heatmap(indices, pmt_letters, title="Efficiencies per string seperated by PMT ring, all PMTs", 
    save_map = "plots/efficiencies/rings", include_mean = True)

# eff_all.plot_heatmap_summarised_ring(indices, pmt_letters, "Efficiencies per floor seperated by PMT ring", 
#     save_map = "plots/efficiencies", string=False, include_mean = True, include_mean_of_mean = True)


heatmap_all_pmts = np.mean(heatmap, axis = 0)

#plot_heatmap_ultra_basic(heatmap_all_pmts, "Efficiencies per DOM", stringlist, floorlist, save_map = "plots/efficiencies", include_mean=True)


print(np.mean(heatmap_all_pmts[6:15, 2]))
print(np.nanmean(heatmap_all_pmts[:, 4]))


gel_mask = ma.masked_greater(heatmap_all_pmts, np.nanmean(heatmap_all_pmts[:, 2]))
gel_mask = ma.filled(gel_mask, np.nan)


#plot_heatmap_ultra_basic(gel_mask, "Efficiencies per suspected gel-issue DOM", stringlist, floorlist, save_map = "plots/gel-issue", include_mean=True)


# Try alternative approach: plot 1d distributions of the efficiencies. 


#for i in range(0, heatmap_all_pmts.shape[1]):
#for i in range(0, 8):
    #plt.hist(heatmap_all_pmts[:, i], bins = 40, range = (0.7, 1.1), label = "string %i" %stringlist[i], histtype='bar', stacked=True)
plt.hist(heatmap_all_pmts[:, :], bins = 40, range = (0.7, 1.1), histtype='bar', stacked=True, label = )
plt.legend()
plt.show()

#TODO correlate these results to rate vs efficiency plots. 