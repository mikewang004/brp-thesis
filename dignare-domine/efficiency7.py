#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 15:54:39 2023

@author: mike

This file should also take into account the two different PMT versions and include a mapping per ring in the DOM. 
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as go

pio.kaleido.scope.mathjax= None
pio.renderers.default='browser'

muon_hit_data_sim = np.load("muon_hit_data-sim-reduced_bins-xx1375x.npy")
#muon_hit_data_real = np.load("muon_hit_data-real-reduced_bins-13754.npy")
muon_hit_data_real = np.load("muon_hit_data-real-reduced_bins-xx1375x.npy")
modid_map = np.loadtxt("map.txt")
eff_list = np.loadtxt("../zee-haarzelf/data-133-144-eff.txt", skiprows = 148, usecols=[1,2,3])
pmt_serial_map = np.loadtxt("../pmt-info/pmt-serials.txt", usecols = 0)
pmt_ring_map = np.loadtxt("../pmt-info/pmt-ring.txt", skiprows = 2, usecols = [0,1,2])
magic_number = 16104 # The major version change happened at serial number 16104 (all PMTs <=16104 are of a certain kind (R12199), all abover are another one (R14374)).
ijk = 0


class map_hit_data():
    """Connects the muon hit data only identified per identifier to a map containing the key to which floor/string it is located in."""
    """Note initial data muon_hit_data is 2d array w/ column headers module-id / pmt number / amount of hits.
    If needed to append more data just append columns next to the last column as it should survive all operations."""

    """Workflow as follows: data loadin is as 2d array with column as above. Then data gets appended, upon which the array gets transformed to 3d 
    with the pmt number being the 3rd dimension. Then data reorganised into heatmap."""
    def __init__(self, muon_hit_data, pmt_id_map, pmt_serial_map, magic_number, floor_str_hit = None):
        self.modid_map = modid_map
        self.eff_list = eff_list    
        self.pmt_serial_map = pmt_serial_map; self.magic_number = magic_number
        if floor_str_hit == None:
            self.muon_hit_data = muon_hit_data
            self.append_eff_data()
            self.append_pmt_serials()
            self.mod_id_to_floor_string()
            self.pmt_no_to_ring_letter()
        else:
            self.floor_str_hit = floor_str_hit

    def append_eff_data(self):
        """Couples the new efficiency data to a DOM-module number."""
        new_muon_hit_data = np.zeros([self.muon_hit_data.shape[0], self.muon_hit_data.shape[1], self.muon_hit_data.shape[2]+1]) #should be ok with other code 
        new_muon_hit_data[:, :, :3] = self.muon_hit_data
        for i in range(0, len(self.muon_hit_data)):
            #Look up which data corresponds in the efficiency list
            for j in range(0, len(self.eff_list)):
                if self.muon_hit_data[i, 0, 0] == self.eff_list[j, 0]:
                    for k in range(0, 31):
                        new_muon_hit_data[i, k, 3] = self.eff_list[j + k, 2]
                    break
        self.muon_hit_data = new_muon_hit_data
        return 0;

    def append_pmt_serials(self):
        """Couples PMT-serials to respective PMTs."""
        a = 0
        new_muon_hit_data = np.zeros([self.muon_hit_data.shape[0], self.muon_hit_data.shape[1], self.muon_hit_data.shape[2]+2])
        new_muon_hit_data[:,:,:4] = self.muon_hit_data
        for i in range(0, len(self.muon_hit_data)):
            for j in range(0, len(self.pmt_serial_map)):
                if self.muon_hit_data[i, 0, 0] == self.pmt_serial_map[j]:
                    for k in range(0, 31):
                        new_muon_hit_data[i, k, 4] = self.pmt_serial_map[j + k+1]
                        if self.pmt_serial_map[j + k + 1] > magic_number:
                            new_muon_hit_data[i, k, 5] = 1
                            a = a + 1
                        else:
                            new_muon_hit_data[i, k, 5] = 0
                    break
        self.muon_hit_data = new_muon_hit_data
        return 0;


    def mod_id_to_floor_string(self):
        """Transforms data from a 2D to a 3D array including the different pmts. 
        Note output is of the form: amount of mod-ids / pmts numbers / [str no; floor no; mod-id; pmt-no; no. hits; eff; pmt_serial; pmt_old/new]
        Also note that mod-id refers to which DOM one is looking at"""
        floor_str_hit = np.zeros([self.muon_hit_data.shape[0], self.muon_hit_data.shape[1], self.muon_hit_data.shape[2]+2])
        mapping = dict(zip(self.modid_map[:,0], range(len(modid_map))))
        for i in range(0, self.muon_hit_data.shape[1]):
            floor_str_hit[:, i, :] = np.hstack((np.array([self.modid_map[mapping[key],1:] for key in self.muon_hit_data[:,i,0]]), self.muon_hit_data[:, i, :]))
        self.floor_str_hit = floor_str_hit
        return 0;

    def pmt_no_to_ring_letter(self):
        """Transforms the number of each PMT to a ring location."""
        print(self.floor_str_hit[0, :, 3])
        for i in range(0, 31):
            self.floor_str_hit[:, i, 3] = pmt_ring_map[i, 1]
        return 0;


    def sum_over_n_pmts(self, indices):
        """Calculates average over any [n] groups of pmts. Group must include starting PMT-no. of group (inclusive) and stopping number (exclusive)"""
        pmt_group_mean = np.zeros([self.floor_str_hit.shape[0], len(indices) - 1, self.floor_str_hit.shape[2]])
        j = 0; k = 0
        print(pmt_group_mean.shape)
        for i in range(1, max(indices)+1):
            if i in indices:
                pmt_group_mean[:, k, 4] = np.nansum(self.floor_str_hit[:, j:i, 4], axis=1)
                pmt_group_mean[:, k, :3] = np.nanmean(self.floor_str_hit[:, j:i, :3], axis = 1)
                pmt_group_mean[:, k, 5:] = np.nanmean(self.floor_str_hit[:, j:i, 5:], axis = 1)
                j = i + 1; k = k + 1
        #np.savetxt("debug-avgpmts12.txt", pmt_group_mean[:, 1, :])
        return pmt_group_mean

    def normalise_over_n_pmts(self, indices):
        """Calculates average over any [n] groups of pmts. Group must include starting PMT-no. of group (inclusive) and stopping number (exclusive)"""
        pmt_group_mean = np.zeros([self.floor_str_hit.shape[0], len(indices) - 1, self.floor_str_hit.shape[2]])
        j = 0; k = 0
        print(pmt_group_mean.shape)
        for i in range(1, max(indices)+1):
            if i in indices:
                pmt_group_mean[:, k, :] = np.nanmean(self.floor_str_hit[:, j:i, :], axis=1)
                j = i + 1; k = k + 1
        #np.savetxt("debug-avgpmts12.txt", pmt_group_mean[:, 1, :])
        return pmt_group_mean


    def heatmap_averages_single_loop(self, pmt_group_pairs):
        #Sort pmt_groups on str and floor
        pmt_group_pairs = pmt_group_pairs[pmt_group_pairs[:, 0].argsort()] 
        #First sort according to string; then according to floor 
        k = 0; l = 0;
        str_floor_length = np.zeros([len(np.unique(pmt_group_pairs[:, 0]))]) #details how many floors in a string; 
        for i in range(1, pmt_group_pairs.shape[0]):
        #for i in range(0, 100):
            if pmt_group_pairs[i, 0] != pmt_group_pairs[i-1, 0] or i == (pmt_group_pairs.shape[0]-1):
                #pmt_group_pairs[k:i, :] = pmt_group_pairs[pmt_group_pairs[k:i, 1].argsort()]
                #print(pmt_group_pairs[k:i, :])
                #Apply new indexing to correct part of array 
                if i == (pmt_group_pairs.shape[0]-1):
                    aux_array = pmt_group_pairs[k:, :]
                    aux_array = aux_array[aux_array[:, 1].argsort()]
                    pmt_group_pairs[k:, :] = aux_array
                else:
                    aux_array = pmt_group_pairs[k:i, :]
                    aux_array = aux_array[aux_array[:, 1].argsort()]
                    pmt_group_pairs[k:i, :] = aux_array
                str_floor_length[l] = i - k
                k = i; l = l + 1
        np.savetxt("debug-pmt-group-pairs-%i.txt" %(ijk), pmt_group_pairs)
        return pmt_group_pairs, str_floor_length

    def heatmap_array_single_group(self, pmt_group_mean_sorted, pmt_group_no, heatmap, floorlist, stringlist, int_rates_or_eff):
        """int_rates_or_eff = 4 for rates; =5 for efficiencies"""
        k = 0; m = 0; x = 0; j = 0; l = 0; #n = 4 for rates; n = 5 for efficiencies 
        for i in range(1, pmt_group_mean_sorted.shape[0]): #first fill in the x/string direction
            #New approach: floorlist index is not the same as pmt group mean sorted index. 
            if pmt_group_mean_sorted[i, 0] != pmt_group_mean_sorted[i-1, 0] or i == pmt_group_mean_sorted.shape[0]-1:
                while j < len(floorlist): #checks if all floor/string combinations are extant
                    if floorlist[j] == pmt_group_mean_sorted[k + l, 1]: #case if combination is extant 
                        heatmap[j, m] = pmt_group_mean_sorted[k + l, int_rates_or_eff]
                        l = l + 1
                    else: #case if combination non-extant 
                        heatmap[j, m] = np.nan
                        x = x + 1
                    j = j + 1
                m = m + 1; k = i; j = 0; l = 0
        return heatmap

    def export_heatmap(self, indices, int_rates_or_eff=4):
        #pmt_group_mean = self.normalise_over_n_pmts(indices)
        pmt_group_mean = self.sum_over_n_pmts(indices)
        pmt_group_mean_sorted = pmt_group_mean
        #pmt_group_pairs = pmt_group_mean[:, 0, :]
        for m in range(0, len(indices)-1):
            pmt_group_pairs = pmt_group_mean[:, m, :]
            pmt_group_mean_sorted[:, m, :], str_floor_length = self.heatmap_averages_single_loop(pmt_group_pairs) #Sorts the thing on string and floor so that they are in sequence. 
        floorlist = np.unique(pmt_group_mean_sorted[:, 0, 1]); stringlist = np.unique(pmt_group_mean_sorted[:, 0, 0])
        heatmap = np.zeros([len(indices)-1, len(floorlist), len(np.unique(stringlist))])
        #Now get heatmap for each group
        for i in range(0, len(indices)-1):
            heatmap[i, :, :] = self.heatmap_array_single_group(pmt_group_mean_sorted[:, i, :], i, heatmap[i, :, :], floorlist, stringlist, int_rates_or_eff)
        np.savetxt("pmt-groups-mean-sorted-debug.txt", pmt_group_mean_sorted[:, 0, :])
        return heatmap, floorlist, stringlist 
    


class heatmap():
    """Class to generate heatmap plots with string and/or floor information. Also useful to perform heatmap operations with."""
    def __init__(self, heatmap, floorlist, stringlist):
        self.heatmap = heatmap
        self.floorlist = floorlist
        self.stringlist = stringlist

    def plot_heatmap(self, indices, pmt_letters, title):
    
        custom_colorscale = [
            [1, 'rgb(255, 255, 0)'],
            [0, 'rgb(0, 0, 255)']
        ]
        for i in range(0, len(indices)-1):
            annotation_text = np.round(self.heatmap[i, :, :], 4)
            layout = go.Layout(
                title = title + ", PMT group %s" %((pmt_letters[i])),
                xaxis = dict(
                    tickmode = "array",
                    tickvals = np.arange(len(self.stringlist)),
                    ticktext = self.stringlist,
                    dtick = 1
                ),
                yaxis = dict(
                    tickmode = "array",
                    tickvals = np.arange(len(self.floorlist)),
                    ticktext = self.floorlist,
                    dtick = 1
                ),
        )
            zmax, zmin = np.nanmax(self.heatmap[i, :, :]), np.nanmin(self.heatmap[i, :, :])
            fig = go.Figure(data = go.Heatmap(z=self.heatmap[i, :, :], text = annotation_text, texttemplate="%{text}", colorscale=custom_colorscale, zmin=zmin, zmax=zmax), layout = layout)
            #print("Average of hits is %f +- %f" %(np.nanmean(self.heatmap[i, :, :]), np.nanstd(self.heatmap[i, :, :])))
            fig.show()
            write_path = str('%s_pmt_%i_%i.pdf' %(title, indices[i], indices[i + 1]))
            pio.write_image(fig, write_path)
        return 0;
    
    def plot_heatmap_matplotlib(self, indices, title):
        for i in range(0, len(indices) - 1):
            annotation_text = np.round(self.heatmap[i, :, :], 4)
            #plt.subplots_adjust(top=0.92, bottom=0.08)
            extent = [-0.5, len(stringlist) - 0.5, -0.5, len(floorlist) - 0.5]
            title_counter = ", PMTs {} - {}".format(indices[i], indices[i + 1])
            fig, ax = plt.subplots()
            im = ax.imshow(self.heatmap[i, :, :], extent=extent)
            ax.set_title(title + title_counter)
            ax.set_xticks(np.arange(len(self.stringlist)))
            ax.set_xticklabels(self.stringlist)
            ax.set_yticks(np.arange(len(self.floorlist)))
            ax.set_yticklabels(self.floorlist)
            ax.set_xlabel('X Label')
            ax.set_ylabel('Y Label')
    
            # Display the z-value on the plot
            for row in range(len(self.floorlist)):
                for col in range(len(self.stringlist)):
                    ax.text(col, row, annotation_text[row, col],
                            ha="center", va="center", color="black")
    
            # Add a colorbar
            cbar = ax.figure.colorbar(im)
            cbar.set_label('Colorbar Label')
    
            write_path = "{}_pmt_{}_{}.pdf".format(title, indices[i], indices[i + 1])
            plt.savefig(write_path)
            
            plt.show()

        return 0;

    

    def compare_upper_lower_pmts_heatmap(self, indices, title):
        if len(indices)-1 != 2:
            Exception("only works if there are only two heat maps to compare!")
        else:
            new_heatmap = self.heatmap[0, :, :]/self.heatmap[1, :, :]
            custom_colorscale = [
                [1, 'rgb(255, 255, 0)'],
                [0, 'rgb(0, 0, 255)']
            ]
            annotation_text = np.round(new_heatmap[:, :], 4)
            layout = go.Layout(
                title = title,
                xaxis = dict(
                    tickmode = "array",
                    tickvals = np.arange(len(self.stringlist)),
                    ticktext = self.stringlist,
                    dtick = 1
                ),
                yaxis = dict(
                    tickmode = "array",
                    tickvals = np.arange(len(self.floorlist)),
                    ticktext = self.floorlist,
                    dtick = 1
                ),
        )
            fig = go.Figure(data = go.Heatmap(z=new_heatmap[:, :], text = annotation_text, texttemplate="%{text}", colorscale=custom_colorscale, zmin=0.9, zmax=1.1), layout = layout)
            write_path = str('comparison_ratio.pdf')
            pio.write_image(fig, write_path)
            fig.show()




def calc_heatmap_ratio(heatmap_real, heatmap_sim):
    return heatmap_sim/heatmap_real
    #return heatmap_real/heatmap_sim




#indices = [0, 12, 30]
indices = [0, 1, 7, 13, 19, 25, 31]
pmt_letters = ["A", "B", "C", "D", "E", "F"]

data_real = map_hit_data(muon_hit_data_real, modid_map, pmt_serial_map, magic_number)
data_sim = map_hit_data(muon_hit_data_sim, modid_map, pmt_serial_map, magic_number)

#data_real.apply_mask_return_floor_str_hit()




real_eff_map, __, __, = data_real.export_heatmap(indices, int_rates_or_eff = 5)
sim_map, __, __ = data_sim.export_heatmap(indices)
sim_eff_map, __, __ = data_sim.export_heatmap(indices, int_rates_or_eff =5)

real_map, floorlist, stringlist = data_real.export_heatmap(indices)







sim_ratio_map = heatmap(calc_heatmap_ratio(real_map, sim_map), floorlist, stringlist)
real_heatmap = heatmap(real_map, floorlist, stringlist)
eff_heatmap = heatmap(sim_eff_map, floorlist, stringlist)

#eff_heatmap.plot_heatmap(indices, "Map of the efficiencies")
#sim_eff_heatmap.plot_heatmap(indices, "Average efficiencies of DOMs, simulated data")

#sim_ratio_map.compare_upper_lower_pmts_heatmap(indices, "Ratio of upper/lower PMTs ratio of simulated/real rates")
#sim_ratio_map.plot_heatmap(indices, "Ratio of simulated vs real rates of PMTs per DOM")
#real_heatmap.plot_heatmap(indices, pmt_letters, "Sum of all hits per DOM for indicated PMT group, real data")
#Create new dataset by laying the efficiency map over the ratio map 

sim_ratio_eff_map = np.zeros([2,sim_eff_map.shape[0], sim_eff_map.shape[1], sim_eff_map.shape[2]])
sim_ratio_eff_map[0, :, :, :] = sim_ratio_map.heatmap[:, :, :]
sim_ratio_eff_map[1, :, :, :] = eff_heatmap.heatmap[:, :, :]


def plot_ratio_eff(sim_ratio_eff_map, indices):
    """This is to plot efficiency vs simulated/real hit ratio for all strings seperately"""
    #for k in range(0, 3):
    for k in range(0, sim_ratio_eff_map.shape[3]):
        plt.figure()
        for l in range(0, len(indices)-1):
    #    for j in range(0, sim_ratio_eff_map.shape[3]):
    #        for i in range(0, sim_ratio_eff_map.shape[2]): #Loops per string over floors
    #            pass
    #            plt.scatter(sim_ratio_eff_map[0,k, :, :], sim_ratio_eff_map[1,k, :, :])
        #plt.scatter(sim_ratio_eff_map[0, k, :, 1], sim_ratio_eff_map[1, k, :, 1], label = "PMTs %i-%i" %(indices[k], indices[k + 1]))
        #plt.scatter(sim_ratio_eff_map[0,:,:,k], sim_ratio_eff_map[1,:,:,k], label="string %i" %(stringlist[k]))
            plt.scatter(sim_ratio_eff_map[0, l, :, k], sim_ratio_eff_map[1, l, :, k], label = "PMTs %i-%i" %(indices[l], indices[l + 1]))
        plt.ylim(0.5, 1.15)
        plt.xlim(1.2, 1.5)
        plt.xlabel("simulated/real hit rates ratio")
        plt.ylabel("efficiency")
        plt.title("DOMs efficiency vs simulated/real hit rate ratio, string %i" %(stringlist[k]))
        plt.legend()
        #plt.savefig("str_plots/t0-ratio_eff-str-%i.pdf" %(stringlist[k]))
        plt.show()

def plot_ratio_eff_one_plot(sim_ratio_eff_map, indices):
    err_colors = ["m", "c"]
    for k in range(0, len(indices)-1):
        plt.scatter(sim_ratio_eff_map[0,k, :, :], sim_ratio_eff_map[1,k, :, :], label = "PMTs %i-%i" %(indices[k], indices[k+1]))
    for k in range(0, len(indices)-1):
        if k == len(indices)-2:
            plt.errorbar(np.nanmean(sim_ratio_eff_map[0,k, :, :]), np.nanmean(sim_ratio_eff_map[1,k, :, :]), 
                        xerr = np.nanstd(sim_ratio_eff_map[0, k, :, :]),
                        yerr = np.nanstd(sim_ratio_eff_map[1, k, :, :]), 
                        fmt = "o", color = err_colors[k],
                        label = "PMTs %i-%i, average" %(indices[k], indices[k+1]))
        else:
            plt.errorbar(np.nanmean(sim_ratio_eff_map[0,k, :, :]), np.nanmean(sim_ratio_eff_map[1,k, :, :]), 
                        xerr = np.nanstd(sim_ratio_eff_map[0, k, :, :]),
                        yerr = np.nanstd(sim_ratio_eff_map[1, k, :, :]), 
                        fmt = "o", color = err_colors[k],
                        label = "PMTs %i-%i, average" %(indices[k], indices[k+1]-1))
    plt.ylim(0.5, 1.15)
    plt.xlim(1.2, 1.5)
    plt.xlabel("simulated/real hit rates ratio")
    plt.ylabel("efficiency")
    plt.title("DOMs efficiency vs simulated/real hit rate ratio, all strings")
    plt.legend()
    plt.savefig("str_plots/t0-ratio_eff-str-fixed.pdf" %(stringlist[k]))
    plt.show()
    for k in range(0, len(indices)-1):
        if k == len(indices)-2:
            print("For PMTs %i-%i the simulated/real hit rates ratio mean is %f +- %f." 
                %(indices[k], indices[k+1],np.nanmean(sim_ratio_eff_map[0,k, :, :]), np.nanstd(sim_ratio_eff_map[0,k, :, :])))
            print("For PMTs %i-%i the efficiency mean is %f +- %f." 
                %(indices[k], indices[k+1],np.nanmean(sim_ratio_eff_map[1,k, :, :]), np.nanstd(sim_ratio_eff_map[1,k, :, :])))
        else:
            print("For PMTs %i-%i the simulated/real hit rates ratio mean is %f +- %f." 
                %(indices[k], indices[k+1]-1,np.nanmean(sim_ratio_eff_map[0,k, :, :]), np.nanstd(sim_ratio_eff_map[0,k, :, :])))
            print("For PMTs %i-%i the efficiency mean is %f +- %f." 
                %(indices[k], indices[k+1]-1,np.nanmean(sim_ratio_eff_map[1,k, :, :]), np.nanstd(sim_ratio_eff_map[1,k, :, :])))

#plot_ratio_eff_one_plot(sim_ratio_eff_map, indices)



