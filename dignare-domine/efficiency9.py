#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 15:54:39 2023

@author: mike

This file should also take into account the two different PMT versions and include a mapping per ring in the DOM. 
Streamlined version of 'efficiency8.py', not backward compatitble.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as go
import plotly.colors as colors

pio.kaleido.scope.mathjax= None
pio.renderers.default='browser'

muon_hit_data_sim = np.load("data/muon_hit_data-sim-reduced_bins-xx1375x.npy")
#muon_hit_data_real = np.load("muon_hit_data-real-reduced_bins-13754.npy")
muon_hit_data_real = np.load("data/muon_hit_data-real-reduced_bins-xx1375x.npy")
modid_map = np.loadtxt("../pmt-info/map.txt")
eff_list = np.loadtxt("data/data-133-144-eff.txt", skiprows = 148, usecols=[1,2,3])
pmt_serial_map = np.loadtxt("../pmt-info/pmt-serials.txt", usecols = 0)
pmt_ring_map = np.loadtxt("../pmt-info/pmt-ring.txt", skiprows = 2, usecols = [0,1,2])
magic_number = 16104 # The major version change happened at serial number 16104 (all PMTs <=16104 are of a certain kind (R12199), all abover are another one (R14374)).

indices = [0, 1, 7, 13, 19, 25, 31]
pmt_letters = ["A", "B", "C", "D", "E", "F"]
floorlist = np.loadtxt("data/floorlist.txt"); stringlist = np.loadtxt("data/stringlist.txt")


class map_hit_data():
    """Connects the muon hit data only identified per identifier to a map containing the key to which floor/string it is located in."""
    """Note initial data muon_hit_data is 2d array w/ column headers module-id / pmt number / amount of hits.
    If needed to append more data just append columns next to the last column as it should survive all operations."""

    """Workflow as follows: data loadin is as 2d array with column as above. Then data gets appended, upon which the array gets transformed to 3d 
    with the pmt number being the 3rd dimension. Then data reorganised into heatmap."""
    def __init__(self, muon_hit_data, pmt_id_map, pmt_serial_map, magic_number, new_versions = None, floor_str_hit = None):
        self.modid_map = modid_map
        self.eff_list = eff_list    
        self.pmt_serial_map = pmt_serial_map; self.magic_number = magic_number
        if floor_str_hit == None:
            self.muon_hit_data = muon_hit_data
            self.append_eff_data()
            self.append_pmt_serials()
            self.mod_id_to_floor_string()
            self.pmt_no_to_ring_letter()
            if new_versions != None: 
                self.apply_pmt_mask(new_versions = new_versions)
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
        for i in range(0, 31):
            self.floor_str_hit[:, i, 3] = pmt_ring_map[i, 1]
        return 0;

    def apply_pmt_mask(self, new_versions=1):
        """Filters for either new pmt version or the old one. Data located in column no 8.
        If new_versions = 1, then only new versions are returned, otherwise only old versions returned."""
        new_floor_str_hit = self.floor_str_hit.copy()
        for i in range(0, 31):
            current_floor_str_hit = self.floor_str_hit[:, i, :]
            mask = self.floor_str_hit[:, i, 7].astype(int)
            if new_versions == 1:
                mask = 1 - mask
            #masked_floor_str_hit = np.ma.masked_array(current_floor_str_hit, np.tile(mask, (8,1)).T)
            #new_floor_str_hit[:, i, :] = masked_floor_str_hit.filled(fill_value= np.nan)
            masked_floor_str_hit = np.ma.masked_array(current_floor_str_hit[:, 4:6], np.tile(mask, (2,1)).T)
            new_floor_str_hit[:, i, 4:6] = masked_floor_str_hit.filled(fill_value= np.nan)
        self.floor_str_hit = new_floor_str_hit
        return 0;

    def sum_over_n_pmts(self, indices):
        """Calculates average over any [n] groups of pmts. Group must include starting PMT-no. of group (inclusive) and stopping number (exclusive)"""
        pmt_group_mean = np.zeros([self.floor_str_hit.shape[0], len(indices) - 1, self.floor_str_hit.shape[2]])
        j = 0; k = 0
        for i in range(1, max(indices)+1):
            if i in indices:
                pmt_group_mean[:, k, 4] = np.nansum(self.floor_str_hit[:, j:i, 4], axis=1)
                pmt_group_mean[:, k, :3] = np.nanmean(self.floor_str_hit[:, j:i, :3], axis = 1)
                pmt_group_mean[:, k, 5:] = np.nanmean(self.floor_str_hit[:, j:i, 5:], axis = 1)
                j = i + 1; k = k + 1
        return pmt_group_mean

    def normalise_over_n_pmts(self, indices):
        """Calculates average over any [n] groups of pmts. Group must include starting PMT-no. of group (inclusive) and stopping number (exclusive)"""
        pmt_group_mean = np.zeros([self.floor_str_hit.shape[0], len(indices) - 1, self.floor_str_hit.shape[2]])
        j = 0; k = 0
        for i in range(1, max(indices)+1):
            if i in indices:
                pmt_group_mean[:, k, :] = np.nanmean(self.floor_str_hit[:, j:i, :], axis=1)
                j = i + 1; k = k + 1
        return pmt_group_mean


    def heatmap_averages_single_loop(self, pmt_group_pairs):
        #Sort pmt_groups on str and floor
        pmt_group_pairs = pmt_group_pairs[pmt_group_pairs[:, 0].argsort()] 
        #First sort according to string; then according to floor 
        k = 0; l = 0;
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
                k = i; l = l + 1
        return pmt_group_pairs

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
            pmt_group_mean_sorted[:, m, :] = self.heatmap_averages_single_loop(pmt_group_pairs) #Sorts the thing on string and floor so that they are in sequence. 
        floorlist = np.unique(pmt_group_mean_sorted[:, 0, 1]).tolist(); stringlist = np.unique(pmt_group_mean_sorted[:, 0, 0]).tolist()
        heatmap = np.zeros([len(indices)-1, len(floorlist), len(np.unique(stringlist))])
        #Now get heatmap for each group
        for i in range(0, len(indices)-1):
            heatmap[i, :, :] = self.heatmap_array_single_group(pmt_group_mean_sorted[:, i, :], i, heatmap[i, :, :], floorlist, stringlist, int_rates_or_eff)
        return heatmap
    

def speedrun_heatmap(muon_hit_data, pmt_id_map, pmt_serial_map, magic_number, new_versions = None, floorlist = floorlist, stringlist = stringlist, int_rates_or_eff = 4):
    hit_runs_arr = map_hit_data(muon_hit_data, pmt_id_map, pmt_serial_map, magic_number, new_versions= new_versions)
    return heatmap(hit_runs_arr.export_heatmap(indices, int_rates_or_eff = int_rates_or_eff))

class heatmap():
    """Class to generate heatmap plots with string and/or floor information. Also useful to perform heatmap operations with."""
    def __init__(self, heatmap, floorlist = floorlist, stringlist = stringlist, x_ax = None):
        self.heatmap = heatmap
        self.floorlist = floorlist.tolist()
        self.stringlist = stringlist.tolist()
        if x_ax != None:
            self.x_ax = x_ax
        

    def append_mean_row_column(self, indices):
        """Appends the mean to all rows and columns in the heatmap. Also appends floorlist and stringlist 
        to reflect this."""
        floorlist = self.floorlist.copy(); stringlist = self.stringlist.copy()
        stringlist.append("mean")
        floorlist.append("mean")
        string_mean = np.nanmean(self.heatmap, axis = 1); floor_mean = np.nanmean(self.heatmap, axis = 2)
        #print(string_mean)
        new_heatmap = np.zeros([self.heatmap.shape[0], self.heatmap.shape[1]+1, self.heatmap.shape[2]+1])
        new_heatmap[:, :-1, :-1] = self.heatmap
        new_heatmap[:, -1, :-1] = string_mean; 
        new_heatmap[:, :-1, -1] = floor_mean
        new_heatmap[:, -1, -1] = np.nan #this corner does not mean anything and should have no data
        return new_heatmap, floorlist, stringlist

    def delete_mean_row_column(self, indices):
        """Deletes the mean of all rows and columns in the heatmap. Reverses above function essentially."""
        old_heatmap = self.heatmap[:, :-1, :-1]
        floorlist = self.floorlist[:-1]; stringlist = self.stringlist[:-1]
        print(old_heatmap.shape)
        return old_heatmap, floorlist, stringlist

    def summarise_per_ring(self, indices, pmt_letters, string = True):
        """Returns mean heatmap of either floor of string along with floor or stringlist."""
        heatmap, __, __ = self.append_mean_row_column(indices)
        if string == True: 
            x_ax = self.stringlist
        else:
            x_ax = self.floorlist
        heatmap_summarised = np.zeros([len(pmt_letters), len(x_ax)])
        
        if string == True:
             heatmap_summarised = heatmap[:, -1, :-1]
        else:
             heatmap_summarised = heatmap[:, :-1, -1]
        return heatmap_summarised, x_ax

    def summarise_per_ring_part(self, indices, pmt_letters, start_index, stop_index, slice_string = True):
        if slice_string == True:
            string_start = start_index; string_stop = stop_index; floor_start, floor_stop = 0, 18
            axis = 2
            #x_ax = self.stringlist[string_start:string_stop]
            x_ax = self.floorlist
        else:
            string_start, string_stop = 0, 20; floor_start, floor_stop = string_start, string_stop
            axis = 1
            #x_ax = self.floorlist[floor_start:floor_stop]
            x_ax = self.stringlist
        mean_heatmap = np.nanmean(self.heatmap[:, floor_start:floor_stop, string_start:string_stop], axis = axis)
        print(mean_heatmap.shape)
        return mean_heatmap, x_ax



    def export_summarised_heatmap(self, indices, pmt_letters, string = True):
        heatmap_summarised, __ =  self.summarise_per_ring(indices, pmt_letters, string)
        return heatmap_summarised
        


    def plot_heatmap(self, indices, pmt_letters, title, save = "Yes", save_map=None, zmax_array = None, zmin_array = None, include_mean = False):
        """Manually set pmt_letters to none if groups are different"""
        colorscale = colors.sequential.Sunset
        colorscale = colorscale[::-1]
        colorscale[0] = '#665679'
        
        if zmax_array is None and zmin_array is None:
            zmin_present, zmax_present = False, False
            zmax_array, zmin_array = np.zeros(len(indices)-1), np.zeros(len(indices)-1)
        else:
            zmin_present, zmax_present = True, True
        if include_mean == True:
            heatmap, floorlist, stringlist = self.append_mean_row_column(indices)
        else:
            heatmap, floorlist, stringlist = self.heatmap, self.floorlist, self.stringlist
        for i in range(0, len(indices)-1):
            if pmt_letters is None:
                title_complete = title
            else:
                title_complete = title + ", PMT group %s" %((pmt_letters[i]))
            heatmap_current = heatmap[i, :, :]
            annotation_text = np.round(heatmap_current, 4)
            layout = go.Layout(
                title = title_complete,
                xaxis = dict(
                    tickmode = "array",
                    tickvals = np.arange(len(stringlist)),
                    ticktext = stringlist,
                    dtick = 1
                ),
                yaxis = dict(
                    tickmode = "array",
                    tickvals = np.arange(len(floorlist)),
                    ticktext = floorlist,
                    dtick = 1
                ),
        )
            if zmax_present == False:
                zmax = np.nanmax(heatmap_current)
                zmax_array[i] = zmax
            else:
                zmax = zmax_array[i]
            if zmin_present == False:
                zmin = np.nanmin(heatmap_current[heatmap_current != 0])
                zmin_array[i] = zmin
            else:
                zmin = zmin_array[i]
            fig = go.Figure(data = go.Heatmap(z=heatmap_current, text = annotation_text, texttemplate="%{text}", colorscale=colorscale, zmin=zmin, zmax=zmax), layout = layout)
            #print("Average of hits is %f +- %f" %(np.nanmean(self.heatmap[i, :, :]), np.nanstd(self.heatmap[i, :, :])))
            #fig.show()
            if save == "Yes":
                if save_map == None:
                    write_path = str('%s_pmt-ring-%s.pdf' %(title.replace(" ", "-"), pmt_letters[i]))
                else:
                    write_path = save_map + str('/%s_pmt-ring-%s.pdf' %(title.replace(" ", "-"), pmt_letters[i]))
                pio.write_image(fig, write_path)
            else:
                pass
        #print(zmax_array, zmin_array)
        #if include_mean == True: 
        #    self.heatmap, self.floorlist, self.stringlist = self.delete_mean_row_column(indices)
        return zmax_array, zmin_array

    def plot_heatmap_summarised_ring(self, indices, pmt_letters, title, save = "Yes", save_map=None, string = True, edit_heatmap = False, x_ax = None):
        """Summarises the mean per string or per floor into a single heatmap. If string = false then floor plot.
        If x_ax != None then already assumsed self.summarise_per_ring executed"""
        colorscale = colors.sequential.Sunset
        colorscale = colorscale[::-1]
        print(colorscale[0])
        colorscale[0] = '#665679'
        if x_ax == None:
            heatmap_summarised, x_ax = self.summarise_per_ring(indices, pmt_letters, string)
        else:
            heatmap_summarised = self.heatmap; x_ax = self.x_ax
        annotation_text = np.round(heatmap_summarised, 4)
        layout = go.Layout(
            title = title,
            xaxis = dict(
                tickmode = "array",
                tickvals = np.arange(len(x_ax)),
                ticktext = x_ax,
                dtick = 1
            ),
            yaxis = dict(
                tickmode = "array",
                tickvals = np.arange(len(pmt_letters)),
                ticktext = pmt_letters,
                dtick = 1
            ),
        )
        zmax = np.nanmax(heatmap_summarised); zmin = np.nanmin(heatmap_summarised)
        fig = go.Figure(data = go.Heatmap(z=heatmap_summarised, text = annotation_text, 
            texttemplate="%{text}", colorscale=colorscale, zmin=zmin, zmax=zmax), layout = layout)
        if save == "Yes":
            if save_map == None:
                write_path = str('%s.pdf' %(title.replace(" ", "-")))
            else:
                write_path = save_map + str('/%s.pdf' %(title.replace(" ", "-")))
                pio.write_image(fig, write_path)
        #self.heatmap, __, __ = self.delete_mean_row_column(indices)
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

    def get_avg_std(self, indices):
        mean = np.nanmean(self.heatmap, axis = (1, 2))
        std = np.nanstd(self.heatmap, axis = (1, 2))
        return mean, std

    

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
    return heatmap_real/heatmap_sim
    #return heatmap_real/heatmap_sim

def summarised_heatmap_ratio(heatmap_num, heatmap_denom, title, indices, pmt_letters, plot_string = True, slice_string = True, save = "Yes", save_map=None,
    start_index = None, stop_index = None):
    """Plots heatmap of a ratio of the numerator map over the denominator map. 
    Try new/better over old/worse maps."""
    if start_index != None and stop_index != None: 
        exportmap_num, x_ax = heatmap_num.summarise_per_ring_part(indices, pmt_letters, start_index, stop_index, slice_string = slice_string)
        exportmap_denom, __ = heatmap_denom.summarise_per_ring_part(indices, pmt_letters, start_index, stop_index, slice_string = slice_string)
    else:
        exportmap_num, x_ax = heatmap_num.summarise_per_ring(indices, pmt_letters, string=plot_string)
        exportmap_denom, __ = heatmap_denom.summarise_per_ring(indices, pmt_letters, string=plot_string)
    heatmap_ratio = heatmap(exportmap_num / exportmap_denom, x_ax = x_ax)
    heatmap_ratio.plot_heatmap_summarised_ring(indices, pmt_letters, title, save_map = save_map, string = plot_string, x_ax = x_ax)

#plot_ratio_eff_one_plot(sim_ratio_eff_map, indices)



