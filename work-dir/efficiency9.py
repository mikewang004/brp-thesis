#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 15:54:39 2023

@author: mike

This file should also take into account the two different PMT versions and include a mapping per ring in the DOM. 
Streamlined version of 'efficiency8.py', not backward compatitble.
"""

from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as go
import plotly.colors as colors
#from extract_root.py import *
pio.kaleido.scope.mathjax= None
pio.renderers.default='browser'

muon_hit_data_sim = np.load("data/muon_hit_data-sim-reduced_bins-xx1375x.npy")
#muon_hit_data_real = np.load("muon_hit_data-real-reduced_bins-13754.npy")
muon_hit_data_real = np.load("data/muon_hit_data-real-reduced_bins-xx1375x.npy")
modid_map = np.loadtxt("../pmt-info/map.txt")
#eff_list = np.loadtxt("data/data-133-144-eff.txt", skiprows = 148, usecols=[1,2,3])
eff_list = np.loadtxt("data/runs-14413-14440-eff.txt")
pmt_serial_map = np.loadtxt("../pmt-info/pmt-serials.txt", usecols = 0)
pmt_ring_map = np.loadtxt("../pmt-info/pmt-ring.txt", skiprows = 2, usecols = [0,1,2])
magic_number = 16104 # The major version change happened at serial number 16104 (all PMTs <=16104 are of a certain kind (R12199), all abover are another one (R14374)).

indices = [0, 1, 7, 13, 19, 25, 31]
pmt_letters = ["A", "B", "C", "D", "E", "F"]
floorlist = np.loadtxt("data/floorlist.txt"); stringlist = np.loadtxt("data/stringlist.txt")

min_std_fac = 2 
max_std_fac = 3


class meanhitrate():
    """Everything that has to do with the efficiency/rate arrays"""
    def __init__(self, mean_hit_rate):
        """Remember mean hit rate array has following structure:
            [data, dom-id / pmt-no / eff / dom / string / pmt-serial / new pmt [yes/no] / rate [khz] ]"""
        self.meanhitrate = mean_hit_rate
        self.block_per_pmt()
        self.filter_data()

    def block_per_pmt(self):
        self.meanhitrate = self.meanhitrate.reshape(int(self.meanhitrate.shape[0]/31), 31, 8) #block per PMT
        return 0; 


        
    def avg_top_bottom_pmts(self, mid_pmt = 12):
        """Averages over the top (0-11) and bottom (12-31) pmts."""
        test = self.meanhitrate.shape
        #self.meanhitrate = self.meanhitrate[~np.isnan(self.meanhitrate).any(axis=0)]

        self.top_avg = np.mean(self.meanhitrate[:, 0:mid_pmt, :], axis=1)
        self.bottom_avg = np.mean(self.meanhitrate[:, mid_pmt:31, :], axis=1)
        #self.pmtavg = np.zeros([self.meanhitrate.shape[0]*2, self.meanhitrate.shape[2]])
        #self.pmtavg[::2, :] = self.top_avg; self.pmtavg[1::2, :] = self.bottom_avg;
        #self.top_length = self.top_avg.shape[0]
        #self.top_mask = self.top_avg[:-1] == self.top_avg[1:]
        return 0;

    def filter_data(self):
        """Removes all outliers greater than say 3 sigma from the average. Default usage is on efficiency and rate"""
        j = 0
        for i in range(1, self.meanhitrate.shape[0]):
        #for i in range(1, 25):
            if self.meanhitrate[i-1, 0, 3] != self.meanhitrate[i,0,  3] or i == 377: #checks for start new string
                avg_eff, std_eff = np.nanmean(self.meanhitrate[j:i, :, 2]), np.nanstd(self.meanhitrate[j:i:, :, 2])
                avg_rate, std_rate = np.nanmean(self.meanhitrate[j:i, :, 7]), np.nanstd(self.meanhitrate[j:i, :, 7])
                #Now create new filter which nans 
                index_eff = np.abs(avg_eff - self.meanhitrate[j:i, :, 2]) < 1 * std_eff
                index_rate = np.abs(avg_rate - self.meanhitrate[j:i, :, 7]) < 1 * std_rate
                index_all = index_eff * index_rate
                self.meanhitrate[j:i, :, 2] = np.where(index_all, self.meanhitrate[j:i, :, 2], np.nan)
                self.meanhitrate[j:i, :, 7] = np.where(index_all, self.meanhitrate[j:i, :, 7], np.nan)
                j = i
        return 0;

class map_hit_data():
    """Connects the muon hit data only identified per identifier to a map containing the key to which floor/string it is located in."""
    """Note initial data muon_hit_data is 2d array w/ column headers module-id / pmt number / amount of hits.
    If needed to append more data just append columns next to the last column as it should survive all operations."""

    """Workflow as follows: data loadin is as 2d array with column as above. Then data gets appended, upon which the array gets transformed to 3d 
    with the pmt number being the 3rd dimension. Then data reorganised into heatmap."""
    def __init__(self, muon_hit_data, pmt_id_map, pmt_serial_map, magic_number, new_versions = None, floor_str_hit = None, apply_shadow_mask = None):
        self.modid_map = modid_map
        self.eff_list = eff_list    
        self.pmt_serial_map = pmt_serial_map; self.magic_number = magic_number
        if floor_str_hit == None:
            self.muon_hit_data = muon_hit_data
            self.append_eff_data()
            self.append_pmt_serials()
            self.append_hit_rates()
            self.mod_id_to_floor_string()
            if apply_shadow_mask != None:
                self.apply_shadow_mask(filter = apply_shadow_mask) #can be True or False
            if new_versions != None: 
                print("test")
                self.apply_pmt_mask(new_versions = new_versions)
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
        new_muon_hit_data = np.zeros([self.muon_hit_data.shape[0], self.muon_hit_data.shape[1], self.muon_hit_data.shape[2]+2])
        new_muon_hit_data[:,:,:4] = self.muon_hit_data
        for i in range(0, len(self.muon_hit_data)):
            for j in range(0, len(self.pmt_serial_map)):
                if self.muon_hit_data[i, 0, 0] == self.pmt_serial_map[j]:
                    for k in range(0, 31):
                        new_muon_hit_data[i, k, 4] = self.pmt_serial_map[j + k+1]
                        if self.pmt_serial_map[j + k + 1] > magic_number:
                            new_muon_hit_data[i, k, 5] = 1
                        else:
                            new_muon_hit_data[i, k, 5] = 0
                    break
        self.muon_hit_data = new_muon_hit_data
        return 0;

    def append_hit_rates(self, single_rates_path = "bins/mean_hit_rate.npy"):
        single_rates = meanhitrate(np.load(single_rates_path))
        new_muon_hit_data = np.zeros([self.muon_hit_data.shape[0], self.muon_hit_data.shape[1], self.muon_hit_data.shape[2]+1])
        new_muon_hit_data[:,:,:-1] = self.muon_hit_data
        for i in range(0, len(self.muon_hit_data)):
            for j in range(0, len(self.muon_hit_data[:, 0, 0])):
                if self.muon_hit_data[i, 0, 0] == single_rates.meanhitrate[j, 0, 0]:
                    for k in range(0, 31):
                        new_muon_hit_data[i, k, 6] = single_rates.meanhitrate[i, k, 7]
                    break
        self.muon_hit_data = new_muon_hit_data
        return 0; 


    def mod_id_to_floor_string(self):
        """Transforms data from a 2D to a 3D array including the different pmts. 
        Note output is of the form: amount of mod-ids / pmts numbers / [str no; floor no; mod-id; pmt-no; no. hits; eff; pmt_serial; pmt_old/new; rates]
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
        # Then sort the thing so that ring A is appearing first, then B, etc. 
        floor_str_hit = np.zeros(self.floor_str_hit.shape)
        floor_str_hit[:, 0, :] = self.floor_str_hit[:, 22, :] #specific case that ring A is DAQ 22
        i = 1
        for l in range(2, 7):
            for k in range(0, 31):
                if self.floor_str_hit[0, k, 3] == l:
                    floor_str_hit[:, i, :] = self.floor_str_hit[:, k, :]
                    i = i + 1
        self.floor_str_hit = floor_str_hit
        return 0;

    def apply_shadow_mask(self, filter = True):
        """Filters out specific PMTs E2 E6 C2 C5 (DAQs 6 11 20 21).
           If filter = False, return only PMTs that are NOT equator-tape-shadowed. Else only return those that are."""
        shadowed_pmts = [6, 11, 20, 21]
        shadowed_mask = np.zeros(self.floor_str_hit.shape[1])
        shadowed_mask[shadowed_pmts] = 1
        new_floor_str_hit = self.floor_str_hit
        if filter != False:
            shadowed_mask = 1 - shadowed_mask
        masked_floor_str_hit = np.ma.masked_array(new_floor_str_hit[:, :, 4], np.tile(shadowed_mask, (1,new_floor_str_hit.shape[0])))
        masked_floor_str_hit_2 = np.ma.masked_array(new_floor_str_hit[:, :, 5], np.tile(shadowed_mask, (1, new_floor_str_hit.shape[0])))
        new_floor_str_hit[:, :, 4] = masked_floor_str_hit.filled(fill_value= np.nan)
        new_floor_str_hit[:, :, 5] = masked_floor_str_hit_2.filled(fill_value = np.nan)
        self.floor_str_hit = new_floor_str_hit
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
            masked_floor_rates_only = np.ma.masked_array(current_floor_str_hit[:, 8], np.tile(mask, (1,1)).T)
            new_floor_str_hit[:, i, 4:6] = masked_floor_str_hit.filled(fill_value= np.nan)
            new_floor_str_hit[:, i, 8] = masked_floor_rates_only.filled(fill_value= np.nan)
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
        """int_rates_or_eff = 4 for rates; =5 for efficiencies; = 9 for single rates"""
        k = 0; m = 0; x = 0; j = 0; l = 0; #n = 4 for rates; n = 5 for efficiencies or 
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

    def export_all_31_layers(self, int_rates_or_eff = 4):
        pmt_group_mean_sorted = self.floor_str_hit.copy()
        for m in range(0, 31):
            pmt_group_pairs = self.floor_str_hit[:, m, :]
            pmt_group_mean_sorted[:, m, :] = self.heatmap_averages_single_loop(pmt_group_pairs) #Sorts the thing on string and floor so that they are in sequence. 
        floorlist = np.unique(pmt_group_mean_sorted[:, 0, 1]).tolist(); stringlist = np.unique(pmt_group_mean_sorted[:, 0, 0]).tolist()
        heatmap = np.zeros([31, len(floorlist), len(np.unique(stringlist))])
        for i in range(0, 31):
            heatmap[i, :, :] = self.heatmap_array_single_group(pmt_group_mean_sorted[:, i, :], i, heatmap[i, :, :], floorlist, stringlist, int_rates_or_eff)
        return heatmap
    

def speedrun_heatmap(muon_hit_data, pmt_id_map, pmt_serial_map, magic_number, new_versions = None, floorlist = floorlist, stringlist = stringlist, int_rates_or_eff = 4, apply_shadow_mask = False):
    hit_runs_arr = map_hit_data(muon_hit_data, pmt_id_map, pmt_serial_map, magic_number, new_versions= new_versions, apply_shadow_mask= apply_shadow_mask)
    return heatmap(hit_runs_arr.export_heatmap(indices, int_rates_or_eff = int_rates_or_eff))

def speedrun_heatmap_31(muon_hit_data, pmt_id_map, pmt_serial_map, magic_number, new_versions = None, floorlist = floorlist, stringlist = stringlist, int_rates_or_eff = 4, apply_shadow_mask = False):
    hit_runs_arr = map_hit_data(muon_hit_data, pmt_id_map, pmt_serial_map, magic_number, new_versions= new_versions, apply_shadow_mask= apply_shadow_mask)
    return heatmap(hit_runs_arr.export_all_31_layers(int_rates_or_eff = int_rates_or_eff))

class heatmap():
    """Class to generate heatmap plots with string and/or floor information. Also useful to perform heatmap operations with."""
    def __init__(self, heatmap, floorlist = floorlist, stringlist = stringlist, x_ax = None):
        self.heatmap = heatmap
        self.floorlist = floorlist.tolist()
        self.stringlist = stringlist.tolist()
        if x_ax != None:
            self.x_ax = x_ax
        

    def append_mean_row_column(self, indices, append_list_1 = None, append_list_2 = None, heatmap = None, include_mean_of_mean = False):
        """Appends the mean to all rows and columns in the heatmap. Also appends floorlist and stringlist 
        to reflect this."""
        if type(append_list_1) == type(None) and type(append_list_2) == type(None):
            floorlist = self.floorlist.copy(); stringlist = self.stringlist.copy()
            stringlist.append("mean")
            floorlist.append("mean")
        else:
            if not isinstance(append_list_1, list):
                append_list_a = append_list_1.tolist()
            else:
                append_list_a = append_list_1.copy()
            if not isinstance(append_list_2, list):
                append_list_b = append_list_2.tolist()
            else:
                append_list_b = append_list_2.copy()
            append_list_a.append("mean")
            append_list_b.append("mean")
        if type(heatmap) == type(None):
            heatmap = self.heatmap
        else:
            heatmap = heatmap
        heatmap_zero_to_nan = heatmap.copy()
        heatmap_zero_to_nan[heatmap_zero_to_nan == 0] = np.nan
        if heatmap.ndim == 2:
            string_mean = np.nanmean(heatmap_zero_to_nan, axis = 0); floor_mean = np.nanmean(heatmap_zero_to_nan, axis = 1)
            new_heatmap = np.zeros([heatmap.shape[0]+1, heatmap.shape[1]+1])
            new_heatmap[:-1, :-1] = heatmap
            new_heatmap[-1, :-1] = string_mean; 
            new_heatmap[:-1, -1] = floor_mean
            if include_mean_of_mean == True:
                new_heatmap[-1, -1] = np.nanmean(heatmap)
            else:
                new_heatmap[-1, -1] = np.nan #this corner does not mean anything and should have no data
        else:
            string_mean = np.nanmean(heatmap_zero_to_nan, axis = 1); floor_mean = np.nanmean(heatmap_zero_to_nan, axis = 2)
        #print(string_mean)
            new_heatmap = np.zeros([heatmap.shape[0], heatmap.shape[1]+1, heatmap.shape[2]+1])
            new_heatmap[:, :-1, :-1] = heatmap
            new_heatmap[:, -1, :-1] = string_mean; 
            new_heatmap[:, :-1, -1] = floor_mean
            if include_mean_of_mean == True:
                for i in range(0, new_heatmap.shape[0]):
                    new_heatmap[i, -1, -1] = np.nanmean(new_heatmap[i, :-1, :-1])
            else:
                new_heatmap[:, -1, -1] = np.nanmean #this corner does not mean anything and should have no data
        if type(append_list_1) == type(None) and type(append_list_2) == type(None):
            return new_heatmap, floorlist, stringlist
        else:
            return new_heatmap, append_list_a, append_list_b

    def delete_mean_row_column(self, indices):
        """Deletes the mean of all rows and columns in the heatmap. Reverses above function essentially."""
        old_heatmap = self.heatmap[:, :-1, :-1]
        floorlist = self.floorlist[:-1]; stringlist = self.stringlist[:-1]
        print(old_heatmap.shape)
        return old_heatmap, floorlist, stringlist

    def summarise_per_ring(self, indices, pmt_letters, string = True):
        """Returns mean heatmap of either floor of string along with floor or stringlist."""
        heatmap, __, __ = self.append_mean_row_column(indices, include_mean_of_mean = True)
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
        
    def export_summarised_heatmap_and_labels(self, indices, pmt_letters, string = True):
        heatmap_summarised, x_ax =  self.summarise_per_ring(indices, pmt_letters, string)
        return heatmap_summarised, x_ax, pmt_letters

    def plot_heatmap(self, indices, pmt_letters, title, save = "Yes", save_map=None, zmax_array = None, zmin_array = None, include_mean = False):
        """Manually set pmt_letters to none if groups are different"""
        global min_std_fac; global max_std_fac
        colorscale = colors.sequential.Sunset
        colorscale = colorscale[::-1]
        colorscale[0] = '#665679'
        
        if zmax_array is None and zmin_array is None:
            zmin_present, zmax_present = False, False
            zmax_array, zmin_array = np.zeros(len(indices)-1), np.zeros(len(indices)-1)
        else:
            zmin_present, zmax_present = True, True
        if include_mean == True:
            heatmap, floorlist, stringlist = self.append_mean_row_column(indices, include_mean_of_mean = True)
        else:
            heatmap, floorlist, stringlist = self.heatmap, self.floorlist, self.stringlist
        zmax2 = np.nanmean(heatmap) + min_std_fac* np.nanstd(heatmap)
        zmin2 = np.nanmean(heatmap) - max_std_fac * np.nanstd(heatmap)
        print(zmax2, zmin2)
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
                zmax = zmax_array
            if zmin_present == False:
                zmin = np.nanmin(heatmap_current[heatmap_current != 0])
                zmin_array[i] = zmin
            else:
                zmin = zmin_array
            if zmax_present == False:
                fig = go.Figure(data = go.Heatmap(z=heatmap_current, text = annotation_text, texttemplate="%{text}", colorscale=colorscale, zmin=zmin2, zmax=zmax2, zauto=False), layout = layout)
            else:
                fig = go.Figure(data = go.Heatmap(z=heatmap_current, text = annotation_text, texttemplate="%{text}", colorscale=colorscale, zmin=zmin_array, zmax=zmax_array, zauto=False), layout = layout)
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
        #return zmax_array, zmin_array
        return self.heatmap, self.floorlist, self.stringlist



    def plot_heatmap_summarised_ring(self, indices, pmt_letters, title, save = "Yes", 
        save_map=None, string = True, x_ax = None, include_mean = False, include_mean_of_mean = False, zmin = None, zmax = None):
        """Summarises the mean per string or per floor into a single heatmap. If string = false then floor plot.
        If x_ax != None then already assumsed self.summarise_per_ring executed"""
        colorscale = colors.sequential.Sunset
        colorscale = colorscale[::-1]
        print(colorscale[0])
        colorscale[0] = '#665679'
        if x_ax == None:
            heatmap_summarised, x_ax = self.summarise_per_ring(indices, pmt_letters, string=string)
        else:
            heatmap_summarised = self.heatmap; x_ax = self.x_ax
        if include_mean == True:
            heatmap_summarised, x_ax, pmt_letters = self.append_mean_row_column(indices, heatmap = heatmap_summarised,
            append_list_1 = x_ax, append_list_2 = pmt_letters, include_mean_of_mean = include_mean_of_mean)
        print(x_ax, pmt_letters)
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
        #zmax = np.nanmax(heatmap_summarised); zmin = np.nanmin(heatmap_summarised)
        if type(zmin) == type(None):
            zmin = np.nanmean(heatmap_summarised) - 2 * np.nanstd(heatmap_summarised)
        if type(zmax) == type(None):
            zmax = np.nanmean(heatmap_summarised) + 2 * np.nanstd(heatmap_summarised)
        fig = go.Figure(data = go.Heatmap(z=heatmap_summarised, text = annotation_text, 
            texttemplate="%{text}", colorscale=colorscale, zmin=zmin, zmax=zmax), layout = layout)
        if save == "Yes":
            if save_map == None:
                write_path = str('%s.pdf' %(title.replace(" ", "-")))
            else:
                write_path = save_map + str('/%s.pdf' %(title.replace(" ", "-")))
                pio.write_image(fig, write_path)
        #self.heatmap, __, __ = self.delete_mean_row_column(indices)
        return heatmap_summarised
    
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


def plot_heatmap_ultra_basic( heatmap_arr, title, x_ax, y_ax, save="Yes", save_map = None, include_mean = False):
    colorscale = colors.sequential.Sunset
    colorscale = colorscale[::-1]
    colorscale[0] = '#665679'
    zmin = np.nanmean(heatmap_arr) - 2 * np.nanstd(heatmap_arr)
    zmax = np.nanmean(heatmap_arr) + 3 * np.nanstd(heatmap_arr)
    if include_mean == True:
        heatmap_cl = heatmap(heatmap_arr)
        heatmap_arr, x_ax, y_ax = heatmap_cl.append_mean_row_column(
            indices, append_list_1 = x_ax, append_list_2 = y_ax, heatmap = heatmap_arr, include_mean_of_mean = True)
    annotation_text = np.round(heatmap_arr, 4)
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
            tickvals = np.arange(len(y_ax)),
            ticktext = y_ax,
            dtick = 1
        ),
)
    fig = go.Figure(data = go.Heatmap(z=heatmap_arr, text = annotation_text, texttemplate="%{text}", colorscale=colorscale,
    zmin = zmin, zmax = zmax), layout = layout)
    #fig.show()
    if save == "Yes":
        if save_map == None:
            write_path = str('%s.pdf' %(title.replace(" ", "-")))
        else:
            write_path = save_map + str('/%s.pdf' %(title.replace(" ", "-")))
        pio.write_image(fig, write_path)
    return 0;

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
    heatmap_ratio_appended, x_ax_new, pmt_letters_appended = heatmap_ratio.append_mean_row_column(indices, append_list_1=x_ax, append_list_2= pmt_letters)
    #plot_heatmap_ultra_basic(heatmap_ratio.heatmap, title, x_ax, pmt_letters, save_map = save_map)
    plot_heatmap_ultra_basic(heatmap_ratio_appended, title, x_ax_new, pmt_letters_appended, save_map = save_map)
    
def get_zminmax(heatmap):
    zmin = np.nanmean(heatmap) - 1 * np.nanstd(heatmap)
    zmax = np.nanmean(heatmap) + 1 * np.nanstd(heatmap)
    return zmin, zmax


class dist_plots():
    """Everything to do with 1d-distributions. Takes an unlabled heatmap as input."""
    def __init__(self, heatmap_array):
        self.heatmap = heatmap_array
        #

    def generate_counts_bins(self, num_bins = 50, range = None, heatmap = None):
        if type(heatmap) != type(None):
            self.heatmap = heatmap
        heatmap = heatmap[~np.isnan(heatmap)]
        #heatmap = heatmap[heatmap > 0.4]
        #heatmap = heatmap[heatmap < 3]
        if range == None:
            counts, bins = np.histogram(heatmap, bins=num_bins)
            mu, sigma = stats.norm.fit(heatmap)
        else:
            counts, bins = np.histogram(heatmap, bins=num_bins, range = range)
            heatmap = heatmap[heatmap > range[0]]
            heatmap = heatmap[heatmap < range[1]]
            mu, sigma = stats.norm.fit(heatmap)
        self.counts = counts; self.bins = bins
        return counts, bins, mu, sigma

    def plot_dist(self, xlabel=None, title=None, num_bins = 50, save_map = None):
        """Plots heatmap into a 1d distributions."""
        counts, bins, mu, sigma = self.generate_counts_bins(num_bins = num_bins)
        plt.figure()
        plt.stairs(counts, bins)
        plt.xlabel(xlabel)
        plt.title(title)
        plt.ylabel("count")
        if save_map != None:
            write_path = save_map + str('/%s.pdf' %(title.replace(" ", "-")))
            plt.savefig(write_path)
        return 0;

    def plot_dist_barebones(self, num_bins = 50, heatmap = None, label = None, range = None):
        counts, bins, mu, sigma = self.generate_counts_bins(num_bins = num_bins, range = range, heatmap = heatmap)
        if type(label) == type(None):
            plt.stairs(counts, bins)
        else:
            label = label + ', $\mu=%.3f, \sigma=%.3f$' %(mu, sigma)
            plt.stairs(counts, bins, label = label)
        return 0;


    def plot_dist_forloop(self, pmt_letters, num_bins = 50):
        "For multiple rings in a row."
        for i in range(0, self.heatmap.shape[0]):
            heatmap = self.heatmap[i, :, :]
            heatmap = heatmap[~np.isnan(heatmap)]
            counts, bins = np.histogram(heatmap)
            plt.stairs(counts, bins, label = "PMT ring %s" %(pmt_letters[i]))
        return 0;

    def plot_dist_save(self, title, save_map):
        plt.title(title)
        write_path = save_map + str('/%s.pdf' %(title.replace(" ", "-")))
        plt.savefig(write_path)
        return 0;

#plot_ratio_eff_one_plot(sim_ratio_eff_map, indices)



def remove_nan(x):
    y = x[~np.isnan(x)]
    return y


