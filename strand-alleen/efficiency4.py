#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 15:54:39 2023

@author: mike
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
import ROOT 
import plotly.io as pio
import plotly.express as px
from collections import Counter

pio.renderers.default='browser'

muon_hit_data = np.load("muon_hit_data2.npy")
modid_map = np.loadtxt("map.txt")

class map_hit_data():
    """Connects the muon hit data only identified per identifier to a map containing the key to which floor/string it is located in."""
    def __init__(self, muon_hit_data, pmt_id_map):
        self.muon_hit_data = muon_hit_data
        self.modid_map = modid_map

    def mod_id_to_floor_string(self):
        """Note output is of the form: amount of mod-ids / pmts numbers / [str no; floor no; mod-id; pmt-no; no. hits]"""
        """Also note that mod-id refers to which DOM one is looking at"""
        floor_str_hit = np.zeros([self.muon_hit_data.shape[0], self.muon_hit_data.shape[1], self.muon_hit_data.shape[2]+2])
        mapping = dict(zip(self.modid_map[:,0], range(len(modid_map))))
        for i in range(0, self.muon_hit_data.shape[1]):
            floor_str_hit[:, i, :] = np.hstack((np.array([self.modid_map[mapping[key],1:] for key in self.muon_hit_data[:,i,0]]), self.muon_hit_data[:, i, :]))
        print(np.unique(floor_str_hit[:, 0, 0], axis = 0))
        self.floor_str_hit = floor_str_hit
        print(floor_str_hit.shape)

    def normalise_over_n_pmts(self, indices):
        """Calculates average over any [n] groups of pmts. Group must include starting PMT-no. of group (inclusive) and stopping number (exclusive)"""
        pmt_group_mean = np.zeros([self.floor_str_hit.shape[0], len(indices) - 1, self.floor_str_hit.shape[2]])
        j = 0; k = 0
        for i in range(1, max(indices)+1):
            if i in indices:
                pmt_group_mean[:, k, :] = np.mean(self.floor_str_hit[:, j:i, :], axis=1)
                j = i + 1; k = k + 1
        return pmt_group_mean


    def heatmap_modid(self):
        #Reshape thing into heatmap w/ x-axis mod-id; y-axis pmt; z-axis hit rate 
        heatmap = np.zeros([self.floor_str_hit.shape[1], self.floor_str_hit.shape[0]])
        for i in range(0, heatmap.shape[0]):
            heatmap[i, :] = self.floor_str_hit[:, i, 4]

        fig = px.imshow(heatmap, labels=dict(x="module-id",y="pmt number", color="number of hits"), 
                        x= self.floor_str_hit[:, 0, 2],y = self.floor_str_hit[0, :, 3],
                        title = "Rates of PMTs per DOM", 
                        text_auto=True, aspect="auto", width=2560, height=1440)
        fig.write_image("mod-id-pmt-test-plot.pdf")
        #fig.show()
        return 0;

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
                new_index_part = pmt_group_pairs[k:i, 1].argsort()
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
        return pmt_group_pairs, str_floor_length

    def heatmap_averages(self, indices):
        pmt_group_mean = self.normalise_over_n_pmts(indices)
        pmt_group_mean_sorted = pmt_group_mean
        print(pmt_group_mean_sorted[:, 0, :])
        #pmt_group_pairs = pmt_group_mean[:, 0, :]
        for m in range(0, len(indices)-1):
            pmt_group_pairs = pmt_group_mean[:, m, :]
            pmt_group_mean_sorted[:, m, :], str_floor_length = self.heatmap_averages_single_loop(pmt_group_pairs)
        floorlist = np.unique(pmt_group_mean_sorted[:, 0, 1]); stringlist = np.unique(pmt_group_mean_sorted[:, 0, 0])
        heatmap = np.zeros([len(floorlist), len(np.unique(stringlist))])
        print(heatmap.shape)
        #TODO fill in heatmap and account for nans.
        k = 0; m = 0; x = 0
        for i in range(0, pmt_group_mean_sorted.shape[0]): #first fill in the x/string direction
            if pmt_group_mean_sorted[i, 0, 0] != pmt_group_mean_sorted[i-1, 0, 0]:
                #Now look if there are any gaps in the data 
                for j in range(0, len(floorlist)):
                    l = 0 
                    if pmt_group_mean_sorted[j+k, 0, 1] == floorlist[j]:
                        heatmap[j, m] = pmt_group_mean_sorted[j+k, 0, 4]
                    else:
                        heatmap[j, m] = np.nan
                        print(stringlist[m], floorlist[j])
                        print(pmt_group_mean_sorted[j+k, 0, :2])
                        #j = j + 1
                        x = x + 1
                        j = j + np.abs(pmt_group_mean_sorted[j+k, 0, 1] - floorlist[j]) #skip ahead to current pmt number
                m = m + 1; k = i
        print(x)
        #print(heatmap)
        






def check_for_unique_pairs(a):
        """Takes an [n, 2] sized array and compares if any combination of the pairs across the 2-dimension 
        occurs again."""
        #ctr = Counter(frozenset(x) for x in a)
        #b = [ctr[frozenset(x)] == 1 for x in a]
        #print(b)
        b = np.unique(a, axis = 0)
        print(b.shape)




indices = [0, 11, 18, 30]
test = map_hit_data(muon_hit_data, modid_map)
test.mod_id_to_floor_string()
test.heatmap_averages(indices)
#test.heatmap_modid()

pmt_group_mean = test.normalise_over_n_pmts(indices)
pmt_group_pairs = pmt_group_mean[:, 0, :2]