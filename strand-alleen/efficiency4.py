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

pio.renderers.default='browser'

muon_hit_data = np.load("muon_hit_data.npy")
modid_map = np.loadtxt("map.txt")

class map_hit_data():
    """Connects the muon hit data only identified per identifier to a map containing the key to which floor/string it is located in."""
    def __init__(self, muon_hit_data, pmt_id_map):
        self.muon_hit_data = muon_hit_data
        self.modid_map = modid_map

    def mod_id_to_floor_string(self):
        """Note output is of the form: no mod-ids / pmts numbers / [str no; floor no; mod-id; pmt-no; no. hits]"""
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
                print(i, j)
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

    def heatmap_averages(self, indices):
        pmt_group_mean = self.normalise_over_n_pmts(indices)
        str_list, floor_list = np.unique(pmt_group_mean[:, 0, 0]), np.unique(pmt_group_mean[:, 0, 1])
        heatmap = np.zeros([len(floor_list), len(str_list)])
        #Sort pmt_groups on str and floor
        print(pmt_group_mean[0, ...])


    def check_for_unique_pairs(a):
        a = np.array(a)
        a.sort(axis=1)
        b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
        _, inv, ct = np.unique(b, return_inverse=True, return_counts=True)
        print(ct[inv] == 1)




test = map_hit_data(muon_hit_data, modid_map)
test.mod_id_to_floor_string()
test.heatmap_averages([0, 11, 18, 30])
#test.heatmap_modid()