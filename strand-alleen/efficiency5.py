#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 15:54:39 2023

@author: mike
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as go

pio.renderers.default='browser'

muon_hit_data_sim = np.load("muon_hit_data-sim.npy")
muon_hit_data_real = np.load("muon_hit_data-real.npy")
modid_map = np.loadtxt("map.txt")

class map_hit_data():
    """Connects the muon hit data only identified per identifier to a map containing the key to which floor/string it is located in."""
    def __init__(self, muon_hit_data, pmt_id_map):
        self.muon_hit_data = muon_hit_data
        self.modid_map = modid_map
        self.mod_id_to_floor_string()

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

    def heatmap_array_single_group(self, pmt_group_mean_sorted, pmt_group_no, heatmap, floorlist, stringlist):
        k = 0; m = 0; x = 0; j = 0; l = 0
        for i in range(1, pmt_group_mean_sorted.shape[0]): #first fill in the x/string direction
        #for i in range(1, 70):
            #New approach: floorlist index is not the same as pmt group mean sorted index. 
            if pmt_group_mean_sorted[i, 0] != pmt_group_mean_sorted[i-1, 0] or i == pmt_group_mean_sorted.shape[0]-1:
                while j < len(floorlist): #checks if all floor/string combinations are extant
                    #print(pmt_group_mean_sorted[k + l, 0, 1])
                    #print(j, l)
                    if floorlist[j] == pmt_group_mean_sorted[k + l, 1]: #case if combination is extant 
                        heatmap[j, m] = pmt_group_mean_sorted[k + l, 4]
                        l = l + 1
                    else: #case if combination non-extant 
                        heatmap[j, m] = np.nan
                        x = x + 1
                    j = j + 1
                #print(pmt_group_mean_sorted[k:i, 0, :2])
                m = m + 1; k = i; j = 0; l = 0
        return heatmap

    def heatmap_equal_bins_single_loop(self, heatmap, indices, group_no, stringlist, floorlist):
        #Define layout 
        annotation_text = np.round(heatmap, 4)
        print(heatmap)
        layout = go.Layout(
            title = "Rates of PMTs per DOM, PMTs %i - %i" %(indices[group_no], indices[group_no + 1]),
            xaxis = dict(
                tickmode = "array",
                tickvals = stringlist,
                ticktext = stringlist,
                dtick = 1
            ),
            yaxis = dict(
                tickmode = "array",
                tickvals = floorlist,
                ticktext = floorlist,
                dtick = 1
            )
        )
        fig = go.Figure(data = go.Heatmap(z=heatmap, text = annotation_text, texttemplate="%{text}"), layout = layout)
        fig.show()
        write_path = str('rates_doms_pmt_%i_%i.pdf' %(indices[group_no], indices[group_no + 1]))
        pio.write_image(fig, write_path)
        return 0;

    def pmt_group_data_to_heatmap(self, indices):
        """Wrapper for the above two functions. Automatically converts the pmt data into groups, normalises over given groups and 
        generates heatmap plots for all of them"""
        pmt_group_mean = self.normalise_over_n_pmts(indices)
        pmt_group_mean_sorted = pmt_group_mean
        #pmt_group_pairs = pmt_group_mean[:, 0, :]
        for m in range(0, len(indices)-1):
            pmt_group_pairs = pmt_group_mean[:, m, :]
            pmt_group_mean_sorted[:, m, :], str_floor_length = self.heatmap_averages_single_loop(pmt_group_pairs) #Sorts the thing on string and floor so that they are in sequence. 
        floorlist = np.unique(pmt_group_mean_sorted[:, 0, 1]); stringlist = np.unique(pmt_group_mean_sorted[:, 0, 0])
        heatmap = np.zeros([len(indices)-1, len(floorlist), len(np.unique(stringlist))])
        #Now get heatmap for each group
        print(len(indices)-1)
        for i in range(0, len(indices)-1):
            heatmap[i, :, :] = self.heatmap_array_single_group(pmt_group_mean_sorted[:, i, :], i, heatmap[i, :, :], floorlist, stringlist)
        #Plot heatmaps
        for i in range(0, len(indices)-1):
        #for i in range(0, 1):
            self.heatmap_equal_bins_single_loop(heatmap[i, :, :], indices, i, stringlist, floorlist)

    def export_heatmap(self, indices):
        pmt_group_mean = self.normalise_over_n_pmts(indices)
        pmt_group_mean_sorted = pmt_group_mean
        #pmt_group_pairs = pmt_group_mean[:, 0, :]
        for m in range(0, len(indices)-1):
            pmt_group_pairs = pmt_group_mean[:, m, :]
            pmt_group_mean_sorted[:, m, :], str_floor_length = self.heatmap_averages_single_loop(pmt_group_pairs) #Sorts the thing on string and floor so that they are in sequence. 
        floorlist = np.unique(pmt_group_mean_sorted[:, 0, 1]); stringlist = np.unique(pmt_group_mean_sorted[:, 0, 0])
        heatmap = np.zeros([len(indices)-1, len(floorlist), len(np.unique(stringlist))])
        #Now get heatmap for each group
        for i in range(0, len(indices)-1):
            heatmap[i, :, :] = self.heatmap_array_single_group(pmt_group_mean_sorted[:, i, :], i, heatmap[i, :, :], floorlist, stringlist)

        return heatmap, floorlist, stringlist 

class heatmap():
    def __init__(self, heatmap, floorlist, stringlist):
        self.heatmap = heatmap
        self.floorlist = floorlist
        self.stringlist = stringlist

    def plot_heatmap(self, indices):

        custom_colorscale = [
            [0.0, 'rgb(0, 0, 255)'],
            [1.0, 'rgb(255, 255, 0)'],
            [2.0, 'rgb(0, 0, 255)']
        ]
        for i in range(0, len(indices)-1):
            annotation_text = np.round(self.heatmap[i, :, :], 4)
            layout = go.Layout(
                title = "Ratio of simulated vs real rates of PMTs per DOM, PMTs %i - %i" %(indices[i], indices[i + 1]),
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
            fig = go.Figure(data = go.Heatmap(z=self.heatmap[i, :, :], text = annotation_text, texttemplate="%{text}"), layout = layout)
            fig.update_traces(
                hovertemplate='x: %{x}<br>y: %{y}<br>z: %{text}',
                textfont=dict(color='black')
                )
            fig.show()
            write_path = str('ratio_rates_doms_pmt_%i_%i.pdf' %(indices[i], indices[i + 1]))
            pio.write_image(fig, write_path)
        return 0;

def calc_heatmap_ratio(heatmap_real, heatmap_sim):
    return heatmap_sim/heatmap_real



indices = [0, 11, 18, 30]

data_real = map_hit_data(muon_hit_data_real, modid_map)
data_sim = map_hit_data(muon_hit_data_sim, modid_map)

real_map, floorlist, stringlist = data_real.export_heatmap(indices)
sim_map, __, __ = data_sim.export_heatmap(indices)

sim_ratio_map = heatmap(calc_heatmap_ratio(real_map, sim_map), floorlist, stringlist)
sim_ratio_map.plot_heatmap(indices)