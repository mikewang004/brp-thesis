#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 13:56:23 2023

@author: mike
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
import ROOT 

def gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def exp_func(x, a, b, c):
    return a*np.exp(b*x) + c

def lin_func(x, a, b):
    return a*x + b

class gaussfit():
    """Routine to get gaussian data fitted and plotted. Only requires
    the x and y data assuming it is gaussian."""
    def __init__(self, x, y):
        self.xdata = x
        self.ydata = y
        
    def gaussfunction(x, a, x0, sigma):
        return a* np.exp(-(x-x0)**2/(2*sigma**2))
        
    def mean(self):
        return np.sum(self.xdata * self.ydata) / np.sum(self.ydata)
        
    def sigma(self):
        return np.sqrt(np.sum(self.ydata * (self.xdata - self.mean())**2) / np.sum(self.ydata))
    
    def gaussfit(self):
        popt, pcov = sp.curve_fit(gauss, self.xdata, self.ydata, p0=[max(self.ydata), self.mean(), self.sigma()], maxfev = 50000)
        return popt, pcov
    
    def gaussplot(self):
        plt.plot()
        popt, pcov = self.gaussfit()
        plt.scatter(self.xdata, self.ydata, label="raw data")
        plt.plot(self.xdata, gauss(self.xdata, *popt), label="gauss fit", color="orange")
        plt.title("Gaussian of a PMT rate distribution.")
        plt.xlabel("Rate [kHz]")
        plt.legend()
        plt.show()
        return popt, pcov
    
    def get_mean_coords(self):
        popt, pcov = self.gaussfit()
        mean = self.mean()
        return mean, gauss(mean, *popt)
    
def fit_bin_size(filename):
    binsizes = np.loadtxt(filename)
    xbins = np.arange(0, len(binsizes[:, 0]))
    #fit now binsize to some sort of relation 
    popt, pcov = sp.curve_fit(exp_func, xbins, binsizes[:,1], p0 = [1, 0.1, 0])
    #plt.scatter(xbins, exp_func(xbins, *popt))
    return popt, pcov
    
    
def get_map_data(eff_map, effs, pmt_per_dom):
    """Substitutes module-id for corresponding efficiency and generates a mapping
    of efficiencies to du and floor number"""
    mapdata = np.zeros([len(effs[:, 0]), 4])
    for i in range(0, len(effs)):
        if i % pmt_per_dom == 0:
            j = 0
            for l_no, line in enumerate(eff_map): #first search for number then use that for next 31 entries.
                if effs[i,0] in line:
                    du, floor = line[1], line[2]
                    #Note PMT-channel is 0 through 30 repetitive 
        mapdata[i, :] = [effs[i, 2],j, du, floor] # efficiency / pmt-channel/ no du/ no floor
        j = j + 1
    return mapdata

def read_mapdata(str_effs, str_effs_map, low_pmt, high_pmt):
    """Reads in data, returns a map of the efficiencies to du and floor and returns 
    a list of dus and floors """
    effs = np.loadtxt(str_effs, skiprows = 149, usecols=[1,2,3])
    eff_map = np.loadtxt(str_effs_map)
    pmt_per_dom = high_pmt - low_pmt
    mapdata = get_map_data(eff_map, effs ,pmt_per_dom)
    return eff_map, mapdata

def get_du_floor_rate(du, floor, hit_data):
    """Loads in data corresponding to a du and floor."""
    domstr = "DU%i" %du; floorstr = "F%i" %floor
    domattr = getattr(hit_data.Detector, domstr)
    try:
        floorattr = getattr(domattr, floorstr)
        domfloordata = floorattr.h_pmt_rate_distributions_Summaryslice
        return domfloordata
    except:
        print(domstr, domattr)
        return None 
    
def get_du_floor_data(domfloordata, pmt_per_dom):
    """Converts rootpy histogram to python-handable data."""
    domfloorhitrate = np.zeros([pmt_per_dom, 100])
    for i in range(0, pmt_per_dom):
        for j in range(0, 100):
            domfloorhitrate[i, j] = domfloordata.Integral(i, i+1, j, j+1)
            
    domfloorhitrate = domfloorhitrate.swapaxes(1, 0)
    return domfloorhitrate

def calc_hit_rate(mapdata):
    """To do: fix loop here"""
    hit_data = ROOT.TFile.Open("jra_133_14307.root")
    mean_hit_rate = np.zeros([mapdata.shape[0], 5]) # eff / pmt / du / floor / mean or something
    mean_hit_rate[:, 0:4] = np.copy(mapdata)
    j = 0
    bin_popt, bin_pcov = fit_bin_size("y-bin_size.txt")
    ytest = np.arange(0, 100)
    #for i in range(0,int(11718/31)):
    for i in range(0,int(mapdata.shape[0]/31)):
        dufloordata = get_du_floor_rate(int(mean_hit_rate[j, 2]), int(mean_hit_rate[j, 3]), hit_data)
        if dufloordata != None:
            dufloorhitrate = get_du_floor_data(dufloordata, 31)
            for k in range(0, 31):
                test = gaussfit(ytest, dufloorhitrate[:, k])
                mean_hit_rate[j + k, 4] = exp_func(test.get_mean_coords()[0], *bin_popt)
        else:
                mean_hit_rate[j:j+31, 4] = np.nan
        j = j + 31
    return mean_hit_rate


def avg_top_bottom_pmts(mean_hit_rate, mid_pmt):
    mean_hit_rate = mean_hit_rate[~np.isnan(mean_hit_rate).any(axis=1)]
    mean_hit_rate = mean_hit_rate.reshape(int(mean_hit_rate.shape[0]/31), 31, 5)
    top_avg = np.mean(mean_hit_rate[:, 0:mid_pmt, :], axis=1)
    bottom_avg = np.mean(mean_hit_rate[:, mid_pmt:31, :], axis=1)
    print(top_avg)
    return top_avg, bottom_avg

def recombine_duf_top_bottom_pmts(top_avg, bottom_avg):
    n_groups = 2
    top_mask = top_avg[:-1] == top_avg[1:]
    no_doms = top_mask.shape[0]-np.count_nonzero(top_mask[:, 2])
    #grouped_doms = np.zeros([no_doms, 18, 2 + n_groups*3])
    j = 0
    for i in range(0, top_avg.shape[0]-1):
        if top_mask[i,2] == False:
            #Concatenate arrays together
            xy = np.zeros([2, np.abs(j-i)*2])
            halfway = np.abs(j-i)
            xy[:, :halfway] = top_avg[j:i, 0], top_avg[j:i, 4]
            xy[:, halfway:halfway*2] = bottom_avg[j:i, 0], bottom_avg[j:i, 4]
            print(xy.shape)
            popt,pcov = sp.curve_fit(lin_func, xy[0, :], xy[1, :])
            plt.plot(xy[0,:], lin_func(xy[0,:], *popt))
            plt.scatter(top_avg[j:i, 0], top_avg[j:i, 4], label="average of pmts 0-11")
            plt.scatter(bottom_avg[j:i, 0], bottom_avg[j:i, 4], label="average of pmts 12-30")
            plt.title("Rate vs efficiency for du no %i" %(top_avg[i, 2]))
            plt.ylim(0, 12.5)
            plt.xlim(0, 1.4)
            plt.xlabel("Efficiency")
            plt.ylabel("Rate [kHz]")
            plt.legend()
            plt.show()
            j = i
    
def main():
    str_effs = "data-133-144-eff.txt"; str_effs_map = "map.txt"
    low_pmt, high_pmt, mid_pmt = 0, 31, 12
    pmt_per_dom = high_pmt - low_pmt
    eff_map, mapdata = read_mapdata(str_effs, str_effs_map, low_pmt, high_pmt)
    mapdata = mapdata[np.lexsort((mapdata[:, 3], mapdata[:, 2]))]
    mean_hit_rate = (calc_hit_rate(mapdata))
    top_avg, bottom_avg = avg_top_bottom_pmts(mean_hit_rate, mid_pmt)
    recombine_duf_top_bottom_pmts(top_avg, bottom_avg)
    
if __name__ == "__main__":
    main()


























