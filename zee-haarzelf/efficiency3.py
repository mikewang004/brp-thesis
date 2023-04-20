#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 13:49:44 2023

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

def lin_func(x, a,b):
    return a*x

def get_mean_std(mean_rate_array, axis=0):
        return np.mean(mean_rate_array, axis=axis), np.std(mean_rate_array, axis =axis)

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
    
class meanhitrate():
    """Everything that has to do with the efficiency/rate arrays"""
    def __init__(self, mean_hit_rate):
        """Remember mean hit rate array has following structure:
            eff / pmt / du / floor / mean or something"""
        self.meanhitrate = mean_hit_rate
        
    def avg_top_bottom_pmts(self, mid_pmt):
        """Averages over the top (0-11) and bottom (12-31) pmts."""
        self.meanhitrate = self.meanhitrate[~np.isnan(self.meanhitrate).any(axis=1)]
        self.meanhitrate = self.meanhitrate.reshape(int(self.meanhitrate.shape[0]/31), 31, 5) #block per DOM 
        self.top_avg = np.mean(self.meanhitrate[:, 0:mid_pmt, :], axis=1)
        self.bottom_avg = np.mean(self.meanhitrate[:, mid_pmt:31, :], axis=1)
        self.pmtavg = np.zeros([self.meanhitrate.shape[0]*2, self.meanhitrate.shape[2]])
        self.pmtavg[::2, :] = self.top_avg; self.pmtavg[1::2, :] = self.bottom_avg;
        self.top_length = self.top_avg.shape[0]
        self.top_mask = self.top_avg[:-1] == self.top_avg[1:]
        return 0;
    
    def plot_top_bottom_pmts(self, fit=True):
        self.top_avg = self.filter_data(self.top_avg)
        self.bottom_avg = self.filter_data(self.bottom_avg)
        j = 0; k = 0
        if fit == True:
            self.fit_top_bottom_pmts()
            xfit = np.linspace(0, 1.4, 2)
        for i in range(0, self.top_length-1):
            if self.top_mask[i,2] == False: #checks if new du starts. 
                #Concatenate arrays together
                if fit == True:
                    #plt.plot(xfit, self.uplinres[k, 1] + self.uplinres[k,0] * xfit, label="top fit")
                    #plt.plot(xfit, self.lowlinres[k, 1] + self.lowlinres[k,0] * xfit, label="bottom fit")
                    
                    plt.plot(xfit, lin_func(xfit, *self.uppopt[k,:]), label="top fit")
                    plt.plot(xfit, lin_func(xfit, *self.lowpopt[k, :]),label="bottom fit")
                    k = k + 1
                plt.scatter(self.top_avg[j:i, 0], self.top_avg[j:i, 4], label="average of pmts 0-11")
                plt.scatter(self.bottom_avg[j:i, 0], self.bottom_avg[j:i, 4], label="average of pmts 12-30")
                plt.scatter(np.mean(self.top_avg[j:i, 0]), np.mean(self.top_avg[j:i, 4]), label="mean point top pmts")
                plt.scatter(np.mean(self.bottom_avg[j:i, 0]), np.mean(self.bottom_avg[j:i, 4]), label="mean point bottom pmts")
                plt.title("Rate vs efficiency for du no %i" %(self.top_avg[i, 2]))
                plt.ylim(0, 12.5)
                plt.xlim(0, 1.4)
                plt.xlabel("Efficiency")
                plt.ylabel("Rate [kHz]")
                plt.legend()
                plt.show()
                j = i  
                
    def plot_top_bottom_performance(self):
        k = 0
        plt.title("Performance of various DUs")
        plt.xlabel("DU number")
        plt.ylabel("Linear fit slope")
        for i in range(0, self.top_length-1):
            if self.top_mask[i,2] == False:
                plt.plot(self.top_avg[i, 2], self.uppopt[k,0], c='red', marker="o", ls="")
                plt.plot(self.bottom_avg[i, 2], self.lowpopt[k,0], c='blue', marker="o", ls="")
                k = k + 1
        plt.scatter(np.mean(self.top_avg[:,2]), np.mean(self.uppopt[:,0]), c='darkred', marker="X", label="average performance top pmts")
        plt.scatter(np.mean(self.bottom_avg[:,2]), np.mean(self.lowpopt[:,0]), c='darkblue', marker="X", label="average performance bottom pmts")
        plt.legend()
        plt.show()
    

    def fit_top_bottom_pmts(self):
        print(self.top_avg, self.bottom_avg)
        top_mask = self.top_avg[:-1] == self.top_avg[1:]
        no_dus = top_mask.shape[0] - top_mask[:, 2].sum()
        self.uppopt, self.uppcov = np.zeros([no_dus,2]), np.zeros([no_dus, 2, 2])
        self.lowpopt, self.lowpcov = np.zeros([no_dus,2]), np.zeros([no_dus, 2, 2])
        j = 0; k = 0
        for i in range(0, self.top_avg.shape[0]-1):
            if top_mask[i,2] == False:
                self.uppopt[k,:], self.uppcov[k, :, :] = sp.curve_fit(lin_func, self.top_avg[j:i, 0], self.top_avg[j:i, 4])
                self.lowpopt[k,:], self.lowpcov[k, :, :] = sp.curve_fit(lin_func, self.bottom_avg[j:i, 0], self.bottom_avg[j:i, 4])
                j = i; k = k + 1
        return 0;


class extract_mean_hit_rate():

    def __init__(self, run_numbers, low_pmt, high_pmt, du_eff_map = "map.txt"):
        self.run_numbers = run_numbers
        self.du_eff_map = np.loadtxt(du_eff_map)
        self.path = "../get-data/"
        self.low_pmt, self.high_pmt = low_pmt, high_pmt

    def get_map_data(self, effs, pmt_per_dom):
        """Substitutes module-id for corresponding efficiency and generates a mapping
        of efficiencies to du and floor number"""
        self.effs_length = len(effs[:, 0])
        mapdata = np.zeros([len(effs[:, 0]), 4])
        for i in range(0, len(effs)):
            if i % pmt_per_dom == 0:
                j = 0
                for l_no, line in enumerate(self.du_eff_map): #first search for number then use that for next 31 entries.
                    if effs[i,0] in line:
                        du, floor = line[1], line[2]
                        #Note PMT-channel is 0 through 30 repetitive 
            mapdata[i, :] = [effs[i, 2],j, du, floor] # efficiency / pmt-channel/ no du/ no floor
            j = j + 1
        return mapdata

    def read_mapdata(self, str_effs):
        """Reads in data, returns a map of the efficiencies to du and floor and returns 
        a list of dus and floors """
        effs = np.loadtxt(str_effs, skiprows = 148, usecols=[1,2,3])
        pmt_per_dom = self.high_pmt - self.low_pmt
        mapdata = self.get_map_data(effs ,pmt_per_dom)
        return mapdata

    def analysis_mul_runs(self):
        for i in range(0, len(self.run_numbers)):
                str_effs = self.path + "KM3NeT_00000133_000%i.v8.0_PMTeff_new.ToT.QE.PMTeff.txt" %(self.run_numbers[i])
                mapdata = self.read_mapdata(str_effs)
                mapdata = mapdata[np.lexsort((mapdata[:, 3], mapdata[:, 2]))]
                if i == 0:
                    effs = np.loadtxt(str_effs, skiprows = 148, usecols=[1,2,3])
                    mapdata_large = np.zeros([len(self.run_numbers), len(effs[:, 0]), 4])
                    del effs
                mapdata_large[i, :, :] = mapdata #runs vs length data vs cols containing efficiency / pmt-channel / no du / no floor
        return mapdata_large


def main():
    run_numbers = np.array([14399, 14400])
    test = extract_mean_hit_rate(run_numbers, 0, 31)
    test2 = test.analysis_mul_runs()
    print(test2.shape)

if __name__ == "__main__":
    main()