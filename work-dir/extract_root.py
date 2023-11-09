

# Test here whether what relation is between new/old pmts via rates vs eff plot. 

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as go
import plotly.colors as colors
import ROOT 

pio.kaleido.scope.mathjax= None
pio.renderers.default='browser'

modid_map = np.loadtxt("../pmt-info/map.txt")
pmt_serial_map = np.loadtxt("../pmt-info/pmt-serials.txt", usecols = 0)
pmt_ring_map = np.loadtxt("../pmt-info/pmt-ring.txt", skiprows = 2, usecols = [0,1,2])
magic_number = 16104 # The major version change happened at serial number 16104 (all PMTs <=16104 are of a certain kind (R12199), all abover are another one (R14374)).

indices = [0, 1, 7, 13, 19, 25, 31]
pmt_letters = ["A", "B", "C", "D", "E", "F"]
floorlist = np.loadtxt("data/floorlist.txt"); stringlist = np.loadtxt("data/stringlist.txt")

def fit_bin_size(filename):
    """Corrects for the fact that the y-bin sizes in the ROOT data are exponential.
    Returns fit parameters for a function y = a*exp(b*x) + c"""
    binsizes = np.loadtxt(filename)
    xbins = np.arange(0, len(binsizes[:, 0]))
    #fit now binsize to some sort of relation 
    popt, pcov = sp.curve_fit(exp_func, xbins, binsizes[:,1], p0 = [1, 0.1, 0])
    #plt.scatter(xbins, exp_func(xbins, *popt))
    return popt, pcov
def exp_func(x, a, b, c):
    return a*np.exp(b*x) + c
def get_du_floor_rate(du, floor, hit_data):
    """Loads in data corresponding to a du and floor."""
    domstr = "DU%i" %du; floorstr = "F%i" %floor
    domattr = getattr(hit_data.Detector, domstr)
    try:
        floorattr = getattr(domattr, floorstr)
        domfloordata = floorattr.h_pmt_rate_distributions_Summaryslice
        return domfloordata
    except:
        print(domattr)
        return None 
def get_du_floor_data(domfloordata, pmt_per_dom):
    """Converts rootpy histogram to python-handable data."""
    domfloorhitrate = np.zeros([pmt_per_dom, 100])
    for i in range(0, pmt_per_dom):
        for j in range(0, 100):
            domfloorhitrate[i, j] = domfloordata.Integral(i, i+1, j, j+1)
            
    domfloorhitrate = domfloorhitrate.swapaxes(1, 0)
    return domfloorhitrate
def gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
def lin_func(x, a):
    return a * x 

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


class extract_mean_hit_rate():

    def __init__(self, run_numbers, du_eff_map = "../pmt-info/map.txt"):
        self.run_numbers = run_numbers
        self.du_eff_map = np.loadtxt(du_eff_map)
        self.path = "../get-data/"
        self.pmt_per_dom = 31
        self.magic_number = 16104
        self.pmt_serial_map = np.loadtxt("../pmt-info/pmt-serials.txt", usecols = 0)

    def get_map_data(self, effs):
        """Substitutes module-id for corresponding efficiency and generates a mapping
        of efficiencies to du and floor number"""
        mapdata = np.zeros([len(effs[:, 0]), 5]) * np.nan
        for i in range(0, len(effs)):
            if i % self.pmt_per_dom == 0:
                for l_no, line in enumerate(self.du_eff_map): #first search for number then use that for next 31 entries.
                    if effs[i, 0] in line:
                        du, floor = line[1], line[2]
                        #Note PMT-channel is 0 through 30 repetitive 
            mapdata[i, :] = [effs[i, 0], effs[i, 1], effs[i, 2],  du, floor] # dom-id  / pmt-channel/ efficiency / no du/ no floor
        mapdata_2 = self.append_pmt_serials(mapdata)
        return mapdata_2

    def append_pmt_serials(self, mapdata):
        mapdata_new = np.zeros([mapdata.shape[0], mapdata.shape[1] + 2])
        mapdata_new[:, :mapdata.shape[1]] = mapdata
        l = 0
        print(mapdata_new.shape); print(self.pmt_serial_map.shape)
        for i in range(0, int(mapdata.shape[0]/31)):
            for j in range(0, len(self.pmt_serial_map)):
                if mapdata_new[l, 0] == self.pmt_serial_map[j]:
                    #print(self.pmt_serial_map[j])
                    for k in range(0, 31):
                        mapdata_new[k+l, 5] = self.pmt_serial_map[j+k+1]
                        if self.pmt_serial_map[j + k + 1] > magic_number:
                            mapdata_new[l + k , 6] = 1
                        else:
                            mapdata_new[l + k , 6] = 0
                    
                    break;
            l = 31 * (i + 1)
            if i >= int(mapdata.shape[0]):
                break
        #np.savetxt("debug.txt", mapdata_new)
        return mapdata_new # dom-id / pmt-channel / efficiency / string / floor / pmt-id / new pmt yes/no

    def read_mapdata(self, str_effs):
        """Reads in data, returns a map of the efficiencies to du and floor and returns 
        a list of dus and floors """
        effs = np.loadtxt(str_effs, skiprows = 148, usecols=[1,2,3])
        mapdata = self.get_map_data(effs)
        return mapdata

    def analysis_mul_runs(self):
        for i in range(0, len(self.run_numbers)):
                str_effs = self.path + "KM3NeT_00000133_000%i.v8.0_PMTeff_new.ToT.QE.PMTeff.txt" %(self.run_numbers[i])
                mapdata = self.read_mapdata(str_effs)
                mapdata = mapdata[np.lexsort((mapdata[:, 4], mapdata[:, 3]))]
                if i == 0:
                    effs = np.loadtxt(str_effs, skiprows = 148, usecols=[1,2,3])
                    mapdata_large = np.zeros([len(self.run_numbers), len(effs[:, 0]), 7])
                    del effs
                mapdata_large[i, :, :] = mapdata #runs vs length data vs cols containing efficiency / pmt-channel / no du / no floor
        self.mapdata_large = mapdata_large
        return mapdata_large

    def get_mean_runs(self, mean_hit_rate):
        return np.mean(mean_hit_rate, axis = 0), np.std(mean_hit_rate, axis = 0)

    def read_root_data(self):
        self.analysis_mul_runs()
        mean_hit_rate = np.zeros([self.mapdata_large.shape[0], self.mapdata_large.shape[1], 8]) * np.nan # eff / pmt / du / floor / mean or something
        #mean_hit_rate.astype(np.float32)
        mean_hit_rate[:, :, 0:7] = np.copy(self.mapdata_large)
        bin_popt, bin_pcov = fit_bin_size("data/y-bin_size.txt") #assume bin sizes constant over all runs
        for l in range(0, len(self.run_numbers)):
            hit_data = ROOT.TFile.Open(self.path + "jra_133_%i.root" %(self.run_numbers[l]))
            j = 0
            ytest = np.arange(0, 100)
            for i in range(0,int(self.mapdata_large.shape[1]/31)):
                dufloordata = get_du_floor_rate(int(mean_hit_rate[l, j, 3]), int(mean_hit_rate[l, j, 4]), hit_data)
                if dufloordata != None:
                    dufloorhitrate = get_du_floor_data(dufloordata, 31)
                    for k in range(0, 31):
                        test = gaussfit(ytest, dufloorhitrate[:, k]) #gets gaussian of hits per rate per pmt 
                        mean_hit_rate[l, j + k, 7] = exp_func(test.get_mean_coords()[0], *bin_popt) 
                #to convert bins to real unit as the y-bin size is log
                else:
                    mean_hit_rate[j:j+31, 7] = np.nan
                j = j + 31
        self.mean_hit_rate = mean_hit_rate
        mean_hit_rate, __ = self.get_mean_runs(mean_hit_rate)
        #np.save("mean_hit_rate.npy", mean_hit_rate)
        print(mean_hit_rate.shape)
        return mean_hit_rate

run_numbers = np.arange(14413,14440, 1)


test2 = extract_mean_hit_rate(run_numbers)
test2.read_root_data()
