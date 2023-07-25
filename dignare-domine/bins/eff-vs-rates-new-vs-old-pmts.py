

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

modid_map = np.loadtxt("../../pmt-info/map.txt")
pmt_serial_map = np.loadtxt("../../pmt-info/pmt-serials.txt", usecols = 0)
pmt_ring_map = np.loadtxt("../../pmt-info/pmt-ring.txt", skiprows = 2, usecols = [0,1,2])
magic_number = 16104 # The major version change happened at serial number 16104 (all PMTs <=16104 are of a certain kind (R12199), all abover are another one (R14374)).

indices = [0, 1, 7, 13, 19, 25, 31]
pmt_letters = ["A", "B", "C", "D", "E", "F"]
floorlist = np.loadtxt("../data/floorlist.txt"); stringlist = np.loadtxt("../data/stringlist.txt")

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

    def __init__(self, run_numbers, du_eff_map = "../../pmt-info/map.txt"):
        self.run_numbers = run_numbers
        self.du_eff_map = np.loadtxt(du_eff_map)
        self.path = "../../get-data/"
        self.pmt_per_dom = 31
        self.magic_number = 16104
        self.pmt_serial_map = np.loadtxt("../../pmt-info/pmt-serials.txt", usecols = 0)

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
        bin_popt, bin_pcov = fit_bin_size("../data/y-bin_size.txt") #assume bin sizes constant over all runs
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
        return mean_hit_rate

class meanhitrate():
#TODO implement outliers filter, check for inconsistensies string 12 and 27. Also implement correct check on new/old pmt. 
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
                # print(index_all.shape)
                # print(np.nonzero(index_eff == False))
                # print()
                # print(np.nonzero(index_rate == False))
                # print()
                # print(np.nonzero(index_all == False))
                #self.meanhitrate[j:i, :, 2] = index_eff; self.meanhitrate[j:i, :, 7] = index_rate
                self.meanhitrate[j:i, :, 2] = np.where(index_all, self.meanhitrate[j:i, :, 2], np.nan)
                self.meanhitrate[j:i, :, 7] = np.where(index_all, self.meanhitrate[j:i, :, 7], np.nan)
                j = i
        return 0;


    #def filter_data(self):


    def plot_all_strings(self):
        "Plots rate vs efficieny for all strings"
        j = 0
        for i in range(0, self.top_avg.shape[0]-1 ):
            if self.top_avg[i, 3] != self.top_avg[i+1, 3]: #checks for start new string
                if int(self.top_avg[i, 6]) == 1:
                    color = "blue"
                else:
                    color = "orange"
                plt.figure()
                plt.title("Rate vs efficiency for du no %i" %(self.top_avg[i, 3]))
                plt.ylim(0, 12.5)
                plt.xlim(0, 1.4)
                plt.xlabel("Efficiency")
                plt.ylabel("Rate [kHz]")
                plt.scatter(self.top_avg[j:i, 2], self.top_avg[j:i, 7], label = "top pmts", c=color, marker="o") 
                plt.scatter(self.bottom_avg[j:i, 2], self.bottom_avg[j:i, 7], label = "bottom pmts", c=color, marker="s")
                plt.legend()
                plt.savefig("plots/rate_eff_string_%i.pdf" %(self.top_avg[i, 3]))
                j = i 
        return 0;


    def apply_pmt_mask(self, meanhitrate, new_versions=1, print_version = False):
        """Filters for either new pmt version or the old one. Data located in column no 7.
        If new_versions = 1, then only new versions are returned, otherwise only old versions returned."""
        newmeanhitrate = meanhitrate.copy()
        mask = newmeanhitrate[:, :, 6].astype(int)
        if new_versions == 1:
            mask = 1 - mask
        masked_floor_str_hit = np.ma.masked_array(newmeanhitrate[:, :, 2], mask)
        masked_floor_str_hit_2 = np.ma.masked_array(newmeanhitrate[:, :, 7], mask)
        newmeanhitrate[:, :, 2] = masked_floor_str_hit.filled(fill_value= np.nan)
        newmeanhitrate[:, :, 7] = masked_floor_str_hit_2.filled(fill_value = np.nan)
        if print_version == True:
            pass
            #print(newmeanhitrate[:, 0, 2])
            #print(np.count_nonzero(np.isnan(newmeanhitrate[:, :, 2])))
            #print(newmeanhitrate[:, :, 2].size)
        if newmeanhitrate[~np.isnan(newmeanhitrate[:, :, 2])].size == 0 or newmeanhitrate[~np.isnan(newmeanhitrate[:, :, 7])].size == 0:
            return None
        elif np.count_nonzero(np.isnan(newmeanhitrate[:, :, 2])) >= int(0.75*newmeanhitrate[:, :, 2].size):
            return None
        return newmeanhitrate

    def fit_lin_func_eff_rate(self, meanhitrate):
        rates, eff = meanhitrate[:, :, 7], meanhitrate[:, :, 2]
        rates, eff = rates[~np.isnan(rates)], eff[~np.isnan(eff)]
        popt, pcov = sp.curve_fit(lin_func, eff, rates)
        return popt

    def plot_all_strings_no_mean(self, pmt_start = 0, pmt_stop = 31, pmt_range = "all", plot = True):
        "Plots rate vs efficieny for all strings. Note pmts are arranged in the second dimension based on rings, so that "
        "ring E-F denoted by index 18-30; pmt_range = {lower, upper, all}"
        j = 0; x_eff = np.linspace(0, 1.4, 100); l = 0
        err = 0.05 # just assume this 
        popt_arr = np.zeros([2, int(self.meanhitrate.shape[0]/18)]) * np.nan
        if pmt_range == "lower":
            pmt_start, pmt_stop = 0, 18
            low_high_str = "lower"
        elif pmt_range == "upper":
            pmt_start, pmt_stop = 18, 31
            low_high_str = "upper"
        for i in range(1, self.meanhitrate.shape[0]): #checks for start new string
        #for i in range(1, 100):
            if self.meanhitrate[i-1, 0, 3] != self.meanhitrate[i,0,  3] or i == 377: #checks for start new string
                meanhitrate_old = self.apply_pmt_mask(self.meanhitrate[j:i, :, :], new_versions = 0)
                meanhitrate_new = self.apply_pmt_mask(self.meanhitrate[j:i, :, :])
                if plot == True:
                    plt.figure()
                    if type(meanhitrate_old) != type(None):
                        popt = self.fit_lin_func_eff_rate(meanhitrate_old[:, pmt_start:pmt_stop, :])
                        xerr, yerr = meanhitrate_old[:, :, 2] * err, meanhitrate_old[:, :, 7] * err
                        for j in range (pmt_start, pmt_stop):
                            plt.errorbar(meanhitrate_old[:, j, 2], meanhitrate_old[:, j, 7], 
                            xerr = xerr[:, j], yerr = yerr[:, j], fmt = ".", color = "blue") # version R12199
                        plt.plot(x_eff, lin_func(x_eff, *popt), color = "blue", label = "old PMT fit, slope %f [kHz]/[eff]" %popt)
                        popt_arr[0, l] = popt
                    if type(meanhitrate_new) != type(None):
                        popt = self.fit_lin_func_eff_rate(meanhitrate_new[:, pmt_start:pmt_stop, :])
                        xerr, yerr = meanhitrate_new[:, :, 2] * err, meanhitrate_new[:, :, 7] * err
                        for j in range(pmt_start, pmt_stop):
                            plt.errorbar(meanhitrate_new[:, j, 2], meanhitrate_new[:, j, 7], 
                            xerr = xerr[:, j], yerr = yerr[:, j], fmt = ".", color = "orange") #version R14374
                        plt.plot(x_eff, lin_func(x_eff, *popt), color = "orange", label = "new PMT fit, slope is %f [kHz]/[eff]" %popt)
                        popt_arr[1, l] = popt
                    plt.ylim(0, 10)
                    plt.xlim(0, 1.3)
                    plt.xlabel("Efficiency")
                    plt.ylabel("Rate [kHz]")
                    plt.legend()
                    if pmt_range == "all":
                        plt.title("Rate vs efficiency for DU %i, all PMTs" %(self.meanhitrate[i-1,0, 3]))
                        plt.savefig("plots/all-data/all_pmts_rate_eff_string_%i.pdf" %(self.meanhitrate[i-1, 0, 3]))
                    else:
                        plt.title("Rate vs efficiency for DU %i, %s PMTs only" %(self.meanhitrate[i-1,0, 3], low_high_str))
                        plt.savefig("plots/all-data/%s_pmts_rate_eff_string_%i.pdf" %(low_high_str, self.meanhitrate[i-1, 0, 3]))
                    j = i 
                    plt.close()
                else: #plot == False
                    if type(meanhitrate_old) != type(None):
                        popt = self.fit_lin_func_eff_rate(meanhitrate_old[:, pmt_start:pmt_stop, :])
                        popt_arr[0, l] = popt
                    if type(meanhitrate_new) != type(None):
                        popt2 = self.fit_lin_func_eff_rate(meanhitrate_new[:, pmt_start:pmt_stop, :])
                        popt_arr[1, l] = popt2
                l = l + 1
            if i == 377:
                return popt_arr

    def plot_all_strings_one_plot(self, pmt_start = 0, pmt_stop = 31, pmt_range = "all"):
        err = 0.05; x_eff = np.linspace(0, 1.4, 100);
        if pmt_range == "lower":
            pmt_start, pmt_stop = 0, 18
            low_high_str = "lower"
        elif pmt_range == "upper":
            pmt_start, pmt_stop = 18, 31
            low_high_str = "upper"
        plt.figure()
        meanhitrate_old = self.apply_pmt_mask(self.meanhitrate[:, :, :], new_versions = 0)
        meanhitrate_new = self.apply_pmt_mask(self.meanhitrate[:, :, :])
        popt = self.fit_lin_func_eff_rate(meanhitrate_old[:, pmt_start:pmt_stop, :])
        xerr, yerr = meanhitrate_old[:, :, 2] * err, meanhitrate_old[:, :, 7] * err
        for j in range (pmt_start, pmt_stop):
            plt.errorbar(meanhitrate_old[:, j, 2], meanhitrate_old[:, j, 7], 
            xerr = xerr[:, j], yerr = yerr[:, j], fmt = ".", color = "blue") # version R12199
        plt.plot(x_eff, lin_func(x_eff, *popt), color = "blue", label = "old PMT fit, slope %f [kHz]/[eff]" %popt)
        popt2 = self.fit_lin_func_eff_rate(meanhitrate_new[:, pmt_start:pmt_stop, :])
        xerr2, yerr2 = meanhitrate_new[:, :, 2] * err, meanhitrate_new[:, :, 7] * err
        for j in range (pmt_start, pmt_stop):
            plt.errorbar(meanhitrate_new[:, j, 2], meanhitrate_new[:, j, 7], 
            xerr = xerr2[:, j], yerr = yerr2[:, j], fmt = ".", color = "orange") # version R12199
        plt.plot(x_eff, lin_func(x_eff, *popt2), color = "orange", label = "new PMT fit, slope %f [kHz]/[eff]" %popt2)
        plt.ylim(0, 10)
        plt.xlim(0, 1.3)
        plt.xlabel("Efficiency")
        plt.ylabel("Rate [kHz]")
        plt.legend()
        if pmt_range == "all":
            plt.title("Rate vs efficiency for all DUs, all PMTs")
            plt.savefig("plots/all-data/all_dus_all_pmts_rate_eff_string.pdf")
        else:
            plt.title("Rate vs efficiency for DU, %s PMTs only" %(low_high_str))
            plt.savefig("plots/all-data/all_dus_%s_pmts_rate_eff_string.pdf" %(low_high_str))
        plt.close()
        return 0;

run_numbers = np.arange(14413,14440, 1)
#run_numbers = np.arange(14413, 14415, 1)
#test = extract_mean_hit_rate(run_numbers)
#test2 = meanhitrate(test.read_root_data())

test2 = meanhitrate(np.load("mean_hit_rate.npy"))

#popt_arr = test2.plot_all_strings_no_mean(pmt_range = "upper")
#test2.plot_all_strings_one_plot(pmt_range = "upper")
pmt_range_list = ["all", "upper", "lower"]
#pmt_range_list = ["all"]

def plot_all_eff_range_str_plots(pmt_range_list = pmt_range_list, plot = True):
    popt_arr = []
    for s in pmt_range_list:
        popt_arr.append(test2.plot_all_strings_no_mean(pmt_range = s, plot = plot))
    return popt_arr

popt_all, popt_upper, popt_lower = plot_all_eff_range_str_plots(plot = True)
print(popt_all)
# For DOM-gel issue:
# String 12 is 3 sigma away from average; affected DOMs in string 10 


#Now plot mean hit rate 


#test.append_pmt_serials(mapdata)

