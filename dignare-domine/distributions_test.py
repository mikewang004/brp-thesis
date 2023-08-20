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


class meanhitrate():
#TODO implement outliers filter, check for inconsistensies string 12 and 27. Also implement correct check on new/old pmt. 
    """Everything that has to do with the efficiency/rate arrays"""
    def __init__(self, mean_hit_rate):
        """Remember mean hit rate array has following structure:
            [data, dom-id / pmt-no / eff / dom / string / pmt-serial / new pmt [yes/no] / rate [khz] ]"""
        self.meanhitrate = mean_hit_rate
        self.block_per_pmt()
        self.filter_data()
        print(self.meanhitrate.shape)

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


run_numbers = np.arange(14413,14440, 1)
#run_numbers = np.arange(14413, 14415, 1)
#test = extract_mean_hit_rate(run_numbers)
#test2 = meanhitrate(test.read_root_data())

test2 = meanhitrate(np.load("bins/mean_hit_rate.npy"))