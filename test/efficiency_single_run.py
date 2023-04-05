import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
import ROOT 

def gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def exp_func(x, a, b, c):
    return a*np.exp(b*x) + c

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
        popt, pcov = sp.curve_fit(gauss, self.xdata, self.ydata, p0=[max(self.ydata), self.mean(), self.sigma()])
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

#data2 = hit_data.Detector.DU10.F10.h_pmt_rate_distributions_Summaryslice

#plt.plot(data2)

#Get average hit rate 

# map map to effficiency rate 

def get_map_data(eff_map, effs, pmt_per_dom):
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

def get_unique_pairs(mapdata, counter, pmt_per_dom):
    pmt_eff = np.zeros([pmt_per_dom, 2])
    for i in range(0, pmt_per_dom):
        pmt_eff = mapdata[pmt_per_dom*counter:pmt_per_dom*(counter+1), 0:2] #only contains data for corresponding domfloor
    dom_floor = mapdata[pmt_per_dom*counter, 2:4]
    counter = counter + 1
    return np.roll(pmt_eff, 1, axis=1), dom_floor, counter



    
    
def get_dom_floor_rate(domfloor, hit_data):
    dom, floor = domfloor[0], domfloor[1]
    domstr = "DU%i" %dom; floorstr = "F%i" %floor 
    domattr = getattr(hit_data.Detector, domstr); floorattr = getattr(domattr, floorstr)
    domfloordata = floorattr.h_pmt_rate_distributions_Summaryslice
    return domfloordata
    
def get_domfloor_data(domfloordata, pmt_per_dom):
    domfloorhitrate = np.zeros([pmt_per_dom, 100])

    for i in range(0, pmt_per_dom):
        for j in range(0, 100):
            domfloorhitrate[i, j] = domfloordata.Integral(i, i+1, j, j+1)
            
    domfloorhitrate = domfloorhitrate.swapaxes(1, 0)
    return domfloorhitrate



#Fit gaussian to the pmt data 
def domfloorgaussian(domfloorhitrate, bin_popt, bin_pcov):
    """Fits gaussian to the data and returns mean of that"""
    """To do: parallise this if routine becomes slow"""
    pmt_mean_gauss_array = np.zeros(31)
    for i in range(0, 31):#range(0, len(domfloorhitrate[0,:])):
        currenthit = domfloorhitrate[:, i]
        print(currenthit)
        ytest = np.arange(0, 100)
        test = gaussfit(ytest, currenthit)
        pmt_mean_gauss_array[i] = test.get_mean_coords()[0]
        #test.gaussplot()
    #Scale pmt data correctly
    pmt_mean_gauss_array = exp_func(pmt_mean_gauss_array, *bin_popt)
    return pmt_mean_gauss_array


def fit_bin_size(filename):
    binsizes = np.loadtxt(filename)
    xbins = np.arange(0, len(binsizes[:, 0]))
    #fit now binsize to some sort of relation 
    popt, pcov = sp.curve_fit(exp_func, xbins, binsizes[:,1], p0 = [1, 0.1, 0])
    #plt.scatter(xbins, exp_func(xbins, *popt))
    return popt, pcov
        

#Link the dom/floor hits up to the efficiencies per dom 
def plot_hit_eff(pmt_eff, domfloormean):
    plt.plot()
    for i in range(0, 11):
        plt.scatter(pmt_eff[:,1], domfloormean, label="pmt %i" %pmt_eff[i,0])
    #plt.legend()
    plt.title('Hit rate vs efficiency; du 10 floor 3')
    plt.xlabel("Efficiency"); plt.ylabel("Rate [kHz]")
    plt.show()

#Doing whole plotting process for multiple doms 

def main():
    work_dir = "/Documents/uni_shit/zee-onzin/test"
    effs = np.loadtxt("data-133-144-eff.txt", skiprows = 149, usecols=[1,2,3])
    eff_map = np.loadtxt("map.txt")
    pmt_per_dom = 31
    hit_data = ROOT.TFile.Open("jra_133_14307.root")
    mapdata = get_map_data(eff_map, effs ,pmt_per_dom)
    counter = 0
    bin_popt, bin_pcov = fit_bin_size("y-bin_size.txt")
    max_runs = 100
    domfloormeanarray = np.zeros([31, max_runs])
    for i in range(0, max_runs):
        pmt_eff, domfloor, counter = get_unique_pairs(mapdata, counter, pmt_per_dom)
    domfloordata = get_dom_floor_rate(domfloor, hit_data)
    domfloorhitrate = get_domfloor_data(domfloordata, pmt_per_dom)
    domfloormean = domfloorgaussian(domfloorhitrate, bin_popt, bin_pcov)
    print(domfloormean)
    #domfloorhitrate = np.where(domfloorhitrate == 0, np.nan, domfloorhitrate)
    plot_hit_eff(pmt_eff, domfloormean)
    
if __name__ == "__main__":
    main()



