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
    """Assume bins are along axis 0."""
    def __init__(self, x, y):
        self.xdata = x
        self.ydata = y
        
    def mean(self):
        return np.sum(self.xdata * self.ydata, axis=0) / np.sum(self.ydata, axis=0)
        
    def sigma(self):
        return np.sqrt(np.sum(self.ydata * (self.xdata - self.mean())**2,axis=0) / np.sum(self.ydata,axis=0))
    
    def gaussfit(self):
        no_pmts = len(self.xdata[0,:])
        lpopt, lpcov = np.zeros(no_pmts), np.zeros(no_pmts)
        for i in range(0, no_pmts):
            #Note: create dict or something as pmt and thus array size will not be larger than 31. 
            popt, pcov = sp.curve_fit(gauss, self.xdata[:,i], self.ydata[:,i], p0=[max(self.ydata[:,i]), self.mean()[i], self.sigma()[i]], maxfev = 50000)
            lpopt[i], lpcov[i] = popt, pcov
        return lpopt, lpcov
    
    def get_mean_coords(self):
        popt, pcov = self.gaussfit()
        mean = self.mean()
        return mean, gauss(mean, *popt)

#Get average hit rate 

# map map to effficiency rate 

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

def get_unique_pairs(mapdata, counter, pmt_per_dom):
    pmt_eff = np.zeros([pmt_per_dom, 2])
    for i in range(0, pmt_per_dom):
        pmt_eff = mapdata[pmt_per_dom*counter:pmt_per_dom*(counter+1), 0:2] #only contains data for corresponding domfloor
    dom_floor = mapdata[pmt_per_dom*counter, 2:4]
    return np.roll(pmt_eff, 1, axis=1), dom_floor


    
    
def get_dom_floor_rate(domfloor, hit_data):
    """Loads in data corresponding to a du and floor."""
    domstr = "DU%i" %domfloor[0]; floorstr = "F%i" %domfloor[1] 
    domattr = getattr(hit_data.Detector, domstr)
    try:
        floorattr = getattr(domattr, floorstr)
        domfloordata = floorattr.h_pmt_rate_distributions_Summaryslice
        return domfloordata
    except:
        print(domstr, domattr)
        return None 
    
def get_domfloor_data(domfloordata, pmt_per_dom):
    """Converts rootpy histogram to python-handable data."""
    domfloorhitrate = np.zeros([pmt_per_dom, 100])

    for i in range(0, pmt_per_dom):
        for j in range(0, 100):
            domfloorhitrate[i, j] = domfloordata.Integral(i, i+1, j, j+1)
            
    domfloorhitrate = domfloorhitrate.swapaxes(1, 0)
    return domfloorhitrate



#Fit gaussian to the pmt data 
def domfloorgaussian(domfloorhitrate, bin_popt, bin_pcov, pmt_per_dom):
    """Fits gaussian to the data and returns mean of that"""
    pmt_mean_gauss_array = np.zeros(pmt_per_dom)
    #print(domfloorhitrate)
    ytest = np.swapaxes(np.tile(np.arange(0,100),(pmt_per_dom,1)),1,0)
    test = gaussfit(ytest, domfloorhitrate)
    l = test.gaussfit()
    #for i in range(0, pmt_per_dom):#rang        print(popt, pcov)e(0, len(domfloorhitrate[0,:])):
        #currenthit = domfloorhitrate[:, i]
        #print(currenthit)
        #ytest = np.arange(0, 100)
        #test = gaussfit(ytest, currenthit)
        #pmt_mean_gauss_array[i] = test.get_mean_coords()[0]
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
    "Deprecated; see main file for new function."
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
    low_pmt, high_pmt = 0, 12
    #low_pmt, high_pmt = 12, 31
    pmt_per_dom = high_pmt - low_pmt
    hit_data = ROOT.TFile.Open("jra_133_14307.root")
    mapdata = get_map_data(eff_map, effs ,pmt_per_dom)
    bin_popt, bin_pcov = fit_bin_size("y-bin_size.txt")
    max_runs = 1
    domfloormeanarray = np.zeros([pmt_per_dom, max_runs]) #pmt x domfloors 
    domfloorlist = eff_map[:, 1:3]
    for i in range(0, max_runs):
        pmt_eff, domfloor = get_unique_pairs(mapdata, i, pmt_per_dom)
        domfloordata = get_dom_floor_rate(domfloor, hit_data)
        if domfloordata == None: #catch if data does not exist
            domfloormeanarray[:, i] = np.nan * np.ones([len(domfloormeanarray[:,i])])
        else:
            domfloorhitrate = get_domfloor_data(domfloordata, pmt_per_dom)
            domfloormean = domfloorgaussian(domfloorhitrate, bin_popt, bin_pcov, pmt_per_dom)
            domfloormeanarray[:, i] = domfloormean
    
    #Plot data now 
    #print(domfloormeanarray.shape) #axis0 pmt numbers, axis1 dom/floor combos
    plt.plot()
    for j in range(0, max_runs):
        pmt_eff, ___ = get_unique_pairs(mapdata, j, pmt_per_dom)
        domfloorrate = domfloormeanarray[:, j]
        plt.scatter(pmt_eff[low_pmt:high_pmt,1], domfloorrate[low_pmt:high_pmt], label="du %i floor %i" %(domfloorlist[j, 0], domfloorlist[j, 1]))
    plt.title('Hit rate vs efficiency; pmts %i through %i' %(low_pmt, high_pmt-1))
    plt.xlabel("Efficiency"); plt.ylabel("Rate [kHz]")
    plt.ylim(0, 12.5)
    plt.xlim(0, 1.4)
    #plt.legend()
    #plt.savefig("plotjes/week-14/hitrate-eff-pmts-%i-through-%i.pdf" %(low_pmt, high_pmt-1))
    plt.show()
        
    
if __name__ == "__main__":
    main()



