import numpy as np
import matplotlib.pyplot as plt
import ROOT 

work_dir = "/Documents/uni_shit/zee-onzin/test"
effs = np.loadtxt("data-133-144-eff.txt", skiprows = 149, usecols=[1,2,3])
eff_map = np.loadtxt("map.txt")
pmt_per_dom = 31
hit_data = ROOT.TFile.Open("jra_133_14307.root")

#data2 = hit_data.Detector.DU10.F10.h_pmt_rate_distributions_Summaryslice

#plt.plot(data2)

#Get average hit rate 

# map map to effficiency rate 

def get_map_data(eff_map):
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

def get_unique_pairs(mapdata, counter):
    pmt_eff = np.zeros([pmt_per_dom, 2])
    for i in range(0, pmt_per_dom):
        pmt_eff = mapdata[pmt_per_dom*counter:pmt_per_dom*(counter+1), 0:2] #only contains data for corresponding domfloor
    dom_floor = mapdata[pmt_per_dom*counter, 2:4]
    counter = counter + 1
    return np.roll(pmt_eff, 1, axis=1), dom_floor, counter

mapdata = get_map_data(eff_map)
counter = 0
for i in range(0, 100):
    pmt_eff, domfloor, counter = get_unique_pairs(mapdata, counter)
    
    
def get_dom_floor_rate(domfloor):
    dom, floor = domfloor[0], domfloor[1]
    domstr = "DU%i" %dom; floorstr = "F%i" %floor 
    domattr = getattr(hit_data.Detector, domstr); floorattr = getattr(domattr, floorstr)
    domfloordata = floorattr.h_pmt_rate_distributions_Summaryslice
    return domfloordata
    
def get_domfloor_data(domfloordata):
    domfloorhitrate = np.zeros([pmt_per_dom, 100])

    for i in range(0, pmt_per_dom):
        for j in range(0, 100):
            domfloorhitrate[i, j] = domfloordata.Integral(i, i+1, j, j+1)
            
    domfloorhitrate = domfloorhitrate.swapaxes(1, 0)
    return domfloorhitrate

domfloordata = get_dom_floor_rate(domfloor)
domfloorhitrate = get_domfloor_data(domfloordata)

domfloorhitrate = np.where(domfloorhitrate == 0, np.nan, domfloorhitrate)

#print(domfloorhitrate)

#Link the dom/floor hits up to the efficiencies per dom 
def plot_hit_eff(pmt_eff, domfloorhitrate):
    plt.plot()
    for i in range(0, 11):
        eff_this_pmt = pmt_eff[i, 1] * np.ones(len(domfloorhitrate[:,1]))
        domfloorbinno = [~np.isnan(domfloorhitrate[:,i])] * np.arange(0, 100)
        plt.scatter(eff_this_pmt, domfloorbinno, label="pmt %i" %pmt_eff[i,0])
    plt.legend()
    plt.title('Hit rate vs efficiency; dom 10 floor 3')
    plt.ylim(20, 70)
    plt.xlabel("Efficiency"); plt.ylabel("Rate [bin number]")
    plt.show()

plot_hit_eff(pmt_eff, domfloorhitrate)



