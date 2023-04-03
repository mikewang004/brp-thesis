import numpy as np
import matplotlib.pyplot as plt
import ROOT 

work_dir = "/Documents/uni_shit/zee-onzin/test"
effs = np.loadtxt("data-133-144-eff.txt", skiprows = 149, usecols=[1,2,3])
eff_map = np.loadtxt("map.txt")
hit_data = ROOT.TFile.Open("jra_133_14307.root")

data2 = hit_data.Detector.DU10.F10.h_pmt_rate_distributions_Summaryslice

plt.plot(data2)

#Get average hit rate 

# map map to effficiency rate 
mapdata = np.zeros([len(effs[:, 0]), 4])
for i in range(0, len(effs)):
    if i % 31 == 0:
        j = 0
        for l_no, line in enumerate(eff_map): #first search for number then use that for next 31 entries.
            if effs[i,0] in line:
                du, floor = line[1], line[2]
    #Note PMT-channel is 0 through 31 repetitive 
    mapdata[i, :] = [effs[i, 2],j, du, floor] # efficiency / no du/ no floor
    print(mapdata[i, :])
    j = j + 1
                
    
    



        



