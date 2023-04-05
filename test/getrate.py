#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 13:50:35 2023

@author: mike
"""

import numpy as np
import matplotlib.pyplot as plt
import ROOT 
max_no_xbins = 31

hit_data = ROOT.TFile.Open("jra_133_14307.root")
data2 = hit_data.Detector.DU10.F10.h_pmt_rate_distributions_Summaryslice # 31 x 100 bins



hit_rate_test = np.zeros([max_no_xbins, 100])

for i in range(0, max_no_xbins):
    for j in range(0, 100):
        hit_rate_test[i, j] = data2.Integral(i, i+1, j, j+1)
        
hit_rate_test = hit_rate_test.swapaxes(1, 0)

yaxis = data2.GetYaxis()
#Get approximation of binwidths
binarray = np.zeros([100, 2])
for i in range(0, 100):
    binarray[i, 0] = yaxis.GetBinWidth(i)
binarray[:, 1] = np.cumsum(binarray[:, 0]) #y-axis bin to rate 


np.savetxt("y-bin_size.txt", binarray)

#Plot array by binsize 





    


