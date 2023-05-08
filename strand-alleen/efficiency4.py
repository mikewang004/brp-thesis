#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 15:54:39 2023

@author: mike
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
import ROOT 
import plotly.io as pio
import plotly.express as px

pio.renderers.default='browser'

muon_hit_data = np.load("muon_hit_data.npy")
pmt_id_map = np.loadtxt("map.txt")
a = np.isin(pmt_id_map[:, 0], muon_hit_data[:, 0, 0])
print(a)