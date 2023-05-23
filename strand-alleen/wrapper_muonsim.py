#!/usr/bin/env python
import aa, ROOT
from ROOT import EventFile, Det, cxx, TH1D, Timer, time_at
import muonsim
from muonsim import residual_map
import numpy as np

#infile = "/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3094.root" #simulated tracks
#infile = "/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/data/KM3NeT_00000133/v8.1.1/reco/datav8.1.1.jchain.aashower.00013754.root" #real data

#infile = ["/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3092.root", 
#"/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3093.root",
#"/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3094.root",
#"/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3095.root",
#"/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3096.root",
#"/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3097.root",
#"/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3098.root",
#"/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3099.root",
#"/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3100.root",
#"/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3101.root",
#"/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3102.root"]

#infile = ["/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3092.root", "/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3093.root"]
#used numbers for above line are 3092-3 
infile = "/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3092.root"
f = EventFile( infile )


rmap = residual_map()

notrack = 0

for evt in f :

    try :
        trk = ROOT.get_best_reconstructed_jppmuon( evt )
    except ROOT.Exception: 
        notrack+=1
        continue

    rmap.fill( evt.hits, trk )


print (notrack,"out of", f.size(),"events had no track")
print()
print ("the following doms have at least one hit ")
print()
print (list(rmap.domids()))
print()

no_pmts = 31
hit_array = np.zeros([len(list(rmap.domids())),no_pmts,  3])
#Note structure as follows and note that pmts without hits are excluded 
#dom-ids / pmt no / no-hits 
domid_map = list(rmap.domids())
print("test if the map is working properly")
lowerbound, upperbound = -50, 60

for j in range(0, len(domid_map)):
    """Structure as follows: [dom number, pmt-number, amount of hits]"""
    for i in range(0, 31):
        hit_array[j, i, 0] = domid_map[j]
        hit_array[j, i, 1] = i
        hit_array[j, i, 2] = rmap.geth(domid_map[j], i).Integral(lowerbound,upperbound)

np.save("muon_hit_data2.npy", hit_array)
