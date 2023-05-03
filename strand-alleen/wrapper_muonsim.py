#!/usr/bin/env python
import aa, ROOT
from ROOT import EventFile, Det, cxx, TH1D, Timer, time_at
import muonsim
from muonsim import residual_map

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

infile = ["/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3092.root", "/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco/mcv8.1.mupage_tuned_100G.sirene.jterbr00013754.jchain.aashower.3093.root"]

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


# plot the resduals for channel 1 of some DOM
rmap.geth( 819047388, 1).Draw()
rmap.geth(819047388, 1).SaveAs("test.png")
print(rmap.geth( 819047388, 1).GetNbinsX())
print(rmap.geth( 819047388, 1).Integral(-1000, 1000))