import numpy as np 
#import wrapper_muonsim


aashower_startno = 137
infile_real = []
for i in range(0, 10):
    infile_real.append("/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/data/KM3NeT_00000133/v8.1.1/reco/datav8.1.1.jchain.aashower.00%i%i.root" %(aashower_startno, 50 + i))
    

sim_events_list = []
sim_path = "/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/mc/atm_muon/KM3NeT_00000133/v8.1/reco"
with open("output.txt") as file:
    for line in file:
        sim_events_list.append(sim_path + line)

for i,n in enumerate(sim_events_list):
    sim_events_list[i] = n.strip()

print(sim_events_list[2])