import numpy as np 
#import wrapper_muonsim


aashower_startno = 137
infile_real = []
for i in range(0, 10):
    infile_real.append("/sps/km3net/repo/data_processing/tag/v8.1/data_processing/prod/data/KM3NeT_00000133/v8.1.1/reco/datav8.1.1.jchain.aashower.00%i%i.root" %(aashower_startno, 50 + i))
print(infile_real[8])