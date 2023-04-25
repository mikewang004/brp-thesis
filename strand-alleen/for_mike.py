#!/usr/bin/env python

import aa, ROOT
from ROOT import EventFile, Det, cxx, TH1D, Timer, time_at
ROOT.gStyle.SetOptStat(0)




# we'll do the histogram filling in c++ for speed.

cxx("""
struct residual_map 
{
    map<int, map< int, TH1D*> > M;

    residual_map() {}

    /*! get histogram by dom and channel id -- create it if needed */
    
    TH1D* geth( int dom_id, int channel_id )
    {
        TH1D*& h = M[ dom_id][ channel_id];

        if (!h) 
        {
            string nam = "res"+str(dom_id)+"_"+str(channel_id);
            h = new TH1D(nam.c_str(),nam.c_str(), 200,-50,150);
        }
        return h;  
    }

    /*! fill the residuals */

    void fill( const vector<Hit>& hits , const Trk& trk ) 
    {
        for( auto& hit : hits )
        {
            double texp = time_at( trk, hit.pos );
            geth( hit.dom_id, hit.channel_id ) -> Fill ( hit.t - texp );
        }
    }

    /*! return all the domids that we know about */
    
    vector<int> domids() const
    {
        return keys( M );
    }
};
""")

from ROOT import residual_map


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

# plot the resduals for channel 1 of some DOM
rmap.geth( 819047388, 1).Draw()
  
    
