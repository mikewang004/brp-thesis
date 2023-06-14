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
        int a = 0 
        if (!h) 
        {
            string nam = "res"+str(dom_id)+"_"+str(channel_id);
            h = new TH1D(nam.c_str(),nam.c_str(), 200,-50,150); // TH1D (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup)
            h->SetDirectory(0);
            a = a + 1 
        }
        if (a == 50) 
        {
            h->Draw(0);
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



  
    
