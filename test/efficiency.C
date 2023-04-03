void efficiency() {
    char name[100];
    char dname[100];
    char fname[100];
    sprintf(fname,"%s","jra_133_14421.root");

    TFile *f = new TFile(fname,"read");

    for (int du=1; du++;du<33){ //loops over du
        for (int floor=1;floor<19;floor++){ //loops over strings 
            sprintf(dname,"%s%i%s%i","Detector/DU",du,"/F",floor);
            sprintf(name,"%s%i%s%i%s","Detector/DU",du,"/F",floor,"/h_pmt_rate_distributions_Summaryslice"); //2d plot; x-axis pmt channel, y rate;
    
}
