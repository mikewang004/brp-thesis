void twodhist(){

  char name[100];
  char dname[100];
  char fname[100];
 // double g1,g2,g3;
  double g2 = 1.0;
  TF1 *gausfit = new TF1("gausfit", "[0]*exp(-0.5*((x-[1])/[2])**2)");


 sprintf(fname,"%s","jra_133_14421.root");

 TFile *f = new TFile(fname,"read");

 for (int du=1; du++;du<33){ //loops over du
  for (int floor=1;floor<19;floor++){ //loops over strings 

  sprintf(dname,"%s%i%s%i","Detector/DU",du,"/F",floor);
  sprintf(name,"%s%i%s%i%s","Detector/DU",du,"/F",floor,"/h_pmt_rate_distributions_Summaryslice"); //2d plot; x-axis pmt channel, y rate;

  if (f->GetDirectory(dname)){

   TH2D* hi = (TH2D*)f->Get(name);

   for (int ch=1;ch<32;ch++){
     TH1D* proj = hi->ProjectionY("proj",ch,ch);

    //g2=proj->GetMean();
    //int n = proj->GetXaxis()->GetNbins();  
    //std::vector<double> x(n);
    //proj->GetXaxis()->GetCenter( &x[0] );
    //const double *y = proj->GetArray(); 

    // gausfit->SetParameter(0,proj->GetMaximum());
    // gausfit->SetParameter(1,proj->GetXaxis()->GetBinCenter(proj->GetMaximumBin()));
    // gausfit->SetParameter(2,0.2);
    
    // g2=proj->GetXaxis()->GetBinCenter(proj->GetMaximumBin());
    // int fitStatus=proj->Fit("gausfit","OQ");


    // if (fitStatus>=0){
     // g1=proj->GetFunction("gausfit")->GetParameter(0);
     // g2=proj->GetFunction("gausfit")->GetParameter(1);
     // g3=proj->GetFunction("gausfit")->GetParameter(2);

//=================  these lines are only to demonstrate plotting the histogram
    if (g2 > 0.5){
     cout <<"HELLO!" <<endl;
     hi->GetXaxis()->SetRangeUser(0.1,30.);
     hi->GetYaxis()->SetRangeUser(0.1,20.);
     hi->GetXaxis()->SetTitle("Detector Unit");
     hi->GetYaxis()->SetTitle("Floors");
     //proj->Draw();
     hi->Draw();
     return;
    }
//=======================================================

   }

    // cout << du << " " << floor << " " << ch << " " << g2 <<  endl;
}
  }
 }
}
//return;
//}
