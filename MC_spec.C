#include "Setting.cc"
const int NCanvas=2;


void MC_spec(){

  gRandom -> SetSeed( time(NULL) ); //seed set by time
  Setting *set = new Setting(); 

  // z=70 (4/11 first est)
  int NSi=700;
  int NBG=80000;

  // z=25 (4/12)
  NSi=2200;
  NBG=8000*4;

  double meanXray =6.18;//keV
  double sigmaXray=0.1;
  double BGmin=5.1;
  double BGmax=6.7;
  h_sum       = new TH1D("h_sum","",100/2,5.7, 6.7);
  h_S         = new TH1D("h_S"  ,"",100/2,5.7, 6.7);
  h_BG        = new TH1D("h_BG" ,"",100/2,5.7, 6.7);
  set->SetTH1(h_sum  ,"","X-ray [keV]" ,"Counts"    ,1,3000,0);
  set->SetTH1(h_S    ,"","X-ray [keV]" ,"Counts"    ,2,3001,2);
  set->SetTH1(h_BG   ,"","X-ray [keV]" ,"Counts"    ,4,3000,0);
  h_sum->SetMinimum(0);
  h_BG ->SetMinimum(0);

  for(int n=0;n<NSi;n++){
    double xray=gRandom->Gaus(meanXray, sigmaXray);
    h_S   ->Fill(xray);
    h_sum ->Fill(xray);
  }
  cout<<"signal done."<<endl;

  for(int n=0;n<NBG;n++){
    double xray=gRandom->Uniform(BGmin, BGmax);
    h_BG  ->Fill(xray);
    h_sum ->Fill(xray);
  }
  cout<<"noise done."<<endl;

  TCanvas *c[NCanvas];
  for(int i=0;i<NCanvas;i++){
    c[i]= new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1200,800 );
  }

  c[0]->cd();
  h_sum->Draw();
  h_S  ->Draw("same");

  c[1]->Divide(2,2);
  c[1]->cd(1);h_sum->Draw();
  c[1]->cd(2);h_S  ->Draw();
  c[1]->cd(3);h_BG ->Draw();
  c[1]->cd(4);h_sum->Draw();h_BG ->Draw("same");h_S  ->Draw("same");


  string pdf_name="pdf/ExpectedSpectrum_0412.pdf";
  c[0]->Print( Form("%s[",pdf_name.c_str()) );
  for(int i=0;i<NCanvas;i++){
    c[i]->Print( Form("%s",pdf_name.c_str()) );
  }
  c[NCanvas-1]->Print( Form("%s]",pdf_name.c_str()) );
}
