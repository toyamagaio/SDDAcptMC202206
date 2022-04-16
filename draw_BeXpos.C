#include "Setting.cc"
const int NCanvas = 4;
const int NFiles=3;
const int NSDD=1;
//
void draw_BeXpos(){
  string rootfile[NFiles];
  rootfile[0] = "root/SDDAcptMC_pmuTest_SDD1_d55_n100000_z25.root";
  rootfile[1] = "root/SDDAcptMC_pmuTest_SDD1_d55_n100000_z25_pmuuni.root";
  rootfile[2] = "root/SDDAcptMC_pmuTest_SDD1_d55_n100000_z25_pmuGau.root";

  Setting *set = new Setting();

  TH2D *h2_yz_pmu_gene[NFiles], *h2_yz_pmu_gene_wBeHit[NFiles];
  TH1D *h_cost_pmu_gene[NFiles], *h_cost_pmu_gene_wBeHit[NFiles];

  TH1D *h_r_muBeXray_gene[NFiles], *h_theta_muBeXray_gene[NFiles];
  TH1D *h_x_muBeXray_gene[NFiles], *h_y_muBeXray_gene[NFiles];
  TH2D *h2_xy_muBeXray_gene[NFiles], *h2_xy_muBeXray_gene_wSDDall[NFiles];
  TH2D *h2_xy_muBeXray_wSDD[NFiles][NSDD], *h2_xy_muBeXray_atSDD[NFiles][NSDD];
  TH1D *h_cost_muBeXray_gene[NFiles], *h_cost_muBeXray_gene_wSDDall[NFiles];

  int col[]={1, 6, 4, 2};
  double max_x_muBe, max_y_muBe;
  max_x_muBe = max_y_muBe = 0.;
  for(int i=0;i<NFiles;i++){
    TFile *file = new TFile(rootfile[i].c_str(),"readonly");
    h2_yz_pmu_gene[i]      = (TH2D*)file->Get("h2_yz_pmu_gene");
    h2_xy_muBeXray_gene[i] = (TH2D*)file->Get("h2_xy_muBeXray_gene");
    h_x_muBeXray_gene[i]   = (TH1D*)file->Get("h_x_muBeXray_gene");
    h_y_muBeXray_gene[i]   = (TH1D*)file->Get("h_y_muBeXray_gene");

    h2_yz_pmu_gene[i]      ->SetMarkerColor(col[i]);
    h2_xy_muBeXray_gene[i] ->SetMarkerColor(col[i]);
    h_x_muBeXray_gene[i]   ->SetLineColor(col[i]);
    h_y_muBeXray_gene[i]   ->SetLineColor(col[i]);
    h_x_muBeXray_gene[i]   ->SetStats(0);
    h_y_muBeXray_gene[i]   ->SetStats(0);
    if(max_x_muBe < h_x_muBeXray_gene[i]->GetBinContent(h_x_muBeXray_gene[i]->GetMaximumBin()))max_x_muBe = h_x_muBeXray_gene[i]->GetBinContent(h_x_muBeXray_gene[i]->GetMaximumBin());
    if(max_y_muBe < h_y_muBeXray_gene[i]->GetBinContent(h_y_muBeXray_gene[i]->GetMaximumBin()))max_y_muBe = h_y_muBeXray_gene[i]->GetBinContent(h_y_muBeXray_gene[i]->GetMaximumBin());
  }
  h_x_muBeXray_gene[0]   ->GetYaxis()->SetRangeUser(0, 1.1*max_x_muBe);
  h_y_muBeXray_gene[0]   ->GetYaxis()->SetRangeUser(0, 1.1*max_y_muBe);

  TCanvas *c[NCanvas];
  for(int i=0;i<NCanvas;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1200, 1000);
  }
  c[0]->Divide(2,2);
  for(int i=0;i<NFiles;i++){
    c[0]->cd(1+i);h2_yz_pmu_gene[i]->Draw("colz");
  }

  c[1]->Divide(2,2);
  for(int i=0;i<NFiles;i++){
    c[1]->cd(1+i);h2_xy_muBeXray_gene[i]->Draw("colz");
  }

  c[2]->cd();
  h_x_muBeXray_gene[0]->Draw();
  for(int i=1;i<NFiles;i++){
    h_x_muBeXray_gene[i]   ->Draw("same");
  }

  c[3]->cd();
  h_y_muBeXray_gene[0]->Draw();
  for(int i=1;i<NFiles;i++){
    h_y_muBeXray_gene[i]   ->Draw("same");
  }

  string pdf_name="pdf/BeXpos.pdf";
  c[0]->Print( Form("%s[",pdf_name.c_str()) );
  for(int i=0;i<NCanvas;i++){
    c[i]->Print( Form("%s",pdf_name.c_str()) );
  }
  c[NCanvas-1]->Print( Form("%s]",pdf_name.c_str()) );
}
