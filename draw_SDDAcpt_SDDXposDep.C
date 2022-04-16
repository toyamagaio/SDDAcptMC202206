#include "Setting.cc"
const int NCanvas = 4;
const int NFiles=7;
const int NSDD=1;
//
void draw_SDDAcpt_SDDXposDep(){
  gStyle->SetOptStat("ie");
  string rootfile[NFiles];
  rootfile[0] = "root/SDDAcptMC_pmuTest_SDD1_d55_n100000_z25_pmuGau_SDDx-20.root";
  rootfile[1] = "root/SDDAcptMC_pmuTest_SDD1_d55_n100000_z25_pmuGau_SDDx-15.root";
  rootfile[2] = "root/SDDAcptMC_pmuTest_SDD1_d55_n100000_z25_pmuGau_SDDx-10.root";
  rootfile[3] = "root/SDDAcptMC_pmuTest_SDD1_d55_n100000_z25_pmuGau_SDDx-5.root";
  rootfile[4] = "root/SDDAcptMC_pmuTest_SDD1_d55_n100000_z25_pmuGau_SDDx0.root";
  rootfile[5] = "root/SDDAcptMC_pmuTest_SDD1_d55_n100000_z25_pmuGau_SDDx5.root";
  rootfile[6] = "root/SDDAcptMC_pmuTest_SDD1_d55_n100000_z25_pmuGau_SDDx10.root";

  Setting *set = new Setting();

  double  sddx[NFiles]={-20, -15, -10, -5, 0, 5, 10};
  double esddx[NFiles]={  0,   0,   0,  0, 0, 0,  0};
  TH2D *h2_yz_pmu_gene[NFiles], *h2_yz_pmu_gene_wBeHit[NFiles];
  TH1D *h_cost_pmu_gene[NFiles], *h_cost_pmu_gene_wBeHit[NFiles];

  TH1D *h_r_muBeXray_gene[NFiles], *h_theta_muBeXray_gene[NFiles];
  TH1D *h_x_muBeXray_gene[NFiles], *h_y_muBeXray_gene[NFiles];
  TH2D *h2_xy_muBeXray_gene[NFiles], *h2_xy_muBeXray_gene_wSDDall[NFiles];
  TH2D *h2_xy_muBeXray_wSDD[NFiles][NSDD], *h2_xy_muBeXray_atSDD[NFiles][NSDD];
  TH1D *h_cost_muBeXray_gene[NFiles], *h_cost_muBeXray_gene_wSDDall[NFiles];

  double NmuBeXray[NFiles], NSDDAcpt[NFiles], SDDAcptRatio[NFiles];
  double eNmuBeXray[NFiles], eNSDDAcpt[NFiles], eSDDAcptRatio[NFiles];

  TGraphErrors *tg_NmuBeXray, *tg_NSDDAcpt, *tg_SDDAcptRatio;
  double max_x_muBe, max_y_muBe;
  max_x_muBe = max_y_muBe = 0.;
  for(int i=0;i<NFiles;i++){
    TFile *file = new TFile(rootfile[i].c_str(),"readonly");
    h2_yz_pmu_gene[i]              = (TH2D*)file->Get("h2_yz_pmu_gene");
    h2_xy_muBeXray_gene[i]         = (TH2D*)file->Get("h2_xy_muBeXray_gene");
    h2_xy_muBeXray_gene_wSDDall[i] = (TH2D*)file->Get("h2_xy_muBeXray_gene_wSDDall");
    h2_xy_muBeXray_atSDD[0][i]     = (TH2D*)file->Get(Form("h2_xy_muBeXray_atSDD%d",1));
    h_x_muBeXray_gene[i]           = (TH1D*)file->Get("h_x_muBeXray_gene");
    h_y_muBeXray_gene[i]           = (TH1D*)file->Get("h_y_muBeXray_gene");

    h2_xy_muBeXray_gene[i]         ->SetStats(1);
    h2_xy_muBeXray_gene_wSDDall[i] ->SetStats(1);
    h2_xy_muBeXray_atSDD[0][i]     ->SetStats(1);
    h_x_muBeXray_gene[i]   ->SetStats(0);
    h_y_muBeXray_gene[i]   ->SetStats(0);
    if(max_x_muBe < h_x_muBeXray_gene[i]->GetBinContent(h_x_muBeXray_gene[i]->GetMaximumBin()))max_x_muBe = h_x_muBeXray_gene[i]->GetBinContent(h_x_muBeXray_gene[i]->GetMaximumBin());
    if(max_y_muBe < h_y_muBeXray_gene[i]->GetBinContent(h_y_muBeXray_gene[i]->GetMaximumBin()))max_y_muBe = h_y_muBeXray_gene[i]->GetBinContent(h_y_muBeXray_gene[i]->GetMaximumBin());

    NmuBeXray[i]    = h2_xy_muBeXray_gene[i]        ->GetEntries();
    NSDDAcpt[i]     = h2_xy_muBeXray_gene_wSDDall[i]->GetEntries();
    eNmuBeXray[i]   = sqrt(NmuBeXray[i]);
    eNSDDAcpt[i]    = sqrt(NSDDAcpt[i] );
    SDDAcptRatio[i] = NSDDAcpt[i]/NmuBeXray[i];
    eSDDAcptRatio[i] = sqrt(SDDAcptRatio[i]*(1.-SDDAcptRatio[i])/NmuBeXray[i]);
    cout<<sddx[i] <<" "<<SDDAcptRatio[i]<<endl;
  }
  h_x_muBeXray_gene[0]   ->GetYaxis()->SetRangeUser(0, 1.1*max_x_muBe);
  h_y_muBeXray_gene[0]   ->GetYaxis()->SetRangeUser(0, 1.1*max_y_muBe);

  tg_NmuBeXray    =new TGraphErrors(NFiles, sddx, NmuBeXray   , esddx, eNmuBeXray   );
  tg_NSDDAcpt     =new TGraphErrors(NFiles, sddx, NSDDAcpt    , esddx, eNSDDAcpt    );
  tg_SDDAcptRatio =new TGraphErrors(NFiles, sddx, SDDAcptRatio, esddx, eSDDAcptRatio);
  set->SetGrErr(tg_NmuBeXray   ,"","","", 1, 1, 22, 1.0);
  set->SetGrErr(tg_NSDDAcpt    ,"","","", 2, 2, 23, 1.0);
  set->SetGrErr(tg_SDDAcptRatio,"","","", 4, 4, 24, 1.0);


  TCanvas *c[NCanvas];
  for(int i=0;i<NCanvas;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1200, 1000);
  }
  c[0]->Divide(3,3);
  for(int i=0;i<NFiles;i++){
    c[0]->cd(1+i);h2_xy_muBeXray_gene_wSDDall[i]->Draw("colz");
  }

  c[1]->Divide(3,3);
  for(int i=0;i<NFiles;i++){
    c[1]->cd(1+i);h2_xy_muBeXray_gene[i]->Draw("colz");
  }

  c[2]->Divide(3,3);
  for(int i=0;i<NFiles;i++){
    c[2]->cd(1+i);h2_xy_muBeXray_atSDD[0][i]->Draw("colz");
  }

  c[3]->Divide(2,2);
  c[3]->cd(1);tg_NmuBeXray    ->Draw("AP");
  c[3]->cd(2);tg_NSDDAcpt     ->Draw("AP");
  c[3]->cd(3);tg_SDDAcptRatio ->Draw("AP");
  //h_frame->Draw();

  string pdf_name="pdf/SDDXposDep.pdf";
  c[0]->Print( Form("%s[",pdf_name.c_str()) );
  for(int i=0;i<NCanvas;i++){
    c[i]->Print( Form("%s",pdf_name.c_str()) );
  }
  c[NCanvas-1]->Print( Form("%s]",pdf_name.c_str()) );
}
