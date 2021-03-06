#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
using namespace std;

#include "TApplication.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "TRandom.h"
#include "Setting.h"

const int NCanvas=5;
const int NSDD=1;//num. of SDD

class SDDAcptMC
{
 public:
         SDDAcptMC();
        ~SDDAcptMC();
  void MakeHist(string ofname);
  void InitSDDsetup();
  void Loop();
  bool Draw(); 
  void SaveCanvas(string PdfFileName); 
  void SetPmuDistType( int type )  { pmu_type = type; }
  void SetMaxEvent( int N )  { ENumMax = N; }
  void SetDiameter( double d )  { D_Be = d; }
  void SetBranch();
  void SetNBeXray(bool flag){fix_NmuBeXray_flag=flag;}
  void SetSDDx(double x){SDDx=x;}
  void   CalcMap();
  void InitializeTreeBr();

  private:
    Setting *set;
    bool fix_NmuBeXray_flag;
    int ENum, ENumMax;
    int pmu_type;
    int BeHitCounter, SDDHitCounter;
    double D_Be;//diameter of Be target [mm]
    double SDDx;//xpos of SDD
    TVector3 GetUniformDistribution(double costMin=-1., double costMax = 1.);
    void ConvertToNKS2Coordinate(TVector3 &dir);
    double GetUniformR(double Rmax);
    int  IsIntersectSDDSurface(TVector3 posX, TVector3 dirX, TVector3 &sectVec);
    bool IsIntersectBeTarget(TVector3 posPmu, TVector3 dirPmu, TVector3 &sectVec);
    TH1D *h_r_pmu_gene, *h_theta_pmu_gene, *h_x_pmu_gene;
    TH2D *h2_yz_pmu_gene, *h2_yz_pmu_gene_wBeHit;
    TH1D *h_cost_pmu_gene, *h_cost_pmu_gene_wBeHit;

    TH1D *h_r_muBeXray_gene, *h_theta_muBeXray_gene;
    TH1D *h_x_muBeXray_gene, *h_y_muBeXray_gene;
    TH2D *h2_xy_muBeXray_gene, *h2_xy_muBeXray_gene_wSDDall;
    TH2D *h2_xy_muBeXray_wSDD[NSDD], *h2_xy_muBeXray_atSDD[NSDD];
    TH1D *h;
    TH1D *h_cost_muBeXray_gene, *h_cost_muBeXray_gene_wSDDall;
    TH2F *h_frame;
    TF1 *f1;

    bool treeout_flag;
    TVector3 dir_muBe_Xray, pos_pmu_stop, pos_at_SDD;
    TVector3 pos_SDD_surf[NSDD];
    double dia_SDD_surf[NSDD], til_SDD_surf[NSDD];

    TFile *ofp;
    TLatex *latex;
    TLine *tline;
    TTree *tree_out;

    TCanvas *c[NCanvas];

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SDDAcptMC::SDDAcptMC()
{

  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");
  gROOT->SetBatch(1);

  gStyle->SetOptDate(0);
  gStyle->SetOptFit(1);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetLineWidth(0);
  gStyle->SetOptDate(0);
  //gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleTextColor(1);
  gStyle->SetGridWidth(0);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetLineWidth(0);
  gStyle->SetNdivisions(510); // tertiary*10000 + secondary*100 + first
  gStyle->SetOptStat("iMen");
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);

  //const Int_t NRGBs = 5;
  //const Int_t NCont = 255;
  //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  //Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  //Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  //Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  //TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  //gStyle->SetNumberContours(NCont);
      
  set = new Setting();

  for(int i=0;i<NCanvas;i++){
    c[i]= new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1200,800 );
  }

  latex = new TLatex();
  tline  = new TLine();
  set->SetTLine(tline,2,2,1);

  gRandom -> SetSeed( time(NULL) ); //seed set by time
  fix_NmuBeXray_flag=false;
  SDDx  = 0.;
}
////////////////////////////////////////////////////////////////////////////
SDDAcptMC::~SDDAcptMC(){

  ofp->Write();
  ofp->Close();
}
////////////////////////////////////////////////////////////////////////////
void SDDAcptMC::MakeHist(string ofname){
  cout<<"MakeHist"<<endl;
  ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  cout<<"newed ofp"<<endl;

  ofp->cd();

  h_cost_pmu_gene         = new TH1D("h_cost_pmu_gene"        ,"", 200,  -1, 1);
  h_cost_pmu_gene_wBeHit  = new TH1D("h_cost_pmu_gene_wBeHit" ,"", 200,  -1, 1);
  h_r_pmu_gene            = new TH1D("h_r_pmu_gene"           ,"", 100,   0, 50);
  h_theta_pmu_gene        = new TH1D("h_theta_pmu_gene"       ,"", 100, -1.*TMath::Pi(), TMath::Pi());
  h_x_pmu_gene            = new TH1D("h_x_pmu_gene"           ,"", 100,-100, 0);
  h2_yz_pmu_gene          = new TH2D("h2_yz_pmu_gene"         ,"", 100, -50, 50, 100, -100, 0); 
  h2_yz_pmu_gene_wBeHit   = new TH2D("h2_yz_pmu_gene_wBeHit"  ,"", 100, -50, 50, 100, -100, 0); 
  set->SetTH1(h_cost_pmu_gene         ,"p#mu cos#theta"                ,"cos#theta"    ,"Counts/0.1"    ,1,3000,0);
  set->SetTH1(h_cost_pmu_gene_wBeHit  ,"p#mu cos#theta (w/ Be hit)"    ,"cos#theta"    ,"Counts/0.1"    ,2,3001,2);
  set->SetTH1(h_r_pmu_gene            ,"p#mu source point (R)"         ,"R [mm]"       ,"Counts"        ,1,3000,0);
  set->SetTH1(h_theta_pmu_gene        ,"p#mu source point (#theta)"    ,"#theta [rad]" ,"Counts"        ,1,3000,0);
  set->SetTH1(h_x_pmu_gene            ,"p#mu source point X"           ,"X [mm]" ,"Counts"        ,1,3000,0);
  set->SetTH2(h2_yz_pmu_gene          ,"p#mu source point"             ,"Y [mm]" ,"Z [mm]" );
  set->SetTH2(h2_yz_pmu_gene_wBeHit   ,"p#mu source point (w/ Be hit)" ,"Y [mm]" ,"Z [mm]" );

  h        = new TH1D("h","",100,-10, 10);
  h_cost_muBeXray_gene         = new TH1D("h_cost_muBeXray_gene"        ,"",200,-1, 1);
  h_cost_muBeXray_gene_wSDDall = new TH1D("h_cost_muBeXray_gene_wSDDall","",200,-1, 1);
  h_r_muBeXray_gene            = new TH1D("h_r_muBeXray_gene"    ,"", 100, 0, 50);
  h_theta_muBeXray_gene        = new TH1D("h_theta_muBeXray_gene","", 100, -1.*TMath::Pi(), TMath::Pi());
  h_x_muBeXray_gene            = new TH1D("h_x_muBeXray_gene"    ,"",  50, -50, 50);
  h_y_muBeXray_gene            = new TH1D("h_y_muBeXray_gene"    ,"",  50, -50, 50);
  h2_xy_muBeXray_gene              = new TH2D("h2_xy_muBeXray_gene"        ,"",       100, -50, 50, 100, -50, 50); 
  h2_xy_muBeXray_gene_wSDDall      = new TH2D("h2_xy_muBeXray_gene_wSDDall","",       100, -50, 50, 100, -50, 50); 
  set->SetTH1(h  ,"","M_{d#pi+} [GeV/#it{c}^{2}]" ,"Counts"    ,1,3001,2);
  set->SetTH1(h_cost_muBeXray_gene          ,"#muBe cos#theta"             ,"cos#theta"    ,"Counts/0.1"    ,1,3000,0);
  set->SetTH1(h_cost_muBeXray_gene_wSDDall  ,"#muBe cos#theta (w/ SDD hit)","cos#theta"    ,"Counts/0.1"    ,2,3001,2);
  set->SetTH1(h_r_muBeXray_gene             , "#muBe source point (R)"     ,"R [mm]"       ,"Counts"        ,1,3000,0);
  set->SetTH1(h_theta_muBeXray_gene         , "#muBe source point (#theta)","#theta [rad]" ,"Counts"        ,1,3000,0);
  set->SetTH1(h_x_muBeXray_gene             , "#muBe source point (X)"     ,"X [mm]"       ,"Counts"        ,1,3000,0);
  set->SetTH1(h_y_muBeXray_gene             , "#muBe source point (Y)"     ,"Y [mm]"       ,"Counts"        ,1,3000,0);
  set->SetTH2(h2_xy_muBeXray_gene            , "#muBe source point"             ,"X [mm]" ,"Y [mm]" );
  set->SetTH2(h2_xy_muBeXray_gene_wSDDall    , "#muBe source point (w/ SDD hit)","X [mm]" ,"Y [mm]" );
  h_cost_muBeXray_gene         ->SetStats(0);
  h_cost_muBeXray_gene_wSDDall ->SetStats(0);

  for(int n=0;n<NSDD;n++){
    h2_xy_muBeXray_wSDD[n]       = new TH2D(Form("h2_xy_muBeXray_wSDD%d" ,n+1),"",       100, -50, 50, 100, -50, 50); 
    h2_xy_muBeXray_atSDD[n]      = new TH2D(Form("h2_xy_muBeXray_atSDD%d",n+1),"",       100, -50, 50, 100, -50, 50); 
  }
}
////////////////////////////////////////////////////////////////////////////
void SDDAcptMC::InitSDDsetup(){

  cout<<"SDDx:"<<SDDx<<endl;
  pos_SDD_surf[0].SetXYZ(SDDx,0.,25.);//z=70: default
  //pos_SDD_surf[0].SetXYZ(0.,0.,25.);//z=70: default
  dia_SDD_surf[0]=20.;
  til_SDD_surf[0]=0.;

}
////////////////////////////////////////////////////////////////////////////
int SDDAcptMC::IsIntersectSDDSurface(TVector3 posX, TVector3 dirX, TVector3 &sectVec){

  // w/o tilt calc.
  for(int n=0;n<NSDD;n++){
    double a = pos_SDD_surf[n].Z()/dirX.Z();
    double x_at_SDDz=posX.X() + a*dirX.X();
    double y_at_SDDz=posX.Y() + a*dirX.Y();
    sectVec.SetXYZ(x_at_SDDz,y_at_SDDz,a*dirX.Z());
    if(( (x_at_SDDz-pos_SDD_surf[n].X())*(x_at_SDDz-pos_SDD_surf[n].X()) + (y_at_SDDz-pos_SDD_surf[n].Y())*(y_at_SDDz-pos_SDD_surf[n].Y()) ) < 0.25*dia_SDD_surf[n]*dia_SDD_surf[n]){return n;}
  }


  return -1;

}
////////////////////////////////////////////////////////////////////////////
bool SDDAcptMC::IsIntersectBeTarget(TVector3 posPmu, TVector3 dirPmu, TVector3 &sectVec){

  double a = -posPmu.Z()/dirPmu.Z();//Be target located at z=0
  double x_at_BeTarget=posPmu.X() + a*dirPmu.X();
  double y_at_BeTarget=posPmu.Y() + a*dirPmu.Y();
  sectVec.SetXYZ(x_at_BeTarget,y_at_BeTarget,0.);
  if(a<0.)return false;
  if((x_at_BeTarget*x_at_BeTarget + y_at_BeTarget*y_at_BeTarget) < 0.25*D_Be*D_Be){return true;}

  return false;

}
////////////////////////////////////////////////////////////////////////////
void SDDAcptMC::Loop(){
  cout<<"Loop"<<endl;

  BeHitCounter = SDDHitCounter = 0;

  for(int n=0;n<ENumMax;n++){
    if(n%10000==0) cout<<n<<" / "<<ENumMax<<endl;
    double r_muBeXray    ;
    double theta_muBeXray;
    double x_muBeXray    ;
    double y_muBeXray    ;

    //////////////////////
    //hot pmu simulation//
    //////////////////////
    
    //source of hot pmu (Solid H2 local coordinate)
    TVector3 pos_pmu, dir_pmu;//pmu position and direction at H2 target surface
    double r_pmu     = 0.;
    double theta_pmu = 0.;
    double x_pmu  = 0.;
    double y_pmu  = 0.;
    if(pmu_type==1){
      r_pmu     = GetUniformR(30.);
      theta_pmu = 2.*TMath::Pi()*(gRandom->Uniform()-0.5);
      x_pmu     = r_pmu*TMath::Cos(theta_pmu);
      y_pmu     = r_pmu*TMath::Sin(theta_pmu);
    }else if(pmu_type==2){
      while(1){
        x_pmu     = gRandom->Gaus(0.,20.);
        y_pmu     = gRandom->Gaus(0.,20.);
        r_pmu     = sqrt(x_pmu*x_pmu + y_pmu*y_pmu);
        theta_pmu = TMath::ATan2(x_pmu, y_pmu);
        if(r_pmu<30.)break;
      }
    }

    pos_pmu.SetXYZ(x_pmu, y_pmu, -20.);
    dir_pmu=GetUniformDistribution(0.);
    double cost_pmu =dir_pmu.Z();

    //convert to global coordinate
    ConvertToNKS2Coordinate(pos_pmu);
    pos_pmu.SetZ(pos_pmu.Z()-50.);//offset
    ConvertToNKS2Coordinate(dir_pmu);
    //cout<<"pos_pmu.X():"<<pos_pmu.X()<<endl;

    h_cost_pmu_gene          ->Fill(cost_pmu);
    h_r_pmu_gene             ->Fill(r_pmu);
    h_theta_pmu_gene         ->Fill(theta_pmu);
    h_x_pmu_gene             ->Fill(pos_pmu.X());
    h2_yz_pmu_gene           ->Fill(pos_pmu.Y(), pos_pmu.Z());
    if(IsIntersectBeTarget(pos_pmu, dir_pmu, pos_pmu_stop)){
      x_muBeXray     = pos_pmu_stop.X();
      y_muBeXray     = pos_pmu_stop.Y();
      r_muBeXray      = sqrt( x_muBeXray*x_muBeXray + y_muBeXray*y_muBeXray );
      theta_muBeXray  = TMath::ATan2(x_muBeXray, y_muBeXray);//[-pi, pi]

      h_cost_pmu_gene_wBeHit   ->Fill(cost_pmu);
      h2_yz_pmu_gene_wBeHit    ->Fill(pos_pmu.Y(), pos_pmu.Z());
 
      if(x_muBeXray<-20.)cout<<Form("pos_pmu.X(): %lf, pos_pmu.Z(): %lf, dir_pmu.X(): %lf, dir_pmu.Z(): %lf, x_muBeXray: %lf, cost_pmu: %lf"
                                    ,     pos_pmu.X(),      pos_pmu.Z(),      dir_pmu.X(),      dir_pmu.Z(),      x_muBeXray,      cost_pmu)<<endl;
      BeHitCounter++;
    }else {if(fix_NmuBeXray_flag)n--;continue;}

    ///////////////////
    //muBe simulation//
    ///////////////////
    
    //if(treeout_flag)InitializeTreeBr();
    //double x_muBeXray = D_Be*(gRandom->Uniform(1)-0.5);
    //double y_muBeXray = D_Be*(gRandom->Uniform(1)-0.5);
    //if( (x_muBeXray*x_muBeXray + y_muBeXray*y_muBeXray) > 0.25*D_Be*D_Be){n--;continue;}
    
    //r_muBeXray     = GetUniformR(0.5*D_Be);
    //theta_muBeXray = 2.*TMath::Pi()*(gRandom->Uniform()-0.5);
    //x_muBeXray     = r_muBeXray*TMath::Cos(theta_muBeXray);
    //y_muBeXray     = r_muBeXray*TMath::Sin(theta_muBeXray);

    h_r_muBeXray_gene     ->Fill(r_muBeXray);
    h_theta_muBeXray_gene ->Fill(theta_muBeXray);
    h2_xy_muBeXray_gene   ->Fill(x_muBeXray, y_muBeXray);
    h_x_muBeXray_gene     ->Fill(x_muBeXray);
    h_y_muBeXray_gene     ->Fill(y_muBeXray);
    //pos_pmu_stop.SetXYZ(x_muBeXray, y_muBeXray, 0.);
    dir_muBe_Xray=GetUniformDistribution(0.);
    h_cost_muBeXray_gene  ->Fill(dir_muBe_Xray.Z());
    //dir_muBe_Xray.SetXYZ(0.,0.,1.);
    int segSDD =IsIntersectSDDSurface( pos_pmu_stop, dir_muBe_Xray, pos_at_SDD);
    if(segSDD>=0){
      h2_xy_muBeXray_gene_wSDDall  ->Fill(x_muBeXray, y_muBeXray);
      h_cost_muBeXray_gene_wSDDall ->Fill(dir_muBe_Xray.Z());
      h2_xy_muBeXray_wSDD[segSDD]  ->Fill(x_muBeXray, y_muBeXray); 
      h2_xy_muBeXray_atSDD[segSDD] ->Fill(pos_at_SDD.X(), pos_at_SDD.Y()); 
      SDDHitCounter++;
    }
  }

  cout<<Form("%d %d %d",ENumMax, BeHitCounter, SDDHitCounter)<<endl;
  cout<<h2_xy_muBeXray_gene->Integral()<<endl;
}
////////////////////////////////////////////////////////////////////////////
void SDDAcptMC::SetBranch(){
  cout<<"make branch"<<endl;
  //ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  
  //ofp->cd();
  //tree_out = new TTree("tree","tree");
  //tree_out->Branch("eg"            , &tr.eg              , "eg/D"               );
  //treeout_flag = true;
}
////////////////////////////////////////////////
void SDDAcptMC::InitializeTreeBr(){
}
//_______________________________________________________________________________________________________________________________

TVector3 SDDAcptMC::GetUniformDistribution(double costMin, double costMax)
{

  // forward peak

  //  G4double uniform_tmp; 
  //  do {
  //    uniform_tmp = G4UniformRand();
  //  } while ( G4UniformRand() >= uniform_tmp) ;
  //  G4double cost = (costMax-costMin)*uniform_tmp + costMin;

  // backford peak

  //  G4double uniform_tmp; 
  //  do {
  //    uniform_tmp = G4UniformRand();
  //  } while( G4UniformRand() <= uniform_tmp );
  //  G4double cost = (costMax-costMin)*uniform_tmp + costMin;
  

  // uniform
  double cost = (costMax-costMin)*gRandom->Uniform(1) + costMin;
  double sint = sqrt(1.-cost*cost);
  double phi  = 2.*TMath::Pi()*(gRandom->Uniform(1)-0.5);
  
  //return TVector3(cost, sint*cos(phi), sint*sin(phi));//NKS2 coordinate (x is the beam axis)
  return TVector3(sint*cos(phi), sint*sin(phi),cost);//normal coordinate (z is the beam axis)
}
//_______________________________________________________________________________________________________________________________
void  SDDAcptMC::ConvertToNKS2Coordinate(TVector3 &dir)
{
  //beam axis z -> x
  //conversion x->ynew, y->znew, z->xnew,
  double xnew=dir.Z();
  double ynew=dir.X();
  double znew=dir.Y();

  dir.SetXYZ(xnew, ynew, znew);
}
//_______________________________________________________________________________________________________________________________

double SDDAcptMC::GetUniformR(double Rmax)
{
  double r, yy;
  
  while(1){
    r = gRandom->Uniform(Rmax);
    yy= gRandom->Uniform(1.);
    if(yy< r/Rmax)break;
  }
  
  return r;
}

//_______________________________________________________________________________________________________________________________

////////////////////////////////////////////////////////////////////////////
bool SDDAcptMC::Draw(){
  int ic=0;
  c[ic]->Clear();c[ic]->Divide(2,2);
  c[ic]->cd(1);gPad->SetLogz(0);  h_cost_pmu_gene          ->Draw("");
  c[ic]->cd(2);gPad->SetLogz(0);  h_r_pmu_gene             ->Draw("");
  c[ic]->cd(3);gPad->SetLogz(0);  h_theta_pmu_gene         ->Draw("");
  c[ic]->cd(4);gPad->SetLogz(0);  h2_yz_pmu_gene           ->Draw("colz");

  ic++;
  if(ic>=NCanvas)return false;
  c[ic]->Clear();c[ic]->Divide(2,2);
  c[ic]->cd(1);gPad->SetLogy(0);h_cost_pmu_gene          ->Draw("");h_cost_pmu_gene_wBeHit          ->Draw("same");
  c[ic]->cd(2);gPad->SetLogy(0);h2_yz_pmu_gene           ->Draw("colz");
  c[ic]->cd(3);gPad->SetLogy(0);h2_yz_pmu_gene_wBeHit    ->Draw("colz");
  c[ic]->cd(4);gPad->SetLogz(0);h_cost_pmu_gene_wBeHit          ->Draw("");//h_x_pmu_gene             ->Draw("");

  ic++;
  if(ic>=NCanvas)return false;
  h2_xy_muBeXray_gene_wSDDall->RebinX(2);
  h2_xy_muBeXray_gene_wSDDall->RebinY(2);
  c[ic]->Clear();c[ic]->Divide(2,2);
  c[ic]->cd(1);gPad->SetLogz(0); h2_xy_muBeXray_gene          ->Draw("colz");
  c[ic]->cd(2);gPad->SetLogz(0); h2_xy_muBeXray_gene_wSDDall  ->Draw("colz");
  c[ic]->cd(3);gPad->SetLogz(0); h2_xy_muBeXray_wSDD[0]       ->Draw("colz");
  c[ic]->cd(4);gPad->SetLogz(0); h2_xy_muBeXray_atSDD[0]      ->Draw("colz");

  ic++;
  if(ic>=NCanvas)return false;
  c[ic]->Clear();c[ic]->Divide(2,2);
  c[ic]->cd(1);gPad->SetLogy(0);h_cost_muBeXray_gene->Draw();
  c[ic]->cd(2);gPad->SetLogy(0);h_cost_muBeXray_gene->Draw();h_cost_muBeXray_gene_wSDDall->Draw("same");
  c[ic]->cd(3);gPad->SetLogy(0);h_cost_muBeXray_gene_wSDDall->Draw();

  ic++;
  if(ic>=NCanvas)return false;
  c[ic]->Clear();c[ic]->Divide(2,2);
  c[ic]->cd(1);gPad->SetLogz(0); h2_xy_muBeXray_gene        ->Draw("colz");
  c[ic]->cd(2);gPad->SetLogz(0); h_y_muBeXray_gene          ->Draw("");
  c[ic]->cd(3);gPad->SetLogz(0); h_x_muBeXray_gene          ->Draw("");
  c[ic]->cd(4);gPad->SetLogz(0); h_r_muBeXray_gene          ->Draw("");
  //c[ic]->cd(3);gPad->SetLogz(0); h_theta_muBeXray_gene      ->Draw("");

  return true;
}
////////////////////////////////////////////////////////////////////////////
void SDDAcptMC::SaveCanvas(string PdfFileName){
  c[0]->Print(Form("%s[",PdfFileName.c_str() ));
  for(int i=0;i<NCanvas;i++)
    c[i]->Print(Form("%s" ,PdfFileName.c_str() ));
  c[NCanvas-1]->Print(Form("%s]",PdfFileName.c_str() ));
  cout<<"saved : "<<Form("%s",PdfFileName.c_str() )<<endl;
}
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ofname = "root/tmp.root";
  string pdfname = "pdf/hoge.pdf";
  string ofname_root;
  //string ofname_pdf = "tmp.pdf";
  bool draw_flag = true;
  double Dinput=80.;
  int ch;
  int MaxNum   = 1000;
  int PmuType = 0;
  int sdd_x = 0;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hw:n:bp:d:t:x:"))!=-1){
    switch(ch){
    case 'w':
      ofname = optarg;
      cout<<"output filename : "<<ofname<<endl;
      break;
    case 'd':
      Dinput = atof(optarg);
      break;
    case 'n':
      MaxNum = atoi(optarg);
      break;
    case 't':
      PmuType = atoi(optarg);
      break;
    case 'x':
      sdd_x = atoi(optarg);
      break;
    case 'b':
      draw_flag=false;
      cout<<"BACH MODE!"<<endl;
      break;
    case 'p':
      pdfname = optarg;
      break;
    case 'h':
      cout<<"-w : output root filename"<<endl;
      cout<<"-p : output pdf filename" <<endl;
      cout<<"-b : batch mode"           <<endl;
      cout<<"-d : diameter of Be target (phi [mm])"     <<endl;
      cout<<"-t : pmu source poit type 0(default): from center only, 1: phi60 mm uniform, 2: gauss (sigma=20mm)"           <<endl;
      cout<<"-n : maximum number of events to be analysed "<<endl;
      return 0;
      break;
    case '?':
      cout<<"unknown option...."<<endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;
    }
  }

  ofname_root = ofname;


  TApplication *theApp = new TApplication("App", &argc, argv);
  SDDAcptMC *calc = new SDDAcptMC();
  calc->SetMaxEvent(MaxNum);
  calc->SetDiameter(Dinput);
  calc->SetPmuDistType(PmuType);
  calc->SetNBeXray(true);//
  calc->SetSDDx((double)sdd_x);//
  calc->InitSDDsetup();
  calc->MakeHist(ofname_root);
  //calc->SetBranch();
  calc->Loop();
  if(draw_flag){
    if(!calc->Draw())cout<<"Number of Canvas is not enough!"<<endl;
    calc->SaveCanvas(pdfname);
  }
  delete calc;
  cout<<"delete calc"<<endl;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}
