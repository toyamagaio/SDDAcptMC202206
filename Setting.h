#ifndef Setting_h
#define Setting_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>

#include "TApplication.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TArrow.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TLegend.h"

class Setting  
{
public:
  Setting();
  ~Setting();
  //virtual ~Setting(){};
  void SetTH1(TH1 *h, TString hname, TString xname, TString yname, int LColor=1, int FStyle=0, int FColor=0);
  void SetTH2(TH2 *h, TString name, TString xname, TString yname, double min=0.8, double MStyle=1, double MSize=1.0);
  void SetTH3(TH3 *h, TString name, TString xname, TString yname, TString zname, double min=0.8, double MStyle=1, double MSize=1.0);
  void SetGr(TGraph *gr, TString hname, TString xname, TString yname, int LColor=1, int MColor=2, int MStyle=21, double MSize=0.8,double Yoffset=0.8);
  void SetGrErr(TGraphErrors *gr, TString hname, TString xname, TString yname, int LColor, int MColor, int MStyle, double MSize = 0.8, double Yoffset=0, double min=0, double max=0);
  void SetGrErr(TGraphAsymmErrors *gr, TString hname, TString xname, TString yname, int LColor, int MColor, int MStyle, double MSize = 0.8, double Yoffset=0, double min=0, double max=0);
  void SetTF1(TF1 *f, int LColor, int LStyle,double LWidth);
  void SetTLatex(TLatex *latex, int TColor, double TSize,int Align);
  void SetTLine(TLine *line, int LColor=2,int LStyle=1, double LWidth=1);
  void SetTArrow(TArrow *arrow, int LColor=2,int LStyle=1, double LWidth=1, float AAngle=45, float ASize=0.01, int FStyle=1001, int FColor=1);
  void SetTLegend(TLegend *legend, int TFont=22, double TSize=0.025, int BSize=1, int FColor=0, int FStyle=1001);

  void FitGaus(TH1 *h, double &gamin, double &gamax, double range=3.0, int itr=10);

};

#endif
