#include "Setting.h"
//____________________________________________________________________________________________
Setting::Setting()
{
}
//____________________________________________________________________________________________
Setting::~Setting()
{
std::cout<<"Setting::deconstructor called"<<std::endl;
}
//____________________________________________________________________________________________
void Setting::SetTH1(TH1 *h, TString hname, TString xname, TString yname, int LColor, int FStyle, int FColor){
  h->SetTitle(hname);
  h->SetLineColor(LColor);
  h->SetLineWidth(1);
  h->SetFillStyle(FStyle);
  h->SetFillColor(FColor);
//  h->SetMinimum(0.8);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(3);
}
//____________________________________________________________________________________________
void Setting::SetTH2(TH2 *h, TString name, TString xname, TString yname, double min, double MStyle, double MSize){
  h->SetTitle(name);
  h->SetMinimum(min);
  h->SetLineWidth(1);
  h->SetTitleSize(0.05,"");
  h->SetMarkerStyle(1);
  h->SetMarkerSize(0.1);
  h->SetMarkerColor(1);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(3);

  h->SetStats(0);
}

//____________________________________________________________________________________________
void Setting::SetTH3(TH3 *h, TString name, TString xname, TString yname, TString zname, double min, double MStyle, double MSize){
  h->SetTitle(name);
  h->SetMinimum(min);
  h->SetLineWidth(1);
  h->SetTitleSize(0.05,"");
  h->SetMarkerStyle(1);
  h->SetMarkerSize(0.1);
  h->SetMarkerColor(1);

  h->GetXaxis()->SetTitle(xname);
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelOffset(0.01);

  h->GetYaxis()->SetTitle(yname);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleOffset(0.7);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelOffset(0.01);
  ((TGaxis*)h->GetYaxis())->SetMaxDigits(3);

  h->SetStats(0);
}

//____________________________________________________________________________________________
void Setting::SetGr(TGraph *gr, TString hname, TString xname, TString yname, int LColor, int MColor, int MStyle, double MSize, double Yoffset){
  gr->SetTitle(hname);
  gr->SetName(hname);
  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleOffset(1.0);
  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->SetTitleOffset(1.0);
  gr->GetYaxis()->CenterTitle();
  gr->SetLineColor(LColor);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
  gr->GetYaxis()->SetTitleOffset(Yoffset);
//  gr->GetYaxis()->SetRangeUser(min,max);
}
//____________________________________________________________________________________________
void Setting::SetGrErr(TGraphErrors *gr, TString hname, TString xname, TString yname, int LColor, int MColor, int MStyle, double MSize, double Yoffset, double min, double max){
  gr->SetTitle(hname);
  gr->SetName(hname);
  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetTitleOffset(1.0);
  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->SetTitleOffset(1.0);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->SetLineColor(LColor);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
  gr->GetYaxis()->SetTitleOffset(Yoffset);
//  gr->GetYaxis()->SetRangeUser(min,max);
}
//____________________________________________________________________________________________
void Setting::SetGrErr(TGraphAsymmErrors *gr, TString hname, TString xname, TString yname, int LColor, int MColor, int MStyle, double MSize, double Yoffset, double min, double max){
  gr->SetTitle(hname);
  gr->SetName(hname);
  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetTitleOffset(1.0);
  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->SetTitleOffset(1.0);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->SetLineColor(LColor);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(MSize);
  gr->GetYaxis()->SetTitleOffset(Yoffset);
//  gr->GetYaxis()->SetRangeUser(min,max);
}
//____________________________________________________________________________________________
void Setting::SetTF1(TF1 *f, int LColor, int LStyle,double LWidth){
  f->SetLineColor(LColor);
  f->SetLineStyle(LStyle);
  f->SetLineWidth(LWidth);
  f->SetNpx(2000);
}
//____________________________________________________________________________________________
void Setting::SetTLatex(TLatex *latex, int TColor, double TSize,int Align){
 latex -> SetTextSize(TSize); 
 latex -> SetTextColor(TColor);
 latex -> SetTextAlign(Align);
 latex -> SetTextFont(22);
} 
//____________________________________________________________________________________________
void Setting::SetTLine(TLine *line, int LColor,int LStyle, double LWidth){
  line -> SetLineColor(LColor);
  line -> SetLineStyle(LStyle);
  line -> SetLineWidth(LWidth);

}
//____________________________________________________________________________________________
void Setting::SetTArrow(TArrow *arrow,  int LColor,int LStyle, double LWidth, float AAngle, float ASize, int FStyle, int FColor){
  arrow -> SetLineColor(LColor);
  arrow -> SetLineStyle(LStyle);
  arrow -> SetLineWidth(LWidth);
  arrow -> SetAngle(AAngle);
  arrow -> SetArrowSize(ASize);
  arrow -> SetFillStyle(FStyle);
  arrow -> SetFillColor(FColor);
     
}
//____________________________________________________________________________________________
void Setting::SetTLegend(TLegend *legend, int TFont, double TSize, int BSize, int FColor, int FStyle){
  legend -> SetBorderSize(BSize);
  legend -> SetFillColor(FColor);
  legend -> SetFillStyle(FStyle);
  legend -> SetTextFont(TFont);
  legend -> SetTextSize(TSize);
  //Font number         X11 Names             Win32/TTF Names
  //    1 :       times-medium-i-normal      "Times New Roman"
  //    2 :       times-bold-r-normal        "Times New Roman"
  //    3 :       times-bold-i-normal        "Times New Roman"
  //    4 :       helvetica-medium-r-normal  "Arial"
  //    5 :       helvetica-medium-o-normal  "Arial"
  //    6 :       helvetica-bold-r-normal    "Arial"
  //    7 :       helvetica-bold-o-normal    "Arial"
  //    8 :       courier-medium-r-normal    "Courier New"
  //    9 :       courier-medium-o-normal    "Courier New"
  //   10 :       courier-bold-r-normal      "Courier New"
  //   11 :       courier-bold-o-normal      "Courier New"
  //   12 :       symbol-medium-r-normal     "Symbol"
  //   13 :       times-medium-r-normal      "Times New Roman"
  //   14 :                                  "Wingdings"
  //   15 :       Symbol italic (derived from Symbol)
}
//____________________________________________________________________________________________
void Setting::FitGaus(TH1 *h, double &gamin, double &gamax, double range, int itr){
  //gamin = h  ->GetBinCenter(h ->GetMaximumBin())-20.0;
  //gamax = h  ->GetBinCenter(h ->GetMaximumBin())+20.0;
  for(Int_t l=0; l<itr; l++){
    TF1 *ga = new TF1("ga","gaus");
    h  ->Fit(ga,"0QR","",gamin,gamax);
    gamin = ga->GetParameter(1) - ga->GetParameter(2)*range;
    gamax = ga->GetParameter(1) + ga->GetParameter(2)*range;
    ga->Clear();
  }
}
