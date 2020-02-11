#include "../include/CentralInclude.h"
#include "TRegexp.h"
#include <sys/stat.h>

using namespace std;

void drawCategory(TString region, TString ptBin);
void setStyle();
TString OutDir = "../Plots/WMass/";
bool scaleQCD=true;
//Rebin mSD
bool RebinMSD=true;
double minMSD=50;
double maxMSD=170;
int N_MSD = 20;
double * new_binning;

// double YMin=0.9*pow(10,1);
// double YMax=pow(10,6);
double YMin=-1;
double YMax=-1;

double xLabelSize=18.0;
double yLabelSize=18.0;
double xTitleSize=20.0;
double yTitleSize=22.0;
double xTitleOffset=4.0;
double yTitleOffset=1.3;


TString InPath=TString(getenv("CMSSW_BASE"))+"/src/UHH2/JetMass/Histograms/W/";
vector<TFile*> BKG_Files={
  new TFile(InPath+"QCD.root","READ"),
  new TFile(InPath+"WMatched.root","READ"),
  new TFile(InPath+"WUnmatched.root","READ")
};
vector<TFile*> SIG_Files={};
vector<TFile*> DAT_Files={new TFile(InPath+"Data.root","READ")};

vector<int> BKG_Colors={867,kOrange-3,kOrange+3};
vector<int> SIG_Colors={};

vector<TString> BKG_Legend={"QCD","W+jets to qq (merged)","W+jets to qq (unmerged)"};
vector<TString> SIG_Legend={};
TString DAT_Legend="Data JetHT 2016";



int main(int argc, char* argv[]){	
  setStyle();
  struct stat buffer;
  if(stat (((std::string) OutDir).c_str(), &buffer) != 0){
    mode_t standard_mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IXGRP | S_IROTH| S_IXOTH;
    mkdir(((std::string) OutDir).c_str(),standard_mode);  
    mkdir(((std::string) (OutDir+"/logY")).c_str(),standard_mode);  
  } 
  double binWidth = (maxMSD-minMSD)/N_MSD;
  vector<double> bins;
  for(int i = 0 ; i<N_MSD ; ++i)bins.push_back(minMSD+i*binWidth);	
  bins.push_back(maxMSD);
  new_binning = &bins[0];
	

  cout << "Plotting Templates from: " << InPath <<endl;
	

  vector<TString> regions = {"_pass","_fail"};
  vector<TString> ptBins = {"","_pt500To550","_pt550To600","_pt600To675","_pt675To800","_pt800To1200"};
  // vector<TString> regions = {"_pass"};
  // vector<TString> ptBins = {""};

  for(auto const region : regions){
    for(auto const ptBin : ptBins){
      drawCategory(region,ptBin);
    }
  }


}


void drawCategory(TString region="pass", TString ptBin=""){
  TString Selection = "N2DDT";
  // TString HistDir = "JetMass_"+Selection+ptBin+"_noJecOnMass"+region;
  TString HistDir = "JetMass_"+Selection+ptBin+region;
  TString HistName = "Mass_central";

  TLegend * legend = new TLegend(0.6,0.70,0.85,0.9);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.03);
  legend->SetMargin(0.4);
  // legend->SetNColumns(2);
  legend->SetColumnSeparation(0.3);
		
  TH1F * dat_Hist = (TH1F*) DAT_Files[0]->Get(HistDir+"/"+HistName);
  for(int i=1;i<DAT_Files.size();++i)dat_Hist->Add((TH1F*) DAT_Files[1]->Get(HistDir+"/"+HistName));

  if(RebinMSD)dat_Hist = (TH1F*) dat_Hist->Rebin(N_MSD,"",new_binning);

  dat_Hist->SetMarkerStyle(8);
  dat_Hist->SetMarkerColor(kBlack);
  dat_Hist->SetLineColor(kBlack);
  legend->AddEntry(dat_Hist,DAT_Legend,"p");

  TH1F * bg_Err_unscaled = (TH1F*) BKG_Files[0]->Get(HistDir+"/"+HistName)->Clone();
  for(int i=1;i<BKG_Files.size();++i)bg_Err_unscaled->Add((TH1F*) BKG_Files[i]->Get(HistDir+"/"+HistName));

  TH1F * qcd_from_data_hist = (TH1F*) dat_Hist->Clone();
  TH1F* qcd_hist = (TH1F*) dat_Hist->Clone();
  qcd_hist->Reset();
  for(int i = 0; i<BKG_Files.size() ; ++i){
    TH1F * bg_Hist =  (TH1F*) BKG_Files[i]->Get(HistDir+"/"+HistName);
    bg_Hist = (TH1F*) bg_Hist->Rebin(N_MSD,"",new_binning);
    if(TString(BKG_Files[i]->GetName()).Contains("QCD")){
      qcd_hist=(TH1F*)bg_Hist->Clone();
      // continue;
    }else{
      qcd_from_data_hist->Add(bg_Hist,-1);
    }
  }
        
  // double QCD_Scale = dat_Hist->Integral()/bg_Err_unscaled->Integral();
  std::cout << "data-w norm: "<<qcd_from_data_hist->Integral() << " qcd norm: " << qcd_hist->Integral() << std::endl;
  double QCD_Scale = qcd_from_data_hist->Integral()/qcd_hist->Integral();
  if(!scaleQCD)QCD_Scale=1.0;

  THStack * bg_Stack = new THStack ("bg_stack","bg_stack");
  TH1F * bg_Err = (TH1F*) BKG_Files[0]->Get(HistDir+"/"+HistName)->Clone();
  if(RebinMSD)bg_Err = (TH1F*) bg_Err->Rebin(N_MSD,"",new_binning);
  bg_Err->Reset();
  for(int i = 0; i<BKG_Files.size() ; ++i){
    TH1F * bg_Hist =  (TH1F*) BKG_Files[i]->Get(HistDir+"/"+HistName);		
    if(RebinMSD)bg_Hist = (TH1F*) bg_Hist->Rebin(N_MSD,"",new_binning);
    if(TString(BKG_Files[i]->GetName()).Contains("QCD")) bg_Hist->Scale(QCD_Scale);
    bg_Hist->SetFillColor(BKG_Colors[i]);
    bg_Hist->SetLineColor(BKG_Colors[i]);
    bg_Stack->Add(bg_Hist);
    legend->AddEntry(bg_Hist,BKG_Legend[i],"f");
    bg_Err->Add(bg_Hist);
  }

			
  bg_Err->SetFillStyle(3204);
  bg_Err->SetFillColor(kGray+2);
  bg_Err->SetLineColor(1);
  legend->AddEntry(bg_Err,"MC stat. Unc.","f");

  bg_Err->GetYaxis()->SetTitle("Events");
  bg_Err->GetYaxis()->SetTitleFont(43);
  bg_Err->GetYaxis()->SetTitleSize(yTitleSize);
  bg_Err->GetYaxis()->SetTitleOffset(yTitleOffset);
  bg_Err->GetYaxis()->SetLabelFont(43);
  bg_Err->GetYaxis()->SetLabelSize(yLabelSize);
  bg_Err->GetXaxis()->SetTitleSize(0.0);
  bg_Err->GetXaxis()->SetLabelSize(0.0);


			
      
  TH1F * ratioHist = (TH1F*) dat_Hist->Clone();
  ratioHist->SetLineColor(kBlack);
  ratioHist->SetStats(0);
  ratioHist->Divide(bg_Err);
  ratioHist->SetMarkerStyle(21);
  ratioHist->SetMarkerSize(0.7);

  ratioHist->GetYaxis()->SetRangeUser(0.45,1.55);
  ratioHist->GetYaxis()->SetTitle("Data/BG");
  ratioHist->GetYaxis()->CenterTitle();
  ratioHist->GetYaxis()->SetTitleFont(43);
  ratioHist->GetYaxis()->SetTitleSize(yTitleSize);
  ratioHist->GetYaxis()->SetTitleOffset(yTitleOffset);
  ratioHist->GetYaxis()->SetLabelFont(43);
  ratioHist->GetYaxis()->SetLabelSize(yLabelSize);
  ratioHist->GetYaxis()->SetNdivisions(506);

  ratioHist->GetXaxis()->SetTitle("m_{SD} [GeV]");
  ratioHist->GetXaxis()->SetTitleFont(43);
  ratioHist->GetXaxis()->SetTitleSize(xTitleSize);
  ratioHist->GetXaxis()->SetTitleOffset(xTitleOffset);
  ratioHist->GetXaxis()->SetLabelFont(43);
  ratioHist->GetXaxis()->SetLabelSize(xLabelSize);
  ratioHist->GetXaxis()->SetTickLength(0.08);
  ratioHist->GetXaxis()->SetNdivisions(506);


  double ratioXMin = ratioHist->GetXaxis()->GetXmin();
  double ratioXMax = ratioHist->GetXaxis()->GetXmax();

  TLine * zeropercent= new TLine(ratioXMin,1,ratioXMax,1);
  TLine * plus10percent= new TLine(ratioXMin,1.1,ratioXMax,1.1);
  TLine * minus10percent= new TLine(ratioXMin,0.9,ratioXMax,0.9);
  plus10percent->SetLineStyle(kDashed);
  minus10percent->SetLineStyle(kDashed);
      
  TCanvas * c1 = new TCanvas("c","c",600,600);
  double yplot=0.7;
  double yratio=0.3;
  double ymax=1.0;
  double xmax=1.0;
  double xmin=0.0;

  int MaxMagnitude=TMath::FloorNint(TMath::Log10(bg_Err->GetMaximum()));
  if(YMin<0){
    YMin=(float) pow(10,MaxMagnitude-3);
    YMin+=(float) pow(10,MaxMagnitude-4);
  }
  if(YMax<0)YMax=(float) pow(10,MaxMagnitude+3);
  // (float) pow(10,MaxMagnitude-5);
  double YMax_lin=bg_Err->GetMaximum()+1*pow(10,MaxMagnitude);
  TPad * plotpad = new TPad("plotpad","Plot",xmin,ymax-yplot,xmax,ymax);
  TPad * ratiopad = new TPad("ratiopad","Ratio",xmin,ymax-yplot-yratio,xmax,ymax-yplot);

  // gROOT->ProcessLine("CMS_lumi(c1);");

  plotpad->SetTopMargin(0.08);
  plotpad->SetLeftMargin(0.1);
  plotpad->SetRightMargin(0.05);
  plotpad->SetTicks();
  plotpad->Draw();

  plotpad->SetBottomMargin(0.016);
  ratiopad->SetTopMargin(0.016);
  ratiopad->SetBottomMargin(0.35);
  ratiopad->SetLeftMargin(0.1);
  ratiopad->SetRightMargin(0.05);
  ratiopad->SetTicks();
  ratiopad->Draw();

  plotpad->cd();
  bg_Err->GetYaxis()->SetRangeUser(0,YMax_lin);
  bg_Err->Draw("E2");
  bg_Stack->Draw("HistSame");
  bg_Err->Draw("E2Same");
  dat_Hist->Draw("PSame");
  plotpad->RedrawAxis();
  TLatex latex = TLatex();
  latex.SetNDC(kTRUE);
  latex.SetTextSize(20);
  latex.DrawLatex(0.69,0.93,TString::Format("%.2f fb^{-1} (13 TeV)",41.86));
  latex.DrawLatex(0.25,0.93,TString::Format("%s N2^{DDT} Selection", region.Contains("p") ? "pass" : "fail" ));

  cout<< "ptBin: " << ptBin <<endl;
  if(ptBin != ""){		
    TString pT_TString_lower=ptBin(TRegexp("[0-9]+To"));
    string pT_string_lower=(string) pT_TString_lower.ReplaceAll("To","");
    double pT_lower = stod(pT_string_lower);
    TString pT_TString_upper=ptBin(TRegexp("To[0-9]+"));
    string pT_string_upper=(string) pT_TString_upper.ReplaceAll("To","");
    double pT_upper = stod(pT_string_upper);
    latex.SetTextSize(25);
    TString ptBin_tex=TString::Format("%.0f GeV #leq p_{T} < %.0f GeV",pT_lower,pT_upper);
    latex.DrawLatex(0.20,0.80,ptBin_tex);
  }


  // latex.DrawLatex(0.1,0.953,"private work");
  legend->Draw("SAME");

	
  ratiopad->cd();
  ratioHist->Draw("EP");
  zeropercent->Draw();
  plus10percent->Draw();
  minus10percent->Draw();

  TString OutName=OutDir+HistName+ptBin+region;				
  c1->SaveAs(OutName+".pdf");

  bg_Err->GetYaxis()->SetRangeUser(YMin,YMax);
  c1->RedrawAxis();
  plotpad->RedrawAxis();
  plotpad->SetLogy();
  c1->SetLogy();
  plotpad->RedrawAxis();
  OutName=OutDir+"logY/"+HistName+ptBin+region;
  c1->SaveAs(OutName+"_logY.pdf");

  delete c1,bg_Stack;
}



void setStyle(){
  //Set TDR styles
  // gROOT->LoadMacro("tdrstyle.C");
  // gROOT->ProcessLine("setTDRStyle();");
	
  // Add CMS text
  // gROOT->LoadMacro("CMS_lumi.C");


  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);

  gStyle->SetTextFont(43);
  //
  gStyle->SetTitleOffset(0.86,"X");
  gStyle->SetTitleOffset(1.6,"Y");
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetMarkerSize(0.5);
  gStyle->SetHistLineWidth(1);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");
  gStyle->SetNdivisions(506, "XYZ");
  gStyle->SetLegendBorderSize(0);
}
