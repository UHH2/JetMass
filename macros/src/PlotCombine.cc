#include "../include/CentralInclude.h"


using namespace std;
using namespace RooFit;

void SetPlotStyle();
// TH1F* GetRatio(TH1F* num, TH1F* den);
void PlotHist(TGraphAsymmErrors* h_dat, TH1F* h_sig, TH1F* h_bkg, TString name);

int main(int argc, char* argv[]){

  TString filename;

  if(argc != 2){
    cout << "specify file!" << endl;
    return 1;
  }
  else filename = argv[1];
  TFile* file = new TFile(filename);





  vector<TString> directories = {"fit_s", "prefit"};
  // vector<TString> channels = {"TopMass200", "TopMass300", "TopMass400"};
  vector<TString> channels = {"TopMass200", "TopMass300", "TopMass400", "WMassPt500pass","WMassPt550pass","WMassPt675pass","WMassPt800pass","WMassPt500fail","WMassPt550fail","WMassPt675fail","WMassPt800fail"};

  for(auto dir: directories){
    for(auto ch: channels){
      TString postfit_dir = "shapes_"+dir+"/"+ch+"/";

      TGraphAsymmErrors* h_dat = (TGraphAsymmErrors*) file->Get(postfit_dir+"data");
      TH1F* h_sig = (TH1F*) file->Get(postfit_dir+"total_signal");
      TH1F* h_bkg = (TH1F*) file->Get(postfit_dir+"total_background");

      TH1F* h_mc = (TH1F*) h_sig->Clone();
      h_mc->Add(h_bkg);

      SetPlotStyle();
      TString plotname = ch+"_"+dir;
      PlotHist(h_dat, h_mc, h_bkg, plotname);
    }
  }
  return 0;
}

void SetPlotStyle(){
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);
  return;
}

// TH1F* GetRatio(TH1F* num, TH1F* den){
//   TH1F* ratio = (TH1F*) num->Clone();
//   ratio->Reset();
//   int Nbins = num->GetXaxis()->GetNbins();
//   for(int bin=1; bin<=Nbins; bin++){
//     double r = 0.0;
//     double e = 0.0;
//     if(den->GetBinContent(bin) != 0.0) r = num->GetBinContent(bin) / den->GetBinContent(bin);
//     if(den->GetBinContent(bin) != 0.0 && num->GetBinContent(bin) != 0.0){
//       e = r * sqrt(pow(num->GetBinError(bin)/num->GetBinContent(bin),2) + pow(den->GetBinError(bin)/den->GetBinContent(bin),2));
//     }
//     ratio->SetBinContent(bin, r);
//     ratio->SetBinError(bin, 2);
//   }
//   return ratio;
// }

void PlotHist(TGraphAsymmErrors* h_dat, TH1F* h_mc, TH1F* h_bkg, TString name){
  TCanvas *c = new TCanvas("c", "c", 600, 600);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.14);
  h_mc->SetTitle(" ");
  h_mc->GetXaxis()->SetTitle("jet mass");
  h_mc->GetYaxis()->SetTitle("Events");
  h_mc->GetYaxis()->SetTitleOffset(1.6);
  h_mc->GetXaxis()->SetTitleOffset(1.1);
  h_mc->GetYaxis()->SetTitleSize(0.05);
  h_mc->GetXaxis()->SetTitleSize(0.05);
  h_mc->GetYaxis()->SetLabelSize(0.05);
  h_mc->GetXaxis()->SetLabelSize(0.05);
  h_mc->GetXaxis()->SetNdivisions(505);
  h_mc->GetYaxis()->SetNdivisions(505);
  h_mc->GetYaxis()->SetRangeUser(0, h_mc->GetMaximum()*1.3);
  h_mc->SetFillColor(kGreen+2);
  h_mc->SetLineColor(kGreen+2);
  if(name.Contains("TopMass")){
    h_mc->SetFillColor(kRed-4);
    h_mc->SetLineColor(kRed-4);
  }
  h_bkg->SetLineColor(kAzure-9);
  h_bkg->SetFillColor(kAzure-9);
  if(name.Contains("TopMass")){
    h_bkg->SetLineColor(13);
    h_bkg->SetFillColor(13);
  }
  h_dat->SetMarkerColor(kBlack);
  h_dat->SetMarkerStyle(8);
  h_dat->SetMarkerSize(1);
  h_mc->Draw("HIST");
  h_bkg->Draw("HIST SAME");
  h_dat->Draw("P SAME");
  gPad->RedrawAxis();
  TLegend *leg = new TLegend(0.61,0.67,0.83,0.87);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h_dat, "pseudo data","pe");
  if(name.Contains("TopMass")) leg->AddEntry(h_mc, "t#bar{t}","f");
  else                         leg->AddEntry(h_mc, "W #rightarrow qq","f");
  leg->AddEntry(h_bkg, "background","f");
  leg->SetTextSize(0.04);
  leg->Draw();
  c->SaveAs("../Plots/"+name+".pdf");
  delete c;
  return;
}
