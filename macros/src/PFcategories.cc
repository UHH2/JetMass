#include "../include/CentralInclude.h"


using namespace std;
void SetPlotStyle();
vector<TH1F*> NormalisePerBin(vector<TH1F*> h_eta_Eweight);

int main(int argc, char* argv[]){

  TFile* file = new TFile("/nfs/dust/cms/user/schwarzd/JetMass/PostSel/uhh2.AnalysisModuleRunner.MC.PseudoData.root");
  vector<TString> categories = {"other", "gamma", "neutralH", "chargedH"};
  vector<TH1F*> h_eta, h_pt, h_eta_Eweight;
  for(auto cat: categories){
    h_eta.push_back((TH1F*) file->Get("PFHists/eta_"+cat));
    h_pt.push_back((TH1F*) file->Get("PFHists/pt_"+cat));
    h_eta_Eweight.push_back((TH1F*) file->Get("PFHists/eta_Eweight_"+cat));
  }

  vector<TH1F*> h_eta_Eweight_norm = NormalisePerBin(h_eta_Eweight);

  THStack *s_eta = new THStack();
  THStack *s_eta_norm = new THStack();
  THStack *s_pt = new THStack();
  Color_t col[] = {13, 798, kAzure+7, kRed-4};
  for(unsigned int i=0; i<categories.size(); i++){
    h_eta[i]->SetFillColor(col[i]);
    s_eta->Add(h_eta[i]);
    h_pt[i]->SetFillColor(col[i]);
    s_pt->Add(h_pt[i]);
    h_eta_Eweight_norm[i]->SetFillColor(col[i]);
    s_eta_norm->Add(h_eta_Eweight_norm[i]);
  }

  SetPlotStyle();
  TCanvas *c_pt = new TCanvas("c_pt","c_pt",600, 600);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.14);
  s_pt->Draw("HIST");
  s_pt->GetXaxis()->SetTitle("p_{T}");
  s_pt->GetYaxis()->SetTitle("PF particles");
  s_pt->GetYaxis()->SetTitleOffset(1.6);
  s_pt->GetXaxis()->SetTitleOffset(1.1);
  s_pt->GetYaxis()->SetTitleSize(0.05);
  s_pt->GetXaxis()->SetTitleSize(0.05);
  s_pt->GetYaxis()->SetLabelSize(0.05);
  s_pt->GetXaxis()->SetLabelSize(0.05);
  s_pt->GetYaxis()->SetNdivisions(505);
  s_pt->GetXaxis()->SetRangeUser(0, 20);
  s_pt->SetMaximum(s_pt->GetMaximum()*1.3);
  TLegend* l_pt = new TLegend(0.64,0.65,0.86,0.85);
  l_pt->SetBorderSize(0);
  l_pt->SetFillStyle(0);
  for(int i=categories.size()-1; i>=0; i--){
    l_pt->AddEntry(h_pt[i], categories[i],"f");
  }
  l_pt->SetTextSize(0.04);
  l_pt->Draw();
  c_pt->SaveAs("../Plots/PFcategories_pt.pdf");

  TCanvas *c_eta = new TCanvas("c_eta","c_eta",600, 600);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.14);
  s_eta->Draw("HIST");
  s_eta->GetXaxis()->SetTitle("#eta");
  s_eta->GetYaxis()->SetTitle("PF particles");
  s_eta->GetYaxis()->SetTitleOffset(1.6);
  s_eta->GetXaxis()->SetTitleOffset(1.1);
  s_eta->GetYaxis()->SetTitleSize(0.05);
  s_eta->GetXaxis()->SetTitleSize(0.05);
  s_eta->GetYaxis()->SetLabelSize(0.05);
  s_eta->GetXaxis()->SetLabelSize(0.05);
  s_eta->GetYaxis()->SetNdivisions(505);
  s_eta->GetXaxis()->SetRangeUser(-3, 3);
  s_eta->SetMaximum(s_eta->GetMaximum()*1.3);
  TLegend* l_eta = new TLegend(0.64,0.65,0.86,0.85);
  l_eta->SetBorderSize(0);
  l_eta->SetFillStyle(0);
  for(int i=categories.size()-1; i>=0; i--){
    l_eta->AddEntry(h_eta[i], categories[i],"f");
  }
  l_eta->Draw();
  c_eta->SaveAs("../Plots/PFcategories_eta.pdf");

  TCanvas *c_eta_norm = new TCanvas("c_eta_norm","c_eta_norm",600, 600);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.14);
  s_eta_norm->Draw("HIST");
  s_eta_norm->GetXaxis()->SetTitle("#eta");
  s_eta_norm->GetYaxis()->SetTitle("energy fraction");
  s_eta_norm->GetYaxis()->SetTitleOffset(1.6);
  s_eta_norm->GetXaxis()->SetTitleOffset(1.1);
  s_eta_norm->GetYaxis()->SetTitleSize(0.05);
  s_eta_norm->GetXaxis()->SetTitleSize(0.05);
  s_eta_norm->GetYaxis()->SetLabelSize(0.05);
  s_eta_norm->GetXaxis()->SetLabelSize(0.05);
  s_eta_norm->GetYaxis()->SetNdivisions(505);
  s_eta_norm->GetXaxis()->SetRangeUser(-3, 3);
  s_eta_norm->SetMaximum(s_eta_norm->GetMaximum());
  TLegend* l_eta_norm = new TLegend(0.64,0.65,0.86,0.85);
  l_eta_norm->SetBorderSize(0);
  l_eta_norm->SetFillStyle(0);
  for(int i=categories.size()-1; i>=0; i--){
    l_eta_norm->AddEntry(h_eta_Eweight_norm[i], categories[i],"f");
  }
  l_eta_norm->Draw();
  c_eta_norm->SaveAs("../Plots/PFcategories_eta_norm.pdf");

  return 0;
}

vector<TH1F*> NormalisePerBin(vector<TH1F*> hists){
  vector<TH1F*> h_norm;
  for(auto h: hists)h_norm.push_back((TH1F*) h->Clone());
  int Nbins = hists[0]->GetXaxis()->GetNbins();
  for(int bin=1; bin<=Nbins; bin++){
    double binTotal = 0;
    for(auto h: hists){
      binTotal += h->GetBinContent(bin);
    }
    for(auto g: h_norm){
      double oldContent = g->GetBinContent(bin);
      if(binTotal == 0) g->SetBinContent(bin, 0);
      else              g->SetBinContent(bin, oldContent/binTotal);
    }
  }
  return h_norm;
}



void SetPlotStyle(){
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);
  return;
}
