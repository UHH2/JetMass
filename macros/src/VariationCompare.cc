#include "../include/CentralInclude.h"


using namespace std;
void SetPlotStyle();
TH1F* GetRatio(TH1F*, TH1F*);

int main(int argc, char* argv[]){

  TFile* file = new TFile("/nfs/dust/cms/user/schwarzd/JetMass/PostSel/uhh2.AnalysisModuleRunner.MC.PseudoData.root");
  vector<TString> dirs = {"JetMass_pt200", "JetMass_pt300", "JetMass_pt400"};
  vector<TString> n_dirs = {"p_{T} > 200", "p_{T} > 300", "p_{T} > 400"};
  vector<TString> vars00 = {"Mass_0_0_chargedH", "Mass_0_0_neutralH", "Mass_0_0_gamma"};
  vector<TString> vars01 = {"Mass_0_1_chargedH", "Mass_0_1_neutralH", "Mass_0_1_gamma"};
  vector<TString> vars02 = {"Mass_0_2_chargedH", "Mass_0_2_neutralH", "Mass_0_2_gamma"};
  vector<TString> vars10 = {"Mass_1_0_chargedH", "Mass_1_0_neutralH", "Mass_1_0_gamma"};
  vector<TString> vars11 = {"Mass_1_1_chargedH", "Mass_1_1_neutralH", "Mass_1_1_gamma"};
  vector<TString> vars12 = {"Mass_1_2_chargedH", "Mass_1_2_neutralH", "Mass_1_2_gamma"};

  vector<TString> vars;
  vars.insert(vars.end(), vars00.begin(), vars00.end());
  vars.insert(vars.end(), vars01.begin(), vars01.end());
  vars.insert(vars.end(), vars02.begin(), vars02.end());
  vars.insert(vars.end(), vars10.begin(), vars10.end());
  vars.insert(vars.end(), vars11.begin(), vars11.end());
  vars.insert(vars.end(), vars12.begin(), vars12.end());


  int Nrebin = 4;

  vector<TH1F*> h_central, r_central;
  vector<vector<TH1F*>> h_up, h_down, r_up, r_down;
  for(auto dir: dirs){
    TH1F* c = (TH1F*) file->Get(dir+"/Mass_central");
    c->Rebin(Nrebin);
    h_central.push_back(c);
    r_central.push_back(GetRatio(c,c));
    vector<TH1F*> h_up_temp, h_down_temp, r_up_temp, r_down_temp;
    for(auto var: vars){
      TH1F* u = (TH1F*) file->Get(dir+"/"+var+"_up");
      TH1F* d = (TH1F*) file->Get(dir+"/"+var+"_down");
      u->Rebin(Nrebin);
      d->Rebin(Nrebin);
      h_up_temp.push_back(u);
      h_down_temp.push_back(d);
      r_up_temp.push_back(GetRatio(u,c));
      r_down_temp.push_back(GetRatio(d,c));
    }
    h_up.push_back(h_up_temp);
    h_down.push_back(h_down_temp);
    r_up.push_back(r_up_temp);
    r_down.push_back(r_down_temp);
  }



  SetPlotStyle();
  for(int i=0; i<dirs.size(); i++){
    for(int j=0; j<vars.size(); j++){
      TCanvas *C = new TCanvas("C","C",600, 600);
      gPad->SetLeftMargin(0.17);
      gPad->SetBottomMargin(0.14);
      r_central[i]->SetLineWidth(3);
      r_central[i]->SetLineColor(13);
      r_central[i]->Draw("HIST");
      r_up[i][j]->SetLineWidth(3);
      r_up[i][j]->SetLineColor(kRed-2);
      r_up[i][j]->Draw("HIST SAME");
      r_down[i][j]->SetLineWidth(3);
      r_down[i][j]->SetLineColor(kAzure+7);
      r_down[i][j]->Draw("HIST SAME");
      TString title = vars[j];
      title.Remove(0,5);
      title.ReplaceAll("_", " ");
      r_central[i]->SetTitle(n_dirs[i]+", "+title);
      r_central[i]->GetXaxis()->SetTitle("m_{jet}");
      r_central[i]->GetYaxis()->SetTitle("#frac{variation}{central}");
      r_central[i]->GetYaxis()->SetTitleOffset(1.6);
      r_central[i]->GetXaxis()->SetTitleOffset(1.1);
      r_central[i]->GetYaxis()->SetTitleSize(0.05);
      r_central[i]->GetXaxis()->SetTitleSize(0.05);
      r_central[i]->GetYaxis()->SetLabelSize(0.05);
      r_central[i]->GetXaxis()->SetLabelSize(0.05);
      r_central[i]->GetXaxis()->SetNdivisions(505);
      r_central[i]->GetYaxis()->SetNdivisions(505);
      r_central[i]->GetYaxis()->SetRangeUser(0.5, 1.5);
      TLegend* leg = new TLegend(0.64,0.65,0.86,0.85);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry(r_central[i], "central", "l");
      leg->AddEntry(r_up[i][j], "up", "l");
      leg->AddEntry(r_down[i][j], "down", "l");
      leg->SetTextSize(0.04);
      leg->Draw();
      C->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/JetMass/VariationComparision/"+dirs[i]+"_"+vars[j]+".pdf");
      delete C;
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

TH1F* GetRatio(TH1F* var, TH1F* central){
  TH1F* h_ratio = (TH1F*) var->Clone();
  int Nbins = var->GetXaxis()->GetNbins();
  for(int bin=1; bin<= Nbins; bin++){
    double denominator = central->GetBinContent(bin);
    double numerator = var->GetBinContent(bin);
    double ratio = 0.0;
    double error = 0.0;
    if(denominator != 0){
      ratio = numerator/denominator;
      error = sqrt(fabs((1-ratio)*ratio/denominator));
    }
    h_ratio->SetBinContent(bin, ratio);
    h_ratio->SetBinError(bin, error);
  }
  return h_ratio;
}
