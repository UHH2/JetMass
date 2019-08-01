#include "../include/CentralInclude.h"


using namespace std;

void SetPlotStyle();

int main(int argc, char* argv[]){

  TString filename;

  if(argc != 2){
    cout << "specify file!" << endl;
    return 1;
  }
  else filename = argv[1];
  TFile* file = new TFile(filename);

  TGraph* input = (TGraph*) file->Get("params_graph");
  TH2F* grid = (TH2F*) file->Get("grid_fit_chargedH");
  TH1F* h_cat = (TH1F*) file->Get("categories");
  vector<TString> categories;
  for(int bin=1; bin<=h_cat->GetXaxis()->GetNbins(); bin++) categories.push_back(h_cat->GetXaxis()->GetBinLabel(bin));
  // get parameter names

  int Npt = grid->GetXaxis()->GetNbins();
  int Neta = grid->GetYaxis()->GetNbins();
  int Ncat = categories.size();
  vector<TString> parname;
  for(int i=0; i<Npt; i++){
    for(int j=0; j<Neta; j++){
      for(int k=0; k<Ncat; k++){
        TString parameter_name = "massScale_";
        parameter_name += "pt" + to_string(i);
        parameter_name += "_eta" + to_string(j);
        parameter_name += "_" + categories[k];
        parname.push_back(parameter_name);
      }
    }
  }
  SetPlotStyle();
  TH1F* frame = new TH1F("frame", " ", parname.size(), 0, parname.size());
  frame->GetYaxis()->SetRangeUser(-2, 2);
  for(unsigned int i=0; i<parname.size(); i++){
    int bin = i+1;
    frame->GetXaxis()->SetBinLabel(bin, parname[i]);
  }

  TCanvas* c = new TCanvas("c", "c", 1800, 600);
  gPad->SetBottomMargin(.2);
  gPad->SetRightMargin(.2);
  frame->SetFillColor(kWhite);
  frame->Draw();
  input->SetMarkerColor(kBlack);
  input->SetMarkerStyle(8);
  input->SetMarkerSize(1);
  input->Draw("P SAME");
  c->SaveAs("../Plots/Nuisance_Input.pdf");

  return 0;
}

void SetPlotStyle(){
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);
  return;
}
