#include "../include/CentralInclude.h"


using namespace std;
using namespace RooFit;

void PlotNuisance(vector<TString> parname, vector<double> par, vector<double> errU, vector<double> errD);
void SetPlotStyle();

int main(int argc, char* argv[]){

  TString filename;

  if(argc != 2){
    cout << "specify file!" << endl;
    return 1;
  }
  else filename = argv[1];
  TFile* file = new TFile(filename);

  // get grid
  TFile* gfile = new TFile("../Histograms/grid.root");
  TH2F* grid = (TH2F*) gfile->Get("grid");
  TH2F* grid_fit = (TH2F*) grid->Clone("grid_fit");
  TH2F* grid_fit_up = (TH2F*) grid->Clone("grid_fit_up");
  TH2F* grid_fit_down = (TH2F*) grid->Clone("grid_fit_down");

  RooFitResult* result = (RooFitResult*) file->Get("fit_s");
  RooArgSet fitargs = result->floatParsFinal();

  int Npt = grid->GetXaxis()->GetNbins();
  int Neta = grid->GetYaxis()->GetNbins();

  vector<TString> parname;
  vector<double> par;
  vector<double> errU;
  vector<double> errD;

  for(int i=0; i<Npt; i++){
    for(int j=0; j<Neta; j++){
      TString parameter_name = "massScale_";
      parameter_name += "pt" + to_string(i);
      parameter_name += "_eta" + to_string(j);
      RooRealVar* params = (RooRealVar*) result->floatParsFinal().find(parameter_name);
      double central = params->getValV();
      double error_up = fabs(params->getErrorHi());
      double error_down = fabs(params->getErrorLo());
      grid_fit->SetBinContent(i+1,j+1,central);
      grid_fit_up->SetBinContent(i+1,j+1,central+error_up);
      grid_fit_down->SetBinContent(i+1,j+1,central-error_down);
      parname.push_back(parameter_name);
      par.push_back(central);
      errU.push_back(error_up);
      errD.push_back(error_down);
    }
  }

  TFile* outputFile=new TFile("../Histograms/grid_fit.root","recreate");
  outputFile->cd();
  grid_fit->Write("grid_fit");
  grid_fit_up->Write("grid_fit_up");
  grid_fit_down->Write("grid_fit_down");
  outputFile->Close();

  SetPlotStyle();
  PlotNuisance(parname, par, errU, errD);


  return 0;
}


void PlotNuisance(vector<TString> parname, vector<double> par, vector<double> errU, vector<double> errD){
  int Npoints = parname.size();
  TGraphAsymmErrors* g = new TGraphAsymmErrors(Npoints);
  for(unsigned int i=0; i<parname.size(); i++){
    g->SetPoint(i, i+0.5, par[i]);
    g->SetPointError(i, 0.0, 0.0, errD[i], errU[i]);
  }

  TH1F* frame = new TH1F("frame", " ", parname.size(), 0, parname.size());
  frame->GetYaxis()->SetRangeUser(-2, 2);
  for(unsigned int i=0; i<parname.size(); i++){
    int bin = i+1;
    frame->GetXaxis()->SetBinLabel(bin, parname[i]);
  }

  TCanvas* c = new TCanvas("c", "c", 600, 600);
  frame->SetFillColor(kWhite);
  frame->Draw();
  g->SetMarkerColor(kBlack);
  g->SetMarkerStyle(8);
  g->SetMarkerSize(1);
  g->Draw("P SAME");
  c->SaveAs("Nuisance.pdf");
  return;
}

void SetPlotStyle(){
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);
  return;
}
