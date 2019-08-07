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

  // size of Variation
  double variation = 0.1;
  cout << "size of variation is assumend to be +- " << variation << endl;

  // get grid
  TFile* gfile = new TFile("../Histograms/grid.root");
  TH2F* grid = (TH2F*) gfile->Get("grid");
  TH1F* h_cat = (TH1F*) gfile->Get("categories");
  vector<TString> categories;
  for(int bin=1; bin<=h_cat->GetXaxis()->GetNbins(); bin++) categories.push_back(h_cat->GetXaxis()->GetBinLabel(bin));

  // create grids for new fit values
  vector<TH2F*> grid_fit, grid_fit_up, grid_fit_down;
  for(auto cat: categories){
    grid_fit.push_back((TH2F*) grid->Clone("grid_fit"+cat));
    grid_fit_up.push_back((TH2F*) grid->Clone("grid_fit_up"+cat));
    grid_fit_down.push_back((TH2F*) grid->Clone("grid_fit_down"+cat));
  }

  RooFitResult* result = (RooFitResult*) file->Get("fit_s");
  RooArgSet fitargs = result->floatParsFinal();

  int Npt = grid->GetXaxis()->GetNbins();
  int Neta = grid->GetYaxis()->GetNbins();
  int Ncat = categories.size();

  vector<TString> parname;
  vector<double> par;
  vector<double> errU;
  vector<double> errD;

  for(int i=0; i<Npt; i++){
    for(int j=0; j<Neta; j++){
      for(int k=0; k<Ncat; k++){
        TString parameter_name = "massScale_";
        parameter_name += "pt" + to_string(i);
        parameter_name += "_eta" + to_string(j);
        parameter_name += "_" + categories[k];
        RooRealVar* params = (RooRealVar*) result->floatParsFinal().find(parameter_name);
        double central = params->getValV();
        double error_up = fabs(params->getErrorHi());
        double error_down = fabs(params->getErrorLo());
        double bincont = 1.0 + central * variation;
        double bincontU = 1.0 + (central+error_up) * variation;
        double bincontD = 1.0 + (central-error_up) * variation;
        grid_fit[k]->SetBinContent(i+1,j+1, bincont);
        grid_fit_up[k]->SetBinContent(i+1,j+1, bincontU);
        grid_fit_down[k]->SetBinContent(i+1,j+1, bincontD);
        parname.push_back(parameter_name);
        par.push_back(central);
        errU.push_back(error_up);
        errD.push_back(error_down);
      }
    }
  }

  TFile* outputFile=new TFile("../Histograms/grid_fit.root","recreate");
  outputFile->cd();
  for(int k=0; k<Ncat; k++){
    grid_fit[k]->Write("grid_fit_"+categories[k]);
    grid_fit_up[k]->Write("grid_fit_up_"+categories[k]);
    grid_fit_down[k]->Write("grid_fit_down_"+categories[k]);
  }
  h_cat->Write("categories");
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

  TCanvas* c = new TCanvas("c", "c", 2000, 600);
  gPad->SetBottomMargin(.2);
  // gPad->SetRightMargin(.2);
  frame->GetYaxis()->SetTitle("#sigma");
  frame->SetFillColor(kWhite);
  frame->Draw();
  g->SetMarkerColor(kBlack);
  g->SetMarkerStyle(8);
  g->SetMarkerSize(1);
  g->Draw("P SAME");
  c->SaveAs("../Plots/Nuisance.pdf");
  return;
}

void SetPlotStyle(){
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);
  return;
}
