#include "../include/CentralInclude.h"


using namespace std;
using namespace RooFit;

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

  for(int i=0; i<Npt; i++){
    for(int j=0; j<Neta; j++){
      TString parameter_name = "massScale_";
      parameter_name += "pt" + to_string(i);
      parameter_name += "_eta" + to_string(j);
      RooRealVar* par = (RooRealVar*) result->floatParsFinal().find(parameter_name);
      double central = par->getValV();
      double error_up = par->getErrorHi();
      double error_down = par->getErrorLo();
      grid_fit->SetBinContent(i+1,j+1,central);
      grid_fit_up->SetBinContent(i+1,j+1,central+error_up);
      grid_fit_down->SetBinContent(i+1,j+1,central-error_down);
    }
  }

  TFile* outputFile=new TFile("../Histograms/grid_fit.root","recreate");
  outputFile->cd();
  grid_fit->Write("grid_fit");
  grid_fit_up->Write("grid_fit_up");
  grid_fit_down->Write("grid_fit_down");
  outputFile->Close();

  return 0;
}
