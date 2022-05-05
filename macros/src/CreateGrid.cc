#include "../include/CentralInclude.h"


using namespace std;

int main(int argc, char* argv[]){

  // vector<double> ptbins = {0, 10, 50};
  // vector<double> etabins = {0, 0.522, 1.305, 5.191};
  vector<double> ptbins = {0, 100000.};
  vector<double> etabins = {0, 9};
  // eta bins from jec: {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.322, 2.65, 3.139, 5.191}
  vector<TString> categories = {"all","chargedH", "neutralH", "gamma", "other"};
  // vector<TString> categories = {"chargedH", "neutralH", "gamma", "other"};
  int Npt  = ptbins.size()-1;
  int Neta = etabins.size()-1;
  TFile* outputFile=new TFile("../Histograms/grid.root","recreate");
  outputFile->cd();
  TH2F* h_grid = new TH2F("grid", "x=pt, y=eta", Npt, &ptbins[0], Neta, &etabins[0]);
  TH1F* h_categories = new TH1F("categories", "categories", categories.size(), 0, categories.size());
  for(unsigned int i=0; i<categories.size();i++){
    int bin = i+1;
    h_categories->GetXaxis()->SetBinLabel(bin, categories[i]);
  }
  h_grid->Write("grid");
  h_categories->Write("categories");
  outputFile->Close();

  return 0;
}
