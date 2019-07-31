#include "../include/CentralInclude.h"


using namespace std;

int main(int argc, char* argv[]){

  vector<double> ptbins = {0, 5, 50};
  vector<double> etabins = {0, 1, 5};
  vector<TString> categories = {"chargedH", "neutralH", "gamma", "other"};
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
