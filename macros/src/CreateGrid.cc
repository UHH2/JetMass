#include "../include/CentralInclude.h"


using namespace std;

int main(int argc, char* argv[]){

  vector<double> ptbins = {0, 5, 50};
  vector<double> etabins = {0, 1, 5};
  int Npt  = ptbins.size()-1;
  int Neta = etabins.size()-1;
  TFile* outputFile=new TFile("../Histograms/grid.root","recreate");
  outputFile->cd();
  TH2F* grid = new TH2F("grid", "x=pt, y=eta", Npt, &ptbins[0], Neta, &etabins[0]);
  grid->Write("grid");
  outputFile->Close();

  return 0;
}
