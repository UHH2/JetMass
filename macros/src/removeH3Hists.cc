#include "../include/CentralInclude.h"
#include "TH3D.h"

using namespace std;

int main(int argc, char* argv[]){
  TString inFileName=argv[1];
  TString inCycleName=inFileName.Copy();
  inCycleName.ReplaceAll(".root","");
  TString outFileName=inCycleName+"_H3_N2_hists.root";

  cout << "Moving 3D Hists from "<< inFileName << " to " << outFileName << endl;
  TFile * h3File = new TFile(outFileName,"RECREATE");

  TFile * oldFile = new TFile(inFileName,"Update");

  for(auto TKey : *oldFile->GetListOfKeys() ){
    TString keyName=TKey->GetName();
    if( keyName.Contains("H3") || keyName.Contains("3D") ){

      TH3D* h3 = (TH3D*) oldFile->Get(keyName+"/N2_v_mSD_v_pT");
      h3->SetName(keyName+"_N2_v_mSD_v_pT");
      h3->SetDirectory(0);
      h3File->cd();
      h3->Write(keyName+"_N2_v_mSD_v_pT");

      TH3D* h3_1 = (TH3D*) oldFile->Get(keyName+"/N2_v_rho_v_pT");
      h3_1->SetName(keyName+"_N2_v_rho_v_pT");
      h3_1->SetDirectory(0);
      h3File->cd();
      h3_1->Write(keyName+"_N2_v_rho_v_pT");
      }
  }
  h3File->Close();

  for(auto TKey : *oldFile->GetListOfKeys() ){
    TString keyName=TKey->GetName();
    if( keyName.Contains("H3") || keyName.Contains("3D") ){
        oldFile->cd(keyName);
        gDirectory->Delete("N2_v_mSD_v_pT;*");
        gDirectory->Delete("N2_v_rho_v_pT;*");
      }
  }
  oldFile->Close();
  return 0;
}
