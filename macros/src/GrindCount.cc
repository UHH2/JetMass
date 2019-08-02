#include "../include/CentralInclude.h"

using namespace std;

int main(int argc, char* argv[]){
  TFile *file = new TFile("../Histograms/top/Data.root");
  TH2F * h_countAll = (TH2F*)file->Get("PFHists/count_all");

  h_countAll->SetTitle(" ");
  h_countAll->GetXaxis()->SetTitle("p_{T}");
  h_countAll->GetYaxis()->SetTitle("#eta");
  h_countAll->GetZaxis()->SetTitle("events");
  h_countAll->GetYaxis()->SetTitleSize(0.06);
  h_countAll->GetXaxis()->SetTitleSize(0.05);
  h_countAll->GetZaxis()->SetTitleSize(0.05);
  h_countAll->GetXaxis()->SetTitleOffset(0.9);
  h_countAll->GetYaxis()->SetTitleOffset(1.1);
  h_countAll->GetZaxis()->SetTitleOffset(0.9);
  h_countAll->GetXaxis()->SetNdivisions(505);
  h_countAll->GetYaxis()->SetNdivisions(505);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);

  TCanvas *b = new TCanvas();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogz();
  h_countAll->Draw("COLZ");
  h_countAll->Draw("BOX SAME");
  b->SaveAs("../Plots/GridCount.pdf");

  //// ---------------------------------------------------------------------------------------------------------------------
  //// ---------------------------------------------------------------------------------------------------------------------


}
