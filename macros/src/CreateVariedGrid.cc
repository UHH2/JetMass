#include "../include/CentralInclude.h"


using namespace std;
using namespace RooFit;

double randomBetween(double, double);
void PlotInputParameters(TGraph*, vector<TString>, TString);
void SetPlotStyle();

int main(int argc, char* argv[]){
  SetPlotStyle();

  // get grid
  TFile* gfile = new TFile("../Histograms/grid.root");
  TH2F* grid = (TH2F*) gfile->Get("grid");
  TH1F* h_cat = (TH1F*) gfile->Get("categories");
  vector<TString> categories;
  for(int bin=1; bin<=h_cat->GetXaxis()->GetNbins(); bin++) categories.push_back(h_cat->GetXaxis()->GetBinLabel(bin));

  vector<TString>  variations = {"chargedHUP_oneBin", "neutralHDOWN_gammaDOWN", "highptUP", "highetaDOWN"};
  for(auto var: variations){
    vector<TH2F*> grid_fit;
    for(auto cat: categories){
      grid_fit.push_back((TH2F*) grid->Clone("grid_fit"+cat));
    }

    int Npt = grid->GetXaxis()->GetNbins();
    int Neta = grid->GetYaxis()->GetNbins();
    int Ncat = categories.size();

    vector<TString> parname;
    vector<double> par;

    double variation = 0.1;

    for(int i=0; i<Npt; i++){
      for(int j=0; j<Neta; j++){
        for(int k=0; k<Ncat; k++){
          TString parameter_name = "massScale_";
          parameter_name += "pt" + to_string(i);
          parameter_name += "_eta" + to_string(j);
          parameter_name += "_" + categories[k];
          double central = 0.0;
          if(var == "chargedHUP_oneBin"){
            if(i==1 && j==0 && categories[k] == "chargedH") central = 1.0;
          }
          else if(var == "neutralHDOWN_gammaDOWN"){
            if(categories[k] == "neutralH") central = -1.5;
            if(categories[k] == "gamma") central = -1.0;
          }
          else if(var == "highptUP"){
            if(i == (Npt-1)) central = 1.2;;
          }
          else if(var == "highetaDOWN"){
            if(j == (Neta-1)) central = -1.2;
          }
          double bincont = 1.0 + central*variation;
          grid_fit[k]->SetBinContent(i+1,j+1,bincont);
          parname.push_back(parameter_name);
          par.push_back(central);
        }
      }
    }
    TFile* outputFile=new TFile("../Histograms/grid_"+var+".root","recreate");
    outputFile->cd();
    for(int k=0; k<Ncat; k++){
      grid_fit[k]->Write("grid_fit_"+categories[k]);
    }
    h_cat->Write("categories");

    // also save variation parameters
    TGraph *g = new TGraph(par.size());
    for(unsigned int i=0; i<parname.size(); i++){
      g->SetPoint(i, i+0.5, par[i]);
    }
    TH1F* frame = new TH1F("frame", " ", parname.size(), 0, parname.size());
    frame->GetYaxis()->SetRangeUser(-2, 2);
    for(unsigned int i=0; i<parname.size(); i++){
      int bin = i+1;
      frame->GetXaxis()->SetBinLabel(bin, parname[i]);
    }
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    gPad->SetBottomMargin(.2);
    gPad->SetRightMargin(.2);
    frame->SetFillColor(kWhite);
    frame->Draw();
    g->SetMarkerColor(kBlack);
    g->SetMarkerStyle(8);
    g->SetMarkerSize(1);
    g->Draw("P SAME");
    //
    g->Write("params_graph");
    c->Write("params");
    outputFile->Close();
    cout << "Created ../Histograms/grid_"+var+".root" << endl;

    // directly plot input parameters
    PlotInputParameters(g, parname, var);
    // delete g;
    delete c;
    // delete frame;
  }

  return 0;
}

double randomBetween(double min, double max){
  double x = ((double) rand() / (RAND_MAX));
  double r = min + ( x * (max - min) );
  return r;
}


void PlotInputParameters(TGraph* input, vector<TString> parname, TString name){
  TH1F* frame2 = new TH1F("frame2", " ", parname.size(), 0, parname.size());
  frame2->GetYaxis()->SetRangeUser(-2, 2);
  for(unsigned int i=0; i<parname.size(); i++){
    int bin = i+1;
    frame2->GetXaxis()->SetBinLabel(bin, parname[i]);
  }

  TCanvas* d = new TCanvas("d", "d", 2000, 600);
  gPad->SetBottomMargin(.2);
  // gPad->SetRightMargin(.15);
  frame2->SetFillColor(kWhite);
  frame2->Draw();
  input->SetMarkerColor(kBlack);
  input->SetMarkerStyle(8);
  input->SetMarkerSize(1);
  input->Draw("P SAME");
  d->SaveAs("../Plots/Nuisance_Input_"+name+".pdf");
  delete d;
  delete frame2;
}

void SetPlotStyle(){
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);
  return;
}
