#include "../include/CentralInclude.h"
#include "TRegexp.h"
#include "TGraph2DErrors.h"
#include "TGraph2D.h"
#include "TF2.h"

using namespace std;
using namespace RooFit;

TString generateBernsteinPolyLinComb(int n_pt,int n_rho);
TString generateBernsteinPoly(int lambda, int n, const char * var="x");
void plotHist(TH1 *hist, TString outName){
  TCanvas *c1  =  new TCanvas();
  hist->Draw("hist");
  c1->SaveAs(outName);
  delete c1;
}


int main(int argc, char* argv[]){
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TString filename;

  if(argc != 2){
    cout << "specify file!" << endl;
    return 1;
  }
  else filename = argv[1];
  TFile* file = new TFile(filename);
  vector<TString> bernsteinParNames;
  vector<double> bernsteinPars;

  // vector<TString> qcdParNames;
  //
  // vector<double> qcdPars,qcdParsPt,qcdParsMSD;
  // set<double> Ptbins,MSDbins;

  RooFitResult* result = (RooFitResult*) file->Get("fit_s");
  RooArgList fitargs = result->floatParsFinal();
  for(int i = 0 ; i < fitargs.getSize() ; ++i){
    TString parName=fitargs.at(i)->GetTitle();

    if(parName.Contains("qcd_pass_ralhpTF")){
      bernsteinParNames.push_back(parName);
			cout << parName << endl;
      RooRealVar* par = (RooRealVar*) result->floatParsFinal().find(parName);
      bernsteinPars.push_back(par->getValV());
    }
    // if(parName.Contains("qcdparam")){
    //   TString pT_TString=parName(TRegexp("Pt[0-9]+"));
    //   string pT_string=(string) pT_TString.ReplaceAll("Pt","");
    //   double pT = stod(pT_string);
    //
    //   TString mSD_TString=parName(TRegexp("msdbin[0-9]+"));
    //   string mSD_string=(string) mSD_TString.ReplaceAll("msdbin","");
    //   double mSD = stod(mSD_string);
    //
    //   RooRealVar* params = (RooRealVar*) result->floatParsFinal().find(parName);
    //
    //   qcdPars.push_back(params->getValV());
    //   qcdParsPt.push_back(pT);
    //   Ptbins.insert(pT);
    //   qcdParsMSD.push_back(mSD);
    //   MSDbins.insert(mSD);
    // }
  }

  // cout << "BernsteinParameter:" << endl;
  // for(auto name : bernsteinParNames)cout << name << endl;

  TString category="JetMass_Tau21DDT_pt";
  vector<TString> FileNames= {"Pseudo.root","WMatched.root","WUnmatched.root"};
  vector<TFile*> Files;
  Files.reserve(FileNames.size());

  for(auto name : FileNames)Files.push_back(new TFile("/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2/CMSSW_10_2_10/src/UHH2/JetMass/Histograms/W/"+name));

  vector<double> PtBins;
  vector<TString> passNames;
  vector<TString> failNames;

  for(auto TKey : *Files[0]->GetListOfKeys()){
    TString dir = (TString) TKey->GetName();
    if(dir.Contains(category) && dir.Contains("pass")){
      TString pT_TString=dir(TRegexp("pt[0-9]+"));
      string pT_string=(string) pT_TString.ReplaceAll("pt","");
      double pT = stod(pT_string);
      PtBins.push_back(pT);
      passNames.push_back(dir);
      failNames.push_back(dir.ReplaceAll("pass","fail"));
    }
  }

  //Rebin mSD
  bool RebinMSD=false;
  double minMSD=50;
  double maxMSD=170;
  int N_MSD = 21;
  double binWidth = (maxMSD-minMSD)/N_MSD;
  vector<double> bins;
  for(int i = 0 ; i<N_MSD ; ++i){
    bins.push_back(minMSD+i*binWidth);
  }
  bins.push_back(maxMSD);
  double * new_binning = &bins[0];


  TH1F*  h_info = (TH1F*) Files[0]->Get(passNames[0]+"/"+"Mass_central");
  if(RebinMSD)h_info=(TH1F*)h_info->Rebin(N_MSD,"",new_binning);

  //TODO Get upperEdge from histDirName
  int N_PtBins = PtBins.size()-1;//This works only if the last PtBin in rootFiles will not be used. (1200ToInf is not used and 1200 is only used as upperEdge of last bin)
  int N_MsdBins = h_info->GetNbinsX();
  int N_RhoBins = N_PtBins*N_MsdBins;
  TGraph2DErrors * R_pf_surf = new TGraph2DErrors(N_PtBins*N_RhoBins);
  vector<double> Ptpts;
  vector<double> Rhopts;
  int lookAtBin =50;

  for(int i = 0 ; i < N_PtBins ; ++i){
    TString histName = "/Mass_central";
    TH1F* h_pass = (TH1F*) Files[0]->Get(passNames[i]+histName);
    if(RebinMSD)h_pass=(TH1F*)h_pass->Rebin(N_MSD,"",new_binning);
    TH1F* h_fail = (TH1F*) Files[0]->Get(failNames[i]+histName);
    if(RebinMSD)h_fail=(TH1F*)h_fail->Rebin(N_MSD,"",new_binning);
    for(int k = 1 ; k < Files.size() ; ++k){
      TH1F* h_pass_add_tmp=(TH1F*) Files[k]->Get(passNames[i]+histName);
      if(RebinMSD)h_pass_add_tmp=(TH1F*)h_pass_add_tmp->Rebin(N_MSD,"",new_binning);
      h_pass->Add(h_pass_add_tmp,-1.);
      TH1F* h_fail_add_tmp=(TH1F*) Files[k]->Get(failNames[i]+histName);
      if(RebinMSD)h_fail_add_tmp=(TH1F*)h_fail_add_tmp->Rebin(N_MSD,"",new_binning);
      h_fail->Add(h_fail_add_tmp,-1.);
      delete h_pass_add_tmp,h_fail_add_tmp;
    }

    TH1F * h_Rpf = (TH1F*) h_pass->Clone();
    h_Rpf->Divide(h_fail);

    for(int j = 0 ; j < N_MsdBins ; ++j){
    // for(int j = 55 ; j < 56 ; ++j){
      double pT = PtBins[i]+0.3*(PtBins[i+1]-PtBins[i]);
      double msd = h_info->GetXaxis()->GetBinLowEdge(j)+0.3*h_info->GetXaxis()->GetBinWidth(j);
      double rho = 2*TMath::Log(msd/pT);

      double R_pf=h_Rpf->GetBinContent(j);
			double R_pf_err=h_Rpf->GetBinError(j);
      if(rho<0 && pT>0 && R_pf >0)R_pf_surf->SetPoint(i*N_MsdBins+j,rho,pT,R_pf);cout << "SetPoint("<<i*N_MsdBins+j<< ")" << endl;
      if(rho<0 && pT>0 && R_pf >0)R_pf_surf->SetPointError(i*N_MsdBins+j,0,0,R_pf_err);
    }
  }


	TF2 * BernsteinPol = new TF2("f1",generateBernsteinPolyLinComb(2,3),-6.0,-1.0,500,900);
  BernsteinPol->SetParameters(&bernsteinPars[0]);
	BernsteinPol->SetLineColor(kOrange);
	cout << "N_fitPar: " << bernsteinPars.size() << endl;
	TCanvas* c = new TCanvas("c", "c", 600, 600);
  R_pf_surf->Draw("P0ERR");
  BernsteinPol->Draw("surfsame");
  // R_pf_surf->Fit(BernsteinPol,"R");
	// R_pf_surf->SetMarkerStyle(24);
	// R_pf_surf->SetMarkerSize(1);

  c->SaveAs("R_pf.root");
  c->SaveAs("R_pf.pdf");



  TFile * OutFile = new TFile("R_pf_Obj.root","RECREATE");
  OutFile->cd();
	BernsteinPol->Write();
  R_pf_surf->Write();
  OutFile->Close();

  // TCanvas* c = new TCanvas("c", "c", 600, 600);
  // gPad->SetBottomMargin(.2);
  // // gPad->SetRightMargin(.2);
  // R_pf_surf->SetMarkerColor(kBlack);
  // R_pf_surf->SetMarkerStyle(8);
  // R_pf_surf->SetMarkerSize(1);
  // R_pf_surf->Draw("P");
  // c->SaveAs("~/www/R_pf.root");


  // TH1F * h_pass = (TH1F*) Files[0]->Get("")

  // cout << "QcdParameter:"<<endl;
  //
  // // TH2D * QCD = new TH2D("QCD","x=m_{SD},y=p_{T}",MSDbins.size(),qcdParsMSD[0],qcdParsMSD[qcdParsMSD.size()-1],Ptbins.size(),qcdParsPt[0],qcdParsPt[qcdParsPt.size()-1]);
  //   const int NbinsX=MSDbins.size();
  //   Double_t BoundariesX[NbinsX] = {};
  //
  //   // TH1F* mass=new TH1F("dijet_mass","M_{jj}",NBINS-1,BOUNDARIES);
  //
  // TH2D * QCD = new TH2D("QCD","x=m_{SD},y=p_{T}",MSDbins.size(),50,170,Ptbins.size(),qcdParsPt[0],qcdParsPt[qcdParsPt.size()-1]);
  // for(int i = 1; i<QCD->GetNbinsX() ; ++i){
  //   cout << "X: BinCenter("<< i << "): " << QCD->GetXaxis()->GetBinLowEdge(i)<< endl;
  // }
  // for(int j = 1; j<QCD->GetNbinsY() ; ++j){
  //   cout << "Y: BinCenter("<< j << "): " << QCD->GetYaxis()->GetBinLowEdge(j)<< endl;
  // }
  return 0;
}


TString generateBernsteinPolyLinComb(int n_pt,int n_rho){
  TString linComb;
  int parCount=0;
  for(int l=0; l<= n_pt; ++l){
    for(int k=0; k<= n_rho;++k){
			// cout << "x,y order" << k << "/" <<n_rho << " , " << l<<"/"<<n_pt << endl;
      linComb+=TString::Format("[%i]*",parCount)+generateBernsteinPoly(k,n_rho,"x")+"*"+generateBernsteinPoly(l,n_pt,"y");
      if(parCount + 1  < (n_pt+1)*(n_rho+1))linComb+="+";
      // linComb+=TString::Format("[%i]*%s*%s",parCount,generateBernsteinPoly(k,n_rho,"x"),generateBernsteinPoly(l,n_pt,"y"));
      ++parCount;
    }
  }
  return linComb;
}

TString generateBernsteinPoly(int lambda, int n, const char * var){
	return TString::Format("TMath::Binomial(%i,%i)*TMath::Power(%s,%i)*TMath::Power((1-%s),%i)",n,lambda,var,lambda,var,n-lambda);
}
