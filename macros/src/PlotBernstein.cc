#include "../include/CentralInclude.h"
#include "TRegexp.h"
#include "TGraph2DErrors.h"
#include "TGraph2D.h"
#include "TF2.h"

using namespace std;
using namespace RooFit;

TString generateBernsteinPolyLinComb(int n_pt,int n_rho);
TString generateBernsteinPoly(int lambda, int n, const char * var="x");
TGraph2DErrors * createRpf(TString filePath,vector<TString> FileNames,TString category,bool useRho=false);
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

	//getting BernsteinParameters from fitDiagnostics
  TString filename;

  if(argc != 2){
    cout << "specify <fitDiagnostics.root> file!" << endl;
    return 1;
  }
  else filename = argv[1];
	TString FitDir=filename.Copy();
	FitDir.ReplaceAll("fitDiagnostics.root","");
	FitDir="";
  TFile* file = new TFile(filename);
  vector<TString> bernsteinParNames;
  vector<double> bernsteinPars;

  RooFitResult* result = (RooFitResult*) file->Get("fit_s");
	if(result){
		RooArgList fitargs = result->floatParsFinal();
		for(int i = 0 ; i < fitargs.getSize() ; ++i){
			TString parName=fitargs.at(i)->GetTitle();

			if(parName.Contains("qcd_pass_ralhpTF")){
				bernsteinParNames.push_back(parName);
				RooRealVar* par = (RooRealVar*) result->floatParsFinal().find(parName);
				bernsteinPars.push_back(par->getValV());
			}
		}
	}else{
		cout << "Ralphs fit seems to have failed! setting BersteinParameter to 0." << endl;
		bernsteinPars={0.};
	}


	//////////////////////////////////////////////////
	/////CONSTRUCT OWN R_PF HIST FROM DATA-WJETSMC////
	//////////////////////////////////////////////////
  TString category="JetMass_N2DDT_pt";
	TString filePath = "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2/CMSSW_10_2_10/src/UHH2/JetMass/Histograms/W/";
	vector<TString> FileNames= {"Data.root","WMatched.root","WUnmatched.root"};

	TGraph2DErrors * R_pf_surf = createRpf(filePath,FileNames,category,false);
	TGraph2DErrors * R_pf_surf_QCD = createRpf(filePath,{"QCD.root"},category,false);
	
	///////////////////////////////////////////
	/////CONSTRUCT OWN R_PF HIST FROM QCDMC////
	///////////////////////////////////////////
	// TFile* File = new TFile("/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2/CMSSW_10_2_10/src/UHH2/JetMass/Histograms/W/QCD.root");
	
  // TGraph2DErrors * R_pf_surf_QCD = new TGraph2DErrors(N_RhoBins);

  // for(int i = 0 ; i < N_PtBins ; ++i){
  //   TString histName = "/Mass_central";
	// 	TString passDir = passNames[i].Copy();
	// 	TString failDir = failNames[i].Copy();
		
  //   TH1F* h_pass = (TH1F*) Files[0]->Get(passDir+histName);
  //   if(RebinMSD)h_pass=(TH1F*)h_pass->Rebin(N_MSD,"",new_binning);
  //   TH1F* h_fail = (TH1F*) Files[0]->Get(failDir+histName);
  //   if(RebinMSD)h_fail=(TH1F*)h_fail->Rebin(N_MSD,"",new_binning);

  //   TH1F * h_Rpf = (TH1F*) h_pass->Clone();
  //   h_Rpf->Divide(h_fail);

  //   for(int j = 0 ; j < N_MsdBins ; ++j){
  //     double pT = PtBins[i]+0.3*(PtBins[i+1]-PtBins[i]);
  //     double msd = h_info->GetXaxis()->GetBinLowEdge(j)+0.3*h_info->GetXaxis()->GetBinWidth(j);
  //     double rho = 2*TMath::Log(msd/pT);

  //     double R_pf=h_Rpf->GetBinContent(j);
	// 		double R_pf_err=h_Rpf->GetBinError(j);
  //     if(rho<0 && pT>0 && R_pf >0)R_pf_surf_QCD->SetPoint(i*N_MsdBins+j,rho,pT,R_pf);
  //     if(rho<0 && pT>0 && R_pf >0)R_pf_surf_QCD->SetPointError(i*N_MsdBins+j,0,0,R_pf_err);
	// 	}
  // }

	///////////////////
	//PLOT EVERYTHING//
	///////////////////
	
	int order_pt,order_rho;
	if(bernsteinPars.size()==1){
	order_pt = 0;
	order_rho = 0;
	}else{
	order_pt = 3;
	order_rho = 3;
	}
	TF2 * BernsteinPolRalph = new TF2("f1",generateBernsteinPolyLinComb(order_pt,order_rho),-6.0,-2.1,500,900);


	BernsteinPolRalph->SetParameters(&bernsteinPars[0]);
	BernsteinPolRalph->SetLineColor(kOrange);

	TF2 * BernsteinPol = new TF2("f2",generateBernsteinPolyLinComb(1,3),-6.0,-2.1,500,900);
	BernsteinPol->SetLineColor(kBlue);
	// BernsteinPol->SetParameters(-1.29947e-05, 4.01311e-07, 1.63559e-05, 1.77193e-08, 1.48458e-08, 3.69993e-08, 1.29183e-05,-4.83653e-07,-1.63950e-05);


	TCanvas* c = new TCanvas("c", "c", 600, 600);

  
	R_pf_surf->Draw("P0ERR");
	R_pf_surf_QCD->SetLineColor(kGreen);
	R_pf_surf_QCD->SetMarkerColor(kGreen);
	R_pf_surf_QCD->SetMarkerStyle(2);
	R_pf_surf_QCD->Draw("P0ERRSAME");
	R_pf_surf->GetXaxis()->SetRangeUser(-6.0,-2.1);
  BernsteinPolRalph->Draw("surfsame");
  // R_pf_surf->Fit(BernsteinPol,"R");
	// BernsteinPol->Draw("surfSAME");

	cout << " Chi2: " << BernsteinPol->GetChisquare() << " NDOF: " << BernsteinPol->GetNDF() << " Chi2/NDOF: " << BernsteinPol->GetChisquare()/BernsteinPol->GetNDF() << " prob: " << BernsteinPol->GetProb() << endl;

	// R_pf_surf->SetMarkerStyle(24);
	// R_pf_surf->SetMarkerSize(1);

  c->SaveAs(FitDir+"R_pf.root");
  c->SaveAs(FitDir+"R_pf.pdf");



  TFile * OutFile = new TFile(FitDir+"R_pf_Obj.root","RECREATE");
  OutFile->cd();
	BernsteinPolRalph->Write();
  R_pf_surf->Write();
  OutFile->Close();

  return 0;
}


TString generateBernsteinPolyLinComb(int n_pt,int n_rho){
  TString linComb;
  int parCount=0;
  for(int l=0; l<= n_pt; ++l){
    for(int k=0; k<= n_rho;++k){
      linComb+=TString::Format("[%i]*",parCount)+generateBernsteinPoly(k,n_rho,"x")+"*"+generateBernsteinPoly(l,n_pt,"y");
      if(parCount + 1  < (n_pt+1)*(n_rho+1))linComb+="+";
      ++parCount;
    }
  }
  return linComb;
}

TString generateBernsteinPoly(int lambda, int n, const char * var){
	return TString::Format("TMath::Binomial(%i,%i)*TMath::Power(%s,%i)*TMath::Power((1-%s),%i)",n,lambda,var,lambda,var,n-lambda);
}




TGraph2DErrors * createRpf(TString filePath,vector<TString> FileNames,TString category,bool useRho){
  //Rebin mSD
  bool RebinMSD=true;
  double minMSD=50;
  double maxMSD=170;
  int N_MSD = 20;
  double binWidth = (maxMSD-minMSD)/N_MSD;
  vector<double> bins;
  for(int i = 0 ; i<N_MSD ; ++i){
    bins.push_back(minMSD+i*binWidth);
  }
  bins.push_back(maxMSD);
  double * new_binning = &bins[0];


	vector<TFile*> Files;
  Files.reserve(FileNames.size());

  for(auto name : FileNames)Files.push_back(new TFile(filePath+name));

  vector<double> PtBins;
  vector<TString> passNames;
  vector<TString> failNames;

	if(TString(Files[0]->GetName()).Contains("Data"))category.ReplaceAll("pt","noConst_pt");
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
	
  TH1F*  h_info = (TH1F*) Files[0]->Get(passNames[0]+"/"+"Mass_central");
  if(RebinMSD)h_info=(TH1F*)h_info->Rebin(N_MSD,"",new_binning);

  //TODO Get upperEdge from histDirName
  int N_PtBins = PtBins.size()-1;//This works only if the last PtBin in rootFiles will not be used. (1200ToInf is not used and 1200 is only used as upperEdge of last bin)
  int N_MsdBins = h_info->GetNbinsX();
	int N_RhoBins = N_PtBins*N_MsdBins;

  TGraph2DErrors * R_pf_surf = new TGraph2DErrors(N_RhoBins);

  for(int i = 0 ; i < N_PtBins ; ++i){
    TString histName = "/Mass_central";
		TString passDir = passNames[i].Copy();
		TString failDir = failNames[i].Copy();
		
    TH1F* h_pass = (TH1F*) Files[0]->Get(passDir+histName);
    if(RebinMSD)h_pass=(TH1F*)h_pass->Rebin(N_MSD,"",new_binning);
    TH1F* h_fail = (TH1F*) Files[0]->Get(failDir+histName);
    if(RebinMSD)h_fail=(TH1F*)h_fail->Rebin(N_MSD,"",new_binning);

		if(TString(Files[0]->GetName()).Contains("Data")){
			passDir.ReplaceAll("_noConst","");
			failDir.ReplaceAll("_noConst","");
		}

		for(int k = 1 ; k < Files.size() ; ++k){
      TH1F* h_pass_add_tmp=(TH1F*) Files[k]->Get(passDir+histName);
      if(RebinMSD)h_pass_add_tmp=(TH1F*)h_pass_add_tmp->Rebin(N_MSD,"",new_binning);
      h_pass->Add(h_pass_add_tmp,-1.);
      TH1F* h_fail_add_tmp=(TH1F*) Files[k]->Get(failDir+histName);
      if(RebinMSD)h_fail_add_tmp=(TH1F*)h_fail_add_tmp->Rebin(N_MSD,"",new_binning);
      h_fail->Add(h_fail_add_tmp,-1.);
      delete h_pass_add_tmp,h_fail_add_tmp;
    }

    TH1F * h_Rpf = (TH1F*) h_pass->Clone();
    h_Rpf->Divide(h_fail);

    for(int j = 0 ; j < N_MsdBins ; ++j){
      double pT = PtBins[i]+0.3*(PtBins[i+1]-PtBins[i]);
      double msd = h_info->GetXaxis()->GetBinLowEdge(j)+0.3*h_info->GetXaxis()->GetBinWidth(j);
      double rho = 2*TMath::Log(msd/pT);

      double R_pf=h_Rpf->GetBinContent(j);
			double R_pf_err=h_Rpf->GetBinError(j);
      if(rho<0 && pT>0 && R_pf >0)R_pf_surf->SetPoint(i*N_MsdBins+j,rho,pT,R_pf);
      if(rho<0 && pT>0 && R_pf >0)R_pf_surf->SetPointError(i*N_MsdBins+j,0,0,R_pf_err);
		}
  }
	return R_pf_surf;
}
