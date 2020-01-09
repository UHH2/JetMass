TH2F* computeDDTMap(double working_point, TH3D* n2_msd_pt );
void setStyle();
void set2DHistStyle(TH2 * hist);
void plotDDTMap(TH2 * hist, TString title,TString outFile);



void getN2Map(){
	setStyle();
// std::string h3_fileName = "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2/CMSSW_10_2_10/src/UHH2/JetMass/Histograms/W/wrongWeights/N2_v_mSD_v_pT.root";
// std::string h3_fileName = "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2/CMSSW_10_2_10/src/UHH2/JetMass/Histograms/Wtmp/QCD_H3_N2_hists.root";
  std::string h3_fileName = "/nfs/dust/cms/user/albrechs/SingleJetTrees/WMassPreSel_TMP/N2_3Dmap.root";
  std::string h3_histName = "N2_3DMap/N2_v_rho_v_pT";
  TFile * mapFile = new TFile(h3_fileName.c_str(),"READ");
  TH3D* N2_v_mSD_v_pT = (TH3D*) mapFile->Get(h3_histName.c_str());
  TH2D * N2_v_rho = new TH2D("N2_v_rho","N2_v_rho",12,-7.0,-1.0,20,0.0,0.5);
  int NBins_msd = N2_v_mSD_v_pT->GetNbinsX();
  int NBins_pt = N2_v_mSD_v_pT->GetNbinsY();
  int NBins_N2 = N2_v_mSD_v_pT->GetNbinsZ();
  std::cout << "Bins(x,y,z): " << NBins_msd << ","<< NBins_pt << ","<< NBins_N2 << std::endl;
  std::cout << "processing ";
  for(int msdBin = 1 ; msdBin <= NBins_msd ; ++msdBin ){
    for(int ptBin = 1 ; ptBin <= NBins_pt ; ++ptBin ){
      if(ptBin % 100 == 0 ) std::cout << ".";
      for(int N2Bin = 1 ; N2Bin <= NBins_N2 ; ++N2Bin ){
        double msd=N2_v_mSD_v_pT->GetXaxis()->GetBinCenter(msdBin);
        double pt=N2_v_mSD_v_pT->GetYaxis()->GetBinCenter(ptBin);
        double n2=N2_v_mSD_v_pT->GetZaxis()->GetBinCenter(N2Bin);
        // double rho= msd;
        double rho= 2* TMath::Log(msd/pt);
        N2_v_rho->Fill(rho,n2,N2_v_mSD_v_pT->GetBinContent(msdBin,ptBin,N2Bin));
      }
    }
  }
  std::cout << "Done!"<< std::endl;
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  // gStyle->SetOptStat(kFALSE);
  // gStyle->SetPadTickY(1);
  // gStyle->SetPadTickX(1);
  // gStyle->SetLegendBorderSize(0);
  // gPad->SetBottomMargin(.2);
  // gPad->SetRightMargin(.2);
  N2_v_rho->GetXaxis()->SetTitle("#rho=log(m^{2}_{SD}/p^{2}_{T})");
  N2_v_rho->GetYaxis()->SetTitle("N_{2}^{#beta=1}");
  N2_v_rho->Draw("colz");
  c->SaveAs("N2_v_rho.pdf");


  TH2F * ddt_0p05= computeDDTMap(0.05,N2_v_mSD_v_pT);
  // gPad->SetRightMargin(0.2);
	plotDDTMap(ddt_0p05,"Simple N2-DDT map","DDT_0p05.pdf");

  TFile * mymapFile=new TFile("myN2Map.root","RECREATE");
  mymapFile->cd();
  ddt_0p05->SetTitle("Rho2D");
  ddt_0p05->Write();
  mymapFile->Close();

  TFile *mitN2Map = new TFile("/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2/CMSSW_10_2_10/src/UHH2/JetMass/Histograms/GridOutput_v13.root","READ");
  TH2F * mit_ddtMap = (TH2F*) mitN2Map->Get("Rho2D");
	plotDDTMap(mit_ddtMap,"MIT N2-DDT map (with knn smearing)","MIT_DDT_0p05.pdf");

	
	// cout << "ddt map N Bins rho: " << ddt_0p05->GetNbinsX() << " pt: "<< ddt_0p05->GetNbinsY() << endl;
  // cout << "MITddt map N Bins rho: " << mit_ddtMap->GetNbinsX() << " pt: "<< mit_ddtMap->GetNbinsY() << endl; 
  // TH2F * ddtMap_Comp = (TH2F*)ddt_0p05->Clone();

	// ddtMap_Comp->Add(mit_ddtMap,-1.);
  
  // TCanvas * Comp_c = new TCanvas("Comp_c","Comp_c",600,600);
  // Comp_c->cd();
  // ddtMap_Comp->GetZaxis()->SetRangeUser(0,.5);
  // ddtMap_Comp->Draw("colz");
  // Comp_c->SaveAs("DDTMap_Comparison.pdf");

}




void plotDDTMap(TH2 * hist, TString title,TString outFile){

	hist->GetXaxis()->SetTitle("#rho");
	hist->GetYaxis()->SetTitle("p_{T} [GeV]");
	hist->GetZaxis()->SetTitle("N2^{DDT} 5/% quantile");

	hist->GetXaxis()->SetRangeUser(-6.0,-2.1);
  hist->GetYaxis()->SetRangeUser(200,1200);
  hist->GetZaxis()->SetRangeUser(0.12,0.3);

	double Font=43;
	double TitleSize=24.0;
	double TitleOffset=1.3;
	double LabelSize=18.0;
	hist->GetYaxis()->SetTitleFont(Font);
	hist->GetYaxis()->SetTitleSize(TitleSize);
	hist->GetYaxis()->SetTitleOffset(TitleOffset);
	hist->GetYaxis()->SetLabelFont(Font);
	hist->GetYaxis()->SetLabelSize(LabelSize);

	hist->GetXaxis()->SetTitleFont(Font);
	hist->GetXaxis()->SetTitleSize(TitleSize);
	hist->GetXaxis()->SetTitleOffset(TitleOffset);
	hist->GetXaxis()->SetLabelFont(Font);
	hist->GetXaxis()->SetLabelSize(LabelSize);

	hist->GetZaxis()->SetTitleFont(Font);
	hist->GetZaxis()->SetTitleSize(TitleSize);
	hist->GetZaxis()->SetTitleOffset(TitleOffset);
	hist->GetZaxis()->SetLabelFont(Font);
	hist->GetZaxis()->SetLabelSize(LabelSize);

	hist->GetZaxis()->SetNdivisions(510);

	hist->SetTitle(title);
	hist->SetTitleFont(43);
	hist->SetTitleSize(20.0);
	
  TCanvas * c1 = new TCanvas("c1","c1",700,600);
  c1->cd();
  hist->Draw("colz");
  c1->SaveAs(outFile);
	delete c1;
	
}

void setStyle(){

	gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  // gStyle->SetOptTitle(0);
  // gStyle->SetTextFont(43);
  //
  gStyle->SetTitleOffset(0.86,"X");
  gStyle->SetTitleOffset(1.6,"Y");

	// gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadRightMargin(0.16);

  gStyle->SetMarkerSize(0.5);
  gStyle->SetHistLineWidth(1);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");
  gStyle->SetNdivisions(506, "XYZ");
  gStyle->SetLegendBorderSize(0);
}

TH2F* computeDDTMap(double working_point, TH3D* n2_msd_pt ){
  TH2F * ddtMap= (TH2F*)n2_msd_pt->Project3D("yx");
  char ddtMapName[100];
  // sprintf(ddtMapName,"%s_ddt_wp%f",n2_msd_pt->GetName(),working_point);
  sprintf(ddtMapName,"Rho2D");
  ddtMap->SetName(ddtMapName);
  ddtMap->Reset();
  ddtMap->SetStats(0);
  ddtMap->SetDirectory(0);
  int nbins_x = ddtMap->GetNbinsX();
  int nbins_y = ddtMap->GetNbinsY();
  for (int xbin = 1; xbin <= nbins_x; ++xbin) {
    for (int ybin = 1; ybin <= nbins_y; ++ybin) {
      sprintf(ddtMapName, "%s_x%d_y%d", n2_msd_pt->GetName(), xbin, ybin);
      TH1D* n2proj = n2_msd_pt->ProjectionZ(ddtMapName, xbin, xbin, ybin, ybin);
      for(int i = 1 ; i < n2proj->GetNbinsX() ; i++) {
        if(n2proj->GetBinCenter(i) < 0){
          n2proj->SetBinContent(i,0);
        }
      }
      if (n2proj->Integral() == 0) {
        double xval = n2_msd_pt->GetXaxis()->GetBinCenter(xbin);
        double yval = n2_msd_pt->GetYaxis()->GetBinCenter(ybin);
        std::cout << "N2 integral = 0 for xbin=" << xbin << " (" << xval << ") / ybin=" << ybin << " (" << yval << ") for hist " << n2_msd_pt->GetName() << ". Setting DDT to 0." << std::endl;
        ddtMap->SetBinContent(xbin, ybin, 0);
        continue;
      }
      double wp_array[1] = {working_point};
      double quantiles[1] = {0.};
      n2proj->GetQuantiles(1, quantiles, wp_array);
      ddtMap->SetBinContent(xbin, ybin, quantiles[0]);
      delete n2proj;
    }
  }
  return ddtMap;
}
