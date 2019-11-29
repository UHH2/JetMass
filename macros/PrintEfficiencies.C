using namespace std;

void PrintEfficiencies(){
	std::vector<TString> Samples = {"QCD", "WUnmatched", "WMatched"};
  TString Hist="Rho_central";

  vector<TString> Dirs={
  // "JetMass_N2DDT_pt450To500",
  "JetMass_N2DDT_pt500To550",
  "JetMass_N2DDT_pt550To600",
  "JetMass_N2DDT_pt600To675",
  "JetMass_N2DDT_pt675To800",
  "JetMass_N2DDT_pt800To1200",
  "JetMass_N2DDT_pt1200ToInf"};

	for(auto sample: Samples){
		cout << "Sample: "<< sample << endl;
		TFile * file = new TFile("../Histograms/W/"+sample+".root");
		double pass=0;
		double fail=0;
		for(auto dir : Dirs){
			TH1F * h_pass = (TH1F*) file->Get(dir+"_pass/"+Hist);
			pass+=h_pass->Integral();
			TH1F * h_fail = (TH1F*) file->Get(dir+"_fail/"+Hist);
			cout << "PtBin: " << dir << " pass/(fail+pass) " << h_pass->Integral() / (h_pass->Integral() + h_fail->Integral()) << endl;
			fail+=h_fail->Integral();
		}
		cout << "All PtBins: pass/(fail+pass) " << pass / (pass + fail) << endl;
	}
}
