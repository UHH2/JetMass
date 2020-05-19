#include "../include/CentralInclude.h"

using namespace std;

void fill_hists_top(TString dir, TString process);
void fill_hists_W(TString dir, TString process, int job_index_, int n_jobs_);
bool passN2ddt(double n2, double pt, double mass);
double derive_kfactor(double pt);

void read_grid(TString gfilename);

TFile* outputFile;
vector<TString> handlenames;
TH2F* ddtmap;
TH1F *h_kfactor, *h_ewcorr;
std::string NLOWeightsDir = "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2/CMSSW_10_2_10/src/UHH2/UHHNtupleConverter/NLOweights";



int main(int argc, char* argv[]){

  bool fill_top = true;
  bool fill_W = true;

  if(argc > 1){
    if(strcmp(argv[1], "top") == 0) fill_W=false;
    else if(strcmp(argv[1], "W") == 0) fill_top=false;
  }


  int job_index = 1;
  int n_jobs = 1;
  if(argc > 2){
    job_index = atoi(argv[2]);
    n_jobs = atoi(argv[3]);
  }

  read_grid("../Histograms/grid.root");

  // read in ddtmap here
  TString ddtMapFile = "../macros/scripts/maps/fixedPU/QCD_2017_smooth_gaus4p00sigma.root";
  TString ddtMapName = "N2_v_pT_v_rho_0p05_smooth_gaus4p00sigma_maps_cleaner";
  TFile * mapFile = new TFile(ddtMapFile);
  ddtmap = (TH2F*)mapFile->Get(ddtMapName);

  TString outname;
  if(n_jobs > 1){
    outname.Form("workdir/Histograms_%i.root",job_index);
    cout << outname << endl;
  }else{
    outname = "Histograms.root";
  }
  outputFile = new TFile(outname,"recreate");

  TString histdir_W = "/nfs/dust/cms/user/albrechs/UHH2/JetMassOutput/WMassTrees/merged/";
  TString histdir_top = "../Histograms/top/";

  vector<TString> processes_W = {"Data", "WMatched", "WUnmatched", "ZMatched", "ZUnmatched", "QCD", "Pseudo"};
  vector<TString> processes_top = {"Data", "TTbar", "SingleTop", "WJets", "other", "Pseudo"};

  vector<TFile*> files_W, files_top;

  if(fill_top){
    for(auto process: processes_top) fill_hists_top(histdir_top, process);
  }
  if(fill_W){
    for(auto process: processes_W) fill_hists_W(histdir_W, process , job_index , n_jobs);
  }


  outputFile->Close();
  return 0;
}


//------------------------------------------------------
//------------------------------------------------------
//------------------------------------------------------

void fill_hists_top(TString dir, TString process){
  cout << "filling hists for " << process << " (top) ..." << endl;
  TFile *file = new TFile(dir+process+".root");
  TTree *tree = (TTree *) file->Get("AnalysisTree");

  // creat a dummy hist for mjet
  TString dummyname = "top_"+process+"__mjet";
  TH1F* h_mjet_dummy = new TH1F(dummyname, "m_{jet} [GeV]", 50, 0, 500);

  // define pt bins
  vector<int> ptbins = {-1, 200, 300, 400, 100000};
  vector<TString> binnames = {"inclusive", "200to300", "300to400", "400toInf"};

  // create vectors of mjet hists
  vector<TH1F*> h_mjet_nominal_pass, h_mjet_nominal_fail;
  vector<vector<vector<TH1F*>>> h_mjet_vars_pass, h_mjet_vars_fail;
  h_mjet_vars_pass.resize(ptbins.size()-1);
  h_mjet_vars_fail.resize(ptbins.size()-1);

  for(int i=0; i<ptbins.size()-1; i++){
    TH1F* h1 = (TH1F*)h_mjet_dummy->Clone(dummyname+"_"+binnames[i]+"_pass");
    TH1F* h2 = (TH1F*)h_mjet_dummy->Clone(dummyname+"_"+binnames[i]+"_fail");
    h_mjet_nominal_pass.push_back(h1);
    h_mjet_nominal_fail.push_back(h2);
  }

  // declare variables
  double taucut = 0.5;
  double weight, mjet, pt, tau32;
  vector<vector<double>*> mjet_variations;
  mjet_variations.resize(handlenames.size());

  // assign to branches
  tree->ResetBranchAddresses();
  tree->SetBranchAddress("weight",&weight);
  tree->SetBranchAddress("mjet",&mjet);
  tree->SetBranchAddress("pt",&pt);
  tree->SetBranchAddress("tau32",&tau32);

  if(process != "Data"){
    for(unsigned int i=0; i<handlenames.size(); i++){
      tree->SetBranchAddress(handlenames[i],&mjet_variations[i]);
      for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
        TH1F* h_up_pass = (TH1F*)h_mjet_nominal_pass[ptbin]->Clone();
        TH1F* h_down_pass = (TH1F*)h_mjet_nominal_pass[ptbin]->Clone();
        TH1F* h_up_fail = (TH1F*)h_mjet_nominal_fail[ptbin]->Clone();
        TH1F* h_down_fail = (TH1F*)h_mjet_nominal_fail[ptbin]->Clone();
        TString oldtitle_pass = h_mjet_nominal_pass[ptbin]->GetName();
        TString newtitle_pass = oldtitle_pass.ReplaceAll("mjet", handlenames[i]);
        TString oldtitle_fail = h_mjet_nominal_fail[ptbin]->GetName();
        TString newtitle_fail = oldtitle_fail.ReplaceAll("mjet", handlenames[i]);
        h_up_pass->SetName(newtitle_pass+"__up");
        h_down_pass->SetName(newtitle_pass+"__down");
        h_up_fail->SetName(newtitle_fail+"__up");
        h_down_fail->SetName(newtitle_fail+"__down");
        h_mjet_vars_pass[ptbin].push_back({h_up_pass, h_down_pass});
        h_mjet_vars_fail[ptbin].push_back({h_up_fail, h_down_fail});
      }
    }
  }

  tree->SetBranchStatus("*",1);

  // loop over tree
  for(Int_t ievent=0; ievent < tree->GetEntriesFast(); ievent++) {
    if(tree->GetEntry(ievent)<=0) break;

    // loop over pt bins
    for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
      if(ptbins[ptbin] == -1 || (pt > ptbins[ptbin] && pt < ptbins[ptbin+1]) ){
        // fill nominal hists
        if(tau32<taucut) h_mjet_nominal_pass[ptbin]->Fill(mjet, weight);
        else             h_mjet_nominal_fail[ptbin]->Fill(mjet, weight);
        // mjet variations
        if(process != "Data"){
          for(int i=0; i<h_mjet_vars_pass[ptbin].size(); i++){
            for(int j=0; j<h_mjet_vars_pass[ptbin][i].size(); j++){
              if(tau32<taucut) h_mjet_vars_pass[ptbin][i][j]->Fill(mjet_variations[i]->at(j), weight);
              else             h_mjet_vars_fail[ptbin][i][j]->Fill(mjet_variations[i]->at(j), weight);
            }
          }
        }
      }
    }
  }

  // write hists
  outputFile->cd();
  for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
    h_mjet_nominal_pass[ptbin]->Write();
    h_mjet_nominal_fail[ptbin]->Write();
    for(int i=0; i<h_mjet_vars_pass[ptbin].size(); i++){
      for(int j=0; j<h_mjet_vars_pass[ptbin][i].size(); j++){
        h_mjet_vars_pass[ptbin][i][j]->Write();
        h_mjet_vars_fail[ptbin][i][j]->Write();
      }
    }
  }


  return;
}

//------------------------------------------------------
//------------------------------------------------------
//------------------------------------------------------

void fill_hists_W(TString dir, TString process, int job_index = 1 , int n_jobs = 1){
  // read root-file with k-factors for VJets samples
  if(process.Contains("W") || process.Contains("Z")){
    std::string NLOWeightsFilename = NLOWeightsDir + (std::string)(process.Contains("W") ? "/WJets" : "/ZJets") + "Corr.root";
    TFile * NLOWeightsFile = new TFile(NLOWeightsFilename.c_str());
    if(h_kfactor) h_kfactor->Reset();
    h_kfactor = (TH1F*) NLOWeightsFile->Get("kfactor");
    if(h_ewcorr) h_ewcorr->Reset();
    h_ewcorr = (TH1F*) NLOWeightsFile->Get("ewcorr");
  }

  TString process_file_name = TString(process);
  process_file_name.ReplaceAll("Matched","").ReplaceAll("Unmatched","");
  cout << "filling hists for " << process << " (W) ..." << endl;
  TFile *file = new TFile(dir+process_file_name+".root");
  TTree *tree = (TTree *) file->Get("AnalysisTree");

  // creat a dummy hist for mjet
  TString dummyname = "W_"+process+"__mjet";
  TH1F* h_mjet_dummy = new TH1F(dummyname, "m_{jet} [GeV]", 50, 0, 500);

  // define pt bins
  vector<int> ptbins = {-1, 500, 550, 600, 675, 800, 1200, 100000};
  vector<TString> binnames = {"inclusive", "500to550", "550to600", "600to675", "675to800", "800to1200", "1200toInf"};

  // create vectors of mjet hists
  vector<TH1F*> h_mjet_nominal_pass, h_mjet_nominal_fail;
  vector<vector<vector<TH1F*>>> h_mjet_vars_pass, h_mjet_vars_fail;
  h_mjet_vars_pass.resize(ptbins.size()-1);
  h_mjet_vars_fail.resize(ptbins.size()-1);

  for(int i=0; i<ptbins.size()-1; i++){
    TH1F* h1 = (TH1F*)h_mjet_dummy->Clone(dummyname+"_"+binnames[i]+"_pass");
    TH1F* h2 = (TH1F*)h_mjet_dummy->Clone(dummyname+"_"+binnames[i]+"_fail");
    h_mjet_nominal_pass.push_back(h1);
    h_mjet_nominal_fail.push_back(h2);
  }

  // declare variables
  double weight, mjet, pt, n2,genjet_V_pt;
  bool matchedV;
  vector<vector<double>*> mjet_variations;
  mjet_variations.resize(handlenames.size());

  // assign to branches
  tree->ResetBranchAddresses();
  tree->SetBranchAddress("weight",&weight);
  tree->SetBranchAddress("mjet",&mjet);
  tree->SetBranchAddress("pt",&pt);
  tree->SetBranchAddress("N2",&n2);
  tree->SetBranchAddress("matchedV",&matchedV);
  tree->SetBranchAddress("genjetpt",&genjet_V_pt);

  if(process != "Data"){
    for(unsigned int i=0; i<handlenames.size(); i++){
      tree->SetBranchAddress(handlenames[i],&mjet_variations[i]);
      for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
        TH1F* h_up_pass = (TH1F*)h_mjet_nominal_pass[ptbin]->Clone();
        TH1F* h_down_pass = (TH1F*)h_mjet_nominal_pass[ptbin]->Clone();
        TH1F* h_up_fail = (TH1F*)h_mjet_nominal_fail[ptbin]->Clone();
        TH1F* h_down_fail = (TH1F*)h_mjet_nominal_fail[ptbin]->Clone();
        TString oldtitle_pass = h_mjet_nominal_pass[ptbin]->GetName();
        TString newtitle_pass = oldtitle_pass.ReplaceAll("mjet", handlenames[i]);
        TString oldtitle_fail = h_mjet_nominal_fail[ptbin]->GetName();
        TString newtitle_fail = oldtitle_fail.ReplaceAll("mjet", handlenames[i]);
        h_up_pass->SetName(newtitle_pass+"__up");
        h_down_pass->SetName(newtitle_pass+"__down");
        h_up_fail->SetName(newtitle_fail+"__up");
        h_down_fail->SetName(newtitle_fail+"__down");
        h_mjet_vars_pass[ptbin].push_back({h_up_pass, h_down_pass});
        h_mjet_vars_fail[ptbin].push_back({h_up_fail, h_down_fail});

      }
    }
  }

  tree->SetBranchStatus("*",1);

  // loop over tree
  int n_events_tree = tree->GetEntriesFast();
  int events_per_job = (int)(n_events_tree / n_jobs);
  int start_event = n_jobs == 1 ? 0 : job_index * events_per_job;
  int end_event = n_jobs == 1 ? n_events_tree : start_event + events_per_job;
  end_event = end_event > n_events_tree ? n_events_tree : end_event;

  for(Int_t ievent=start_event; ievent < end_event; ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    if(n_jobs==1 && ievent % 10000 == 0) cout << "\r processing Event ("<< ievent << "/" << end_event << ")"<< flush;
    // loop over pt bins
    for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
      if(ptbins[ptbin] == -1 || (pt > ptbins[ptbin] && pt < ptbins[ptbin+1]) ){
        bool pass_N2ddt = passN2ddt(n2, pt, mjet);
        if(process.Contains("Matched") && (! matchedV)) continue;
        if(process.Contains("Unmatched") && (matchedV)) continue;
        if(process.Contains("W") || process.Contains("Z")){
          double kfactor = derive_kfactor(genjet_V_pt);
          weight *= kfactor;
        }
        // fill nominal hists
        if(pass_N2ddt) h_mjet_nominal_pass[ptbin]->Fill(mjet, weight);
        else                        h_mjet_nominal_fail[ptbin]->Fill(mjet, weight);

        // mjet variations
        if(process != "Data"){
          for(int i=0; i<h_mjet_vars_pass[ptbin].size(); i++){
            for(int j=0; j<h_mjet_vars_pass[ptbin][i].size(); j++){
              if(pass_N2ddt) h_mjet_vars_pass[ptbin][i][j]->Fill(mjet_variations[i]->at(j), weight);
              else                        h_mjet_vars_fail[ptbin][i][j]->Fill(mjet_variations[i]->at(j), weight);
            }
          }
        }
      }
    }
  }
  cout << endl;

  // write hists
  outputFile->cd();
  for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
    h_mjet_nominal_pass[ptbin]->Write();
    h_mjet_nominal_fail[ptbin]->Write();

    for(int i=0; i<h_mjet_vars_pass[ptbin].size(); i++){
      for(int j=0; j<h_mjet_vars_pass[ptbin][i].size(); j++){
        h_mjet_vars_pass[ptbin][i][j]->Write();
        h_mjet_vars_fail[ptbin][i][j]->Write();
      }
    }
  }
  return;
}

//------------------------------------------------------
//------------------------------------------------------
//------------------------------------------------------

void read_grid(TString gfilename){
  // read configuration from root file
  TFile* gfile = new TFile(gfilename);
  TH2F* grid = (TH2F*) gfile->Get("grid");
  TH1F* h_cat = (TH1F*) gfile->Get("categories");
  vector<TString> categories;
  for(int bin=1; bin<=h_cat->GetXaxis()->GetNbins(); bin++) categories.push_back(h_cat->GetXaxis()->GetBinLabel(bin));

  // get binnings and size of variation
  int Nbins_pt =  grid->GetXaxis()->GetNbins();
  int Nbins_eta = grid->GetYaxis()->GetNbins();
  int Nbins_cat = categories.size();

  // name and create all handles for mjet variations
  for(int i=0; i<Nbins_pt; i++){
    for(int j=0; j<Nbins_eta; j++){
      for(int k=0; k<Nbins_cat; k++){
        TString bin_name = "_" + to_string(i) + "_" + to_string(j) + "_" + categories[k];
        TString handlename = "mjet"+bin_name;
        handlenames.push_back(handlename);
      }
    }
  }
  return;
}

//------------------------------------------------------
//------------------------------------------------------
//------------------------------------------------------

bool passN2ddt(double n2, double pt, double mass){
  bool pass = false;

// deriving bin for pt and rho
  int pt_bin = ddtmap->GetYaxis()->FindFixBin(pt);
  if(pt_bin > ddtmap->GetYaxis()->GetNbins()){
    pt_bin = ddtmap->GetYaxis()->GetNbins();
  }else if(pt_bin <=0){
    pt_bin = 1;
  }

  double rho = 2 * TMath::Log(mass/pt);
  int rho_bin = ddtmap->GetXaxis()->FindFixBin(rho);
  if(rho_bin > ddtmap->GetXaxis()->GetNbins()){
    rho_bin = ddtmap->GetXaxis()->GetNbins();
  }else if(rho_bin <= 0){
    rho_bin = 1;
  }

  double N2ddt_map_value = ddtmap->GetBinContent(rho_bin,pt_bin);

  if(n2 < N2ddt_map_value) pass = true;

  return pass;
}

double derive_kfactor(double gen_pt){

  double kfactor_pt = gen_pt;
  double ewk_pt = gen_pt;

  if( kfactor_pt > 3000 ) kfactor_pt = 2800;
  if( kfactor_pt < 200 ) kfactor_pt = 205;

  float kfactor_bin = h_kfactor->GetXaxis()->FindBin(ewk_pt);

  float w= h_kfactor->GetBinContent(kfactor_bin);

  if( ewk_pt > 1205 ) ewk_pt = 1205;
  if( ewk_pt < 160 ) ewk_pt = 165;

  float ewk_bin = h_ewcorr->GetXaxis()->FindBin(ewk_pt);

  float w_ew= h_ewcorr->GetBinContent(ewk_bin);

  return w * w_ew;
}
