#include "../include/CentralInclude.h"

using namespace std;

void fill_hists_top(TString dir, TString process, int job_index_, int n_jobs_);
void fill_hists_WZ(TString dir, TString process, int job_index_, int n_jobs_, TString selection);
void fill_hists_WfromTop(TString dir, TString process, int job_index_, int n_jobs_);
float getN2ddt(double n2, double pt, double mass);
double derive_kfactor(double pt);

void read_grid(TString gfilename);

TFile* outputFile;
vector<TString> handlenames;
TH2F* ddtmap;
TH1F *h_kfactor, *h_ewcorr;
std::string NLOWeightsDir = "..//NLOweights";
bool apply_jec = true;
bool apply_jer = false;

// TString output_directory="../Histograms/";
TString output_directory="./";
std::map<TString,TString> skimmed_Trees = {};
// std::map<TString,TString> skimmed_Trees = {{"ZJets","/skimmedTrees/ZJets_ZToBB.root"}};
// std::map<TString,TString> skimmed_Trees = {{"ZJets","/skimmedTrees/ZJets_ZToLightQuarks.root"}};


int main(int argc, char* argv[]){

  bool fill_top = true;
  bool fill_W = true;
  bool fill_Zbb = true;
  if(argc > 1){
    if(strcmp(argv[1], "top") == 0){
      fill_W=false;
      fill_Zbb=false;
    }else if(strcmp(argv[1], "W") == 0){
      fill_top=false;
      fill_Zbb=false;
    }else if(strcmp(argv[1], "Zbb") == 0){
      fill_top=false;
      fill_W=false;
    }
  }


  int job_index = 1;
  int n_jobs = 1;
  if(argc > 2){
    job_index = atoi(argv[2]);
    n_jobs = atoi(argv[3]);
  }

  TString scaleStudy = "AllPFFlavours";  
  // TString scaleStudy = "PF_flavours";  
  read_grid("../Histograms/grid_"+scaleStudy+".root");
  // read_grid("../Histograms/grid.root");
  TString year = "2017";
  
  // read in ddtmap here
  // TString ddtMapFile = "../Histograms/ddtmaps/QCD_2017_PFMass_smooth_gaus4p00sigma.root";
  // TString ddtMapName = "N2_v_pT_v_rho_0p1_smooth_gaus4p00sigma_maps_cleaner_PFMass";
  TString ddtMapFile = "../Histograms/ddtmaps.root";
  TString ddtMapName = "n2ddt_map_"+year+"_smooth_4_0p05";

  TFile * mapFile = new TFile(ddtMapFile);
  ddtmap = (TH2F*)mapFile->Get(ddtMapName);

  TString outname;
  if(n_jobs > 1){
    outname.Form(output_directory+"/workdir/Histograms_%i.root",job_index);
    cout << outname << endl;
  }else{
    outname = output_directory+"/Histograms.root";
  }
  outputFile = new TFile(outname,"recreate");

  //this macro will attempt to create a TChain using the following histdirs+processes as patterns to any matching root-file. (no need to hadd root files resulting from sframe-batch)
  // TString histdir_W = "../Histograms/VJets/scaleStudy/"+scaleStudy+"/workdir_Vqq_scaleStudy/";
  // TString histdir_top = "../Histograms/TTBar/scaleStudy/"+scaleStudy+"/workdir_Top_scaleStudy/";
  TString histdir_W = "../Output/vjetsTrees/workdir_vjets_"+year+"/";
  // TString histdir_W = "../Histograms/VJets/noExtraBTagCounter/workdir_Vqq_jets_PreSel_0p005Variation/";
  TString histdir_top = "../Output/ttbarTrees/workdir_ttbar_"+year+"/";

  vector<TString> processes_W = {"Data", "WJetsMatched", "WJetsUnmatched", "ZJetsMatched", "ZJetsUnmatched", "TTToHadronic", "TTToSemiLeptonic", "ST_tW_top", "ST_tW_antitop", "QCD"};
  // vector<TString> processes_top = {"Data", "WJets", "DYJets", "TTbar_mergedTop", "TTbar_mergedW", "TTbar_mergedQB", "TTbar_semiMergedTop", "TTbar_notMerged", "ST_tch_top", "ST_tch_antitop","ST_tch","ST_tWch_top", "ST_tWch_antitop","ST_tWch", "ST_sch", "QCD"};
  vector<TString> processes_top = {"Data", "WJets", "DYJets", "TTbar_had", "TTbar_semilep_mergedTop", "TTbar_semilep_mergedW", "TTbar_semilep_mergedQB", "TTbar_semilep_semiMergedTop", "TTbar_semilep_notMerged", "TTbar_dilep", "ST_tch_top", "ST_tch_antitop","ST_tch","ST_tWch_top", "ST_tWch_antitop","ST_tWch", "ST_sch", "QCD"};

  // vector<TFile*> files_W, files_top;

  if(fill_top){
    for(auto process: processes_top) fill_hists_top(histdir_top, process , job_index , n_jobs);
  }
  if(fill_W){
    for(auto process: processes_W) fill_hists_WZ(histdir_W, process , job_index , n_jobs, "W");
  }
  if(fill_Zbb){
    for(auto process: processes_W) fill_hists_WZ(histdir_W, process , job_index , n_jobs, "Zbb");
  }
  
  outputFile->Close();
  return 0;
}


//------------------------------------------------------
//------------------------------------------------------
//------------------------------------------------------

void fill_hists_top(TString dir, TString process, int job_index = 1 , int n_jobs = 1){
  TString process_file_name = TString(process);
  // if(process_file_name.Contains("TTbar")) process_file_name = "TT";
  if(process_file_name.Contains("TTbar_semilep")) process_file_name = "TTbar_semilep";
  // if(process_file_name.Contains("TT_Mtt")) process_file_name = "TT_Mtt";
  process_file_name.ReplaceAll("Data","DATA");
  cout << "filling hists for " << process << " (top) ..." << endl;
  TChain *tree = new TChain("AnalysisTree");
  // std::cout << "adding files: "<< (dir + "*" + process_file_name + "*.root") << std::endl;
  tree->Add(dir + "*" + process_file_name + "*.root");
  
  // creat a dummy hist for mjet
  TString dummyname = "top_"+process+"__mjet";
  TH1F* h_mjet_dummy = new TH1F(dummyname, "m_{jet} [GeV]", 500, 0, 500);

  // define pt bins
  //use -2 and Placeholder as shown here to use multiple sets of pt bins
  // vector<int> ptbins = {-1, 200, 300, 400, 500, 600,-2,300,500,-2,200,250,300,350,400,450,500,550,600,100000,-2,300,100000,-2,400,100000,-2,500,100000};
  // vector<TString> binnames = {"inclusive", "200to300", "300to400", "400to500","500to600","PlaceHolder1","PlaceHolder2","300to500","PlaceHolder3","PlaceHolder4","200to250","250to300","300to350","350to400","400to450","450to500","500to550","550to600","600toInf","PlaceHolder5","PlaceHolder6","300toInf","PlaceHolder7","PlaceHolder8","400toInf","PlaceHolder9","PlaceHolder10","500toInf"};
  vector<int> ptbins = {-1, 200, 300, 400, 500, 650,100000};
  vector<TString> binnames = {"inclusive", "200to300", "300to400", "400to500","500to650","650toInf"};

  // create vectors of mjet hists
  vector<TH1F*> h_mjet_nominal_pass, h_mjet_nominal_passW, h_mjet_nominal_fail;
  vector<vector<vector<TH1F*>>> h_mjet_vars_pass, h_mjet_vars_passW, h_mjet_vars_fail;
  h_mjet_vars_pass.resize(ptbins.size()-1);
  h_mjet_vars_passW.resize(ptbins.size()-1);
  h_mjet_vars_fail.resize(ptbins.size()-1);

  for(int i=0; i<ptbins.size()-1; i++){
    TH1F* h1 = (TH1F*)h_mjet_dummy->Clone(dummyname+"_"+binnames[i]+"_pass");
    TH1F* h2 = (TH1F*)h_mjet_dummy->Clone(dummyname+"_"+binnames[i]+"_passW");
    TH1F* h3 = (TH1F*)h_mjet_dummy->Clone(dummyname+"_"+binnames[i]+"_fail");
    h_mjet_nominal_pass.push_back(h1);
    h_mjet_nominal_passW.push_back(h2);
    h_mjet_nominal_fail.push_back(h3);
  }

  vector<vector<TH1F*>> h_mjet_jer_var_pass, h_mjet_jer_var_passW, h_mjet_jer_var_fail;
  h_mjet_jer_var_pass.resize(ptbins.size()-1);
  h_mjet_jer_var_passW.resize(ptbins.size()-1);
  h_mjet_jer_var_fail.resize(ptbins.size()-1);

  // declare variables
  double tau32cut = 0.5;
  double tau21cut = 0.45;
  double weight, mjet, pt, tau32, tau21, jecfactor;
  bool IsMergedTop, IsMergedQB, IsMergedWZ, IsNotMerged;
  
  vector<vector<double>*> mjet_variations;
  mjet_variations.resize(handlenames.size());

  double jerfactor_nominal,jerfactor_up,jerfactor_down;
  // assign to branches
  tree->ResetBranchAddresses();
  tree->SetBranchAddress("weight",&weight);
  // tree->SetBranchAddress("weight_pre_ttbar_reweight",&weight);

  tree->SetBranchAddress("mjet",&mjet);
  tree->SetBranchAddress("pt",&pt);
  tree->SetBranchAddress("tau32",&tau32);
  tree->SetBranchAddress("tau21",&tau21);
  tree->SetBranchAddress("jecfactor_SD",&jecfactor);
  tree->SetBranchAddress("IsMergedTop",&IsMergedTop);
  tree->SetBranchAddress("IsMergedQB",&IsMergedQB);
  tree->SetBranchAddress("IsMergedWZ",&IsMergedWZ);
  tree->SetBranchAddress("IsNotMerged",&IsNotMerged);

  if(process != "Data"){
    //only stored for MC for now (17.03.21)
    tree->SetBranchAddress("jerfactor_SD_JEC_nominal",&jerfactor_nominal);
    tree->SetBranchAddress("jerfactor_SD_JEC_up",&jerfactor_up);
    tree->SetBranchAddress("jerfactor_SD_JEC_down",&jerfactor_down);
    
    for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
      TH1F* h_jer_up_pass = (TH1F*)h_mjet_nominal_pass[ptbin]->Clone();
      TH1F* h_jer_down_pass = (TH1F*)h_mjet_nominal_pass[ptbin]->Clone();
      TH1F* h_jer_up_passW = (TH1F*)h_mjet_nominal_passW[ptbin]->Clone();
      TH1F* h_jer_down_passW = (TH1F*)h_mjet_nominal_passW[ptbin]->Clone();
      TH1F* h_jer_up_fail = (TH1F*)h_mjet_nominal_fail[ptbin]->Clone();
      TH1F* h_jer_down_fail = (TH1F*)h_mjet_nominal_fail[ptbin]->Clone();
      TString oldtitle_pass = h_mjet_nominal_pass[ptbin]->GetName();
      TString oldtitle_passW = h_mjet_nominal_passW[ptbin]->GetName();
      TString oldtitle_fail = h_mjet_nominal_fail[ptbin]->GetName();
      
      h_jer_up_pass->SetName(oldtitle_pass+"_jer_up");
      h_jer_down_pass->SetName(oldtitle_pass+"_jer_down");
      h_jer_up_passW->SetName(oldtitle_passW+"_jer_up");
      h_jer_down_passW->SetName(oldtitle_passW+"_jer_down");
      h_jer_up_fail->SetName(oldtitle_fail+"_jer_up");
      h_jer_down_fail->SetName(oldtitle_fail+"_jer_down");
      h_mjet_jer_var_pass[ptbin] = {h_jer_up_pass, h_jer_down_pass};
      h_mjet_jer_var_passW[ptbin] = {h_jer_up_passW, h_jer_down_passW};
      h_mjet_jer_var_fail[ptbin] = {h_jer_up_fail, h_jer_down_fail};
    }
    
    
    
    // if(process != "Data"){
    for(unsigned int i=0; i<handlenames.size(); i++){
      tree->SetBranchAddress(handlenames[i],&mjet_variations[i]);
      for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
        TH1F* h_up_pass = (TH1F*)h_mjet_nominal_pass[ptbin]->Clone();
        TH1F* h_down_pass = (TH1F*)h_mjet_nominal_pass[ptbin]->Clone();
        TH1F* h_up_passW = (TH1F*)h_mjet_nominal_passW[ptbin]->Clone();
        TH1F* h_down_passW = (TH1F*)h_mjet_nominal_passW[ptbin]->Clone();
        TH1F* h_up_fail = (TH1F*)h_mjet_nominal_fail[ptbin]->Clone();
        TH1F* h_down_fail = (TH1F*)h_mjet_nominal_fail[ptbin]->Clone();
        TString oldtitle_pass = h_mjet_nominal_pass[ptbin]->GetName();
        TString newtitle_pass = oldtitle_pass.ReplaceAll("mjet", handlenames[i]);
        TString oldtitle_passW = h_mjet_nominal_passW[ptbin]->GetName();
        TString newtitle_passW = oldtitle_passW.ReplaceAll("mjet", handlenames[i]);
        TString oldtitle_fail = h_mjet_nominal_fail[ptbin]->GetName();
        TString newtitle_fail = oldtitle_fail.ReplaceAll("mjet", handlenames[i]);
        h_up_pass->SetName(newtitle_pass+"__up");
        h_down_pass->SetName(newtitle_pass+"__down");
        h_up_passW->SetName(newtitle_passW+"__up");
        h_down_passW->SetName(newtitle_passW+"__down");
        h_up_fail->SetName(newtitle_fail+"__up");
        h_down_fail->SetName(newtitle_fail+"__down");
        h_mjet_vars_pass[ptbin].push_back({h_up_pass, h_down_pass});
        h_mjet_vars_passW[ptbin].push_back({h_up_passW, h_down_passW});
        h_mjet_vars_fail[ptbin].push_back({h_up_fail, h_down_fail});
      }
    }
  }

  tree->SetBranchStatus("*",1);

  // loop over tree
  int n_events_tree = tree->GetEntries();
  int events_per_job = (int)(n_events_tree / n_jobs);
  int start_event = n_jobs == 1 ? 0 : job_index * events_per_job;
  int end_event = n_jobs == 1 ? n_events_tree : start_event + events_per_job;
  end_event = end_event > n_events_tree ? n_events_tree : end_event;

  for(Int_t ievent=start_event; ievent < end_event; ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    if(n_jobs==1 && ievent % 10000 == 0) cout << "\r processing Event ("<< ievent << "/" << end_event << ")"<< flush;

    vector<double> jer_variations = {1.,1.};
    if(process != "Data"){
      if(jerfactor_nominal >0){
        jer_variations = {jerfactor_up/jerfactor_nominal,jerfactor_down/jerfactor_nominal};
      }
    }
    

    if(process.Contains("mergedTop") && (!IsMergedTop)) continue;
    if(process.Contains("mergedW") && (!IsMergedWZ)) continue;
    if(process.Contains("mergedQB") && (!IsMergedQB)) continue;
    if(process.Contains("semiMergedTop") && ( !(IsMergedWZ || IsMergedQB) ) ) continue;
    // if(process.Contains("notMerged") && (!IsNotMerged) && !(!IsMergedTop && !IsMergedWZ && !IsMergedQB && !IsNotMerged) ) continue;
    if(process.Contains("notMerged") && (!IsNotMerged)) continue;

    double correction = 1.0;    
    if(apply_jec){
      correction *= jecfactor;
    }else{//if we don't want jec applied we have to undo this for the jet pt here:
      pt /= jecfactor;
    }
    if(!apply_jer)jerfactor_nominal = 1.0;
    if(process != "Data") correction *= jerfactor_nominal;
    
    // loop over pt bins
    for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
      if(ptbins[ptbin] == -1 || (pt > ptbins[ptbin] && pt < ptbins[ptbin+1]) ){
        // fill nominal hists
        if(tau32<tau32cut) h_mjet_nominal_pass[ptbin]->Fill(mjet*correction, weight);
        else if(tau21 < tau21cut) h_mjet_nominal_passW[ptbin]->Fill(mjet*correction, weight);
        else             h_mjet_nominal_fail[ptbin]->Fill(mjet*correction, weight);
        // mjet variations
        if(process != "Data"){
          //jer smearing/variations
          for(int i=0; i<jer_variations.size(); i++){
            if(tau32<tau32cut) h_mjet_jer_var_pass[ptbin][i]->Fill( mjet*jer_variations[i]*correction, weight);
            else if(tau21<tau21cut) h_mjet_jer_var_passW[ptbin][i]->Fill(mjet*jer_variations[i]*correction, weight);
            else             h_mjet_jer_var_fail[ptbin][i]->Fill(mjet*jer_variations[i]*correction, weight);
          }
          //constituent scale variations
          for(int i=0; i<h_mjet_vars_pass[ptbin].size(); i++){
            for(int j=0; j<h_mjet_vars_pass[ptbin][i].size(); j++){
              if(tau32<tau32cut) h_mjet_vars_pass[ptbin][i][j]->Fill(mjet_variations[i]->at(j)*correction, weight);
              else if(tau21<tau21cut) h_mjet_vars_passW[ptbin][i][j]->Fill(mjet_variations[i]->at(j)*correction, weight);
              else             h_mjet_vars_fail[ptbin][i][j]->Fill(mjet_variations[i]->at(j)*correction, weight);
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
    if(ptbins[ptbin] == -2 || ptbins[ptbin+1] == -2) continue;
    h_mjet_nominal_pass[ptbin]->Write();
    h_mjet_nominal_passW[ptbin]->Write();
    h_mjet_nominal_fail[ptbin]->Write();
    for(int i=0; i<h_mjet_vars_pass[ptbin].size(); i++){
      for(int j=0; j<h_mjet_vars_pass[ptbin][i].size(); j++){
        h_mjet_vars_pass[ptbin][i][j]->Write();
        h_mjet_vars_passW[ptbin][i][j]->Write();
        h_mjet_vars_fail[ptbin][i][j]->Write();
      }
    }
    if(process != "Data"){
      for(int i=0; i<h_mjet_jer_var_pass[ptbin].size(); i++){
        h_mjet_jer_var_pass[ptbin][i]->Write();
        h_mjet_jer_var_passW[ptbin][i]->Write();
        h_mjet_jer_var_fail[ptbin][i]->Write();
      }
    }
  }
  
  return;
}
  
//------------------------------------------------------
//------------------------------------------------------
//------------------------------------------------------

void fill_hists_WZ(TString dir, TString process, int job_index = 1 , int n_jobs = 1, TString selection="W"){
  // read root-file with k-factors for VJets samples
  if(process.Contains("WJets") || process.Contains("ZJets")){
    std::string NLOWeightsFilename = NLOWeightsDir + (std::string)(process.Contains("W") ? "/WJets" : "/ZJets") + "Corr.root";
    TFile * NLOWeightsFile = new TFile(NLOWeightsFilename.c_str());
    if(h_kfactor) h_kfactor->Reset();
    h_kfactor = (TH1F*) NLOWeightsFile->Get("kfactor");
    if(h_ewcorr) h_ewcorr->Reset();
    h_ewcorr = (TH1F*) NLOWeightsFile->Get("ewcorr");
  }

  TString process_file_name = TString(process);
  TString process_file_pattern = dir;

  process_file_name.ReplaceAll("Matched","").ReplaceAll("Unmatched","");//.ReplaceAll("Data","DATA");
  // if(process_file_pattern.Contains("2017"))process_file_name.ReplaceAll("Data","DATA");
  
  if( skimmed_Trees.count( process_file_name ) > 0 ){
    process_file_pattern += skimmed_Trees[process_file_name];
  }else{
    process_file_pattern += + "*" + process_file_name + "*.root";
  }

  cout << "filling hists for " << process << " ("<<selection<<") ..." << endl;
  TChain *tree = new TChain("AnalysisTree");
  cout << "adding files to TChain with pattern: " << process_file_pattern << endl;
  tree->Add(process_file_pattern);
  // creat a dummy hist for mjet
  TString dummyname = selection+"_"+process+"__mjet";
  TH1F* h_mjet_dummy = new TH1F(dummyname, "m_{jet} [GeV]", 500, 0, 500);

  // define pt bins
  // vector<int> ptbins = {-1, 500, 550, 600, 675, 800, 1200, 100000,-2,1200,1600,2000,3000,4000,5000,6000,100000};
  // vector<TString> binnames = {"inclusive", "500to550", "550to600", "600to675", "675to800", "800to1200", "1200toInf","Placeholder1","Placeholder2","1200to1600","1600to2000","2000to3000","3000to4000","4000to5000","5000to6000","6000toInf"};
  vector<int> ptbins = {-1, 500, 650, 800, 1200, 100000};
  vector<TString> binnames = {"inclusive", "500to650", "650to800", "800to1200", "1200toInf"};


  vector<double> nbjet_bins = {-1,0,1,100};
  vector<TString> nbjet_binnames = {"inclusive","eq0","gt0"};
  


  // create vectors of mjet hists
  vector<vector<TH1F*>> h_mjet_nominal_pass, h_mjet_nominal_fail;
  vector<vector<vector<vector<TH1F*>>>> h_mjet_vars_pass, h_mjet_vars_fail;
  h_mjet_nominal_pass.resize(ptbins.size()-1);
  h_mjet_nominal_fail.resize(ptbins.size()-1);
  h_mjet_vars_pass.resize(ptbins.size()-1);
  h_mjet_vars_fail.resize(ptbins.size()-1);


  
  // vector<vector<TH1F*>> h_mjet_jer_var_pass, h_mjet_jer_var_fail;
  vector<vector<vector<TH1F*>>> h_mjet_jer_var_pass, h_mjet_jer_var_fail;
  h_mjet_jer_var_pass.resize(ptbins.size()-1);
  h_mjet_jer_var_fail.resize(ptbins.size()-1);


  
  for(int i=0; i<ptbins.size()-1; i++){
    h_mjet_jer_var_pass[i].resize(nbjet_bins.size()-1);
    h_mjet_jer_var_fail[i].resize(nbjet_bins.size()-1);
    // h_mjet_nominal_pass[i].resize(nbjet_bins.size()-1);
    // h_mjet_nominal_fail[i].resize(nbjet_bins.size()-1);
    h_mjet_vars_pass[i].resize(nbjet_bins.size()-1);
    h_mjet_vars_fail[i].resize(nbjet_bins.size()-1);
    for(int nbjet_bin=0; nbjet_bin < nbjet_bins.size()-1; nbjet_bin++){
      TString nbjet_binname = nbjet_bin>0 ? ("_Nbjet"+nbjet_binnames[nbjet_bin]):"";
      TH1F* h1 = (TH1F*)h_mjet_dummy->Clone(dummyname+"_"+binnames[i]+nbjet_binname+"_pass");
      TH1F* h2 = (TH1F*)h_mjet_dummy->Clone(dummyname+"_"+binnames[i]+nbjet_binname+"_fail");
      h_mjet_nominal_pass[i].push_back(h1);
      h_mjet_nominal_fail[i].push_back(h2);
    }
  }

  //old setup with chf and chf bins
  // //define energy fraction bins
  // vector<double> chf_bins = {-1,0.0,0.6,1.0};
  // vector<TString> chf_binnames = {"inclusive","0to0p6","0p6to1p0"};
  // vector<double> nhf_bins = {-1,0.0,0.6,1.0};
  // vector<TString> nhf_binnames = {"inclusive","0to0p6","0p6to1p0"};

  // // create vectors of mjet hists
  // vector<vector<vector<TH1F*>>> h_mjet_nominal_pass, h_mjet_nominal_fail;
  // vector<vector<vector<vector<vector<TH1F*>>>>> h_mjet_vars_pass, h_mjet_vars_fail;
  // h_mjet_nominal_pass.resize(ptbins.size()-1);
  // h_mjet_nominal_fail.resize(ptbins.size()-1);
  // h_mjet_vars_pass.resize(ptbins.size()-1);
  // h_mjet_vars_fail.resize(ptbins.size()-1);
  // for(int i=0; i<ptbins.size()-1; i++){
  //   h_mjet_nominal_pass[i].resize(chf_bins.size()-1);
  //   h_mjet_nominal_fail[i].resize(chf_bins.size()-1);
  //   h_mjet_vars_pass[i].resize(chf_bins.size()-1);
  //   h_mjet_vars_fail[i].resize(chf_bins.size()-1);
  //   for(int chfbin=0; chfbin < chf_bins.size()-1; chfbin++){
  //     // h_mjet_nominal_pass[i][chfbin].resize(nhf_bins.size()-1);
  //     // h_mjet_nominal_fail[i][chfbin].resize(nhf_bins.size()-1);
  //     h_mjet_vars_pass[i][chfbin].resize(nhf_bins.size()-1);
  //     h_mjet_vars_fail[i][chfbin].resize(nhf_bins.size()-1);
  //     for(int nhfbin=0; nhfbin < nhf_bins.size()-1; nhfbin++){
  //       TString energy_fraction_binname = (chfbin>0?("_CHF"+chf_binnames[chfbin]):"")+(nhfbin>0?("_NHF"+nhf_binnames[nhfbin]):"");
  //       TH1F* h1 = (TH1F*)h_mjet_dummy->Clone(dummyname+"_"+binnames[i]+energy_fraction_binname+"_pass");
  //       TH1F* h2 = (TH1F*)h_mjet_dummy->Clone(dummyname+"_"+binnames[i]+energy_fraction_binname+"_fail");
  //       h_mjet_nominal_pass[i][chfbin].push_back(h1);
  //       h_mjet_nominal_fail[i][chfbin].push_back(h2);
  //     }
  //   }
  // }

  // declare variables
  double weight, mjet, pt, n2, genjet_V_pt, jecfactor;//, CHF, NHF;
  int nbjets;
  double zbb_tagger;
  // double ddb_hbb,ddb_qcd;
  // double zbb_tagger_wp = 0.70;
  // double zbb_tagger_wp = 0.86;
  double zbb_tagger_wp = 0.91;
  // +-----------------------------------+
  // |        MIDeepDoubleBHbbprob       |
  // +------+------+-------+------+------+
  // |  wp  | ZQQ  |  ZBB  | WQQ  | QCD  |
  // +------+------+-------+------+------+
  // | 0.70 | 6.62 | 17.62 | 3.85 | 2.43 |
  // | 0.86 | 3.12 | 10.73 | 1.12 | 0.63 |
  // | 0.89 | 2.20 |  8.01 | 0.68 | 0.35 |
  // | 0.91 | 1.54 |  5.86 | 0.43 | 0.20 |
  // | 0.92 | 1.21 |  4.69 | 0.32 | 0.14 |
  // +------+------+-------+------+------+
  bool matchedV;
  vector<vector<double>*> mjet_variations;
  mjet_variations.resize(handlenames.size());

  double jerfactor_nominal,jerfactor_up,jerfactor_down;

  // assign to branches
  tree->ResetBranchAddresses();
  tree->SetBranchAddress("weight",&weight);
  // tree->SetBranchAddress("weight_pre_ttbar_reweight",&weight);
  tree->SetBranchAddress("mjet",&mjet);
  tree->SetBranchAddress("pt",&pt);
  tree->SetBranchAddress("N2",&n2);
  // tree->SetBranchAddress("MDDeepBoostZbbQCD", &zbb_tagger);
  tree->SetBranchAddress("MIDeepDoubleBHbbprob", &zbb_tagger);
  // tree->SetBranchAddress("MIDeepDoubleBQCDprob", &ddb_qcd);
  // tree->SetBranchAddress("matchedV",&matchedV);
  tree->SetBranchAddress("IsMergedWZ",&matchedV);
  tree->SetBranchAddress("genjetpt",&genjet_V_pt);
  tree->SetBranchAddress("jecfactor_SD",&jecfactor);

  // tree->SetBranchAddress("CHF",&CHF);
  // tree->SetBranchAddress("NHF",&NHF);
  tree->SetBranchAddress("NextraMBtagDR0p8",&nbjets);
  
  if(process != "Data"){
    //only stored for MC for now (17.03.21)
    tree->SetBranchAddress("jerfactor_SD_JEC_nominal",&jerfactor_nominal);
    tree->SetBranchAddress("jerfactor_SD_JEC_up",&jerfactor_up);
    tree->SetBranchAddress("jerfactor_SD_JEC_down",&jerfactor_down);

    for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
      for(int nbjet_bin=0; nbjet_bin < nbjet_bins.size()-1; nbjet_bin++){
      
      TH1F* h_jer_up_pass = (TH1F*)h_mjet_nominal_pass[ptbin][nbjet_bin]->Clone();
      TH1F* h_jer_down_pass = (TH1F*)h_mjet_nominal_pass[ptbin][nbjet_bin]->Clone();
      TH1F* h_jer_up_fail = (TH1F*)h_mjet_nominal_fail[ptbin][nbjet_bin]->Clone();
      TH1F* h_jer_down_fail = (TH1F*)h_mjet_nominal_fail[ptbin][nbjet_bin]->Clone();      
      TString oldtitle_pass = h_mjet_nominal_pass[ptbin][nbjet_bin]->GetName();
      TString oldtitle_fail = h_mjet_nominal_fail[ptbin][nbjet_bin]->GetName();
      
      h_jer_up_pass->SetName(oldtitle_pass+"_jer_up");
      h_jer_down_pass->SetName(oldtitle_pass+"_jer_down");
      h_jer_up_fail->SetName(oldtitle_fail+"_jer_up");
      h_jer_down_fail->SetName(oldtitle_fail+"_jer_down");
      h_mjet_jer_var_pass[ptbin][nbjet_bin] = {h_jer_up_pass, h_jer_down_pass};
      h_mjet_jer_var_fail[ptbin][nbjet_bin] = {h_jer_up_fail, h_jer_down_fail};
      }
    }

    for(unsigned int i=0; i<handlenames.size(); i++){
      tree->SetBranchAddress(handlenames[i],&mjet_variations[i]);
      for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
        for(int nbjet_bin=0; nbjet_bin < nbjet_bins.size()-1; nbjet_bin++){
          TH1F* h_up_pass = (TH1F*)h_mjet_nominal_pass[ptbin][nbjet_bin]->Clone();
          TH1F* h_down_pass = (TH1F*)h_mjet_nominal_pass[ptbin][nbjet_bin]->Clone();
          TH1F* h_up_fail = (TH1F*)h_mjet_nominal_fail[ptbin][nbjet_bin]->Clone();
          TH1F* h_down_fail = (TH1F*)h_mjet_nominal_fail[ptbin][nbjet_bin]->Clone();
          TString oldtitle_pass = h_mjet_nominal_pass[ptbin][nbjet_bin]->GetName();
          TString newtitle_pass = oldtitle_pass.ReplaceAll("mjet", handlenames[i]);
          TString oldtitle_fail = h_mjet_nominal_fail[ptbin][nbjet_bin]->GetName();
          TString newtitle_fail = oldtitle_fail.ReplaceAll("mjet", handlenames[i]);
          h_up_pass->SetName(newtitle_pass+"__up");
          h_down_pass->SetName(newtitle_pass+"__down");
          h_up_fail->SetName(newtitle_fail+"__up");
          h_down_fail->SetName(newtitle_fail+"__down");
          h_mjet_vars_pass[ptbin][nbjet_bin].push_back({h_up_pass, h_down_pass});
          h_mjet_vars_fail[ptbin][nbjet_bin].push_back({h_up_fail, h_down_fail});
        }
      }
    }
    
    // old setup with chf and nhf bins
    // for(unsigned int i=0; i<handlenames.size(); i++){
    //   tree->SetBranchAddress(handlenames[i],&mjet_variations[i]);
    //   for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
    //     for(int chfbin=0; chfbin < chf_bins.size()-1; chfbin++){
    //       for(int nhfbin=0; nhfbin < nhf_bins.size()-1; nhfbin++){
    //         // TString energy_fraction_binname = (chfbin>0?("_CHF"+chf_binnames[chfbin]):"")+(nhfbin>0?("_NHF"+nhf_binnames[nhfbin]):"");
    //         TH1F* h_up_pass = (TH1F*)h_mjet_nominal_pass[ptbin][chfbin][nhfbin]->Clone();
    //         TH1F* h_down_pass = (TH1F*)h_mjet_nominal_pass[ptbin][chfbin][nhfbin]->Clone();
    //         TH1F* h_up_fail = (TH1F*)h_mjet_nominal_fail[ptbin][chfbin][nhfbin]->Clone();
    //         TH1F* h_down_fail = (TH1F*)h_mjet_nominal_fail[ptbin][chfbin][nhfbin]->Clone();
    //         TString oldtitle_pass = h_mjet_nominal_pass[ptbin][chfbin][nhfbin]->GetName();
    //         TString newtitle_pass = oldtitle_pass.ReplaceAll("mjet", handlenames[i]);
    //         TString oldtitle_fail = h_mjet_nominal_fail[ptbin][chfbin][nhfbin]->GetName();
    //         TString newtitle_fail = oldtitle_fail.ReplaceAll("mjet", handlenames[i]);
    //         h_up_pass->SetName(newtitle_pass+"__up");
    //         h_down_pass->SetName(newtitle_pass+"__down");
    //         h_up_fail->SetName(newtitle_fail+"__up");
    //         h_down_fail->SetName(newtitle_fail+"__down");
    //         h_mjet_vars_pass[ptbin][chfbin][nhfbin].push_back({h_up_pass, h_down_pass});
    //         h_mjet_vars_fail[ptbin][chfbin][nhfbin].push_back({h_up_fail, h_down_fail});
    //       }
    //     }
    //   }
    // }
  }
  
  tree->SetBranchStatus("*",1);

  // loop over tree
  int n_events_tree = tree->GetEntries();
  int events_per_job = (int)(n_events_tree / n_jobs);
  int start_event = n_jobs == 1 ? 0 : job_index * events_per_job;
  int end_event = n_jobs == 1 ? n_events_tree : start_event + events_per_job;
  end_event = end_event > n_events_tree ? n_events_tree : end_event;
  for(Int_t ievent=start_event; ievent < end_event; ievent++) {
    if(tree->GetEntry(ievent)<=0) break;
    if(n_jobs==1 && ievent % 10000 == 0) cout << "\r processing Event ("<< ievent << "/" << end_event << ")"<< flush;
    // loop over pt bins

    vector<double> jer_variations = {1.,1.};
    if(process != "Data"){
      if(jerfactor_nominal >0){
        jer_variations = {jerfactor_up/jerfactor_nominal,jerfactor_down/jerfactor_nominal};
      }
    }
    
    double correction = 1.0;    
    if(apply_jec){
      correction *= jecfactor;
      pt *= jecfactor;
    // }else{//if we don't want jec applied we have to undo this for the jet pt here:
      // pt /= jecfactor;
    }
    if((process.Contains("W") || process.Contains("Z")) && !process.Contains("tW")){
      double kfactor = derive_kfactor(genjet_V_pt);
      weight *= kfactor;
    }
    
    if(!apply_jer)jerfactor_nominal = 1.0;
    if(process != "Data") correction *= jerfactor_nominal;
    
    for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){      
      for(int nbjet_bin=0; nbjet_bin < nbjet_bins.size()-1; nbjet_bin++){
      // for(int chfbin=0; chfbin < chf_bins.size()-1; chfbin++){
      //   for(int nhfbin=0; nhfbin < nhf_bins.size()-1; nhfbin++){

        bool in_pt_bin = ptbins[ptbin] == -1 || (pt > ptbins[ptbin] && pt < ptbins[ptbin+1]);
        // bool in_chf_bin = chf_bins[chfbin] == -1 || (CHF >= chf_bins[chfbin] && CHF < chf_bins[chfbin+1]);
        // bool in_nhf_bin = nhf_bins[nhfbin] == -1 || (NHF >= nhf_bins[nhfbin] && NHF < nhf_bins[nhfbin+1]);          
        // if(in_pt_bin && in_chf_bin && in_nhf_bin){
        bool in_nbjet_bin = nbjet_bins[nbjet_bin] == -1 || (nbjets >= nbjet_bins[nbjet_bin] && nbjets < nbjet_bins[nbjet_bin+1]);
        if(in_pt_bin && in_nbjet_bin){
          
          bool pass_tagger = false;
          // zbb_tagger = ddb_hbb/(ddb_hbb+ddb_qcd);
          if(selection == "W") pass_tagger = getN2ddt(n2, pt, mjet) < 0;
          else if(selection == "Zbb") pass_tagger = zbb_tagger >= zbb_tagger_wp;
          
          double rho = 2 * TMath::Log(correction*mjet/pt);
          if(rho < -6.0 || rho > -2.1) continue;
          
          if(process.Contains("Matched") && (! matchedV)) continue;
          if(process.Contains("Unmatched") && (matchedV)) continue;
          // fill nominal hists
          // if(pass_tagger) h_mjet_nominal_pass[ptbin][chfbin][nhfbin]->Fill(mjet*correction, weight);
          // else h_mjet_nominal_fail[ptbin][chfbin][nhfbin]->Fill(mjet*correction, weight);
          // cout << "mjet:weight = " << mjet << ":"<< weight << endl;
          if(pass_tagger) h_mjet_nominal_pass[ptbin][nbjet_bin]->Fill(mjet*correction, weight);
          else h_mjet_nominal_fail[ptbin][nbjet_bin]->Fill(mjet*correction, weight);
          
          // mjet variations
          if(process != "Data"){
            // if(nbjet_bins[nbjet_bin] == -1){
              for(int i=0; i<jer_variations.size(); i++){
                if(pass_tagger) h_mjet_jer_var_pass[ptbin][nbjet_bin][i]->Fill( mjet*jer_variations[i]*correction, weight);
                else h_mjet_jer_var_fail[ptbin][nbjet_bin][i]->Fill(mjet*jer_variations[i]*correction, weight);
              }
            // }
            
            for(int i=0; i<h_mjet_vars_pass[ptbin][nbjet_bin].size(); i++){
              for(int j=0; j<h_mjet_vars_pass[ptbin][nbjet_bin][i].size(); j++){
                if(pass_tagger) h_mjet_vars_pass[ptbin][nbjet_bin][i][j]->Fill(mjet_variations[i]->at(j)*correction, weight);
                else h_mjet_vars_fail[ptbin][nbjet_bin][i][j]->Fill(mjet_variations[i]->at(j)*correction, weight);
              }
            }
            
            
            //fill jer variation hists. do this only for the inclusive CHF/NHF bin since we don't use those at the moment
            // if(chf_bins[chfbin] == -1 && nhf_bins[nhfbin] == -1){
            //   for(int i=0; i<jer_variations.size(); i++){
            //     if(pass_tagger) h_mjet_jer_var_pass[ptbin][i]->Fill( mjet*jer_variations[i]*correction, weight);
            //     else h_mjet_jer_var_fail[ptbin][i]->Fill(mjet*jer_variations[i]*correction, weight);
            //   }
            // }
            
            // for(int i=0; i<h_mjet_vars_pass[ptbin][chfbin][nhfbin].size(); i++){
            //   for(int j=0; j<h_mjet_vars_pass[ptbin][chfbin][nhfbin][i].size(); j++){
            //     if(pass_tagger) h_mjet_vars_pass[ptbin][chfbin][nhfbin][i][j]->Fill(mjet_variations[i]->at(j)*correction, weight);
            //     else h_mjet_vars_fail[ptbin][chfbin][nhfbin][i][j]->Fill(mjet_variations[i]->at(j)*correction, weight);
            //   }
            // }
          }
        }
      }
    }
  }
  
  cout << endl;
  // write hists
  outputFile->cd();
  cout << "writing hists"<<endl;

  for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
    for(int nbjet_bin=0; nbjet_bin < nbjet_bins.size()-1; nbjet_bin++){
      if(ptbins[ptbin] == -2 || ptbins[ptbin+1] == -2) continue;
      h_mjet_nominal_pass[ptbin][nbjet_bin]->Write();
      h_mjet_nominal_fail[ptbin][nbjet_bin]->Write();
      
      for(int i=0; i<h_mjet_vars_pass[ptbin][nbjet_bin].size(); i++){
        for(int j=0; j<h_mjet_vars_pass[ptbin][nbjet_bin][i].size(); j++){
          h_mjet_vars_pass[ptbin][nbjet_bin][i][j]->Write();
          h_mjet_vars_fail[ptbin][nbjet_bin][i][j]->Write();
        }
      }
    // }

  // for(int ptbin=0; ptbin<ptbins.size()-1; ptbin++){
  //   for(int chfbin=0; chfbin < chf_bins.size()-1; chfbin++){
  //     for(int nhfbin=0; nhfbin < nhf_bins.size()-1; nhfbin++){
  //       if(ptbins[ptbin] == -2 || ptbins[ptbin+1] == -2) continue;
  //       h_mjet_nominal_pass[ptbin][chfbin][nhfbin]->Write();
  //       h_mjet_nominal_fail[ptbin][chfbin][nhfbin]->Write();
        
  //       for(int i=0; i<h_mjet_vars_pass[ptbin][chfbin][nhfbin].size(); i++){
  //         for(int j=0; j<h_mjet_vars_pass[ptbin][chfbin][nhfbin][i].size(); j++){
  //           h_mjet_vars_pass[ptbin][chfbin][nhfbin][i][j]->Write();
  //           h_mjet_vars_fail[ptbin][chfbin][nhfbin][i][j]->Write();
  //         }
  //       }
  //     }
  //   }
    if(process != "Data"){ // again do this only for MC for now since JERis not used in Data
      for(int i=0; i<h_mjet_jer_var_pass[ptbin][nbjet_bin].size(); i++){        
        h_mjet_jer_var_pass[ptbin][nbjet_bin][i]->Write();
        h_mjet_jer_var_fail[ptbin][nbjet_bin][i]->Write();
      }
    }
    }}
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

float getN2ddt(double n2, double pt, double mass){
  bool pass = false;
  if(n2<0)return n2;
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
  
  return n2-N2ddt_map_value;
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
