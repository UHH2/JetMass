#include "UHH2/JetMass/include/JetMassHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

JetMassHists::JetMassHists(uhh2::Context & ctx, const std::string & dirname, double variation_, TString mode): Hists(ctx, dirname){
  JetMassHists(ctx, dirname, mode);
  variation = variation_;
}


JetMassHists::JetMassHists(Context & ctx, const string & dirname, TString mode ): Hists(ctx, dirname){
  auto dataset_type = ctx.get("dataset_type");
  isMC = dataset_type == "MC";

  // read configuration from root file
  TString gfilename = ctx.get("GridFile");
  TFile* gfile = new TFile(gfilename);
  grid = (TH2F*) gfile->Get("grid");
  TH1F* h_cat = (TH1F*) gfile->Get("categories");
  for(int bin=1; bin<=h_cat->GetXaxis()->GetNbins(); bin++) categories.push_back(h_cat->GetXaxis()->GetBinLabel(bin));
  ConstructOtherIDs();
  if(categories.size() > all_cat.size()+1) throw runtime_error("you gave too many categories");
  if(categories.size() + otherIDs.size() -1 != 8) throw runtime_error("categories and size of 'other' does not match");

  // should SD be used
  use_SD = false;
  if(mode.Contains("SD")) use_SD = true;

  // use constituents?
  use_constituents = true;
  if(mode.Contains("NC")) use_constituents = false;

  // do variations?
  do_variations = false;
  if(mode.Contains("VAR")) do_variations = true;
  if(!use_constituents) do_variations = false;

  JEC_on_Mass = false;
  if(mode.Contains("JEC")) JEC_on_Mass=true;
  
  variation = 0.1;

  // get binnings and size of variation
  Nbins_pt =  grid->GetXaxis()->GetNbins();
  Nbins_eta = grid->GetYaxis()->GetNbins();
  Nbins_cat = categories.size();

  // setup jetmass hists
  TString xtitleMass = "jet mass";
  double minMass = 50;
  double maxMass = 250;
  int nMass = 100;
  TString xtitleRho = "rho";
  double minRho = -10;
  double maxRho = 0;
  int nRho = 150;

  //for 2D pt_v_mass hists
  TString titlePt = "p_{T}";
  double minPt = 500;
  double maxPt = 1200;
  int nPt = 100;

  // resize vectors of histograms
  h_mass_UP.resize(Nbins_pt);
  h_mass_DOWN.resize(Nbins_pt);
  h_pt_v_mass_UP.resize(Nbins_pt);
  h_pt_v_mass_DOWN.resize(Nbins_pt);
  for (int i = 0; i < Nbins_pt; i++){
    h_mass_UP[i].resize(Nbins_eta);
    h_mass_DOWN[i].resize(Nbins_eta);
    h_pt_v_mass_UP[i].resize(Nbins_eta);
    h_pt_v_mass_DOWN[i].resize(Nbins_eta);
    for (int j = 0; j < Nbins_eta; j++){
      h_mass_UP[i][j].resize(Nbins_cat);
      h_mass_DOWN[i][j].resize(Nbins_cat);
      h_pt_v_mass_UP[i][j].resize(Nbins_cat);
      h_pt_v_mass_DOWN[i][j].resize(Nbins_cat);
    }
  }

  // book some hists without variations
  h_mass = book<TH1F>("Mass_central", xtitleMass, nMass, minMass, maxMass);
  h_mass_jet = book<TH1F>("Mass_central_jet", xtitleMass, nMass, minMass, maxMass);
  h_mass_raw_subjets = book<TH1F>("Mass_central_raw_subjets", xtitleMass, nMass, minMass, maxMass);
  h_rho = book<TH1F>("Rho_central", xtitleRho, nRho, minRho, maxRho);
  h_rho_jet = book<TH1F>("Rho_central_jet", xtitleRho, nRho, minRho, maxRho);
  h_particle_pt = book<TH1F>("particle_pt", "p_T", 50, 0, 50);
  h_particle_eta = book<TH1F>("particle_eta", "eta", 50, -5, 5);
  h_weights = book<TH1F>("weights", "weight", 100, -1, 2);
  h_pt_v_mass = book<TH2F>("pt_v_mass","x="+xtitleMass+",y="+titlePt, nMass, minMass, maxMass, nPt, minPt, maxPt);

  // book hists for variations
  for(int i=0; i<Nbins_pt; i++){
    for(int j=0; j<Nbins_eta; j++){
      for(int k=0; k<Nbins_cat; k++){
        TString mass_name = "Mass";
        TString bin_name = "_" + to_string(i) + "_" + to_string(j) + "_" + categories[k];
        h_mass_UP[i][j][k] = book<TH1F>(mass_name+bin_name+"_up", xtitleMass, nMass, minMass, maxMass);
        h_mass_DOWN[i][j][k] = book<TH1F>(mass_name+bin_name+"_down", xtitleMass, nMass, minMass, maxMass);
        h_pt_v_mass_UP[i][j][k] = book<TH2F>("Pt_v_"+mass_name+bin_name+"_up", "x="+xtitleMass+",y="+titlePt, nMass, minMass, maxMass, nPt, minPt, maxPt);
        h_pt_v_mass_DOWN[i][j][k] = book<TH2F>("Pt_v_"+mass_name+bin_name+"_down", "x="+xtitleMass+",y="+titlePt, nMass, minMass, maxMass, nPt, minPt, maxPt);
      }
    }
  }

}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// fill Hists
void JetMassHists::fill(const Event & event){
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  vector<TopJet>* topjets = event.topjets;
  if(topjets->size() < 1) return;
	vector<PFParticle>* allparticles;
	if(use_constituents)allparticles = event.pfparticles;
	else allparticles = new std::vector<PFParticle>();

  h_weights->Fill(weight);

  // get the PFParticles inside the first topjet
  // if mode is soft drop, get particles from subjets
  vector<PFParticle> particles;
  LorentzVector sum_subjets_v4;  
  float JEC_factor_raw = topjets->at(0).JEC_factor_raw();
  if(allparticles->size() != 0 && use_constituents){
    if(use_SD){
      vector<Jet> subjets = topjets->at(0).subjets();
      for(auto subjet: subjets){
        LorentzVectorXYZE old_v4XYZE = toXYZ(subjet.v4());
        sum_subjets_v4 += toPtEtaPhi(old_v4XYZE * JEC_factor_raw);
        for(const auto candInd : subjet.pfcand_indexs()){
          particles.push_back(allparticles->at(candInd));
        }
      }
    }
    else{
      for(const auto candInd : topjets->at(0).pfcand_indexs()){
        particles.push_back(allparticles->at(candInd));
      }
    }
  }
  // fill central histograms
  double mjet;
  if(use_SD) mjet = topjets->at(0).softdropmass();
  else       mjet = topjets->at(0).v4().M();
  double ptjet = topjets->at(0).v4().Pt();
  double rhojet = 2*TMath::Log(mjet/ptjet);
  h_mass_jet->Fill(mjet, weight);
  h_rho_jet->Fill(rhojet, weight);
  // fill again but get values from PF particles
  double mass, rho;
  if(use_constituents){
    float JEC_factor = JEC_on_Mass ? 1./topjets->at(0).JEC_factor_raw() : 1.0;
    mass = CalculateMJet(particles)*JEC_factor;
    rho = 2*TMath::Log(mass/ptjet);
  }
  else{
    mass = mjet;
    rho = rhojet;
  }
  h_mass->Fill(mass, weight);
  h_rho->Fill(rho, weight);
  h_pt_v_mass->Fill(mass, ptjet, weight);

  h_mass_raw_subjets->Fill(sum_subjets_v4.M(), weight);

  // fill some particle histograms
  if(use_constituents){
    for(auto p: particles){
      h_particle_pt->Fill(p.pt(), weight);
      h_particle_eta->Fill(p.eta(), weight);
    }
  }

  // loop over every bin in pt and eta
  // in every bin, vary the grid, apply to pf particles, compute jetmas
  // and fill up/down histograms
  if(do_variations){
    for(int i=0; i<Nbins_pt; i++){
      for(int j=0; j<Nbins_eta; j++){
        for(int k=0; k<Nbins_cat; k++){
          int ptbin = i+1;
          int etabin = j+1;
          if(isMC){
            vector<double> mjet_var = CalculateMJetVariation(particles, ptbin, etabin, categories[k]);
            h_mass_UP[i][j][k]->Fill(mjet_var[0], weight);
            h_mass_DOWN[i][j][k]->Fill(mjet_var[1], weight);
            h_pt_v_mass_UP[i][j][k]->Fill(mjet_var[0], ptjet, weight);
            h_pt_v_mass_DOWN[i][j][k]->Fill(mjet_var[1], ptjet, weight);
          }
          else{
            h_mass_UP[i][j][k]->Fill(mass, weight);
            h_mass_DOWN[i][j][k]->Fill(mass, weight);
            h_pt_v_mass_UP[i][j][k]->Fill(mass, ptjet, weight);
            h_pt_v_mass_DOWN[i][j][k]->Fill(mass, ptjet, weight);
          }
        }
      }
    }
  }
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// this applies the factor to PF particles according to the grid
vector<PFParticle> JetMassHists::VaryParticles(vector<PFParticle> oldParticles, int i, int j, TString cat, TString var){
  vector<PFParticle> newParticles;
  for(auto p:oldParticles){
    double pt = p.pt();
    double eta = fabs(p.eta());
    int ptbin = grid->GetXaxis()->FindBin(pt);
    int etabin = grid->GetYaxis()->FindBin(eta);
    // if bin is overflow, set manually to last bin
    if(ptbin > Nbins_pt) ptbin = Nbins_pt;
    if(etabin > Nbins_eta) etabin = Nbins_eta;
    double factor = 1.0;
    if(isMC && ptbin == i && etabin == j && inCategory(p,cat)){
      if(var == "up") factor += variation;
      else            factor -= variation;
    }
    PFParticle newP;
    newP.set_v4(factor * p.v4());
    newParticles.push_back(newP);
  }
  return newParticles;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// return varied mjet as {mjet_up, mjet_down}
vector<double> JetMassHists::CalculateMJetVariation(vector<PFParticle>particles, int i, int j, TString cat){
  LorentzVector up, down;
  for(auto p:particles){
    double pt = p.pt();
    double eta = fabs(p.eta());
    int ptbin = grid->GetXaxis()->FindBin(pt);
    int etabin = grid->GetYaxis()->FindBin(eta);
    // if bin is overflow, set manually to last bin
    if(ptbin > Nbins_pt) ptbin = Nbins_pt;
    if(etabin > Nbins_eta) etabin = Nbins_eta;
    double factor = 1.0;
    if(isMC && ptbin == i && etabin == j && inCategory(p,cat)){
      up += (factor+variation) * p.v4();
      down += (factor-variation) * p.v4();
    }
    else{
      up += p.v4();
      down += p.v4();
    }
  }
  return (vector<double>) {up.M(), down.M()};
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// check if particle is in category
bool JetMassHists::inCategory(PFParticle p, TString cat){
  bool inCat = false;
  int id = p.particleID();
  for(unsigned int i=0; i<all_cat.size(); i++){
    if(cat == all_cat[i] && abs(id) == i) inCat = true;
  }
  if(cat = "other"){
    for(auto oID: otherIDs) if(id == oID) inCat = true;
  }
  return inCat;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// this function creates a vector containing all IDs that are not covered
void JetMassHists::ConstructOtherIDs(){
  vector<bool> AddIDToOther(8, true);
  for(unsigned int i=0; i<all_cat.size(); i++){
    for(auto cat: categories){
      if(cat == all_cat[i]) AddIDToOther[i] = false;
    }
  }
  for(unsigned int id=0; id<AddIDToOther.size(); id++){
    if(AddIDToOther[id]) otherIDs.push_back(id);
  }
  return;
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// calculate jet mass from vector of PF particles
double JetMassHists::CalculateMJet(vector<PFParticle> Particles){
  LorentzVector jet_v4;
  for(auto p:Particles){
    jet_v4 += p.v4();
  }
  double mjet = jet_v4.M();
  return mjet;
}


JetMassHists::~JetMassHists(){}
