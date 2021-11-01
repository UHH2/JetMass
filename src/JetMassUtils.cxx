#include "UHH2/JetMass/include/JetMassUtils.h"

#include <TFile.h>
#include <TH1F.h>


// const GenParticle * hardestParticle(uhh2::Event & event, const GenTopJet & jet){
//   double pt_max(-1.0);
//   GenParticle * hardest(NULL);
//   for(auto subjet: jet.subjets()){
//     for(auto igenpart: subjet.genparticles_indices()){
//       GenParticle *p = &event.genparticles->at(igenpart);
//       double pt(p->pt());
//       if(pt>pt_max){
//         pt_max = pt;
//         hardest = p;
//         std::cout << igenpart << " " << pt << " " << p->pdgId()  << " " << p->status() << std::endl;
//       }
//     }
//   }
//   return hardest;
// }

double calculateGenTopPtSum(uhh2::Event & event){
  double genht(-1.0);
  if(event.gentopjets->size() == 0)return genht;
  for(auto jet : *event.gentopjets){
    genht += jet.pt();
  }
  return genht;
}

double calculateGenJetPtSum(uhh2::Event & event){
  double genht(-1.0);
  assert(event.genjets);
  if(event.genjets->size() == 0)return genht;
  for(auto jet : *event.genjets){
    genht += jet.pt();
  }
  return genht;
}

int get_V_index(std::vector<GenParticle> particles){
  int v_index(-1);
  for(unsigned int i=0; i<particles.size();i++){
    GenParticle p = particles.at(i);
    if( (abs(p.pdgId()) == 23) || (abs(p.pdgId()) == 24) ) v_index = i ;
  }
  return v_index;
}


template<>
bool JetSelector<TopJet>::process(uhh2::Event & event){
  assert(event.topjets);
  if(event.topjets->size() == 0) return false;
  const TopJet * selected_jet = &(event.topjets->at(0));
  // float n2_min = 1000.;
  // for(unsigned int i=0; i<event.topjets->size();i++){
  //   const TopJet* jet = &(event.topjets->at(i));
  //   float this_n2 = jet->ecfN2_beta1();
  //   if(this_n2<0) continue;
  //   if(this_n2 < n2_min){
  //     n2_min = this_n2;
  //     selected_jet = jet;
  //   }
  // }
  event.set(jet_handle,selected_jet);
  return true;
}
template<>
bool JetSelector<GenTopJet>::process(uhh2::Event & event){
  assert(event.gentopjets);
  if(event.gentopjets->size() == 0) return false;
  const GenTopJet * selected_jet = &(event.gentopjets->at(0));
  // float n2_min = 1000.;
  // for(unsigned int i=0; i<event.gentopjets->size();i++){
  //   const GenTopJet* jet = &(event.gentopjets->at(i));
  //   float this_n2 = jet->ecfN2_beta1();
  //   if(this_n2<0) continue;
  //   if(this_n2 < n2_min){
  //     n2_min = this_n2;
  //     selected_jet = jet;
  //   }
  // }
  event.set(jet_handle,selected_jet);
  return true;
}


NLOWeights::NLOWeights(uhh2::Context & ctx, const std::string & genjet_handlename){
  TString version = ctx.get("dataset_version");
  TString selection = ctx.get("selection");
  apply_nloweights = false;
  if(selection == "vjets" && (version.Contains("WJets") || version.Contains("ZJets")) ){
    genjet_handle = ctx.get_handle<const GenTopJet*>(genjet_handlename);
    // recojet_handle = ctx.get_handle<const TopJet*>(recojet_handlename);
    
    std::string NLOWeightsFilename =  "JetMass/NLOweights" + (std::string)(version.Contains("W") ? "/WJets" : "/ZJets") + "Corr.root";
    
    TFile * NLOWeightsFile = new TFile(locate_file(NLOWeightsFilename).c_str());
    h_kfactor = (TH1F*) NLOWeightsFile->Get("kfactor");
    h_ewcorr = (TH1F*) NLOWeightsFile->Get("ewcorr");
      
    apply_nloweights = true;
  }  
}

bool NLOWeights::process(uhh2::Event & event){
  if(!apply_nloweights) return false;

  if(!event.is_valid(genjet_handle))throw std::runtime_error("NLOWeights: genjet handle is invalid!");
  const GenTopJet *genjet(NULL);
  genjet = event.get(genjet_handle);
  if(genjet == nullptr)throw std::runtime_error("NLOWeights: genjet from handle is nullptr! ");

  float genjetpt = genjet->pt();
  double kfactor_pt = genjetpt;
  double ewk_pt = genjetpt;
  
  if( kfactor_pt > 3000 ) kfactor_pt = 2800;
  if( kfactor_pt < 200 ) kfactor_pt = 205;
  
  float kfactor_bin = h_kfactor->GetXaxis()->FindBin(kfactor_pt);
  
  float w= h_kfactor->GetBinContent(kfactor_bin);
      
  if( ewk_pt > 1205 ) ewk_pt = 1205;
  if( ewk_pt < 160 ) ewk_pt = 165;
  
  float ewk_bin = h_ewcorr->GetXaxis()->FindBin(ewk_pt);
  
  float w_ew= h_ewcorr->GetBinContent(ewk_bin);
  float nlo_weight = w * w_ew;
  
  event.weight *= nlo_weight;

  return true;
}


template<typename GenericTopJet>
N2DDTComputer<GenericTopJet>::N2DDTComputer(const std::string & map_file_path, const std::string map_name){
  TFile * map_file = new TFile(locate_file(map_file_path).c_str());
  h_n2ddt_map.reset((TH2D*) map_file->Get(map_name.c_str()));
}

template<typename GenericTopJet>
double N2DDTComputer<GenericTopJet>::computeDDTValue(const GenericTopJet * jet){
  if(jet == nullptr) throw std::runtime_error("N2DDTComputer: passed jet is nullptr!");
  double n2 = jet->ecfN2_beta1();
  double pt = jet->pt();
  double mass = jet->softdropmass();
  
  // deriving bin for pt and rho
  int pt_bin = h_n2ddt_map->GetYaxis()->FindFixBin(pt);
  if(pt_bin > h_n2ddt_map->GetYaxis()->GetNbins()){
    pt_bin = h_n2ddt_map->GetYaxis()->GetNbins();
  }else if(pt_bin <=0){
    pt_bin = 1;
  }

  double rho = 2 * TMath::Log(mass/pt);
  int rho_bin = h_n2ddt_map->GetXaxis()->FindFixBin(rho);
  if(rho_bin > h_n2ddt_map->GetXaxis()->GetNbins()){
    rho_bin = h_n2ddt_map->GetXaxis()->GetNbins();
  }else if(rho_bin <= 0){
    rho_bin = 1;
  }

  double N2ddt_map_value = h_n2ddt_map->GetBinContent(rho_bin,pt_bin);
  
  return n2-N2ddt_map_value;
}
