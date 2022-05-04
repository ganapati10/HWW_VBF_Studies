#include "LatinoAnalysis/MultiDraw/interface/TTreeFunction.h"
#include "LatinoAnalysis/MultiDraw/interface/FunctionLibrary.h"
#include "TSystem.h"
#include "iostream"
#include "vector"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TSystem.h"
#include <map>
#include "TString.h"
#include "momemta/ConfigurationReader.h"
#include "momemta/MoMEMta.h"
#include "momemta/Types.h"

using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;

using namespace momemta;

void normalizeInput(LorentzVector& p4) {
  if (p4.M() > 0) return;

  // Increase the energy until M is positive                                                                                                                                                              
  p4.SetE(p4.P());
  while (p4.M2() < 0) {
    double delta = p4.E() * 1e-5;
    p4.SetE(p4.E() + delta);
  };
}


float recoMoMEMta(higgspx, higgspy, higgspz, higgse, jet1pt, jet2pt, jet1phi, jet2phi, jet1eta, jet2eta, path){
  
  logging::set_level(logging::level::off);
  
  //Initializing 4-vectors
  TLorentzVector J1(0.,0.,0.,0.);
  TLorentzVector J2(0.,0.,0.,0.);
  
  J1.SetPtEtaPhiM(jet1pt, jet1eta, jet1phi, 0.0);
  J2.SetPtEtaPhiM(jet2pt, jet2eta, jet2phi, 0.0);
  
  momemta::Particle higgs { "higgs", LorentzVector(higgspx, higgspy, higgspz, higgse), 25 }; // Higgs                                                                                                    
  momemta::Particle Z { "Z", LorentzVector(higgspx, higgspy, higgspz, higgse), 23 }; // Z, same 4 vector as Higgs                                                                             
  momemta::Particle jet1 { "jet1", LorentzVector(J1.Px(), J1.Py(), J1.Pz(), J1.E()), 1 };
  momemta::Particle jet2 { "jet2", LorentzVector(J2.Px(), J2.Py(), J2.Pz(), J2.E()), -1 };

  normalizeInput(higgs.p4);
  normalizeInput(Z.p4);
  normalizeInput(jet1.p4);
  normalizeInput(jet2.p4);
  
  logging::set_level(logging::level::off);

  // Higgs                                                                                                                                                                                                

  ConfigurationReader configuration_VBF("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/qqH_hww_ME/higgs_jets.lua");
  MoMEMta weight_VBF(configuration_VBF.freeze());

  // DY                                                                                                                                                                                                   

  ConfigurationReader configuration_DY("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/DY_ME/DY_ME.lua");
  MoMEMta weight_DY(configuration_DY.freeze());

  ParameterSet lua_parameters;
  lua_parameters.set("USE_TF", true);
  lua_parameters.set("USE_PERM", true);

  std::vector<std::pair<double, double>> weights_VBF = weight_VBF.computeWeights({higgs, jet1, jet2});
  std::vector<std::pair<double, double>> weights_DY = weight_DY.computeWeights({Z, jet1, jet2});

  double vbf = (double)weights_VBF.back().first;
  double dy = (double)weights_DY.back().first;

  return 150 * abs(vbf) / (150 * abs(vbf) + abs(dy));
  
  
}
