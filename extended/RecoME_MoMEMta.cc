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
  if (p4.M() > 0)
    return;

  // Increase the energy until M is positive                                                                                                                                                                
  p4.SetE(p4.P());
  while (p4.M2() < 0) {
    double delta = p4.E() * 1e-5;
    p4.SetE(p4.E() + delta);
  };
}

class RecoME : public multidraw::TTreeFunction {
public:
  //Class Constructor 
  RecoME(char const* name);
  //Class Destructor 
  ~RecoME() {
  }
  //Functions from Multidraw namespace (TTreeFunction class)
  char const* getName() const override {return "RecoME"; }
  TTreeFunction* clone() const override {return new RecoME(name_.c_str());}
  unsigned getNdata() override {return 1; }
  //This function will return the required value
  double evaluate(unsigned) override;

protected:
  void bindTree_(multidraw::FunctionLibrary&) override;

  //name of the required ME
  std::string name_;

  //Needed variables to select the events
  UIntValueReader*  nCleanJet{};
  FloatArrayReader* CleanJet_pt{};
  FloatArrayReader* CleanJet_eta{};
  FloatArrayReader* CleanJet_phi{};

  IntArrayReader* Lepton_pdgId{};
  UIntValueReader*  nLepton{};
  FloatArrayReader* Lepton_pt{};
  FloatArrayReader* Lepton_eta{};
  FloatArrayReader* Lepton_phi{};

  FloatValueReader* MET_pt{};
  FloatValueReader* PuppiMET_pt{};
  FloatValueReader* PuppiMET_phi{};

private:

  Double_t LHCsqrts_= 13., mh_= 125.;
  
};

RecoME::RecoME(char const* name):
  TTreeFunction()
{
  name_ = name;
}

double
RecoME::evaluate(unsigned)
{
	
  //Initializing 4-vectors
  TLorentzVector L1(0.,0.,0.,0.);
  TLorentzVector L2(0.,0.,0.,0.);
  TLorentzVector LL(0.,0.,0.,0.);
  TLorentzVector NuNu(0.,0.,0.,0.);
  TLorentzVector Nu1(0.,0.,0.,0.);
  TLorentzVector Nu2(0.,0.,0.,0.);
  TLorentzVector W1(0.,0.,0.,0.);
  TLorentzVector W2(0.,0.,0.,0.);
  TLorentzVector Higgs(0.,0.,0.,0.);
  TLorentzVector J1(0.,0.,0.,0.);
  TLorentzVector J2(0.,0.,0.,0.);

  //Getting some values to select the events
  unsigned ncleanjet{*nCleanJet->Get()};
  unsigned nlep{*nLepton->Get()};
  float Pmet_pt{*PuppiMET_pt->Get()};
  float Pmet_phi{*PuppiMET_phi->Get()};

  //Conditions to select the event
  if(ncleanjet>=2 && nlep>1){
	 
    //STEP-1
    //4-vectors of the leptons
    //Select one electron and one muon
    int muons = 0;
    int electrons = 0;
    int lep1 = 0;
    int lep2 = 0;
	  
    // Loop over muons and electrons
    for (unsigned int ilep = 0; ilep<nlep; ilep++){
     if (abs(Lepton_pdgId->At(ilep)) == 13){
    	++muons;
    	if (muons == 1 && Lepton_pt->At(ilep) > 13){
    	  L1.SetPtEtaPhiM(Lepton_pt->At(ilep), Lepton_eta->At(ilep), Lepton_phi->At(ilep), 0.0); //Muon
        lep1 = Lepton_pdgId->At(ilep);
    	}
      }
      if (abs(Lepton_pdgId->At(ilep)) == 11){
    	++electrons;
    	if (electrons == 1 && Lepton_pt->At(ilep) > 13){
    	  L2.SetPtEtaPhiM(Lepton_pt->At(ilep), Lepton_eta->At(ilep), Lepton_phi->At(ilep), 0.0); //Electron
        lep2 = Lepton_pdgId->At(ilep);
    	}
     }
    }

    if (muons<1 || electrons<1){
      return -9999; //If there is not an electron and a muon, return -9999
    }
    
    LL = L1 + L2;
	  
    //Reconstructing Higgs 4 vector with MET
    double nunu_px = Pmet_pt*cos(Pmet_phi);
    double nunu_py = Pmet_pt*sin(Pmet_phi);
    double nunu_pz = LL.Pz();
    double nunu_m = 30.0; //Why 30? --> https://indico.cern.ch/event/850505/contributions/3593915/

    double nunu_e = sqrt(nunu_px*nunu_px + nunu_py*nunu_py + nunu_pz*nunu_pz + nunu_m*nunu_m);
    NuNu.SetPxPyPzE(nunu_px, nunu_py, nunu_pz, nunu_e);
    Higgs = LL + NuNu;

	  
    // Reconstructing the W 4 vector using an aproximation
    Nu1.SetPxPyPzE(-L2.Px(), -L2.Py(), 0, L2.Px()*L2.Px()+L2.Py()*L2.Py());
    Nu2.SetPxPyPzE(-L1.Px(), -L1.Py(), 0, L1.Px()*L1.Px()+L1.Py()*L1.Py());

    W1 = L1 + Nu1;  // if lep1 < 0 -> W1 is W+ pdgID=24                                                                                                                                                     
    W2 = L2 + Nu2;
	  
    //Selection for 2 jets
    int jetn = 0;
    for (unsigned int ijet = 0; ijet<ncleanjet; ijet++){
      if (CleanJet_pt->At(ijet)>30){ //Jet pt condition
	      ++jetn;
	      if (jetn==1) J1.SetPtEtaPhiM(CleanJet_pt->At(0), CleanJet_eta->At(0), CleanJet_phi->At(0), 0.0);
	      if (jetn==2) J2.SetPtEtaPhiM(CleanJet_pt->At(1), CleanJet_eta->At(1), CleanJet_phi->At(1), 0.0);
      }

    }

    if (jetn < 2) return -9999; //low number of jets

	  
    if(name_=="top"){
    	
	ConfigurationReader configuration("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/top_leptonic/top_leptonic.lua");
    	configuration.getGlobalParameters().set("top_mass", 173.);  
	MoMEMta weight(configuration.freeze());

    	logging::set_level(logging::level::off);

    	ParameterSet lua_parameters;
    	lua_parameters.set("USE_TF", true);
    	lua_parameters.set("USE_PERM", true);
	
    	momemta::Particle wminus { "wminus", LorentzVector(W1.Px(), W1.Py(), W1.Pz(), W1.E()), -24 }; // W-                                                                                                     
    	momemta::Particle wplus { "wplus", LorentzVector(W2.Px(), W2.Py(), W2.Pz(), W2.E()), 24 }; // W+                                                                                                        
	
    	if (lep1<0){
    	  momemta::Particle wminus { "wminus", LorentzVector(W2.Px(), W2.Py(), W2.Pz(), W2.E()), -24 }; // W-                                                                                                   
    	  momemta::Particle wplus { "wplus", LorentzVector(W1.Px(), W1.Py(), W1.Pz(), W1.E()), 24 }; // W+                                                                                                      
    	}
	
    	momemta::Particle bjet1 { "bjet1", LorentzVector(J1.Px(), J1.Py(), J1.Pz(), J1.E()), 5 }; // Not necessary a bjet, but passed to MoMEMta as if it is                                                    
    	momemta::Particle bjet2 { "bjet2", LorentzVector(J2.Px(), J2.Py(), J2.Pz(), J2.E()), -5 };
	
    	// normalize input for numerical estability                                                                                                                                                             
    	normalizeInput(wminus.p4);
    	normalizeInput(wplus.p4);
    	normalizeInput(bjet1.p4);
    	normalizeInput(bjet2.p4);
	
    	std::vector<std::pair<double, double>> weights = weight.computeWeights({wplus, bjet1, wminus, bjet2});
	
    	return (double)weights.back().first;
    
    }else if(name_=="top_leptons"){
    
    	ConfigurationReader configuration("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/TTbar_FullyLeptonic/TTbar_FullyLeptonic.lua");

    	if (lep1 < 0){
    	  configuration.getGlobalParameters().set("top_mass", 173.);
    	}else{
    	  ConfigurationReader configuration("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/TTbar_FullyLeptonic/TTbar_FullyLeptonic_mue.lua");
    	  configuration.getGlobalParameters().set("top_mass", 173.);
    	}
	    
	MoMEMta weight(configuration.freeze());

    	logging::set_level(logging::level::off);
	
    	ParameterSet lua_parameters;
    	lua_parameters.set("USE_TF", true);
    	lua_parameters.set("USE_PERM", true);
	
    	momemta::Particle lepton1 { "lepton1", LorentzVector(L1.Px(), L1.Py(), L1.Pz(), L1.E()), lep1 }; // muon                                                                                                
    	momemta::Particle lepton2 { "lepton2", LorentzVector(L2.Px(), L2.Py(), L2.Pz(), L2.E()), lep2 }; // electron                                                                                            
    	momemta::Particle bjet1 { "bjet1", LorentzVector(J1.Px(), J1.Py(), J1.Pz(), J1.E()), 5 }; // Not necessary a bjet, but passed to MoMEMta as if it is                                                    
    	momemta::Particle bjet2 { "bjet2", LorentzVector(J2.Px(), J2.Py(), J2.Pz(), J2.E()), -5 };
	
    	// normalize input for numerical estability                                                                                                                                                             
    	normalizeInput(lepton1.p4);
    	normalizeInput(lepton2.p4);
    	normalizeInput(bjet1.p4);
    	normalizeInput(bjet2.p4);
	
    	LorentzVector met_p4 {NuNu.Px(), NuNu.Py(), NuNu.Pz(), NuNu.E()};
	
    	std::vector<std::pair<double, double>> weights = weight.computeWeights({lepton1, bjet1, lepton2, bjet2}, met_p4);
	
    	return (double)weights.back().first;
    
    }else if(name_=="VBF"){
    
    	ConfigurationReader configuration("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/qqH_hww_ME/higgs_jets.lua");
    	MoMEMta weight(configuration.freeze());
	
    	logging::set_level(logging::level::off);
	
    	ParameterSet lua_parameters;
    	lua_parameters.set("USE_TF", true);
    	lua_parameters.set("USE_PERM", true);
	                                                                                         
    	momemta::Particle higgs { "higgs", LorentzVector(Higgs.Px(), Higgs.Py(), Higgs.Pz(), Higgs.E()), 25 }; // Higgs                                                                                         
    	momemta::Particle jet1 { "jet1", LorentzVector(J1.Px(), J1.Py(), J1.Pz(), J1.E()), 1 };                                                     
    	momemta::Particle jet2 { "jet2", LorentzVector(J2.Px(), J2.Py(), J2.Pz(), J2.E()), -1 };

	normalizeInput(higgs.p4);
    	normalizeInput(jet1.p4);
    	normalizeInput(jet2.p4);

    	std::vector<std::pair<double, double>> weights = weight.computeWeights({higgs, jet1, jet2});
	
    	return (double)weights.back().first;
    
    }else if(name_=="DY"){
    
    	ConfigurationReader configuration("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/DY_ME/DY_ME.lua");
    	MoMEMta weight(configuration.freeze());
	
    	logging::set_level(logging::level::off);
	
    	ParameterSet lua_parameters;
    	lua_parameters.set("USE_TF", true);
    	lua_parameters.set("USE_PERM", true);
	                                                                                         
    	momemta::Particle Z { "Z", LorentzVector(Higgs.Px(), Higgs.Py(), Higgs.Pz(), Higgs.E()), 23 }; // Z, same 4 vector as Higgs                                                                                                 
    	momemta::Particle jet1 { "jet1", LorentzVector(J1.Px(), J1.Py(), J1.Pz(), J1.E()), 1 };                                                      
    	momemta::Particle jet2 { "jet2", LorentzVector(J2.Px(), J2.Py(), J2.Pz(), J2.E()), -1 };
	
    	normalizeInput(Z.p4);
    	normalizeInput(jet1.p4);
    	normalizeInput(jet2.p4);
	
    	std::vector<std::pair<double, double>> weights = weight.computeWeights({Z, jet1, jet2});
	
    	return (double)weights.back().first;
    
    }else if(name_=="WW"){
    
    	ConfigurationReader configuration("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/WW_ME/ww_df.lua");

    	MoMEMta weight(configuration.freeze());
	
    	logging::set_level(logging::level::off);
	
    	ParameterSet lua_parameters;
    	lua_parameters.set("USE_TF", true);
    	lua_parameters.set("USE_PERM", true);
	
    	momemta::Particle wminus { "wminus", LorentzVector(W1.Px(), W1.Py(), W1.Pz(), W1.E()), -24 }; // W-                                                                                                     
    	momemta::Particle wplus { "wplus", LorentzVector(W2.Px(), W2.Py(), W2.Pz(), W2.E()), 24 }; // W+                                                                                                        
	
    	if (lep1<0){
    	  momemta::Particle wminus { "wminus", LorentzVector(W2.Px(), W2.Py(), W2.Pz(), W2.E()), -24 }; // W-                                                                                                   
    	  momemta::Particle wplus { "wplus", LorentzVector(W1.Px(), W1.Py(), W1.Pz(), W1.E()), 24 }; // W+                                                                                                      
    	}
	
    	momemta::Particle jet1 { "jet1", LorentzVector(J1.Px(), J1.Py(), J1.Pz(), J1.E()), 1 };                                                       
    	momemta::Particle jet2 { "jet2", LorentzVector(J2.Px(), J2.Py(), J2.Pz(), J2.E()), -1 };
	
    	// normalize input for numerical estability                                                                                                                                                             
    	normalizeInput(wminus.p4);
    	normalizeInput(wplus.p4);
    	normalizeInput(jet1.p4);
    	normalizeInput(jet2.p4);
	
    	std::vector<std::pair<double, double>> weights = weight.computeWeights({wplus, jet1, wminus, jet2});
	
    	return (double)weights.back().first;
    
    }else if(name_=="WW_leptons"){
    
    	ConfigurationReader configuration("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/WW_leptonic_ME/WW_leptonic.lua");

    	if (lep1 < 0){
    	  logging::set_level(logging::level::off);
    	}else{
    	  ConfigurationReader configuration("/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/WW_leptonic_ME/WW_leptonic_mue.lua");
    	}
	    
	MoMEMta weight(configuration.freeze());

    	logging::set_level(logging::level::off);
	
    	ParameterSet lua_parameters;
    	lua_parameters.set("USE_TF", true);
    	lua_parameters.set("USE_PERM", true);
	
    	momemta::Particle lepton1 { "lepton1", LorentzVector(L1.Px(), L1.Py(), L1.Pz(), L1.E()), lep1 }; // muon                                                                                                
    	momemta::Particle lepton2 { "lepton2", LorentzVector(L2.Px(), L2.Py(), L2.Pz(), L2.E()), lep2 }; // electron                                                                                            
    	momemta::Particle jet1 { "jet1", LorentzVector(J1.Px(), J1.Py(), J1.Pz(), J1.E()), 1 };                                                     
    	momemta::Particle jet2 { "jet2", LorentzVector(J2.Px(), J2.Py(), J2.Pz(), J2.E()), -1 };
	
    	// normalize input for numerical estability                                                                                                                                                             
    	normalizeInput(lepton1.p4);
    	normalizeInput(lepton2.p4);
    	normalizeInput(jet1.p4);
    	normalizeInput(jet2.p4);
	
    	LorentzVector met_p4 {NuNu.Px(), NuNu.Py(), NuNu.Pz(), NuNu.E()};
	
    	std::vector<std::pair<double, double>> weights = weight.computeWeights({lepton1, jet1, lepton2, jet2}, met_p4);
	
    	return (double)weights.back().first;
    
    }
    
  }
  //End if(nCleanJet>=2 && nLepton>1)
  else return -9999; 
}
void
RecoME::bindTree_(multidraw::FunctionLibrary& _library)
{
  //CleanJets
  _library.bindBranch(nCleanJet, "nCleanJet");
  _library.bindBranch(CleanJet_pt, "CleanJet_pt");
  _library.bindBranch(CleanJet_eta, "CleanJet_eta");
  _library.bindBranch(CleanJet_phi, "CleanJet_phi");
  //Leptons
  _library.bindBranch(Lepton_pdgId, "Lepton_pdgId");
  _library.bindBranch(nLepton, "nLepton");
  _library.bindBranch(Lepton_pt, "Lepton_pt");
  _library.bindBranch(Lepton_eta, "Lepton_eta");
  _library.bindBranch(Lepton_phi, "Lepton_phi");
  //MET
  _library.bindBranch(MET_pt, "MET_pt");
  _library.bindBranch(PuppiMET_pt, "PuppiMET_pt");
  _library.bindBranch(PuppiMET_phi, "PuppiMET_phi");
}
