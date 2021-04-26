#include "LatinoAnalysis/MultiDraw/interface/TTreeFunction.h"
#include "LatinoAnalysis/MultiDraw/interface/FunctionLibrary.h"
#include "TSystem.h"
#include "iostream"
#include "vector"
#include "TLorentzVector.h"
#include "TMath.h"
//#include "ZZMatrixElement/MELA/interface/Mela.h" //Use ZZMatrixElement or JHUGen(recommended)
#include "JHUGenMELA/MELA/interface/Mela.h"
#include "TSystem.h"
#include <map>
#include "TString.h"

class RecoLevelME : public multidraw::TTreeFunction {
public:
  //Class Constructor 
  RecoLevelME(char const* name);
  //Class Destructor 
  ~RecoLevelME() {
  }
  //Functions from Multidraw namespace (TTreeFunction class)
  char const* getName() const override {return "RecoLevelME"; }
  TTreeFunction* clone() const override {return new RecoLevelME(name_.c_str());}
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

  UIntValueReader*  nSubJet{};
  FloatArrayReader* SubJet_pt{};
  FloatArrayReader* SubJet_eta{};
  FloatArrayReader* SubJet_mass{};
  FloatArrayReader* SubJet_phi{};

private:

  Double_t LHCsqrts_= 13., mh_= 125.;
  TVar::VerbosityLevel verbosity_ = TVar::SILENT;
  
  static Mela* mela;

};
Mela* RecoLevelME :: mela = 0;

RecoLevelME::RecoLevelME(char const* name):
  TTreeFunction()
{
  name_ = name;
  if(mela == 0){
    mela = new Mela(LHCsqrts_, mh_, verbosity_);
  }
}

double
RecoLevelME::evaluate(unsigned)
{
  //Map to store the ME
  std::map<TString, float> MatrixElementsMap;

  //Initializing 4-vectors
  TLorentzVector L1(0.,0.,0.,0.);
  TLorentzVector L2(0.,0.,0.,0.);
  TLorentzVector LL(0.,0.,0.,0.);
  TLorentzVector NuNu(0.,0.,0.,0.);
  TLorentzVector Higgs(0.,0.,0.,0.);
  TLorentzVector J1(0.,0.,0.,0.);
  TLorentzVector J2(0.,0.,0.,0.);
  TLorentzVector J3(0.,0.,0.,0.);

  //Getting some values to select the events
  unsigned ncleanjet{*nCleanJet->Get()};
  unsigned nlep{*nLepton->Get()};
  unsigned nsubjet{*nSubJet->Get()};
  float Pmet_pt{*PuppiMET_pt->Get()};
  float Pmet_phi{*PuppiMET_phi->Get()};

  //Conditions to select the event
  if(ncleanjet>=2 && nlep>1){
	 
    //STEP-1
    //4-vectors of the leptons
    //Select one electron and one muon
    int muons = 0;
    int electrons = 0;
    for (unsigned int ilep = 0; ilep<nlep; ilep++){
     if (abs(Lepton_pdgId->At(ilep)) == 13){
    	++muons;
    	if (muons == 1 && Lepton_pt->At(ilep) > 13){
    	  L1.SetPtEtaPhiM(Lepton_pt->At(ilep), Lepton_eta->At(ilep), Lepton_phi->At(ilep), 0.0); //Muon
    	}
      }
      if (abs(Lepton_pdgId->At(ilep)) == 11){
    	++electrons;
    	if (electrons == 1 && Lepton_pt->At(ilep) > 13){
    	  L2.SetPtEtaPhiM(Lepton_pt->At(ilep), Lepton_eta->At(ilep), Lepton_phi->At(ilep), 0.0); //Electron
    	}
     }
    }

    if (muons<1 || electrons<1){
      return -9999; //If there is not an electron and a muon
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

    //Selection for 2 jets
    int jetn = 0;
    bool use3jet = false;
    for (unsigned int ijet = 0; ijet<ncleanjet; ijet++){

      if (CleanJet_pt->At(ijet)>30){ //Jet pt condition
	++jetn;
	if (jetn==1) J1.SetPtEtaPhiM(CleanJet_pt->At(0), CleanJet_eta->At(0), CleanJet_phi->At(0), 0.0);
	if (jetn==2) J2.SetPtEtaPhiM(CleanJet_pt->At(1), CleanJet_eta->At(1), CleanJet_phi->At(1), 0.0);
	//if (jetn==3) J3.SetPtEtaPhiM(CleanJet_pt->At(2), CleanJet_eta->At(2), CleanJet_phi->At(2), 0.0);
	//	if (jetn==3) use3jet = true;
      }

    }

    if (jetn < 2) return -9999; //low number of jets

    SimpleParticleCollection_t daughter;
    SimpleParticleCollection_t associated;
    SimpleParticleCollection_t mother;

    daughter.push_back(SimpleParticle_t(25, Higgs)); //If studing productions it's only necessary the Higgs 4-vector

    //daughter.push_back(SimpleParticle_t(13, L1));
    //daughter.push_back(SimpleParticle_t(11, L2));
    //daughter.emplace_back(12, NuNu);
    associated.push_back(SimpleParticle_t(0,J1));
    associated.push_back(SimpleParticle_t(0,J2));

    //if (use3jet){
      //associated.push_back(SimpleParticle_t(0,J3));
      //}

    if (Higgs.Pt() == 0 || Higgs.M()==0){
      return -9999;
    }

    //MELA MATRIX ELEMENTS CALCULATION (STEP-2)
    //mela->setCandidateDecayMode(TVar::CandidateDecay_Stable);
    mela->setCandidateDecayMode(TVar::CandidateDecay_WW); //Decay to WW
    mela->setInputEvent(&daughter, &associated, 0, false);
    //mela->setInputEvent(&daughter_coll, &associated_coll, 0, 0);
    mela->setCurrentCandidateFromIndex(0);

    //->VBF Processes
    float RecoLevel_me_VBF_hsm = 0.;
    float RecoLevel_me_VBF_hm = 0.;
    float RecoLevel_me_VBF_hp = 0.;
    float RecoLevel_me_VBF_hl = 0.;
	
    //VBF Angles
    float Q2V1 = 0.;
    float Q2V2 = 0.;

    float costheta1 = 0.;
    float costheta2 = 0.;
    float costhetastar = 0.;
    
    float phi = 0.;
    float phi1 = 0.;

    //Compute and save VBF angles
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
    mela->computeVBFAngles(Q2V1, Q2V2, costheta1, costheta2, phi, costhetastar, phi1);
    mela->computeProdP(RecoLevel_me_VBF_hsm, true);	  
    MatrixElementsMap.insert({"Q2V1", Q2V1});
    MatrixElementsMap.insert({"Q2V2", Q2V2});
    MatrixElementsMap.insert({"costheta1", costheta1});
    MatrixElementsMap.insert({"costheta2", costheta2});
    MatrixElementsMap.insert({"costhetastar", costhetastar});
    MatrixElementsMap.insert({"phi", phi});
    MatrixElementsMap.insert({"phi1", phi1});
    MatrixElementsMap.insert({"RecoLevel_me_VBF_hsm", RecoLevel_me_VBF_hsm});

    mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(RecoLevel_me_VBF_hm, true);
    MatrixElementsMap.insert({"RecoLevel_me_VBF_hm", RecoLevel_me_VBF_hm});

    mela->setProcess(TVar::H0hplus, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(RecoLevel_me_VBF_hp, true);
    MatrixElementsMap.insert({"RecoLevel_me_VBF_hp", RecoLevel_me_VBF_hp});

    mela->setProcess(TVar::H0_g1prime2, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(RecoLevel_me_VBF_hl, true);
    MatrixElementsMap.insert({"RecoLevel_me_VBF_hl", RecoLevel_me_VBF_hl});

    //->QCD Processes                                                                                                                                                                                      
    float RecoLevel_me_QCD_hsm = 0.;
    float RecoLevel_me_QCD_hm = 0.;
    float RecoLevel_me_QCD_hp = 0.;
    float RecoLevel_me_QCD_hl = 0.;

    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
    mela->computeProdP(RecoLevel_me_QCD_hsm, true);
    MatrixElementsMap.insert({"RecoLevel_me_QCD_hsm", RecoLevel_me_QCD_hsm});

    mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::JJQCD);
    mela->computeProdP(RecoLevel_me_QCD_hm, true);
    MatrixElementsMap.insert({"RecoLevel_me_QCD_hm", RecoLevel_me_QCD_hm});

    mela->setProcess(TVar::H0hplus, TVar::JHUGen, TVar::JJQCD);
    mela->computeProdP(RecoLevel_me_QCD_hp, true);
    MatrixElementsMap.insert({"RecoLevel_me_QCD_hp", RecoLevel_me_QCD_hp});

    mela->setProcess(TVar::H0_g1prime2, TVar::JHUGen, TVar::JJQCD);
    mela->computeProdP(RecoLevel_me_QCD_hl, true);
    MatrixElementsMap.insert({"RecoLevel_me_QCD_hl", RecoLevel_me_QCD_hl});
	  
    //->WW/ZZ/Zgamma non-Higgs Backgrounds                                                                                                                                                                                      
    float RecoLevel_me_WW_bkg = 0.;
    float RecoLevel_me_ZZ_bkg = 0.;
    float RecoLevel_me_Zgamma_bkg = 0.;
    float RecoLevel_me_WWZZ_bkg = 0.;
    float RecoLevel_me_ZJets_bkg = 0.;

    mela->setProcess(TVar::bkgWW, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(RecoLevel_me_WW_bkg, true);
    MatrixElementsMap.insert({"RecoLevel_me_WW_bkg", RecoLevel_me_WW_bkg});
	  
    mela->setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(RecoLevel_me_ZZ_bkg, true);
    MatrixElementsMap.insert({"RecoLevel_me_ZZ_bkg", RecoLevel_me_ZZ_bkg});

    mela->setProcess(TVar::bkgZGamma, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(RecoLevel_me_Zgamma_bkg, true);
    MatrixElementsMap.insert({"RecoLevel_me_Zgamma_bkg", RecoLevel_me_Zgamma_bkg});
	  
    mela->setProcess(TVar::bkgWWZZ, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(RecoLevel_me_WWZZ_bkg, true);
    MatrixElementsMap.insert({"RecoLevel_me_WWZZ_bkg", RecoLevel_me_WWZZ_bkg});
	  
    mela->setProcess(TVar::bkgZJets, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(RecoLevel_me_ZJets_bkg, true);
    MatrixElementsMap.insert({"RecoLevel_me_ZJets_bkg", RecoLevel_me_ZJets_bkg});
	  

    //Reset Event and return results
    mela->resetInputEvent(); 
    
    float required_matrixelement = MatrixElementsMap.find(name_)->second;

    return (double)required_matrixelement;
    
  }
  //End if(nCleanJet>=2 && nLepton>1)
  else return -9999; 
}
void
RecoLevelME::bindTree_(multidraw::FunctionLibrary& _library)
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
  //Subjets
  _library.bindBranch(nSubJet, "nSubJet");
  _library.bindBranch(SubJet_pt, "SubJet_pt");
  _library.bindBranch(SubJet_eta, "SubJet_eta");
  _library.bindBranch(SubJet_phi, "SubJet_phi");
  _library.bindBranch(SubJet_mass, "SubJet_mass");
}
