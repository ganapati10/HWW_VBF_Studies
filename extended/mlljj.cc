#include "LatinoAnalysis/MultiDraw/interface/TTreeFunction.h"
#include "LatinoAnalysis/MultiDraw/interface/FunctionLibrary.h"

#include <vector>
#include "TLorentzVector.h"

#include "TVector2.h"
#include "Math/Vector4Dfwd.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

#include <iostream>

class mlljj : public multidraw::TTreeFunction {
public:
  mlljj();

  char const* getName() const override { return "mlljj"; }
  TTreeFunction* clone() const override { return new mlljj(); }

  unsigned getNdata() override { return 1; }
  double evaluate(unsigned) override;

protected:
  void bindTree_(multidraw::FunctionLibrary&) override;
  
  FloatValueReader* MET_pt{};
  FloatValueReader* MET_phi{};
  
  //IntArrayReader* Lepton_pdgId{};
  //BoolArrayReader* Lepton_isPrompt{};
  FloatArrayReader* Lepton_pt{};
  FloatArrayReader* Lepton_phi{};
  FloatArrayReader* Lepton_eta{};
  //FloatArrayReader* Lepton_mass{};
  UIntValueReader*  nLepton{};
  
  FloatArrayReader* CleanJet_pt{};
  FloatArrayReader* CleanJet_phi{};
  FloatArrayReader* CleanJet_eta{};
  //FloatArrayReader* CleanJet_mass{};
  UIntValueReader*  nCleanJet{};
};

mlljj::mlljj() :
  TTreeFunction()
{
}

// The actual function begins here, above is C++ bloat <3

double mlljj::evaluate(unsigned)
{
  unsigned nJ{*nCleanJet->Get()};
  
  unsigned nL{*nLepton->Get()};

  // Create vector of prompt gen leptons

  std::vector<unsigned> iPromptL{};
  iPromptL.reserve(nL);

  for (unsigned iL{0}; iL != nL; ++iL) {
    iPromptL.push_back(iL);
  }

  // If there are no prompt leptons don't bother checking and just count the jets w/ pt > 30
  // Why not returning default value?

  if (iPromptL.size() == 0) {
    unsigned n{0};
    for (unsigned iJ{0}; iJ != nJ; ++iJ) {
      if (CleanJet_pt->At(iJ) > 30.){
        ++n;
      }
    }
    //return n;
    return -9999;
  }

  // Create prompt gen leptons 4-vectors

  std::vector<ROOT::Math::PtEtaPhiMVector> dressedLeptons{};
  for (unsigned iL : iPromptL) {
    dressedLeptons.emplace_back(
				Lepton_pt->At(iL),
				Lepton_eta->At(iL),
				Lepton_phi->At(iL),
				0.0);
				//Lepton_mass->At(iL));
  }


  // Do the actual cleaning, hacked to stop at 2

  // If there are less than 2 jets, return underflow value

  if(nJ < 2){ return -9999;}

  TLorentzVector j1{};
  TLorentzVector j2{};

  // Loop opver prompt gen jets and check overlap with prompt gen leptons
  // If there is no overlap, keep the jet (up to two jets kept)

  unsigned n{0};
  for (unsigned iJ{0}; iJ != nJ; ++iJ) {
    if (CleanJet_pt->At(iJ) < 30.){
      continue;
    }

    bool overlap{false};
    for (auto& p4 : dressedLeptons) {
      if (p4.pt() < 10.){
        continue;
      }

      double dEta{p4.eta() - CleanJet_eta->At(iJ)};
      double dPhi{TVector2::Phi_mpi_pi(p4.phi() - CleanJet_phi->At(iJ))};
      if (dEta * dEta + dPhi * dPhi < 0.016) {
        overlap = true;
        break;
      }
    }
    if (!overlap) {
      if (n == 0) { j1.SetPtEtaPhiM(CleanJet_pt->At(iJ), CleanJet_eta->At(iJ), CleanJet_phi->At(iJ), 0.0); }
      else if (n == 1){ j2.SetPtEtaPhiM(CleanJet_pt->At(iJ), CleanJet_eta->At(iJ), CleanJet_phi->At(iJ), 0.0); }
      ++n;
    }
  }
  if ( (n < 2) || (*nLepton->Get() < 2) ){ 
    return -9999;
  }
  
  double pt3{-1};
  if (dressedLeptons.size() > 2){ 
    pt3 = dressedLeptons[2].pt();
  }
  

  // Define additional useful variables (for fiducial region definition)

  TLorentzVector MET{};
  MET.SetPtEtaPhiM(*MET_pt->Get(), 0., *MET_phi->Get(), 0.);
  TLorentzVector l1{};
  l1.SetPtEtaPhiM(dressedLeptons[0].pt(), dressedLeptons[0].eta(), dressedLeptons[0].phi(), dressedLeptons[0].M());
  TLorentzVector l2{};
  l2.SetPtEtaPhiM(dressedLeptons[1].pt(), dressedLeptons[1].eta(), dressedLeptons[1].phi(), dressedLeptons[1].M());
  
  double mass = (l1 + l2 + j1 + j2).M();
  return mass;
}

void mlljj::bindTree_(multidraw::FunctionLibrary& _library)
{
  _library.bindBranch(nLepton, "nLepton");
  //_library.bindBranch(Lepton_isPrompt, "Lepton_isPrompt");
  // _library.bindBranch(DressedLepton_pdgId, "DressedLepton_pdgId");
  // _library.bindBranch(DressedLepton_pt, "DressedLepton_pt");
  // _library.bindBranch(DressedLepton_eta, "DressedLepton_eta");
  // _library.bindBranch(DressedLepton_phi, "DressedLepton_phi");
  // _library.bindBranch(DressedLepton_mass, "DressedLepton_mass");
  //_library.bindBranch(Lepton_pdgId, "Lepton_pdgId");
  _library.bindBranch(Lepton_pt, "Lepton_pt");
  _library.bindBranch(Lepton_eta, "Lepton_eta");
  _library.bindBranch(Lepton_phi, "Lepton_phi");
  //_library.bindBranch(Lepton_mass, "Lepton_mass");
  _library.bindBranch(nCleanJet, "nCleanJet");
  _library.bindBranch(CleanJet_pt, "CleanJet_pt");
  _library.bindBranch(CleanJet_eta, "CleanJet_eta");
  _library.bindBranch(CleanJet_phi, "CleanJet_phi");
  //_library.bindBranch(CleanJet_mass, "CleanJet_mass");
  _library.bindBranch(MET_pt, "MET_pt");
  _library.bindBranch(MET_phi, "MET_phi");
}
