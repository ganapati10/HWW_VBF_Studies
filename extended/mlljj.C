#include <iostream>
#include <TLorentzVector.h>
#include "LatinoAnalysis/MultiDraw/interface/TTreeFunction.h"
#include "LatinoAnalysis/MultiDraw/interface/FunctionLibrary.h"
#include <vector>
#include "TVector2.h"
#include "Math/Vector4Dfwd.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include <TMath.h>
#include <math.h>

using namespace std;


float mlljj(float lep1pt,float lep1eta, float lep1phi, float lep2pt,float lep2eta, float lep2phi, float jet1pt,float jet1eta, float jet1phi, float jet2pt,float jet2eta, float jet2phi){

  TLorentzVector l1;
  TLorentzVector l2;
  TLorentzVector j1;
  TLorentzVector j2;

  l1.SetPtEtaPhiM(lep1pt, lep1eta, lep1phi, 0.);
  l2.SetPtEtaPhiM(lep2pt, lep2eta, lep2phi, 0.);
  j1.SetPtEtaPhiM(jet1pt, jet1eta, jet1phi, 0.);
  j2.SetPtEtaPhiM(jet2pt, jet2eta, jet2phi, 0.);

  return (l1+l2+j1+j2).M();

}


###############################################################################################

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

  UIntValueReader* nLeptonGen{};
  BoolArrayReader* LeptonGen_isPrompt{};
  // IntArrayReader* DressedLepton_pdgId{};
  // FloatArrayReader* DressedLepton_pt{};
  // FloatArrayReader* DressedLepton_eta{};
  // FloatArrayReader* DressedLepton_phi{};
  // FloatArrayReader* DressedLepton_mass{};
  IntArrayReader* LeptonGen_pdgId{};
  FloatArrayReader* LeptonGen_pt{};
  FloatArrayReader* LeptonGen_eta{};
  FloatArrayReader* LeptonGen_phi{};
  FloatArrayReader* LeptonGen_mass{};
  UIntValueReader* nPhotonGen{};
  FloatArrayReader* PhotonGen_pt{};
  FloatArrayReader* PhotonGen_eta{};
  FloatArrayReader* PhotonGen_phi{};
  FloatArrayReader* PhotonGen_mass{};
  UIntValueReader* nGenJet{};
  FloatArrayReader* GenJet_pt{};
  FloatArrayReader* GenJet_eta{};
  FloatArrayReader* GenJet_phi{};
  FloatArrayReader* GenJet_mass{};
  FloatValueReader* GenMET_pt{};
  FloatValueReader* GenMET_phi{};
};

mlljj::mlljj() :
  TTreeFunction()
{
}

// The actual function begins here, above is C++ bloat <3

double mlljj::evaluate(unsigned)
{
  unsigned nJ{*nGenJet->Get()};
  
  unsigned nL{*nLeptonGen->Get()};

  // Create vector of prompt gen leptons

  std::vector<unsigned> iPromptL{};
  iPromptL.reserve(nL);

  for (unsigned iL{0}; iL != nL; ++iL) {
    if (!LeptonGen_isPrompt->At(iL))
      continue;

    unsigned absId{static_cast<unsigned>(std::abs(LeptonGen_pdgId->At(iL)))};
    if (absId != 11 && absId != 13)
      continue;

    iPromptL.push_back(iL);
  }

  // If there are no prompt leptons don't bother checking and just count the jets w/ pt > 30
  // Why not returning default value?

  if (iPromptL.size() == 0) {
    unsigned n{0};
    for (unsigned iJ{0}; iJ != nJ; ++iJ) {
      if (GenJet_pt->At(iJ) > 30.)
        ++n;
    }
    //return n;
    return -9999;
  }

  // Create prompt gen leptons 4-vectors

  std::vector<ROOT::Math::PtEtaPhiMVector> dressedLeptons{};
  for (unsigned iL : iPromptL) {
    dressedLeptons.emplace_back(
      LeptonGen_pt->At(iL),
      LeptonGen_eta->At(iL),
      LeptonGen_phi->At(iL),
      LeptonGen_mass->At(iL));
  }

  // Add to prompt gen leptons the closest gen photon 4-vector 
  // if close enough (dR < 0.09) - FSR

  unsigned nP{*nPhotonGen->Get()};

  for (unsigned iP{0}; iP != nP; ++iP) {
    double minDR2{1000.};
    int iDMin{-1};
    for (unsigned iD{0}; iD != iPromptL.size(); ++iD) {
      unsigned iL{iPromptL[iD]};
      double dEta{LeptonGen_eta->At(iL) - PhotonGen_eta->At(iP)};
      double dPhi{TVector2::Phi_mpi_pi(LeptonGen_phi->At(iL) - PhotonGen_phi->At(iP))};
      double dR2{dEta * dEta + dPhi * dPhi};
      if (dR2 < minDR2) {
        minDR2 = dR2;
        iDMin = iD;
      }
    }

    if (minDR2 < 0.09)
      dressedLeptons[iDMin] += ROOT::Math::PtEtaPhiMVector(
        PhotonGen_pt->At(iP),
        PhotonGen_eta->At(iP),
        PhotonGen_phi->At(iP),
        PhotonGen_mass->At(iP));
  }

  // Do the actual cleaning, hacked to stop at 2

  // If there are less than 2 jets, return underflow value

  if(nJ < 2) return -9999;

  TLorentzVector j1{};
  TLorentzVector j2{};

  // Loop opver prompt gen jets and check overlap with prompt gen leptons
  // If there is no overlap, keep the jet (up to two jets kept)

  unsigned n{0};
  for (unsigned iJ{0}; iJ != nJ; ++iJ) {
    if (GenJet_pt->At(iJ) < 30.)
      continue;

    bool overlap{false};
    for (auto& p4 : dressedLeptons) {
      if (p4.pt() < 10.)
        continue;

      double dEta{p4.eta() - GenJet_eta->At(iJ)};
      double dPhi{TVector2::Phi_mpi_pi(p4.phi() - GenJet_phi->At(iJ))};
      if (dEta * dEta + dPhi * dPhi < 0.016) {
        overlap = true;
        break;
      }
    }
    if (!overlap) {
      if (n == 0) { j1.SetPtEtaPhiM(GenJet_pt->At(iJ), GenJet_eta->At(iJ), GenJet_phi->At(iJ), GenJet_mass->At(iJ)); }
      else if (n == 1){ j2.SetPtEtaPhiM(GenJet_pt->At(iJ), GenJet_eta->At(iJ), GenJet_phi->At(iJ), GenJet_mass->At(iJ)); }
      ++n;
    }
  }
  if ( (n < 2) || (*nLeptonGen->Get() < 2) ) return -9999;
  else if ( LeptonGen_pdgId->At(0) * LeptonGen_pdgId->At(1) != -11*13 ) return -9999;
  
  double pt3{-1};
  if (dressedLeptons.size() > 2) pt3 = dressedLeptons[2].pt();
  

  // Define additional useful variables (for fiducial region definition)

  TLorentzVector MET{};
  MET.SetPtEtaPhiM(*GenMET_pt->Get(), 0., *GenMET_phi->Get(), 0.);
  TLorentzVector l1{};
  l1.SetPtEtaPhiM(dressedLeptons[0].pt(), dressedLeptons[0].eta(), dressedLeptons[0].phi(), dressedLeptons[0].M());
  TLorentzVector l2{};
  l2.SetPtEtaPhiM(dressedLeptons[1].pt(), dressedLeptons[1].eta(), dressedLeptons[1].phi(), dressedLeptons[1].M());
  
  double mass = (l1 + l2 + j1 + j2).M();
  return mass;
}

void mlljj::bindTree_(multidraw::FunctionLibrary& _library)
{
  _library.bindBranch(nLeptonGen, "nLeptonGen");
  _library.bindBranch(LeptonGen_isPrompt, "LeptonGen_isPrompt");
  // _library.bindBranch(DressedLepton_pdgId, "DressedLepton_pdgId");
  // _library.bindBranch(DressedLepton_pt, "DressedLepton_pt");
  // _library.bindBranch(DressedLepton_eta, "DressedLepton_eta");
  // _library.bindBranch(DressedLepton_phi, "DressedLepton_phi");
  // _library.bindBranch(DressedLepton_mass, "DressedLepton_mass");
  _library.bindBranch(LeptonGen_pdgId, "LeptonGen_pdgId");
  _library.bindBranch(LeptonGen_pt, "LeptonGen_pt");
  _library.bindBranch(LeptonGen_eta, "LeptonGen_eta");
  _library.bindBranch(LeptonGen_phi, "LeptonGen_phi");
  _library.bindBranch(LeptonGen_mass, "LeptonGen_mass");
  _library.bindBranch(nPhotonGen, "nPhotonGen");
  _library.bindBranch(PhotonGen_pt, "PhotonGen_pt");
  _library.bindBranch(PhotonGen_eta, "PhotonGen_eta");
  _library.bindBranch(PhotonGen_phi, "PhotonGen_phi");
  _library.bindBranch(PhotonGen_mass, "PhotonGen_mass");
  _library.bindBranch(nGenJet, "nGenJet");
  _library.bindBranch(GenJet_pt, "GenJet_pt");
  _library.bindBranch(GenJet_eta, "GenJet_eta");
  _library.bindBranch(GenJet_phi, "GenJet_phi");
  _library.bindBranch(GenJet_mass, "GenJet_mass");
  _library.bindBranch(GenMET_pt, "GenMET_pt");
  _library.bindBranch(GenMET_phi, "GenMET_phi");
}











