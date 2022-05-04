#include "TLorentzVector.h"

float mlj(float CleanJet_pt, float CleanJet_phi, float CleanJet_eta, float Lepton_pt, float Lepton_phi, float Lepton_eta){

  TLorentzVector j;
  TLorentzVector l;

  j.SetPtEtaPhiM(CleanJet_pt, CleanJet_eta, CleanJet_phi, 0.0);
  l.SetPtEtaPhiM(Lepton_pt, Lepton_eta, Lepton_phi, 0.0);

  double mass = (j + l).M();
  return mass;
}
