#include "TLorentzVector.h"

float mlljj(float nCleanJet, float nLepton, float CleanJet_pt0, float CleanJet_pt1, float CleanJet_phi0, float CleanJet_phi1, float CleanJet_eta0, float CleanJet_eta1, float Lepton_pt0, float Lepton_pt1,float Lepton_phi0, float Lepton_phi1, float Lepton_eta0, float Lepton_eta1){

  if (nCleanJet < 2 || nLepton < 2) return -9999.;


  TLorentzVector j1;
  TLorentzVector j2;
  TLorentzVector l1;
  TLorentzVector l2;

  j1.SetPtEtaPhiM(CleanJet_pt0, CleanJet_eta0, CleanJet_phi0, 0.0);
  j2.SetPtEtaPhiM(CleanJet_pt1, CleanJet_eta1, CleanJet_phi1, 0.0);
  l1.SetPtEtaPhiM(Lepton_pt0, Lepton_eta0, Lepton_phi0, 0.0);
  l2.SetPtEtaPhiM(Lepton_pt1, Lepton_eta1, Lepton_phi1, 0.0);

  if (CleanJet_pt0 < 30 || CleanJet_pt1 < 30 || Lepton_pt0 < 10 || Lepton_pt1 < 10) return -9999.;

  double mass = (j1 + j2 + l1 + l2).M();
  return mass;
}

