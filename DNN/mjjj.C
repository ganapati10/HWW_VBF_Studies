#include "TLorentzVector.h"

float mjjj(float nCleanJet, float CleanJet_pt0, float CleanJet_pt1, float CleanJet_pt2, float CleanJet_phi0, float CleanJet_phi1, float CleanJet_phi2, float CleanJet_eta0, float CleanJet_eta1, float CleanJet_eta2){

  
  if (nCleanJet < 3) return -9999.;

  TLorentzVector j1;
  TLorentzVector j2;
  TLorentzVector j3;

  j1.SetPtEtaPhiM(CleanJet_pt0, CleanJet_eta0, CleanJet_phi0, 0.0); 
  j2.SetPtEtaPhiM(CleanJet_pt1, CleanJet_eta1, CleanJet_phi1, 0.0);
  j3.SetPtEtaPhiM(CleanJet_pt2, CleanJet_eta2, CleanJet_phi2, 0.0);

  if (CleanJet_pt0 < 30 || CleanJet_pt1 < 30 || CleanJet_pt2 < 30) return -9999.;

  double mass = (j1 + j2 + j3).M();
  return mass;

}

