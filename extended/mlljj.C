#include <iostream>
#include <TMath.h>
#include <TLorentzVector.h>

using namespace std;


float mlljj(float lep1pt,float lep1eta, float lep1phi, float lep2pt,float lep2eta, float lep2phi, float jet1pt,float jet1eta, float jet1phi, float jet1mass, float jet2pt,float jet2eta, float jet2phi, float jet2mass){

  TLorentzVector L1,L2;
  TLorentzVector J1,J2;
  
  L1.SetPtEtaPhiM(lep1pt, lep1eta, lep1phi, 0.);
  L2.SetPtEtaPhiM(lep2pt, lep2eta, lep2phi, 0.);
  J1.SetPtEtaPhiM(jet1pt, jet1eta, jet1phi, jet1mass);
  J2.SetPtEtaPhiM(jet2pt, jet2eta, jet2phi, jet2mass);
 
  return (L1 + L2 + J1 + J2).M();
}
