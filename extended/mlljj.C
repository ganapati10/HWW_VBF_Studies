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
