#include "LatinoAnalysis/MultiDraw/interface/TTreeFunction.h"
#include "LatinoAnalysis/MultiDraw/interface/FunctionLibrary.h"

#include <vector>
#include "TLorentzVector.h"

#include "TVector2.h"
#include "Math/Vector4Dfwd.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

#include <iostream>

class mjjj : public multidraw::TTreeFunction {
public:
  mjjj();

  char const* getName() const override { return "mjjj"; }
  TTreeFunction* clone() const override { return new mjjj(); }

  unsigned getNdata() override { return 1; }
  double evaluate(unsigned) override;

protected:
  void bindTree_(multidraw::FunctionLibrary&) override;
  
  
  FloatArrayReader* CleanJet_pt{};
  FloatArrayReader* CleanJet_phi{};
  FloatArrayReader* CleanJet_eta{};
  //FloatArrayReader* CleanJet_mass{};
  UIntValueReader*  nCleanJet{};
};

mjjj::mjjj() :
  TTreeFunction()
{
}

// The actual function begins here, above is C++ bloat <3

double mjjj::evaluate(unsigned)
{
  unsigned nJ{*nCleanJet->Get()};




  if(nJ < 3){ return -9999;}

  TLorentzVector j1{};
  TLorentzVector j2{};
  TLorentzVector j3{};

  // Loop opver prompt gen jets and check overlap with prompt gen leptons
  // If there is no overlap, keep the jet (up to two jets kept)

  unsigned n{0};
  for (unsigned iJ{0}; iJ != nJ; ++iJ) {
    if (CleanJet_pt->At(iJ) < 30.){
      continue;
    }

    if (n == 0) { j1.SetPtEtaPhiM(CleanJet_pt->At(iJ), CleanJet_eta->At(iJ), CleanJet_phi->At(iJ), 0.0); }
    else if (n == 1){ j2.SetPtEtaPhiM(CleanJet_pt->At(iJ), CleanJet_eta->At(iJ), CleanJet_phi->At(iJ), 0.0); }
    else if (n == 2){j3.SetPtEtaPhiM(CleanJet_pt->At(iJ), CleanJet_eta->At(iJ), CleanJet_phi->At(iJ), 0.0); }
    ++n;
    
  }
  if ( (n < 3)){ 
    return -9999;
  }
 
  double mass = (j1 + j2 + j3).M();
  return mass;
}

void mjjj::bindTree_(multidraw::FunctionLibrary& _library)
{
  _library.bindBranch(nCleanJet, "nCleanJet");
  _library.bindBranch(CleanJet_pt, "CleanJet_pt");
  _library.bindBranch(CleanJet_eta, "CleanJet_eta");
  _library.bindBranch(CleanJet_phi, "CleanJet_phi");
  //_library.bindBranch(CleanJet_mass, "CleanJet_mass");
}
