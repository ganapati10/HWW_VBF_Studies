#include "LatinoAnalysis/MultiDraw/interface/TTreeFunction.h"
#include "LatinoAnalysis/MultiDraw/interface/FunctionLibrary.h"

#include <vector>
#include "TLorentzVector.h"

#include "TVector2.h"
#include "Math/Vector4Dfwd.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

#include <iostream>

class GetJetDis : public multidraw::TTreeFunction {
public:
  GetJetDis();

  char const* getName() const override { return "GetJetDis"; }
  TTreeFunction* clone() const override { return new GetJetDis(); }

  unsigned getNdata() override { return 1; }
  double evaluate(unsigned) override;

protected:
  void bindTree_(multidraw::FunctionLibrary&) override;
  
  
  FloatArrayReader* CleanJet_pt{};
  //FloatArrayReader* CleanJet_phi{};
  FloatArrayReader* CleanJet_eta{};
  //FloatArrayReader* CleanJet_mass{};
  UIntValueReader*  nCleanJet{};
};

GetJetDis::GetJetDis() :
  TTreeFunction()
{
}

// The actual function begins here, above is C++ bloat <3

double GetJetDis::evaluate(unsigned)
{
  unsigned nJ{*nCleanJet->Get()};

  if(nJ < 3){ return 9999;}

  double eta1 = -9999;
  double eta2 = -9999;
  double eta3 = -9999;

  // Loop opver prompt gen jets and check overlap with prompt gen leptons
  // If there is no overlap, keep the jet (up to two jets kept)

  unsigned n{0};
  for (unsigned iJ{0}; iJ != nJ; ++iJ) {
    if (CleanJet_pt->At(iJ) < 30.){
      continue;
    }

    if (n == 0) {  eta1 = CleanJet_eta->At(iJ); }
    else if (n == 1){  eta2 = CleanJet_eta->At(iJ); }
    else if (n == 2){ eta3 = CleanJet_eta->At(iJ); }
    ++n;
    
  }
  if ( (n < 2) || eta1 == -9999 || eta2 == -9999 || eta3 == -9999){ 
    return 9999;
  }
 
  double output = eta3 - (eta1 + eta2)/2;
  return output;
}

void GetJetDis::bindTree_(multidraw::FunctionLibrary& _library)
{
  _library.bindBranch(nCleanJet, "nCleanJet");
  _library.bindBranch(CleanJet_pt, "CleanJet_pt");
  _library.bindBranch(CleanJet_eta, "CleanJet_eta");
  //_library.bindBranch(CleanJet_phi, "CleanJet_phi");
  //_library.bindBranch(CleanJet_mass, "CleanJet_mass");
}
