#include <iostream>


float drjj(float jet1_eta,
		float jet2_eta,
		float jet1_phi,
		float jet2_hi){
    
  return TMath::sqrt(jet1_eta * jet2_eta + jet1_phi * jet2_phi)   
    
}  
