#!/usr/bin/env python

import tensorflow.keras
from keras.utils import np_utils
import tensorflow.keras.callbacks as cb
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation
from tensorflow.keras import regularizers
from tensorflow.keras import backend as K
from tensorflow.keras import optimizers
from tensorflow.keras import callbacks

from tensorflow.keras.models import load_model

import numpy as np
import sys

argv = sys.argv
sys.argv = argv[:1]

def compute_prediction(mjj,
    Ctot,
    Jet_qgl0,
    Jet_qgl1,
    detajj,
    detall,
    drjj,
    dphill,
    dphijjmet,
    dphilljetjet,
    drll,
    Lepton_eta0,
    Lepton_eta1,
    Lepton_pt0,
    Lepton_pt1,
    Lepton_phi0,
    Lepton_phi1,  
    CleanJet_eta0,
    CleanJet_eta1,
    CleanJet_phi0,
    CleanJet_phi1,
    CleanJet_pt0,
    CleanJet_pt1,
    PuppiMET_pt,
    PuppiMET_phi,  
    mth,
    mTi,
    mtw2,
    detal1j1,
    detal1j2,
    detal2j1,
    detal2j2,
    ptll,
    mlljj,
    PTotal, 
    mll,
    RecoMELA_VBF,
    RecoMELA_Phi,
    RecoMELA_CT
    ):
  
  print("Loading DNN model from .h5 file")
  
  model = load_model('Models/model_Higgs_Complex.h5')
  
  print("Model already load and compile")
  
  input = np.array([mjj,
    Ctot,
    Jet_qgl0,
    Jet_qgl1,
    detajj,
    detall,
    drjj,
    dphill,
    dphijjmet,
    dphilljetjet,
    drll,
    Lepton_eta0,
    Lepton_eta1,
    Lepton_pt0,
    Lepton_pt1,
    Lepton_phi0,
    Lepton_phi1,  
    CleanJet_eta0,
    CleanJet_eta1,
    CleanJet_phi0,
    CleanJet_phi1,
    CleanJet_pt0,
    CleanJet_pt1,
    PuppiMET_pt,
    PuppiMET_phi,  
    mth,
    mTi,
    mtw2,
    detal1j1,
    detal1j2,
    detal2j1,
    detal2j2,
    ptll,
    mlljj,
    PTotal, 
    mll,
    RecoMELA_VBF,
    RecoMELA_Phi,
    RecoMELA_CT])
  
  result = model.predict(input)
  
  return result[0]

if __name__=='__main__':
    
    input_var = sys.argv
    
    input = 0

    if(input_var[0]=='Load_DNN.py'):
        input = input_var[1:]
    elif(len(input_var==39)):
        input = input_var
      
    if (input==0):
        print("Error loading input variables")
        return -999
    
    return compute_prediction(input)
    
    
    
