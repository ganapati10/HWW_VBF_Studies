#
# Prepare samples
# 
#

from __future__ import print_function
import pickle
import ROOT 
import numpy as np  
import pandas as pd
import sys
import os
import root_numpy



dir16 = '/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Summer16_102X_nAODv7_Full2016v7/MCl1loose2016v7__MCCorr2016v7__l2loose__l2tightOR2016v7/'
dir17 = '/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Fall2017_102X_nAODv7_Full2017v7/MCl1loose2017v7__MCCorr2017v7__l2loose__l2tightOR2017v7/'
dir18 = '/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Autumn18_102X_nAODv7_Full2018v7/MCl1loose2018v7__MCCorr2018v7__l2loose__l2tightOR2018v7/'


#
#
# VBF Monte Carlo
#
#


def load_dataset_vbf ( max_entries = -1 ):
  _branches = [
    "mjj",
    "log((abs(2*Lepton_eta[0]-CleanJet_eta[0]-CleanJet_eta[1])+abs(2*Lepton_eta[1]-CleanJet_eta[0]-CleanJet_eta[1]))/detajj)",
    "Jet_qgl[0]",
    "Jet_qgl[1]",
    "detajj",
    "Lepton_eta[0]-Lepton_eta[1]",
    "sqrt(CleanJet_eta[0]*CleanJet_eta[1] + CleanJet_phi[0]*CleanJet_phi[1])",
    "dphill",
    "dphijjmet",
    "dphilljetjet",
    "drll",
    "Lepton_eta[0]",
    "Lepton_eta[1]",
    "Lepton_pt[0]",
    "Lepton_pt[1]",
    "Lepton_phi[0]",
    "Lepton_phi[1]",  
    "CleanJet_eta[0]",
    "CleanJet_eta[1]",
    "CleanJet_eta[2]",
    "CleanJet_phi[0]",
    "CleanJet_phi[1]",
    "CleanJet_phi[2]",
    "abs(CleanJet_eta[2]-(CleanJet_eta[0]+CleanJet_eta[1])/2)*(CleanJet_pt[2]>30)",
    "CleanJet_pt[0]",
    "CleanJet_pt[1]",
    "MET_pt",
    "mth",
    "ptll",
    "mlljj(Sum$(CleanJet_pt>30), nLepton, CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1], Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1])",
    "mjjj(Sum$(CleanJet_pt>30), CleanJet_pt[0], CleanJet_pt[1], CleanJet_pt[2], CleanJet_phi[0], CleanJet_phi[1], CleanJet_phi[2], CleanJet_eta[0], CleanJet_eta[1], CleanJet_eta[2])",
    "Lepton_pt[0] + Lepton_pt[1] + CleanJet_pt[0] + CleanJet_pt[1] + MET_pt", 
    "mll",
    "RecoMELA_VBF(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1], Lepton_pdgId[0], Lepton_pdgId[1])",
    "RecoMELA_Q2V1(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",
    #"RecoMELA_Q2V2(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",
    "RecoMELA_Phi(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])"
    #"RecoMELA_CT1(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",
    #"RecoMELA_CT2(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",  
    ]
 
  ROOT.gROOT.ProcessLineSync(".L /eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/mlljj.C+")
  ROOT.gROOT.ProcessLineSync(".L /eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/mjjj.C+")
  
  
  ROOT.gROOT.ProcessLineSync('gSystem->Load("/eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/JHUGenMELA/MELA/data/slc7_amd64_gcc700/libmcfm_707.so","", kTRUE);')
  ROOT.gROOT.ProcessLineSync('gSystem->Load("/eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/lib/slc7_amd64_gcc700/libJHUGenMELAMELA.so","", kTRUE);')  

  ROOT.gROOT.ProcessLineSync(".L /eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/RecoMELA_VBF.C+")
  ROOT.gROOT.ProcessLineSync(".L /eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/RecoMELA_Q2V1.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_Q2V2.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_CT1.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_CT2.C+")
  ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_Phi.C+")

  chain = ROOT.TChain('Events')

  ### Full 2016 ###
    
  chain.Add(dir16+'nanoLatino_VBFHToWWTo2L2Nu_M125__part0.root')
  chain.Add(dir16+'nanoLatino_VBFHToWWTo2L2Nu_M125__part1.root')

  ### Full 2017 ###

  #chain.Add(dir17+'nanoLatino_VBFHToWWTo2L2Nu_M125__part0.root')
  
  for i in range(0, 22, 1):
        chain.Add(dir17+'nanoLatino_VBFHToWWTo2L2Nu_M125__part'+str(i)+'.root')
  
  ### Full 2018 ###  

  for i in range(0, 12, 1):
        chain.Add(dir18+'nanoLatino_VBFHToWWTo2L2Nu_M125__part'+str(i)+'.root')
    
    
  print(chain.GetEntries())
  
  _dataset = root_numpy.tree2array (chain,
      branches = _branches,
      selection = 'Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13 && Lepton_pt[0] > 25. && Lepton_pt[1] > 13. && (abs(Lepton_pdgId[1]) == 13 || Lepton_pt[1] > 13.) && (nLepton >= 2 && Alt$(Lepton_pt[2], 0) < 10.) && Alt$(CleanJet_pt[1], 0) > 30. && abs(CleanJet_eta[0]) < 4.7 && abs(CleanJet_eta[1]) < 4.7 && mth > 40 && mth < 125 && Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217) == 0 && Sum$(CleanJet_pt>30) >= 2 && Sum$(CleanJet_pt>30) <= 3', 
      stop = max_entries
     )

  return { b : _dataset[b] for b in _branches }


#
#
# GGH Monte Carlo
#
#


def load_dataset_ggh ( max_entries = -1 ):
  _branches = [
    "mjj",
    "log((abs(2*Lepton_eta[0]-CleanJet_eta[0]-CleanJet_eta[1])+abs(2*Lepton_eta[1]-CleanJet_eta[0]-CleanJet_eta[1]))/detajj)",
    "Jet_qgl[0]",
    "Jet_qgl[1]",
    "detajj",
    "Lepton_eta[0]-Lepton_eta[1]",
    "sqrt(CleanJet_eta[0]*CleanJet_eta[1] + CleanJet_phi[0]*CleanJet_phi[1])",
    "dphill",
    "dphijjmet",
    "dphilljetjet",
    "drll",
    "Lepton_eta[0]",
    "Lepton_eta[1]",
    "Lepton_pt[0]",
    "Lepton_pt[1]",
    "Lepton_phi[0]",
    "Lepton_phi[1]",  
    "CleanJet_eta[0]",
    "CleanJet_eta[1]",
    "CleanJet_eta[2]",
    "CleanJet_phi[0]",
    "CleanJet_phi[1]",
    "CleanJet_phi[2]",
    "abs(CleanJet_eta[2]-(CleanJet_eta[0]+CleanJet_eta[1])/2)*(CleanJet_pt[2]>30)",
    "CleanJet_pt[0]",
    "CleanJet_pt[1]",
    "MET_pt",
    "mth",
    "ptll",
    "mlljj(Sum$(CleanJet_pt>30), nLepton, CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1], Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1])",
    "mjjj(Sum$(CleanJet_pt>30), CleanJet_pt[0], CleanJet_pt[1], CleanJet_pt[2], CleanJet_phi[0], CleanJet_phi[1], CleanJet_phi[2], CleanJet_eta[0], CleanJet_eta[1], CleanJet_eta[2])",
    "Lepton_pt[0] + Lepton_pt[1] + CleanJet_pt[0] + CleanJet_pt[1] + MET_pt", 
    "mll",
    "RecoMELA_VBF(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1], Lepton_pdgId[0], Lepton_pdgId[1])",
    "RecoMELA_Q2V1(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",
    #"RecoMELA_Q2V2(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",
    "RecoMELA_Phi(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])"
    #"RecoMELA_CT1(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",
    #"RecoMELA_CT2(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",  
    ]
 

  ROOT.gROOT.ProcessLineSync(".L ./mlljj.C+")
  ROOT.gROOT.ProcessLineSync(".L ./mjjj.C+")

   
  ROOT.gROOT.ProcessLineSync('gSystem->Load("/eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/JHUGenMELA/MELA/data/slc7_amd64_gcc700/libmcfm_707.so","", kTRUE);')
  ROOT.gROOT.ProcessLineSync('gSystem->Load("/eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/lib/slc7_amd64_gcc700/libJHUGenMELAMELA.so","", kTRUE);')  

  ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_VBF.C+")
  ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_Q2V1.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_Q2V2.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_CT1.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_CT2.C+")
  ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_Phi.C+")


  chain = ROOT.TChain('Events')  
    

  ### Full 2016 ###
    
  #chain.Add(dir16+'nanoLatino_GluGluHToWWTo2L2Nu_alternative_M125__part0.root')

  for i in range(0, 21, 1):
        chain.Add(dir16+'nanoLatino_GGHjjToWWTo2L2Nu_minloHJJ_M125__part'+str(i)+'.root')

  ### Full 2017 ###

  #chain.Add(dir17+'nanoLatino_GluGluHToWWTo2L2Nu_M125__part0.root')

  for i in range(0, 37, 1):
        chain.Add(dir17+'nanoLatino_GluGluHToWWTo2L2Nu_M125__part'+str(i)+'.root')
        
  for i in range(0, 44, 1):
        chain.Add(dir17+'nanoLatino_GGHjjToWWTo2L2Nu_minloHJJ_M125__part'+str(i)+'.root')      

  ### Full 2018 ###
    
  chain.Add(dir18+'nanoLatino_GluGluHToWWTo2L2Nu_M125__part0.root')
  chain.Add(dir18+'nanoLatino_GluGluHToWWTo2L2Nu_M125__part1.root')
    
  for i in range(0, 57, 1):
    chain.Add(dir18+'nanoLatino_GGHjjToWWTo2L2Nu_minloHJJ_M125__part'+str(i)+'.root')


  
  print(chain.GetEntries())

  _dataset = root_numpy.tree2array (chain, 
      selection = 'Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13 && Lepton_pt[0] > 25. && Lepton_pt[1] > 13. && (abs(Lepton_pdgId[1]) == 13 || Lepton_pt[1] > 13.) && (nLepton >= 2 && Alt$(Lepton_pt[2], 0) < 10.) && Alt$(CleanJet_pt[1], 0) > 30. && abs(CleanJet_eta[0]) < 4.7 && abs(CleanJet_eta[1]) < 4.7 && mth > 40 && mth < 125 && Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217) == 0 && Sum$(CleanJet_pt>30) >= 2 && Sum$(CleanJet_pt>30) <= 3',
      branches = _branches,
      stop = max_entries
     )

  return { b : _dataset[b] for b in _branches }



#
#
# TOP Monte Carlo
#
#



def load_dataset_top ( max_entries = -1 ):
  _branches = [
    "mjj",
    "log((abs(2*Lepton_eta[0]-CleanJet_eta[0]-CleanJet_eta[1])+abs(2*Lepton_eta[1]-CleanJet_eta[0]-CleanJet_eta[1]))/detajj)",
    "Jet_qgl[0]",
    "Jet_qgl[1]",
    "detajj",
    "Lepton_eta[0]-Lepton_eta[1]",
    "sqrt(CleanJet_eta[0]*CleanJet_eta[1] + CleanJet_phi[0]*CleanJet_phi[1])",
    "dphill",
    "dphijjmet",
    "dphilljetjet",
    "drll",
    "Lepton_eta[0]",
    "Lepton_eta[1]",
    "Lepton_pt[0]",
    "Lepton_pt[1]",
    "Lepton_phi[0]",
    "Lepton_phi[1]",  
    "CleanJet_eta[0]",
    "CleanJet_eta[1]",
    "CleanJet_eta[2]",
    "CleanJet_phi[0]",
    "CleanJet_phi[1]",
    "CleanJet_phi[2]",
    "abs(CleanJet_eta[2]-(CleanJet_eta[0]+CleanJet_eta[1])/2)*(CleanJet_pt[2]>30)",
    "CleanJet_pt[0]",
    "CleanJet_pt[1]",
    "MET_pt",
    "mth",
    "ptll",
    "mlljj(Sum$(CleanJet_pt>30), nLepton, CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1], Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1])",
    "mjjj(Sum$(CleanJet_pt>30), CleanJet_pt[0], CleanJet_pt[1], CleanJet_pt[2], CleanJet_phi[0], CleanJet_phi[1], CleanJet_phi[2], CleanJet_eta[0], CleanJet_eta[1], CleanJet_eta[2])",
    "Lepton_pt[0] + Lepton_pt[1] + CleanJet_pt[0] + CleanJet_pt[1] + MET_pt", 
    "mll",
    "RecoMELA_VBF(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1], Lepton_pdgId[0], Lepton_pdgId[1])",
    "RecoMELA_Q2V1(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",
    #"RecoMELA_Q2V2(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",
    "RecoMELA_Phi(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])"
    #"RecoMELA_CT1(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",
    #"RecoMELA_CT2(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",  
    ]
 

  ROOT.gROOT.ProcessLineSync(".L ./mlljj.C+")
  ROOT.gROOT.ProcessLineSync(".L ./mjjj.C+")

   
  ROOT.gROOT.ProcessLineSync('gSystem->Load("/eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/JHUGenMELA/MELA/data/slc7_amd64_gcc700/libmcfm_707.so","", kTRUE);')
  ROOT.gROOT.ProcessLineSync('gSystem->Load("/eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/lib/slc7_amd64_gcc700/libJHUGenMELAMELA.so","", kTRUE);')  

  ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_VBF.C+")
  ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_Q2V1.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_Q2V2.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_CT1.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_CT2.C+")
  ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_Phi.C+")


  chain = ROOT.TChain('Events')

  
  ### Full 2016 ###

  #chain.Add(dir16+'nanoLatino_TTTo2L2Nu__part0.root')  
    
  for i in range(0, 72, 1):
        chain.Add(dir16+'nanoLatino_TTTo2L2Nu__part'+str(i)+'.root')
        
  chain.Add(dir16+'nanoLatino_ST_s-channel__part0.root')
  chain.Add(dir16+'nanoLatino_ST_s-channel__part1.root')
  chain.Add(dir16+'nanoLatino_ST_s-channel__part2.root')
  
  for i in range(0, 35, 1):
        chain.Add(dir16+'nanoLatino_ST_t-channel_antitop__part'+str(i)+'.root')  
        
  for i in range(0, 143, 1):
        chain.Add(dir16+'nanoLatino_ST_t-channel_top__part'+str(i)+'.root')     
        
  for i in range(0, 9, 1):
        chain.Add(dir16+'nanoLatino_ST_tW_antitop__part'+str(i)+'.root') 
  
  chain.Add(dir16+'nanoLatino_ST_tW_top__part0.root')
  chain.Add(dir16+'nanoLatino_ST_tW_top__part1.root')
  chain.Add(dir16+'nanoLatino_ST_tW_top__part2.root')
  chain.Add(dir16+'nanoLatino_ST_tW_top__part3.root')
  chain.Add(dir16+'nanoLatino_ST_tW_top__part4.root')
  chain.Add(dir16+'nanoLatino_ST_tW_top__part5.root')
  chain.Add(dir16+'nanoLatino_ST_tW_top__part6.root')
        
    
  ### Full 2017 ### 
    
  #chain.Add(dir17+'nanoLatino_TTTo2L2Nu_PSWeights__part0.root')  
    
  for i in range(0, 75, 1):
        chain.Add(dir17+'nanoLatino_TTTo2L2Nu_PSWeights__part'+str(i)+'.root')  
        
  for i in range(0, 10, 1):
        chain.Add(dir17+'nanoLatino_ST_s-channel__part'+str(i)+'.root')
        
  for i in range(0, 17, 1):
        chain.Add(dir17+'nanoLatino_ST_t-channel_antitop__part'+str(i)+'.root') 
        
  for i in range(0, 7, 1):
        chain.Add(dir17+'nanoLatino_ST_t-channel_top__part'+str(i)+'.root')
        
  for i in range(0, 11, 1):
        chain.Add(dir17+'nanoLatino_ST_tW_antitop__part'+str(i)+'.root')
        
  for i in range(0, 7, 1):
        chain.Add(dir17+'nanoLatino_ST_tW_top__part'+str(i)+'.root')
        
    
  ### Full 2018 ###

  #chain.Add(dir18+'nanoLatino_TTTo2L2Nu__part0.root')

  for i in range(0, 59, 1):
        chain.Add(dir18+'nanoLatino_TTTo2L2Nu__part'+str(i)+'.root')
  
  for i in range(0, 19, 1):
        chain.Add(dir18+'nanoLatino_ST_s-channel_ext1__part'+str(i)+'.root')
        
  for i in range(0, 60, 1):
        chain.Add(dir18+'nanoLatino_ST_t-channel_antitop__part'+str(i)+'.root')
        
  for i in range(0, 102, 1):
        chain.Add(dir18+'nanoLatino_ST_t-channel_top__part'+str(i)+'.root')
  
  for i in range(0, 8, 1):
        chain.Add(dir18+'nanoLatino_ST_tW_antitop_ext1__part'+str(i)+'.root')  
        
  for i in range(0, 9, 1):
        chain.Add(dir18+'nanoLatino_ST_tW_top_ext1__part'+str(i)+'.root') 

    
  
  print(chain.GetEntries())
  
  _dataset = root_numpy.tree2array (chain, 
      selection = 'Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13 && Lepton_pt[0] > 25. && Lepton_pt[1] > 13. && (abs(Lepton_pdgId[1]) == 13 || Lepton_pt[1] > 13.) && (nLepton >= 2 && Alt$(Lepton_pt[2], 0) < 10.) && Alt$(CleanJet_pt[1], 0) > 30. && abs(CleanJet_eta[0]) < 4.7 && abs(CleanJet_eta[1]) < 4.7 && mth > 40 && mth < 125 && Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217) == 0 && Sum$(CleanJet_pt>30) >= 2 && Sum$(CleanJet_pt>30) <= 3',
      branches = _branches,
      stop = max_entries
     )

  return { b : _dataset[b] for b in _branches }


#
#
# WW Monte Carlo
#
#


def load_dataset_ww ( max_entries = -1 ):
  _branches = [
    "mjj",
    "log((abs(2*Lepton_eta[0]-CleanJet_eta[0]-CleanJet_eta[1])+abs(2*Lepton_eta[1]-CleanJet_eta[0]-CleanJet_eta[1]))/detajj)",
    "Jet_qgl[0]",
    "Jet_qgl[1]",
    "detajj",
    "Lepton_eta[0]-Lepton_eta[1]",
    "sqrt(CleanJet_eta[0]*CleanJet_eta[1] + CleanJet_phi[0]*CleanJet_phi[1])",
    "dphill",
    "dphijjmet",
    "dphilljetjet",
    "drll",
    "Lepton_eta[0]",
    "Lepton_eta[1]",
    "Lepton_pt[0]",
    "Lepton_pt[1]",
    "Lepton_phi[0]",
    "Lepton_phi[1]",  
    "CleanJet_eta[0]",
    "CleanJet_eta[1]",
    "CleanJet_eta[2]",
    "CleanJet_phi[0]",
    "CleanJet_phi[1]",
    "CleanJet_phi[2]",
    "abs(CleanJet_eta[2]-(CleanJet_eta[0]+CleanJet_eta[1])/2)*(CleanJet_pt[2]>30)",
    "CleanJet_pt[0]",
    "CleanJet_pt[1]",
    "MET_pt",
    "mth",
    "ptll",
    "mlljj(Sum$(CleanJet_pt>30), nLepton, CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1], Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1])",
    "mjjj(Sum$(CleanJet_pt>30), CleanJet_pt[0], CleanJet_pt[1], CleanJet_pt[2], CleanJet_phi[0], CleanJet_phi[1], CleanJet_phi[2], CleanJet_eta[0], CleanJet_eta[1], CleanJet_eta[2])",
    "Lepton_pt[0] + Lepton_pt[1] + CleanJet_pt[0] + CleanJet_pt[1] + MET_pt", 
    "mll",
    "RecoMELA_VBF(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1], Lepton_pdgId[0], Lepton_pdgId[1])",
    "RecoMELA_Q2V1(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",
    #"RecoMELA_Q2V2(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",
    "RecoMELA_Phi(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])"
    #"RecoMELA_CT1(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",
    #"RecoMELA_CT2(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",  
    ]
 

  ROOT.gROOT.ProcessLineSync(".L /eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/mlljj.C+")
  ROOT.gROOT.ProcessLineSync(".L /eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/mjjj.C+")

   
  ROOT.gROOT.ProcessLineSync('gSystem->Load("/eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/JHUGenMELA/MELA/data/slc7_amd64_gcc700/libmcfm_707.so","", kTRUE);')
  ROOT.gROOT.ProcessLineSync('gSystem->Load("/eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/lib/slc7_amd64_gcc700/libJHUGenMELAMELA.so","", kTRUE);')  

  ROOT.gROOT.ProcessLineSync(".L /eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/RecoMELA_VBF.C+")
  ROOT.gROOT.ProcessLineSync(".L /eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/RecoMELA_Q2V1.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_Q2V2.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_CT1.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_CT2.C+")
  ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_Phi.C+")


  chain = ROOT.TChain('Events')
  
  ### Full 2016 ###
    
  chain.Add(dir16+'nanoLatino_WWTo2L2Nu__part0.root')
  chain.Add(dir16+'nanoLatino_WWTo2L2Nu__part1.root')
  chain.Add(dir16+'nanoLatino_WWTo2L2Nu__part2.root')
  chain.Add(dir16+'nanoLatino_WWTo2L2Nu__part3.root')
  chain.Add(dir16+'nanoLatino_WWTo2L2Nu__part4.root')

  chain.Add(dir16+'nanoLatino_WpWmJJ_EWK__part0.root')
  chain.Add(dir16+'nanoLatino_WpWmJJ_EWK__part1.root')

  ### Full 2017 ###
  
  chain.Add(dir17+'nanoLatino_WWTo2L2Nu__part0.root')
  chain.Add(dir17+'nanoLatino_WWTo2L2Nu__part1.root')

  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part0.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part1.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part2.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part3.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part4.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part5.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part6.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part7.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part8.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part9.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part10.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part11.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part12.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part13.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part14.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part15.root')
  chain.Add(dir17+'nanoLatino_WpWmJJ_EWK__part16.root')
    
  ### Full 2018 ###
  
  chain.Add(dir18+'nanoLatino_WWTo2L2Nu__part0.root')
  chain.Add(dir18+'nanoLatino_WWTo2L2Nu__part1.root')
  chain.Add(dir18+'nanoLatino_WWTo2L2Nu__part2.root')
  chain.Add(dir18+'nanoLatino_WWTo2L2Nu__part3.root')
  chain.Add(dir18+'nanoLatino_WWTo2L2Nu__part4.root')
  chain.Add(dir18+'nanoLatino_WWTo2L2Nu__part5.root')
  chain.Add(dir18+'nanoLatino_WWTo2L2Nu__part6.root')
  chain.Add(dir18+'nanoLatino_WWTo2L2Nu__part7.root')
  chain.Add(dir18+'nanoLatino_WWTo2L2Nu__part8.root')
  chain.Add(dir18+'nanoLatino_WWTo2L2Nu__part9.root')
  
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part0.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part1.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part2.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part3.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part4.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part5.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part6.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part7.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part8.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part9.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part10.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part11.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part12.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part13.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part14.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part15.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part16.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part17.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part18.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part19.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part20.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part21.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part22.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part23.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part24.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part25.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part26.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part27.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part28.root')
  chain.Add(dir18+'nanoLatino_WpWmJJ_EWK__part29.root')

    
  print(chain.GetEntries())

  _dataset = root_numpy.tree2array (chain, 
      selection = 'Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13 && Lepton_pt[0] > 25. && Lepton_pt[1] > 13. && (abs(Lepton_pdgId[1]) == 13 || Lepton_pt[1] > 13.) && (nLepton >= 2 && Alt$(Lepton_pt[2], 0) < 10.) && Alt$(CleanJet_pt[1], 0) > 30. && abs(CleanJet_eta[0]) < 4.7 && abs(CleanJet_eta[1]) < 4.7 && mth > 40 && mth < 125 && Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217) == 0 && Sum$(CleanJet_pt>30) >= 2 && Sum$(CleanJet_pt>30) <= 3',
      branches = _branches,
      stop = max_entries
     )

  return { b : _dataset[b] for b in _branches }



#
#
# DY Monte Carlo
#
#


def load_dataset_dy ( max_entries = -1 ):
  _branches = [
    "mjj",
    "log((abs(2*Lepton_eta[0]-CleanJet_eta[0]-CleanJet_eta[1])+abs(2*Lepton_eta[1]-CleanJet_eta[0]-CleanJet_eta[1]))/detajj)",
    "Jet_qgl[0]",
    "Jet_qgl[1]",
    "detajj",
    "Lepton_eta[0]-Lepton_eta[1]",
    "sqrt(CleanJet_eta[0]*CleanJet_eta[1] + CleanJet_phi[0]*CleanJet_phi[1])",
    "dphill",
    "dphijjmet",
    "dphilljetjet",
    "drll",
    "Lepton_eta[0]",
    "Lepton_eta[1]",
    "Lepton_pt[0]",
    "Lepton_pt[1]",
    "Lepton_phi[0]",
    "Lepton_phi[1]",  
    "CleanJet_eta[0]",
    "CleanJet_eta[1]",
    "CleanJet_eta[2]",
    "CleanJet_phi[0]",
    "CleanJet_phi[1]",
    "CleanJet_phi[2]",
    "abs(CleanJet_eta[2]-(CleanJet_eta[0]+CleanJet_eta[1])/2)*(CleanJet_pt[2]>30)",
    "CleanJet_pt[0]",
    "CleanJet_pt[1]",
    "MET_pt",
    "mth",
    "ptll",
    "mlljj(Sum$(CleanJet_pt>30), nLepton, CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1], Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1])",
    "mjjj(Sum$(CleanJet_pt>30), CleanJet_pt[0], CleanJet_pt[1], CleanJet_pt[2], CleanJet_phi[0], CleanJet_phi[1], CleanJet_phi[2], CleanJet_eta[0], CleanJet_eta[1], CleanJet_eta[2])",
    "Lepton_pt[0] + Lepton_pt[1] + CleanJet_pt[0] + CleanJet_pt[1] + MET_pt", 
    "mll",
    "RecoMELA_VBF(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1], Lepton_pdgId[0], Lepton_pdgId[1])",
    "RecoMELA_Q2V1(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1], 0)",
    #"RecoMELA_Q2V1(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1], 1)",
    #"RecoMELA_Q2V2(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",
    #"RecoMELA_Phi(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])"
    #"RecoMELA_CT1(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",
    #"RecoMELA_CT2(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])",  
    ]
 

  ROOT.gROOT.ProcessLineSync(".L /eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/mlljj.C+")
  ROOT.gROOT.ProcessLineSync(".L /eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/mjjj.C+")

   
  ROOT.gROOT.ProcessLineSync('gSystem->Load("/eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/JHUGenMELA/MELA/data/slc7_amd64_gcc700/libmcfm_707.so","", kTRUE);')
  ROOT.gROOT.ProcessLineSync('gSystem->Load("/eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/lib/slc7_amd64_gcc700/libJHUGenMELAMELA.so","", kTRUE);')  

  ROOT.gROOT.ProcessLineSync(".L /eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/RecoMELA_VBF.C+")
  ROOT.gROOT.ProcessLineSync(".L /eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/RecoMELA_Q2V1.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_Q2V2.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_CT1.C+")
  #ROOT.gROOT.ProcessLineSync(".L ./RecoMELA_CT2.C+")
  #ROOT.gROOT.ProcessLineSync(".L /eos/user/s/sblancof/SWAN_projects/HWW_VBF/CMSSW_10_6_4/src/RecoMELA_Phi.C+")


  chain = ROOT.TChain('Events')
  
  ### Full 2016 ###
    
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part0.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part1.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part2.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part3.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part4.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part5.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part6.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part7.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part8.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part9.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part10.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part11.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part12.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part13.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part14.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part15.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part16.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part17.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part18.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part19.root')
  chain.Add(dir16+'nanoLatino_DYJetsToTT_MuEle_M-50_ext1__part20.root')

  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part0.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part1.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part2.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part3.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part4.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part5.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part6.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part7.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part8.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part9.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part10.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part11.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part12.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part13.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part14.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part15.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part16.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part17.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part18.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part19.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part20.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part21.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part22.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part23.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part24.root')
  chain.Add(dir16+'nanoLatino_DYJetsToLL_M-10to50-LO__part25.root')
    
  ### Full 2017 ###

  #chain.Add(dir17+'nanoLatino_DYJetsToTT_MuEle_M-50__part0.root') # For test
  
  for i in range(0, 73, 1):
    chain.Add(dir17+'nanoLatino_DYJetsToTT_MuEle_M-50__part'+str(i)+'.root')
    
    
  #chain.Add(dir17+'nanoLatino_DYJetsToLL_M-10to50-LO_ext1__part0.root') # For test
    
  for i in range(0, 76, 1):
       chain.Add(dir17+'nanoLatino_DYJetsToLL_M-10to50-LO_ext1__part'+str(i)+'.root')
    

  ### Full 2018 ###
    
  #chain.Add(dir18+'nanoLatino_DYJetsToTT_MuEle_M-50__part0.root') # For test 
    
  for i in range(0, 101, 1):
    chain.Add(dir18+'nanoLatino_DYJetsToTT_MuEle_M-50__part'+str(i)+'.root')
    
  
  #chain.Add(dir18+'nanoLatino_DYJetsToLL_M-10to50-LO_ext1__part0.root') # For test
    
  for i in range(0, 78, 1):
        chain.Add(dir18+'nanoLatino_DYJetsToLL_M-10to50-LO_ext1__part'+str(i)+'.root')
        
        
    
  print(chain.GetEntries())

  _dataset = root_numpy.tree2array (chain, 
      selection = 'Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13 && Lepton_pt[0] > 25. && Lepton_pt[1] > 13. && (abs(Lepton_pdgId[1]) == 13 || Lepton_pt[1] > 13.) && (nLepton >= 2 && Alt$(Lepton_pt[2], 0) < 10.) && Alt$(CleanJet_pt[1], 0) > 30. && abs(CleanJet_eta[0]) < 4.7 && abs(CleanJet_eta[1]) < 4.7 && mth > 40 && mth < 125 && Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217) == 0 && Sum$(CleanJet_pt>30) >= 2 && Sum$(CleanJet_pt>30) <= 3',
      branches = _branches,
      stop = max_entries
     )

  return { b : _dataset[b] for b in _branches }


if __name__ == '__main__':
    print(load_dataset_vbf(10))
    #print(load_dataset_ggh(10))
    #print(load_dataset_top(10))
    #print(load_dataset_ww(10))
    #print(load_dataset_dy(10))
    
    VARS = ["mjj", "log((abs(2*Lepton_eta[0]-CleanJet_eta[0]-CleanJet_eta[1])+abs(2*Lepton_eta[1]-CleanJet_eta[0]-CleanJet_eta[1]))/detajj)", 
            "Jet_qgl[0]", "Jet_qgl[1]", "detajj", "Lepton_eta[0]-Lepton_eta[1]", "sqrt(CleanJet_eta[0]*CleanJet_eta[1] + CleanJet_phi[0]*CleanJet_phi[1])",
            "dphill", "dphijjmet", "dphilljetjet", "drll", "Lepton_eta[0]", "Lepton_eta[1]", "Lepton_pt[0]", "Lepton_pt[1]", "Lepton_phi[0]", "Lepton_phi[1]", 
            "CleanJet_eta[0]", "CleanJet_eta[1]", "CleanJet_eta[2]", "CleanJet_phi[0]", "CleanJet_phi[1]", "CleanJet_phi[2]", "abs(CleanJet_eta[2]-(CleanJet_eta[0]+CleanJet_eta[1])/2)*(CleanJet_pt[2]>30)",
            "CleanJet_pt[0]", "CleanJet_pt[1]", "MET_pt", "mth", "ptll", 
            "mlljj(Sum$(CleanJet_pt>30), nLepton, CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1], Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1])",
            "mjjj(Sum$(CleanJet_pt>30), CleanJet_pt[0], CleanJet_pt[1], CleanJet_pt[2], CleanJet_phi[0], CleanJet_phi[1], CleanJet_phi[2], CleanJet_eta[0], CleanJet_eta[1], CleanJet_eta[2])", 
            "Lepton_pt[0] + Lepton_pt[1] + CleanJet_pt[0] + CleanJet_pt[1] + MET_pt",  "mll", 
            "RecoMELA_VBF(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1], Lepton_pdgId[0], Lepton_pdgId[1])", 
            "RecoMELA_Q2V1(Sum$(CleanJet_pt>30), nLepton, PuppiMET_pt, PuppiMET_phi, Lepton_pt[0], Lepton_pt[1], Lepton_phi[0], Lepton_phi[1], Lepton_eta[0], Lepton_eta[1], CleanJet_pt[0], CleanJet_pt[1], CleanJet_phi[0], CleanJet_phi[1], CleanJet_eta[0], CleanJet_eta[1],  Lepton_pdgId[0], Lepton_pdgId[1])"
    ]

    NDIM = len(VARS)

    '''
    print("Starting VBF daatsets \n")
    
    dataset_vbf = load_dataset_vbf(-1)
    
    print("Closing VBF datasets \n")
    
    print("Starting ggH datasets \n")
    
    dataset_ggh = load_dataset_ggh(-1)
    
    print("Closing ggH datasets \n")
    
    print("Starting Top datasets \n")
    
    
    dataset_top = load_dataset_top(-1)
    
    print("Closing Top datasets \n")
    '''
    print("Starting WW datasets \n")
    
    dataset_ww = load_dataset_ww(-1)
    
    print("Closing WW datasets \n")
    '''
    print("Starting DY datasets \n")
    
    dataset_dy = load_dataset_dy(-1)
    
    print("Closing DY datasets \n")
    '''

    #df_vbf = pd.DataFrame(dataset_vbf,columns=VARS)
    #df_ggh = pd.DataFrame(dataset_ggh,columns=VARS)
    #df_top = pd.DataFrame(dataset_top,columns=VARS)
    df_ww = pd.DataFrame(dataset_ww,columns=VARS)
    #df_dy = pd.DataFrame(dataset_dy,columns=VARS)
    
    '''
    print("Moving datasets to pickle files \n")
    
    df_vbf.to_pickle('dataset_vbf.pkl')
    
    print("First done") 
    
    df_ggh.to_pickle('dataset_ggh.pkl')
    
    print("Second done") 
    
    df_top.to_pickle('dataset_top.pkl')
    '''
    print("Third done") 
    
    df_ww.to_pickle('dataset_ww.pkl')
    '''
    print("Fourth done") 
    
    df_dy.to_pickle('dataset_dy.pkl')
    
    print("All done") 
    ''' 
    
