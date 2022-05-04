import ROOT
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class HWW_VBF_DNN(Module):
  
    def __init__(self, sample):
        print '####################', sample
        self.sample = sample
        self.cmssw_base = os.getenv('CMSSW_BASE')
        self.cmssw_arch = os.getenv('SCRAM_ARCH')

        ROOT.gSystem.AddIncludePath("-I"+self.cmssw_base+"/interface/")
        ROOT.gSystem.AddIncludePath("-I"+self.cmssw_base+"/src/")
        ROOT.gSystem.Load(self.cmssw_base+"/lib/libmomemta.so.1.0.0")
        ROOT.gSystem.Load(self.cmssw_base+"/lib/libmomemta.so.1")
        ROOT.gSystem.Load(self.cmssw_base+"/lib/libmomemta.so")
        
        ROOT.gSystem.Load(self.cmssw_base+"/src/JHUGenMELA/MELA/data/slc7_amd64_gcc820/libmcfm_707.so")
        ROOT.gSystem.Load("libJHUGenMELAMELA.so")
        

        try:
            ROOT.gROOT.LoadMacro(self.cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/recoMoMEMta.C+g')
            ROOT.gROOT.LoadMacro(self.cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/recoMELA_VBF.C+g')
            ROOT.gROOT.LoadMacro(self.cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/recoMELA_VBF_VH.C+g')
            ROOT.gROOT.LoadMacro(self.cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/recoMELA_QCD_VH.C+g')
            ROOT.gROOT.LoadMacro(self.cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/mlj.C+g')
        except RuntimeError:
            ROOT.gROOT.LoadMacro(self.cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/recoMoMEMta.C++g')
            ROOT.gROOT.LoadMacro(self.cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/recoMELA_VBF.C++g')
            ROOT.gROOT.LoadMacro(self.cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/recoMELA_VBF_VH.C++g')
            ROOT.gROOT.LoadMacro(self.cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/recoMELA_QCD_VH.C++g')
            ROOT.gROOT.LoadMacro(self.cmssw_base+'/src/LatinoAnalysis/Gardener/python/variables/mlj.C++g')

        self.bVetoCut = 0.2217
        self.metpt = 'event.MET_pt'
        self.metphi = 'event.MET_phi'
        self.nlepton = 'event.nLepton'
        self.njets = 'event.nCleanJet'
        self.lep1pt = 'event.Lepton_pt[0]'
        self.lep2pt = 'event.Lepton_pt[1]'
        self.lep1eta = 'event.Lepton_eta[0]'
        self.lep2eta = 'event.Lepton_eta[1]'
        self.lep1phi = 'event.Lepton_phi[0]'
        self.lep2phi = 'event.Lepton_phi[1]'
        self.jet1pt = 'event.CleanJet_pt[0]'
        self.jet2pt = 'event.CleanJet_pt[1]'
        self.jet1eta = 'event.CleanJet_eta[0]'
        self.jet2eta = 'event.CleanJet_eta[1]'
        self.jet1phi = 'event.CleanJet_phi[0]'
        self.jet2phi = 'event.CleanJet_phi[1]'
        
   def beginJob(self):
       pass
   def endJob(self):
       pass
   def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        
        self.inputFile = inputFile

        self.out = wrappedOutputTree        
        self.out.branch('mjj','F')
        self.out.branch('Ctot','F')
        self.out.branch('detajj','F')
        self.out.branch('drll','F')
        self.out.branch('jet1eta','F')
        self.out.branch('jet2eta','F')
        self.out.branch('PuppiMET_pt','F')
        self.out.branch('PuppiMET_phi','F')
        self.out.branch('mth','F')
        self.out.branch('ptll','F')
        self.out.branch('mlj_00','F')
        self.out.branch('mlj_01','F')
        self.out.branch('mlj_10','F')
        self.out.branch('mlj_11','F')
        self.out.branch('mll','F')
        self.out.branch('btagDeepB','F')
        self.out.branch('D_VBF_QCD','F')
        self.out.branch('D_VBF_VH','F')
        self.out.branch('D_QCD_VH','F')
        self.out.branch('D_VBF_DY', 'F')

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        CleanJet = Collection(event, 'CleanJet')
        Lepton = Collection(event, 'Lepton')
        njets = CleanJet._len
        nlepton = Lepton._len
        
        D_VBF_QCD = 0.0
        D_VBF_VH = 0.0
        D_QCD_VH = 0.0
        D_VBF_DY = 0.0
        
        if (nlepton>=2 or njets>=2):
          
          D_VBF_DY = ROOT.recoMoMEMta(njets, nlepton, event.MET_pt, event.MET_phi, Leptons[0].pt, Leptons[1].pt, Leptons[0].phi, Leptons[1].phi, Leptons[0].eta, Leptons[1].eta, CleanJets[0].pt, CleanJets[1].pt, CleanJets[0].phi, CleanJets[1].phi, CleanJets[0].eta, CleanJets[1].eta, Leptons[0].pdgId, Leptons[1].pdgId)
          D_VBF_QCD = ROOT.recoMELA_VBF(njets, nlepton, event.MET_pt, event.MET_phi, Leptons[0].pt, Leptons[1].pt, Leptons[0].phi, Leptons[1].phi, Leptons[0].eta, Leptons[1].eta, CleanJets[0].pt, CleanJets[1].pt, CleanJets[0].phi, CleanJets[1].phi, CleanJets[0].eta, CleanJets[1].eta, Leptons[0].pdgId, Leptons[1].pdgId)
          D_VBF_VH = ROOT.recoMELA_VBF_VH(njets, nlepton, event.MET_pt, event.MET_phi, Leptons[0].pt, Leptons[1].pt, Leptons[0].phi, Leptons[1].phi, Leptons[0].eta, Leptons[1].eta, CleanJets[0].pt, CleanJets[1].pt, CleanJets[0].phi, CleanJets[1].phi, CleanJets[0].eta, CleanJets[1].eta, Leptons[0].pdgId, Leptons[1].pdgId)
          D_QCD_VH = ROOT.recoMELA_QCD_VH(njets, nlepton, event.MET_pt, event.MET_phi, Leptons[0].pt, Leptons[1].pt, Leptons[0].phi, Leptons[1].phi, Leptons[0].eta, Leptons[1].eta, CleanJets[0].pt, CleanJets[1].pt, CleanJets[0].phi, CleanJets[1].phi, CleanJets[0].eta, CleanJets[1].eta, Leptons[0].pdgId, Leptons[1].pdgId)
          
          mlj_00 = mlj(CleanJet[0].pt, CleanJet[0].phi, CleanJet[0].eta, Lepton[0].pt,  Lepton[0].phi, Lepton[0].eta)
          mlj_01 = mlj(CleanJet[0].pt, CleanJet[0].phi, CleanJet[0].eta, Lepton[1].pt,  Lepton[1].phi, Lepton[1].eta)
          mlj_10 = mlj(CleanJet[1].pt, CleanJet[1].phi, CleanJet[1].eta, Lepton[0].pt,  Lepton[0].phi, Lepton[0].eta)
          mlj_11 = mlj(CleanJet[1].pt, CleanJet[1].phi, CleanJet[1].eta, Lepton[1].pt,  Lepton[1].phi, Lepton[1].eta)
          
        else:
          D_VBF_QCD = -999.9
          D_VBF_VH = -999.9
          D_QCD_VH = -999.9
          D_VBF_DY = -999.9
          
          
        if (njets==0):
          Jet_btagDeepB_CleanJet_jetIdx_0_ = -2
        elif (njets == 1):
          jetIdx0 = CleanJet[0].jetIdx
          Jet_btagDeepB_CleanJet_jetIdx_0_ = event.Jet_btagDeepB[jetIdx0]
        else:
          jetIdx0 = CleanJet[0].jetIdx
          Jet_btagDeepB_CleanJet_jetIdx_0_ = event.Jet_btagDeepB[jetIdx0]
        

          
        ctot = np.log((abs(2*Lepton[0].eta-CleanJet[0].eta-CleanJet[1].eta)+abs(2*Lepton[1].eta-CleanJet[0].eta-CleanJet[1].eta))/event.detajj)
          
        self.out.fillBranch('mjj', event.mjj)
        self.out.fillBranch('Ctot', ctot)
        self.out.fillBranch('detajj', event.detajj)
        self.out.fillBranch('drll', event.drll)
        self.out.fillBranch('jet1eta', self.jet1eta)
        self.out.fillBranch('jet2eta', self.jet2eta)
        self.out.fillBranch('PuppiMET_pt', event.PuppiMet_pt)
        self.out.fillBranch('PuppiMET_phi', event.PuppiMet_phi)
        self.out.fillBranch('mth', event.mth)
        self.out.fillBranch('ptll', event.ptll)
        self.out.fillBranch('mlj_00', mlj_00)
        self.out.fillBranch('mlj_01', mlj_01)
        self.out.fillBranch('mlj_10', mlj_10)
        self.out.fillBranch('mlj_11', mlj_11)
        self.out.fillBranch('mll', event.mll)
        self.out.fillBranch('btagDeepB', Jet_btagDeepB_CleanJet_jetIdx_0_)
        self.out.fillBranch('D_VBF_QCD', D_VBF_QCD)
        self.out.fillBranch('D_VBF_VH', D_VBF_VH)
        self.out.fillBranch('D_QCD_VH', D_QCD_VH)
        self.out.fillBranch('D_VBF_DY', D_VBF_DY)

        return True
        

  
  
