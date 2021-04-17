#VBF cuts                                                                                                                                                                                               


#-------------------------------------------------------------------------------                                                                                                                           
# supercut                                                                                                                                                                                                 
#-------------------------------------------------------------------------------                                                                                                                           
_tmp = [
     'Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13',
     'Lepton_pt[0] > 25.',
     'Lepton_pt[1] > 13.',
     '(abs(Lepton_pdgId[1]) == 13 || Lepton_pt[1] > 13.)',
     '(nLepton >= 2 && Alt$(Lepton_pt[2], 0) < 10.)'
     ]

supercut = ' && '.join(_tmp)


def addcut(name, exprs):
    cuts[name] = ' && '.join(exprs)

#-------------------------------------------------------------------------------                                                                                                                           
# VBF                                                                                                                                                                                                       
#-------------------------------------------------------------------------------                                                                                                                           
_tmp = [
     'Alt$(CleanJet_pt[1], 0) > 30.',
     'abs(CleanJet_eta[0]) < 4.5',
     'abs(CleanJet_eta[1]) < 4.5',
     'mth > 40 && mth < 125',
     'drll < 2.5',
     'mjj > 200',
     'Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217) == 0',
     'Sum$(CleanJet_pt>30) >= 2 && Sum$(CleanJet_pt>30) <= 3',
     'mll < 75',
     'abs(CleanJet_eta[0]) > 1.5 && abs(CleanJet_eta[1]) > 1.5',
     'abs(Lepton_eta[0]) < 1.5 && abs(Lepton_eta[1]) < 1.5',
     'log((abs(2*Lepton_eta[0]-CleanJet_eta[0]-CleanJet_eta[1])+abs(2*Lepton_eta[1]-CleanJet_eta[0]-CleanJet_eta[1]))/detajj) < 1.'
     ]

addcut('VBF', _tmp)


#-------------------------------------------------------------------------------                                                                                                                           
# VBF_Bkg DY                                                                                                                                                                                                
#-------------------------------------------------------------------------------                                                                                                                           
_tmp = [
     'Alt$(CleanJet_pt[1], 0) > 30.',
     'abs(CleanJet_eta[0]) < 4.5',
     'abs(CleanJet_eta[1]) < 4.5',
     'mth < 40',
     'Sum$(CleanJet_pt>20) >= 2',
     'Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217) == 0',
     'detall < 1.6',
     'mll < 80'
     ]

addcut('DY', _tmp)


#-------------------------------------------------------------------------------                                                                                                                           
# VBF_Bkg top                                                                                                                                                                                              
#-------------------------------------------------------------------------------                                                                                                                           
_tmp = [
     'Alt$(CleanJet_pt[1], 0) > 30.',
     'abs(CleanJet_eta[0]) < 4.5',
     'abs(CleanJet_eta[1]) < 4.5',
     'Sum$(CleanJet_pt>30) >= 1',
     'Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217) >= 1',
     'mll > 80'
     ]

addcut('top', _tmp)
