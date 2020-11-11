 # VH2j cuts


#-------------------------------------------------------------------------------
# supercut
#-------------------------------------------------------------------------------
_tmp = [
     'Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13',
     'Lepton_pt[0] > 25.',
     'Lepton_pt[1] > 13.',
     '(abs(Lepton_pdgId[1]) == 13 || Lepton_pt[1] > 13.)',
     '(nLepton >= 2 && Alt$(Lepton_pt[2], 0) < 10.)',
     'mll > 12.',
     'ptll > 30.',
     'PuppiMET_pt > 20.', 
     'mjj > 200',
     ]

supercut = ' && '.join(_tmp)


def addcut(name, exprs):
    cuts[name] = ' && '.join(exprs)


#-------------------------------------------------------------------------------
# VH_2j_em
#-------------------------------------------------------------------------------
_tmp = [
     'Alt$(CleanJet_pt[1], 0) > 20.', 
     'abs(CleanJet_eta[0]) < 2.5',
     'abs(CleanJet_eta[1]) < 2.5',
     'mth > 60.',
     'mth < 125.',
     'drll < 2.',
     #'detajj < 3.5',
     'bVeto',
     ]

addcut('VH_2j_emu', _tmp)


#-------------------------------------------------------------------------------
# VH_2j_topemu
#-------------------------------------------------------------------------------
_tmp = [
     'Alt$(CleanJet_pt[1], 0) > 20.',  
     'abs(CleanJet_eta[0]) < 2.5',
     'abs(CleanJet_eta[1]) < 2.5',
     #'detajj < 3.5', 
     'bReq',
     'mll > 50.',
     ]

addcut('VH_2j_topemu', _tmp)


#-------------------------------------------------------------------------------
# VH_2j_DYtautau
#-------------------------------------------------------------------------------
_tmp = [
     'Alt$(CleanJet_pt[1], 0) > 20.',  
     'abs(CleanJet_eta[0]) < 2.5',
     'abs(CleanJet_eta[1]) < 2.5',
     'mth < 60.',
     'drll < 2.',
     'bVetoDY',
     #'detajj < 3.5',
     'mll > 40.',
     'mll < 80.',
     ]

addcut('VH_2j_DYtautau', _tmp)


#-------------------------------------------------------------------------------
# Test
#-------------------------------------------------------------------------------
### _tmp = [
###      'VH2j_TMVAReader(Entry$) > 0.1',
###      ]
### 
### addcut('VH_2j_test', _tmp);
