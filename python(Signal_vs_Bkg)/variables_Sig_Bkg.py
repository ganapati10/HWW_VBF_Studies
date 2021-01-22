variables['mjj'] = {'name'  : 'mjj',
                    'range' : (16, 0., 400.),
                    'xaxis' : 'm_{jj} [GeV]',
                    'fold'  : 0}

variables['mll'] = {      'name'  : 'mll',
                          'range' : (8, 0., 200.),
                          'xaxis' : 'm_{ll} [GeV]',
                          'fold'  : 0}

variables['mth'] = {      'name'  : 'mth',
                          'range' : (10, 0., 200.),
                          'xaxis' : 'm_{T}^{H} [GeV]',
                          'fold'  : 0}


variables['njet'] = {     'name'  : 'Sum$(CleanJet_pt>30)',     
                          'range' : (5, 0, 5),   
                          'xaxis' : 'number of jets',
                          'fold'  : 2}

variables['bjet'] = {     'name'  : 'Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217)',     
                          'range' : (5, 0, 5),   
                          'xaxis' : 'number of jets',
                          'fold'  : 2}

variables['Ctot'] =  {
                          'name': 'log((abs(2*Lepton_eta[0]-CleanJet_eta[0]-CleanJet_eta[1])+abs(2*Lepton_eta[1]-CleanJet_eta[0]-CleanJet_eta[1]))/detajj)',
                          'range' : (10, -4, 6),
                          'xaxis' : 'Ctot',
                          'fold'  : 3,}

variables['detajj'] = {   'name'  : 'detajj',
                          'range' : (10, 0., 4.),
                          'xaxis' : '|#Delta#eta_{jj}|',
                          'fold'  : 0}

variables['detal1j1'] = { 'name'  : 'abs(Lepton_eta[0]-CleanJet_eta[0])',
                          'range' : (10, 0., 4.),
                          'xaxis' : '#Delta#eta_{l1j1}',
                          'fold'  : 0}

variables['detall'] = {   'name'  : 'abs(Lepton_eta[1]-Lepton_eta[0])',
                          'range' : (10, 0., 4.),
                          'xaxis' : '#Delta#eta_{ll}',
                          'fold'  : 0}

variables['drll'] = {     'name'  : 'drll',
                          'range' : (10, 0., 5.),
                          'xaxis' : '#DeltaR_{ll}',
                          'fold'  : 2}

variables['pt1'] = {      'name'  : 'Lepton_pt[0]',     
                          'range' : (20, 0., 200.),   
                          'xaxis' : 'p_{T} 1st lepton [GeV]',
                          'fold'  : 3}

variables['jetpt1'] = {   'name'  : 'CleanJet_pt[0]*(CleanJet_pt[0]>30)',     
                          'range' : (20, 0., 200.),   
                          'xaxis' : 'p_{T} 1st jet',
                          'fold'  : 2}

variables['ptll'] = {     'name'  : 'ptll',
                          'range' : (10, 30., 200.),
                          'xaxis' : 'p_{T}^{ll} [GeV]',
                          'fold'  : 0}

variables['dphill'] = {   'name'  : 'abs(dphill)',     
                          'range' : (10, 0., 3.2),   
                          'xaxis' : '#Delta#phi_{ll}',
                          'fold'  : 3}

variables['dphillj'] = { 'name'  : 'abs(dphilljet)',     
                          'range' : (10, 0., 3.2),   
                          'xaxis' : '#Delta#phi_{llj}',
                          'fold'  : 3}

variables['dphilljj'] = { 'name'  : 'abs(dphilljetjet)',     
                          'range' : (10, 0., 3.2),   
                          'xaxis' : '#Delta#phi_{lljj}',
                          'fold'  : 3}

variables['dphijjmet'] = { 'name'  : 'abs(dphijjmet)',     
                          'range' : (10, 0., 3.2),   
                          'xaxis' : '#Delta#phi_{jjmet}',
                          'fold'  : 3}

variables['qgl0'] = {     'name'  : 'Jet_qgl[0]',
                          'range' : (5, 0., 1.),
                          'xaxis' : 'Quark Gluon likelihood (1^{st} Jet)',
                          'fold'  : 3}


variables['qgl1'] = {     'name'  : 'Jet_qgl[1]',
                          'range' : (5, 0., 1.),
                          'xaxis' : 'Quark Gluon likelihood (2^{nd} Jet)',
                          'fold'  : 3}


variables['jeteta1'] = {  'name'  : 'abs(CleanJet_eta[0])',
                          'range' : (10, 0., 4.),
                          'xaxis' : '#eta 1st jet',
                          'fold'  : 0}

variables['eta1'] = {     'name'  : 'abs(Lepton_eta[0])',
                          'range' : (10, 0., 4.),
                          'xaxis' : '#eta 1st lepton',
                          'fold'  : 0}

variables['mpmet'] = {    'name'  : 'mpmet',      
                          'range' : (10, 0., 150.),  
                          'xaxis' : 'min. (proj. tk. E_{T}^{miss}, proj. E_{T}^{miss}) [GeV]', 
                          'fold'  : 3}

variables['pfmet'] = {    'name'  : 'MET_pt',     
                          'range' : (10, 0., 150.),   
                          'xaxis' : 'PF MET [GeV]',
                          'fold'  : 3}

variables['puppimet'] = { 'name'  : 'PuppiMET_pt',
                          'range' : (10, 0., 150.),
                          'xaxis' : 'puppi MET [GeV]',
                          'fold'  : 3}

variables['rawmet'] = {   'name'  : 'RawMET_pt',     
                          'range' : (10, 0., 150.),   
                          'xaxis' : 'raw MET [GeV]',
                          'fold'  : 3}

variables['TkMET'] = {    'name'  : 'TkMET_pt',     
                          'range' : (10, 0., 150.),   
                          'xaxis' : 'tracker MET [GeV]',
                          'fold'  : 3}

variables['btagDeepB'] = { 'name'  : 'Jet_btagDeepB[CleanJet_jetIdx]',
                          'range' : (10, 0., 1.),
                          'xaxis' : 'Deep B discriminator',
                          'fold'  : 3}

variables['btagCSVv2'] = { 'name'  : 'Jet_btagCSVV2[CleanJet_jetIdx]',
                          'range' : (10, 0., 1.),
                          'xaxis' : 'CSVv2 discriminator',
                          'fold'  : 3}

variables['btagCMVA'] = { 'name'  : 'Jet_btagCMVA[CleanJet_jetIdx]',
                          'range' : (10, -1., 1.),
                          'xaxis' : 'CMVA discriminator',
                          'fold'  : 3}

variables['btagDeepFlavB'] = { 'name'  : 'Jet_btagDeepFlavB[CleanJet_jetIdx]',
                          'range' : (10, 0., 1.),
                          'xaxis' : 'Deep Flavour B discriminator',
                          'fold'  : 3}

variables['ptlljjmet'] = { 'name'  : 'Lepton_pt[0] + Lepton_pt[1] + CleanJet_pt[0] + CleanJet_pt[1] + MET_pt',
                          'range' : (16, 0., 400.),
                          'xaxis' : 'p_{T} leptons+jets+MET',
                          'fold'  : 3}

variables['mlljj'] = {    'name'       : 'mlljj(Lepton_pt[0], Lepton_eta[0], Lepton_phi[0], Lepton_pt[1], Lepton_eta[1], Lepton_phi[1], CleanJet_pt[0], CleanJet_eta[0], CleanJet_phi[0], CleanJet_mass[0], CleanJet_pt[1], CleanJet_eta[1], CleanJet_phi[1], CleanJet_mass[1])',
                          'range'      : (16, 0., 400.),
                          'xaxis'      : 'm_{lljj}',
                          'fold'       : 3,
                          'linesToAdd' : ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_4/src/PlotsConfigurations/Configurations/WW/Full2016_v6/extended/mlljj.C++']}













