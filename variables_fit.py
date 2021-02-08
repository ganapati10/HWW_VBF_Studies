# Varaibles for fit

variables['Ctot'] = {     'name': 'log((abs(2*Lepton_eta[0]-CleanJet_eta[0]-CleanJet_eta[1])+abs(2*Lepton_eta[1]-CleanJet_eta[0]-CleanJet_eta[1]))/detajj)',
                          'range' : (20, -4, 6),
                          'xaxis' : 'Ctot',
                          'fold'  : 3}

variables['dphilljj'] = { 'name'  : 'abs(dphilljetjet)',     
                          'range' : (20, 0., 3.2),   
                          'xaxis' : '#Delta#phi_{lljj}',
                          'fold'  : 3}

variables['events'] = {   'name'  : '1',
                          'range' : (1, 0, 2),
                          'xaxis' : 'events',
                          'fold'  : 3}

variables['ptll'] = {     'name'  : 'ptll',
                          'range' : (40, 30., 200.),
                          'xaxis' : 'p_{T}^{ll} [GeV]',
                          'fold'  : 0}

variables['mll'] = {      'name'  : 'mll',
                          'range' : (8, 0., 200.),
                          'xaxis' : 'm_{ll} [GeV]',
                          'fold'  : 0}

variables['drjj'] = {     'name'  : 'sqrt(CleanJet_eta[0]*CleanJet_eta[1] + CleanJet_phi[0]*CleanJet_phi[1])',
                         'range' : (30, 0., 5.),
                         'xaxis' : '#DeltaR_{jj}',
                         #'linesToAdd': ['.L $CMSSW_BASE/src/PlotsConfigurations/Configurations/WW/Full2016_v6/extended/drjj.C+'], #if want to use a script
                         'fold'  : 2}

variables['ptlljjmet'] = { 'name'  : 'Lepton_pt[0] + Lepton_pt[1] + CleanJet_pt[0] + CleanJet_pt[1] + MET_pt',
                          'range' : (16, 0., 400.),
                          'xaxis' : 'p_{T} leptons+jets+MET',
                          'fold'  : 3}





