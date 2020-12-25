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

variables['dphilljj'] = { 'name'  : 'abs(dphilljetjet)',     
                          'range' : (10, 0., 3.2),   
                          'xaxis' : '#Delta#phi_{lljj}',
                          'fold'  : 3}

variables['qgl0'] = {     'name'  : 'Jet_qgl[0]',
                          'range' : (5, 0., 1.),
                          'xaxis' : 'Quark Gluon likelihood (1^{st} Jet)',
                          'fold'  : 3}

variables['qgl01'] = {    'name'  : 'Jet_qgl[0]',
                          'range' : (5, 0., 1.),
                          'xaxis' : 'Quark Gluon likelihood (1^{st} Jet)',
                          'fold'  : 3}

variables['qgl1'] = {     'name'  : 'Jet_qgl[1]',
                          'range' : (5, 0., 1.),
                          'xaxis' : 'Quark Gluon likelihood (2^{nd} Jet)',
                          'fold'  : 3}

variables['qgl11'] = {    'name'  : 'Jet_qgl[1]',
                          'range' : (5, 0., 1.),
                          'xaxis' : 'Quark Gluon likelihood (2^{nd} Jet)',
                          'fold'  : 3}




