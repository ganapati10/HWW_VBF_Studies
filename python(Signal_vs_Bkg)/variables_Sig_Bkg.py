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




