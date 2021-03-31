# Varaibles for fit

variables['Ctot'] = {     'name': 'log((abs(2*Lepton_eta[0]-CleanJet_eta[0]-CleanJet_eta[1])+abs(2*Lepton_eta[1]-CleanJet_eta[0]-CleanJet_eta[1]))/detajj)',
                          #'range' : (10, -4, 2),
                          'range' : ([-4., -2.34102564, -1.82820513, -1.37179487, -0.6405, -0.08974359, 1.],),
                          'xaxis' : 'Ctot',
                          'fold'  : 3}

variables['Ctot_assym'] = {     'name': 'log((abs(2*Lepton_eta[0]-CleanJet_eta[0]-CleanJet_eta[1])+abs(2*Lepton_eta[1]-CleanJet_eta[0]-CleanJet_eta[1]))/detajj)',
                                'range' : ([-4., -2., -1.5, -1.0, -0.5, 0., 1.],),
                                'xaxis' : 'Ctot',
                                'fold'  : 3}

variables['dphilljj'] = { 'name'  : 'abs(dphilljetjet)',     
                          #'range' : (10, 0., 3.2),   
                          'range' : ([0., 1.18974359, 1.76410256, 2.09230769, 2.42051282, 2.58461538, 2.91282051, 3.2],),
                          'xaxis' : '#Delta#phi_{lljj}',
                          'fold'  : 3}

variables['events'] = {   'name'  : '1',
                          'range' : (1, 0, 2),
                          'xaxis' : 'events',
                          'fold'  : 3}

variables['ptll'] = {     'name'  : 'ptll',
                          #'range' : (10, 30., 200.),
                          'range' : ([ 30., 36.53846154, 49.61538462, 71.41025641, 93.20512821, 128.07692308, 200.],),
                          'xaxis' : 'p_{T}^{ll} [GeV]',
                          'fold'  : 3}

variables['ptll_assym'] = {     'name'  : 'ptll',
                                #'range' : ([30., 47., 64., 81., 98., 115., 132., 149., 200.],),
                                'range' : ([30., 47., 64., 81., 98., 115., 132., 149., 166., 200.],),
                                'xaxis' : 'p_{T}^{ll} [GeV]',
                                'fold'  : 3}

variables['mll'] = {      'name'  : 'mll',
                          #'range' : (10, 0., 80.),
                          'range': ([ 0., 13.46153846, 26.02564103, 52.94871795, 61.92307692, 70., 80.],),
                          'xaxis' : 'm_{ll} [GeV]',
                          'fold'  : 3}

variables['mll_assym'] = {      'name'  : 'mll',
                                'range' : ([0., 16., 24., 32., 40., 48., 56., 64., 72., 80.],),
                                'xaxis' : 'm_{ll} [GeV]',
                                'fold'  : 3}

variables['mll_top'] = {      'name'  : 'mll',
                              'range' : (15, 80., 200.),
                              'xaxis' : 'm_{ll} [GeV]',
                              'fold'  : 3}

variables['mjj'] = {      'name'  : 'mjj',
                          #'range' : (15, 200., 4000.),
                          'range' : ([400., 630.76923077, 907.6921, 1284.61538462, 1646.15384615, 4000.],),
                          'xaxis' : 'm_{jj} [GeV]',
                          'fold'  : 3}

variables['drjj'] = {     'name'  : 'sqrt(CleanJet_eta[0]*CleanJet_eta[1] + CleanJet_phi[0]*CleanJet_phi[1])',
                          #'range' : (10, 0., 5.),
                          'range': ([0., 0.96153846, 1.73076923, 2.37179487, 3.52564103, 3.91025641, 5.],),
                          'xaxis' : '#DeltaR_{jj}',
                         #'linesToAdd': ['.L $CMSSW_BASE/src/PlotsConfigurations/Configurations/WW/Full2016_v6/extended/drjj.C+'], #if want to use a script
                          'fold'  : 3}

variables['drjj_assymetric'] = {    'name'  : 'sqrt(CleanJet_eta[0]*CleanJet_eta[1] + CleanJet_phi[0]*CleanJet_phi[1])',
                                    'range' : ([0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0],),
                                    'xaxis' : '#DeltaR_{jj}',
                                    #'linesToAdd': ['.L $CMSSW_BASE/src/PlotsConfigurations/Configurations/WW/Full2016_v6/extended/drjj.C+'], #if want to use a script
                                    'fold'  : 3}

variables['ptlljjmet'] = { 'name'  : 'Lepton_pt[0] + Lepton_pt[1] + CleanJet_pt[0] + CleanJet_pt[1] + MET_pt',
                           #'range' : (16, 137.5, 800.),
                           'range': ([137.5, 179.96794872, 213.94230769, 298.8782051, 374.83974359, 451.76282051, 556.69871795, 800.],),
                           'xaxis' : 'p_{T} leptons+jets+MET',
                           'fold'  : 3}
