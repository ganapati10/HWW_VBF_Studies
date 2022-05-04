  'HWW_VBF_DNN' : {
                  'isChain'  : True ,
                  'do4MC'    : True  ,
                  'do4Data'  : True ,
                  'selection': '"(nLepton >= 2 && \                                                                                                                                                                                           
                               Alt$(Lepton_pt[2],0) < 10. && \                                                                                                                                                                                
                               Lepton_pt[0] > 25 && \                                                                                                                                                                                         
                               Lepton_pt[1] > 15 && \                                                                                                                                                                                         
                               Lepton_pdgId[0]*Lepton_pdgId[1] == -11*13 && \
                               Sum$(CleanJet_pt>30) >= 2)"',
                  'subTargets' : ['HWW_VBF_DNN'],
                  'outputbranchsel': os.getenv('CMSSW_BASE') + '/src/LatinoAnalysis/NanoGardener/python/data/HWW_VBF_DNN_branches.txt',
               },
