# VH2j configuration file                                                                                                                                                                                   

treeName = 'Events'

tag = 'WW_2016'

# used by mkShape to define output directory for root files                                                                                                                                                 
outputDir = 'rootFile_Sig_Bkg'

# file with TTree aliases                                                                                                                                                                                   
aliasesFile = 'aliases.py'

# file with list of variables                                                                                                                                                                               
variablesFile = 'variables_Sig_Bkg.py'

# file with list of cuts                                                                                                                                                                                    
cutsFile = 'cuts_Sig_Bkg.py'

# file with list of samples                                                                                                                                                                                 
samplesFile = 'samples_ONE.py'

# file with plot configuration                                                                                                                                                                              
plotFile = 'Plot_Sig_Bkg.py'

# luminosity to normalize to (in 1/fb)                                                                                                                                                                      
lumi = 35.867

# used by mkPlot to define output directory for plots                                                                                                                                                       
# different from "outputDir" to do things more tidy                                                                                                                                                         
outputDirPlots = 'plotWW_2016_Sig_Bkg'

# used by mkDatacards to define output directory for datacards                                                                                                                                              
outputDirDatacard = 'datacards'

# structure file for datacard                                                                                                                                                                               
structureFile = 'structure.py'

# nuisances file for mkDatacards and for mkShape                                                                                                                                                            
nuisancesFile = 'nuisances.py'

