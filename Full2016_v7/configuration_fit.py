# VBF Fit configuration file                                                                                                                                                                               

treeName = 'Events'

tag = 'WW_2016'

# used by mkShape to define output directory for root files                                                                                                                                                
outputDir = 'rootFile'

# file with TTree aliases                                                                                                                                                                                  
aliasesFile = 'aliases.py'

# file with list of variables                                                                                                                                                                              
variablesFile = 'variables_fit.py'

# file with list of cuts                                                                                                                                                                                   
cutsFile = 'cuts.py'

# file with list of samples                                                                                                                                                                                
samplesFile = 'samples.py'

# file with plot configuration                                                                                                                                                                             
plotFile = 'plot_fit.py'

# luminosity to normalize to (in 1/fb)                                                                                                                                                                     
lumi = 35.867

# used by mkPlot to define output directory for plots                                                                                                                                                      
# different from "outputDir" to do things more tidy                                                                                                                                                        
outputDirPlots = 'plotWW_2016'

# used by mkDatacards to define output directory for datacards                                                                                                                                             
outputDirDatacard = 'datacards'

# structure file for datacard                                                                                                                                                                              
structureFile = 'structure.py'

# nuisances file for mkDatacards and for mkShape                                                                                                                                                           
nuisancesFile = 'nuisances_ALL.py'

