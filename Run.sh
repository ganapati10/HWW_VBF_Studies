#!/bin/bash
#SBATCH --mail-user=sergio.blancof@alumnos.unican.es  # email address                                                                                                                                       
#SBATCH --mail-type=ALL                    # Alerts sent when job begins, ends, or aborts  

# Additional information can be found in the Altamira Users Guide
# https://confluence.ifca.es/display/IC/Altamira+Users+Guide#AltamiraUsersGuide-Batchsystem


#These scripts are going to be runned:
mkShapesMulti.py --pycfg=configuration.py --inputDir=/gpfs/projects/cms/data/LatinosSkims/nanoAOD/ --treeName=Events

mkShapesMulti.py --doHadd=1 --doNotCleanup --nThreads=8

mkPlot.py --pycfg=configuration.py --inputFile=rootFile/plots_WW_2016.root --minLogC=0.01 --minLogCratio=0.01 --maxLogC=1000 --maxLogCratio=1000 --showIntegralLegend=1






