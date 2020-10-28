# Master

## Analysis of VBF Higgs 

First of all, let's start building the programming enviroment. The next part explain the instructions to execute the code. 

First, log in gridui:

```ssh -Y ***@gridui.ifca.es -o ServerAliveInterval=240```


Set CMS enviroment for the first time:

```
bash -l

source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH=slc7_amd64_gcc530

cmsrel CMSSW_10_6_4_patch1
```


Enter CMS eviroment on gridui:

```
bash -l

source /cvmfs/cms.cern.ch/cmsset_default.sh

cd CMSSW_10_6_4_patch1/src

cmsenv
```

Compile code:

```scram b -j 8```


Get some code from gitHub:

```git clone https://github.com/BlancoFS/...```


## Run code:

submit job

```
mkShapesMulti.py --pycfg=configuration.py --inputDir=/gpfs/projects/cms/data/LatinosSkims/nanoAOD/ --treeName=Events --batchSplit=Samples,Files --doBatch=True
```

On slurm, to send a job Run.sh to the queue:

```
sbatch -o logfile.log -e errofile.err --qos=gridui_sort --partition=cloudcms run.sh
```


