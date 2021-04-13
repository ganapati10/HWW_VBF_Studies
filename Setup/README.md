## Setup changes:

All the links to Github has to be changed in order to replace the:

```git@github.com:*```

 for full links 
 
 ```https://github.com/*```


## Setup process:

```
cmsrel CMSSW_10_6_4

cd CMSSW_10_6_4/src/

source /cvmfs/cms.cern.ch/cmsset_default.sh

cmsenv

git clone --branch 13TeV https://github.com/latinos/setup LatinosSetup

rm -rf LatinoAnalysis
```

Now, open LatinosSetup/SetupShapeOnly.sh and change a few links. The comment lines has been changed.
The new file is the one in this folder.

```
source LatinosSetup/SetupShapeOnly.sh

mv LatinoAnalysis/MultiDraw ./

git clone https://github.com/scodella/LatinoAnalysis

mv MultiDraw ./LatinoAnalysis/

git clone --branch ww-inicial https://github.com/calderona/PlotsConfigurations

cd LatinoAnalysis/Tools/python/

cp userConfig_TEMPLATE.py userConfig.py

cd ../../..

mv NanoAODTools PhysicsTools/

```
The userConfig.py (From src folder: cd LatinoAnalysis/Tools/python/) must be changed in the next way:

```
basedir = '/afs/cern.ch/user/x/xjanssen/cms/HWW2015/'
```
To
```
basedir = '/gpfs/users/blancoser/CMSSW_10_6_4/src/PlotsConfigurations/Configurations/WW/Full2026_v6/'
```

Then:

``` 
cd ../../..

scram b -j 4

```

## Setup Combine framework

The next lines are the instructions to be followed for installing the Combine framework. Please, start from your src path:

```
mkdir combine
cd combine

cmsrel CMSSW_10_2_13

cd CMSSW_10_2_13/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester

scramv1 b clean
scramv1 b -j 8





