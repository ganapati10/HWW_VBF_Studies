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

source LatinosSetup/SetupShapeOnly.sh

mv LatinoAnalysis/MultiDraw ./

git clone https://github.com/scodella/LatinoAnalysis

mv MultiDraw ./LatinoAnalysis/

git clone --branch ww-inicial https://github.com/calderona/PlotsConfigurations

cd LatinoAnalysis/Tools/python/

cp userConfig_TEMPLATE.py userConfig.py

mv NanoAODTools PhysicsTools

```
The userConfig.py must be changed in the next way:

```
basedir = '/afs/cern.ch/user/x/xjanssen/cms/HWW2015/'
```
To
```
basedir = '/gpfs/users/blancoser/CMSSW_10_6_4'
```

Then:

``` 
cd ../../..

scram b -j 4

```
