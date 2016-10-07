# L1TauStudy
```
#!/bin/bash
scramv1 project CMSSW CMSSW_6_2_0_SLHC12_patch1
cd CMSSW_6_2_0_SLHC12_patch1/src
eval `scramv1 runtime -sh`

git cms-init
git cms-addpkg SLHCUpgradeSimulations/L1TrackTrigger
git cms-merge-topic EmanuelPerez:TTI_62X_TrackTriggerObjects

nohup scramv1 b -j 8
```
