#!/bin/sh                                                                                                                              
#NeutrinoGun.txt
voms-proxy-init --voms cms --valid 100:00                                                                                              

cat test-runTaus.py > SUBPhase2.py
cat submit.py >> SUBPhase2.py

rm -r /nfs_scratch/ojalvo/NeutrinoGun-$1-SUBPhase2/
mkdir /nfs_scratch/ojalvo/NeutrinoGun-$1-SUBPhase2/
##make dag dir
mkdir -p /nfs_scratch/ojalvo/NeutrinoGun-$1-SUBPhase2/dags
mkdir -p /nfs_scratch/ojalvo/NeutrinoGun-$1-SUBPhase2/dags/daginputs

## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/$1/server?SFN=/hdfs/store/user/ojalvo/NeutrinoGun-$1-SUBPhase2/
#Matching E and H activity
#farmoutAnalysisJobs  --input-file-list=RAW-Run2015D.txt --no-shared-fs  --submit-dir=/nfs_scratch/ojalvo/NeutrinoGun-$1-SUBPhase2/submit --output-dag-file=/nfs_scratch/ojalvo/NeutrinoGun-$1-SUBPhase2/dags/dag --output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/$1/server?SFN=/hdfs/store/user/ojalvo/NeutrinoGun-$1-SUBPhase2/ NeutrinoGun-$1   $CMSSW_BASE  $CMSSW_BASE/src/L1Trigger/L1TNtuplizer/test/SUBPhase2.py 

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=NeutrinoGun.txt \
--submit-dir=/nfs_scratch/ojalvo/NeutrinoGun-$1-SUBPhase2/submit \
--output-dag-file=/nfs_scratch/ojalvo/NeutrinoGun-$1-SUBPhase2/dags/dag \
NeutrinoGun-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/SLHCUpgradeSimulations/L1TauStudy/test/SUBPhase2.py $2
#SLHCUpgradeSimulations/L1TauStudy/test
rm SUBPhase2.py

