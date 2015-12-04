#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc491
cd /afs/cern.ch/user/k/ktos/BoostedDiTau/CMSSW_7_4_1_patch1/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp /afs/cern.ch/user/k/ktos/BoostedDiTau/CMSSW_7_4_1_patch1/src/BoostedTauAnalysis/TauAnalyzer/BSUB/DIRNAME/PRODUCER.py . 
cmsRun PRODUCER.py
cmsStage -f DIRNAME.root /store/user/ktos/DIRNAME
rm PRODUCER.py 
exit 0
