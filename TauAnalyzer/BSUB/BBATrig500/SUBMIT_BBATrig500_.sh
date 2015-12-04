#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc491
cd /afs/cern.ch/user/k/ktos/BoostedDiTau/CMSSW_7_4_1_patch1/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp /afs/cern.ch/user/k/ktos/BoostedDiTau/CMSSW_7_4_1_patch1/src/BoostedTauAnalysis/TauAnalyzer/BSUB/BBATrig500/BBA_trig_cfg_500_BBATrig500_.py . 
cmsRun BBA_trig_cfg_500_BBATrig500_.py
cmsStage -f BBATrig500.root /store/user/ktos/BBATrig500
rm BBA_trig_cfg_500_BBATrig500_.py 
exit 0
