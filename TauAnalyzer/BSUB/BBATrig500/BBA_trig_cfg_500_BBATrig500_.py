import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile('/afs/cern.ch/user/k/ktos/BoostedDiTau/CMSSW_7_4_1_patch1/src/BoostedTauAnalysis/TauAnalyzer/FileLists/inFileList_BBA_a30_500.txt')

process = cms.Process("RecoTrigProducer")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_71_V1::All', '')

####################
# Message Logger
####################
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
## switch to uncheduled mode
#process.options.allowUnscheduled = cms.untracked.bool(True)

####################
# Input File List
####################
#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/user/ktos/RECO_Step2_BBA_a30/RECO_Step2_BBA_a30_.root'),
#    skipEvents = cms.untracked.uint32(0)
#    )
readFiles = cms.untracked.vstring(*mylist)
process.source = cms.Source("PoolSource",
    fileNames = readFiles,
    skipEvents = cms.untracked.uint32(0)
    )


#this will produce a ref to the original muon collection
process.muonsRef = cms.EDFilter('MuonRefSelector',
                                         src = cms.InputTag('muons'),
                                         cut = cms.string('pt > .001'),
                                         filter = cms.bool(True)
                                         )
##this will produce a ref to the original muon collection
#process.jetRef = cms.EDFilter('JetRefSelector',
#                                         src = cms.InputTag('ak4PFJetsRefs'),
#                                         cut = cms.string('pt > .001'),
#                                         filter = cms.bool(True)
#                                         )

#muon trigger object filter
process.bbaTriggerFilter = cms.EDFilter(
    'BBATriggerObjectFilter',
    mRecoObjTag = cms.InputTag("muonsRef"),
    bRecoObjTag = cms.InputTag("ak4PFJetsRefs"),
    outFileName = cms.string("/afs/cern.ch/user/k/ktos/BoostedDiTau/CMSSW_7_4_1_patch1/src/BoostedTauAnalysis/TauAnalyzer/BSUB/BBATrig500/BBATrig500.root"),
    genParticleTag = cms.InputTag("genParticles"),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    triggerDelRMatch = cms.untracked.double(0.1),
    hltTags = cms.VInputTag(cms.InputTag("HLT_Mu10_CentralPFJet30_BTagCSV0p5PF_v1", "", "HLT"),
                            cms.InputTag("HLT_Mu10_CentralPFJet30_BTagCSV0p5PF_v2", "", "HLT")
                            ),
    theRightHLTTag         = cms.InputTag("HLT_Mu10_CentralPFJet30_BTagCSV0p5PF_v1"),
    mSubFilter = cms.InputTag("hltL3fL1sMu0L1f0L2f3QL3Filtered10Q"),#hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3cr"),
    bSubFilter        = cms.InputTag("hltPFJetForBtag"),
    minNumObjsToPassFilter = cms.uint32(2)
    )

# talk to output module
process.out = cms.OutputModule("PoolOutputModule",
	fileName = cms.untracked.string("file:BBATrig500.root") 
)

process.p = cms.Path(process.muonsRef
#	*process.jetRef
	*process.bbaTriggerFilter )
process.outpath = cms.EndPath(process.out)

