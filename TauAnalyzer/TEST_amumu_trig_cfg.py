import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile('/afs/cern.ch/user/k/ktos/BoostedDiTau/CMSSW_7_4_1_patch1/src/BoostedTauAnalysis/TauAnalyzer/FileLists/inFileList_amumu_a5.txt')

process = cms.Process("RecoTrigProducer")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_71_V1::All', '')
process.load("Configuration.StandardSequences.MagneticField_cff")

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
readFiles = cms.untracked.vstring(*mylist)
process.source = cms.Source("PoolSource",
    fileNames = readFiles,
    skipEvents = cms.untracked.uint32(0)
    )


#muon trigger object filter
process.muonTriggerObjectFilter = cms.EDFilter(
    'MuonTriggerObjectFilter',
    mRecoObjTag = cms.InputTag("muons"),
    bRecoObjTag = cms.InputTag("ak4PFJets"),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    triggerDelRMatch = cms.untracked.double(0.1),
    hltTags = cms.VInputTag(cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v1", "", "HLT"),
                            cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v2", "", "HLT"),
                            cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v3", "", "HLT"),
                            cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v4", "", "HLT"),
                            cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v5", "", "HLT"),
                            cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v6", "", "HLT"),
                            cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v7", "", "HLT"),
                            cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v8", "", "HLT"),
                            cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v9", "", "HLT"),
                            cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v10", "", "HLT"),
                            cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v11", "", "HLT"),
                            cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v12", "", "HLT"),
                            cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v13", "", "HLT"),
                            cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v14", "", "HLT"),
                            cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia_v15", "", "HLT")
                            ),
    theRightHLTTag         = cms.InputTag("HLT_Mu16_TkMu0_dEta18_Onia"),
    globalTrkMuonSubFilter = cms.InputTag("hltGlbTrkMuons"),#hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3cr"),
    diMuonSubFilter        = cms.InputTag("hltDiMuonGlbFiltered16TrkFiltered0"),
    mu16SubFilter	   = cms.InputTag("hltL3fL1sMu14erorMu16L1f0L2f0L3Filtered16"),
    minNumObjsToPassFilter = cms.uint32(2)
    )

# talk to output module
process.out = cms.OutputModule("PoolOutputModule",
	fileName = cms.untracked.string("file:RecoTrigProd_TEST.root") 
)

process.p = cms.Path(process.muonTriggerObjectFilter )
process.outpath = cms.EndPath(process.out)

