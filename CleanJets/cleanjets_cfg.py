import FWCore.ParameterSet.Config as cms

process = cms.Process("CLEANJETS")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load('BoostedTauAnalysis/CleanJets/RecoPFTauTag_cff_TEST')

###############################
# Fixing several problems
############################### 
#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
process.load('RecoBTag.Configuration.RecoBTag_cff')
process.load('RecoJets.Configuration.RecoJetAssociations_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.GlobalTag.globaltag = cms.string('74X_dataRun2_Prompt_v0')

###########
#  Input
###########
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
	 fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/user/mshi/gg2H2aa2mumutautau_STEP_2_9GeV/gg2H2aa2mumutautau_STEP_2_9GeV_977.root') 
)

#########################################################
# this will produce a ref to the original muon collection
#########################################################
process.muonsRef = cms.EDFilter('MuonRefSelector',
	src = cms.InputTag('muons'),
        cut = cms.string('pt > .001'),
        filter = cms.bool(True)
)

##############################
# Clean Jets Definition
##############################
process.CleanJets = cms.EDProducer(
    'CleanJets',
    jetSrc = cms.InputTag("ak4PFJets"),
    muonSrc = cms.InputTag("muonsRef"),
    PFCandSrc = cms.InputTag("pfIsolatedMuonsEI"),
    outFileName = cms.string('file:OUTPUT/CleanJets_Plots.root'),
    momPDGID = cms.uint32(36),
    genMuTauPTMin = cms.double(0.0), #GeV
    genMuPTMin = cms.double(20.0), #GeV                           
    cutOnGenMatches = cms.bool(True),
    thisTag = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')
    )

#######################################
# HPS Tau Reconstruction alterations 
#######################################
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.ak4PFJetTracksAssociatorAtVertex.jets = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')
process.recoTauAK4PFJets08Region.src = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')

process.ak4PFJetsLegacyHPSPiZeros.jetSrc = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')
process.ak4PFJetsRecoTauChargedHadrons.jetSrc = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')
process.combinatoricRecoTaus.jetSrc = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'CLEANJETS')


#####################
# output Definition
####################
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string("file:/afs/cern.ch/user/k/ktos/BoostedDiTau/CMSSW_7_4_1_patch1/src/BoostedTauAnalysis/CleanJets/OUTPUT/CleanJets_edmOuput.root") )
process.p = cms.Path(process.muonsRef*process.CleanJets*process.PFTau)
#process.p = cms.Path(process.PFTau)
process.outpath = cms.EndPath(process.out)
