import FWCore.ParameterSet.Config as cms

CleanJets = cms.EDProducer(
    'CleanJets',
    jetSrc = cms.InputTag("ak4PFJets"),
    muonSrc = cms.InputTag("muons"),
    PFCandSrc = cms.InputTag("pfIsolatedMuonsEI"),
    outFileName = cms.string('file:CleanJets_TEST.root'),
    momPDGID = cms.uint32(36),
    genMuTauPTMin = cms.double(0.0), #GeV
    genMuPTMin = cms.double(20.0), #GeV                           
    cutOnGenMatches = cms.bool(True),
    thisTag = cms.InputTag('CleanJets', 'ak4PFJetsNoMu', 'ANALYSIS')
    )
