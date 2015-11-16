import FWCore.ParameterSet.Config as cms


####################################################
# Initial setup
# Loading config files
# Et cetera

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("RecoTauTag.RecoTau.RecoTauPiZeroProducer_cfi")

process.load("BoostedTauAnalysis.CleanJets.cleanjets_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'root://eoscms//eos/cms/store/user/ktos/' 
    )
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START52_V9::All')

#mylist = FileUtils.loadListFromFile ('/afs/cern.ch/user/f/friccita/CMSSW_5_2_5/src/BoostedTauAnalysis/NSJAnalyzer/WToMuNu_skimfiles.txt')
#readFiles = cms.untracked.vstring( *mylist)
#process.source = cms.Source("PoolSource",
#                            fileNames = readFiles
##                           skipEvents = cms.untracked.uint32(6763)
#                           )


####################################################
# Modify the producers in RecoPFTauTag_cff.py
# to run on ak5PFJetsNoMu objects

process.recoTauAK5PFJets08Region.src = cms.InputTag("CleanJets", "ak4PFJetsNoMu", "OWNPARTICLES")
process.ak4PFJetsRecoTauPiZeros.jetSrc = cms.InputTag("CleanJets", "ak4PFJetsNoMu", "OWNPARTICLES")
process.combinatoricRecoTaus.jetSrc = cms.InputTag("CleanJets", "ak4PFJetsNoMu", "OWNPARTICLES")
process.ak4PFJetTracksAssociatorAtVertex.jets = cms.InputTag("CleanJets", "ak4PFJetsNoMu", "OWNPARTICLES")
process.ak4PFJetsLegacyHPSPiZeros.jetSrc = cms.InputTag("CleanJets", "ak4PFJetsNoMu", "OWNPARTICLES")

####################################################
# Rebuild the paths
process.recoTauCommonSequence = cms.Sequence(
            process.CleanJets *
            process.ak4PFJetTracksAssociatorAtVertex *
            process.recoTauAK5PFJets08Region *
            process.recoTauPileUpVertices *
            process.pfRecoTauTagInfoProducer
)

process.PFTau = cms.Sequence(
            # Jet production
            process.recoTauCommonSequence *
            # Make classic HPS taus
            process.recoTauClassicHPSSequence
)


####################################################
# Produce the output

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:/data1/friccita/NMSSMHiggs_gg_alltaus_PFTausNoMu.root')
)
  
process.p = cms.Path(process.PFTau)

process.e = cms.EndPath(process.out)

#processDumpFile = open('PFTausNoMu.dump', 'w'
#)

#print >> processDumpFile, process.dumpPython()
