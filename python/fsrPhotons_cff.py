import FWCore.ParameterSet.Config as cms

## PF Photons
fsrPhotons = cms.EDProducer("FSRPhotonProducer",
    srcCands = cms.InputTag("packedPFCandidates"),
    muons = cms.InputTag("slimmedMuons"), 
    ptThresh = cms.double(1.0), ## will tighten to 2 at analysis level
    extractMuonFSR = cms.bool(False),
)

import PhysicsTools.PatAlgos.producersLayer1.pfParticleProducer_cfi 
boostedFsrPhotons = PhysicsTools.PatAlgos.producersLayer1.pfParticleProducer_cfi.patPFParticles.clone(
    pfCandidateSource = 'fsrPhotons'
)

fsrPhotonSequence = cms.Sequence(
    fsrPhotons +
    boostedFsrPhotons
)
