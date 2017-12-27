import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsMCGenParticleListDrawer = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint  = cms.untracked.int32(1)
)

