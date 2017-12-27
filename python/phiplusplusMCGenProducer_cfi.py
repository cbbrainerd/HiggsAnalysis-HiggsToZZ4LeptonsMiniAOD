import FWCore.ParameterSet.Config as cms

phiplusplusMCGenProducer = cms.EDProducer("PhiPlusPlusMCGenProducer",
    src = cms.InputTag("genParticles")
)


