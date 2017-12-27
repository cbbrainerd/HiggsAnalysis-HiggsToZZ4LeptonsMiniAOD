import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsConstraintFitProducer = cms.EDProducer("HZZ4LeptonsConstraintFitProducer",
      VertexLabel  = cms.InputTag("offlineSlimmedPrimaryVertices"),
      RECOcollName = cms.InputTag("hTozzTo4leptonsBestCandidateProducer:hToZZTo4LeptonsBestCandidateMother"),
      nParticles   = cms.uint32(4)
)
