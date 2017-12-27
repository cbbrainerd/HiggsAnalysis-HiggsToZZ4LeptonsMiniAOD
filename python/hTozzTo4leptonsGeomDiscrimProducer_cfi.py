import FWCore.ParameterSet.Config as cms


hTozzTo4leptonsGeomDiscrimProducer = cms.EDProducer("HZZ4LeptonsGeomDiscrimProducer",
    decaychannel   = cms.string('2e2mu'),
    BeamSpotLabel  = cms.InputTag("offlineBeamSpot"),
    RECOcollName   = cms.InputTag("hTozzTo4leptonsBestCandidateProducer:hToZZTo4LeptonsBestCandidateMother")
)


