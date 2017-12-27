import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsHiggsFrame = cms.EDProducer("HZZ4LeptonsHiggsFrame",
    decayChain = cms.string('hToZZTo4LeptonsHiggsRestFrame'),
    prodinst = cms.string('hTozzTo4leptonsBestCandidateProducer'),
    RECOcollName = cms.VInputTag(cms.InputTag("hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hToZZTo4LeptonsBestCandidateBoson1"))
)


