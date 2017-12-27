import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsBestCandidateProducer = cms.EDProducer("HZZ4LeptonsBestCandidate",
    decaychannel = cms.string('2e2mu'),
    decayChain = cms.string('hToZZTo4LeptonsBestCandidate'),
    RECOcollName = cms.VInputTag(cms.InputTag("hTozzTo4leptonsLooseIsol"))
)


