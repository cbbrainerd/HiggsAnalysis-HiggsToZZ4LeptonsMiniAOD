import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsCP = cms.EDProducer("HZZ4LeptonsCP",
    decayChain = cms.string('hToZZTo4LeptonsCP'),
    RECOcollName= cms.InputTag("hTozzTo4leptonsHiggsFrame"),
    debug = cms.untracked.bool(False)
)


