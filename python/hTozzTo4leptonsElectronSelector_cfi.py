import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsElectronSelector = cms.EDProducer("HZZ4LeptonsElectronSelector",
    electronCollection = cms.InputTag("slimmedElectrons"),
    electronEtaMax     = cms.double(2.5),
    electronPtMin      = cms.double(5.0)
)


