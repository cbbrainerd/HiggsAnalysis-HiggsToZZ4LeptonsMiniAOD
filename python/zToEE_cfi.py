import FWCore.ParameterSet.Config as cms

zToEE = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('hTozzTo4leptonsElectronSelector@+ hTozzTo4leptonsElectronSelector@-'),
    cut = cms.string('0.0 < mass < 20000.0'),
    moduleLabel = cms.untracked.string('zToEE')
)

