import FWCore.ParameterSet.Config as cms

zToMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('hTozzTo4leptonsMuonSelector@+ hTozzTo4leptonsMuonSelector@-'),
    cut = cms.string('0.0 < mass < 20000.0'),
    moduleLabel = cms.untracked.string('zToMuMu')
)

