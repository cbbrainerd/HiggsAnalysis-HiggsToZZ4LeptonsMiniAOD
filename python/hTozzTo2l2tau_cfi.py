import FWCore.ParameterSet.Config as cms

hTozzTo2l2tau = cms.EDProducer("CandViewShallowCloneCombiner",
                               decay = cms.string('zToEE zToMuMu'),
                               cut = cms.string('0.0 < mass < 20000.0'),
                               checkCharge = cms.bool(False)
                               )

