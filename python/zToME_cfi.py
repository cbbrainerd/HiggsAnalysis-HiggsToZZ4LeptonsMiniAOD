import FWCore.ParameterSet.Config as cms

zToME = cms.EDProducer("CandViewShallowCloneCombiner",
	decay = cms.string("hTozzTo4leptonsMuonSelector@- hTozzTo4leptonsElectronSelector@+"),
	cut = cms.string('0.0 < mass < 20000.0')

)


