import FWCore.ParameterSet.Config as cms

zToCrossLeptons = cms.EDProducer("CandViewShallowCloneCombiner",
	decay       = cms.string("hTozzTo4leptonsMuonSelector hTozzTo4leptonsElectronSelector"),
        checkCharge = cms.bool(False),
	cut         = cms.string('mass > 0')

)


