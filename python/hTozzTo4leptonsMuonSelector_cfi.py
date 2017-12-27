import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsMuonSelector = cms.EDProducer("HZZ4LeptonsMuonSelector",
    isGlobalMuon            = cms.bool(True),
    isTrackerMuon           = cms.bool(False),
    muonCollection = cms.InputTag("slimmedMuons"),  
    muonPtMin      = cms.double(5.0),
    muonEtaMax = cms.double(2.4)
)


