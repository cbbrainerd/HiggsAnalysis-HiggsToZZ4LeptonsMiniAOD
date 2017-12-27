import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsMuonCalibrator = cms.EDProducer("HZZ4LeptonsMuonCalibrator",
    muonCollection = cms.InputTag("slimmedMuons"),
    isData         = cms.bool(False)                               
)


