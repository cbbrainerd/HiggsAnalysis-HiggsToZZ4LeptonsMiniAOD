import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsMuonRochesterCalibrator = cms.EDProducer("HZZ4LeptonsMuonRochesterCalibrator",
    muonCollection = cms.InputTag("slimmedMuons"),
    isData         = cms.bool(False),  
    goodMuonMCMatch      = cms.InputTag("goodMuonMCMatchCalib"),
    myMuons              = cms.InputTag("myMuonsCalib"),
    MCTruth    = cms.bool(True)                            
)


