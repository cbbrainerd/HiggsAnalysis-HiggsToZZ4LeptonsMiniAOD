import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsMuonIsolationProducerMu = cms.EDProducer("HZZ4LeptonsMuonIsolationProducerMu",
                                 MuonsLabel     = cms.InputTag("hTozzTo4leptonsMuonSelector"),
                                 muTkIso = cms.InputTag("muIsoFromDepsTkOptimized"),
                                 muEcalIso = cms.InputTag("muIsoFromDepsEcalOptimized"),
                                 muHcalIso = cms.InputTag("muIsoFromDepsHcalOptimized"),
                                 threshold = cms.double(0.7), 
                                 coeffTk = cms.double(1), 
                                 coeffEcal = cms.double(1), 
                                 coeffHcal = cms.double(1),
                                 useRelativeIso = cms.bool(True)
                                 )

