import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsElectronIsolationProducerEgamma = cms.EDProducer("HZZ4LeptonsElectronIsolationProducerEgamma",
                                 ElectronsLabel     = cms.InputTag("hTozzTo4leptonsElectronSelector"),
                                 eleTkIso = cms.InputTag("eleIsoFromDepsTkOptimized"),
                                 eleEcalIso = cms.InputTag("eleIsoFromDepsEcalFromHitsByCrystalOptimized"),
                                 eleHcalIso = cms.InputTag("eleIsoFromDepsHcalFromTowersOptimized"),
                                 threshold = cms.double(0.7), 
                                 coeffTk = cms.double(1), 
                                 coeffEcal = cms.double(1), 
                                 coeffHcal = cms.double(1),
                                 useRelativeIso = cms.bool(True)
                                 )

