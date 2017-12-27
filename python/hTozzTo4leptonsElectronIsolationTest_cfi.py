import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsElectronIsolationTest = cms.EDProducer("HZZ4LeptonsElectronIsolationTest",
    electronCollection = cms.InputTag("gedGsfElectrons"),
    electronEtaMax     = cms.double(2.5),
    electronPtMin      = cms.double(5.0),
    ElectronPFIsoValueChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPFBRECO"),
    ElectronPFIsoValueCharged    = cms.InputTag("elPFIsoValueCharged03PFIdPFBRECO"),
    ElectronPFIsoValueNeutral    = cms.InputTag("elPFIsoValueNeutral03PFIdPFBRECO"),
    ElectronPFIsoValueGamma      = cms.InputTag("elPFIsoValueGamma03PFIdPFBRECO"),
    ElectronPFIsoValuePU         = cms.InputTag("elPFIsoValuePU03PFIdPFBRECO")
)


