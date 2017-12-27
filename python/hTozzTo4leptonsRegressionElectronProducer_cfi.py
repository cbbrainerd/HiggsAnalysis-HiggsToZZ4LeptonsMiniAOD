import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsRegressionElectronProducer = cms.EDProducer("RegressionElectronProducer",
    eleCollection = cms.InputTag("gsfElectrons"),
    # electron regression
    eleRegressionEnergyErrorLabel  = cms.InputTag("eleRegressionEnergy:eneErrorRegForGsfEle"),
    eleRegressionEnergyLabel       = cms.InputTag("eleRegressionEnergy:eneRegForGsfEle")
)


