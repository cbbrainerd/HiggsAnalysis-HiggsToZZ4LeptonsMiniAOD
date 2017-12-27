import FWCore.ParameterSet.Config as cms

ConvValueMapProd = cms.EDProducer('ConvValueMapProd',
    gsfLabel = cms.untracked.InputTag("hTozzTo4leptonsElectronSelector"),
    tkLabel = cms.untracked.InputTag("generalTracks")
)

