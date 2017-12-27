import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsHLTAnalysisFilter = cms.EDFilter("HZZ4LeptonsHLTAnalysisFilter",
    HLTInfoFired = cms.InputTag("hTozzTo4leptonsHLTInfo"),                                           
)
