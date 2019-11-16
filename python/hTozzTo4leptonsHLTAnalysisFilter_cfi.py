import FWCore.ParameterSet.Config as cms

from triggerMenus_cfi import get_triggers

hTozzTo4leptonsHLTAnalysisFilter = cms.EDFilter("HZZ4LeptonsHLTAnalysisFilter",
    HLTInfoFired = cms.InputTag("hTozzTo4leptonsHLTInfo"),                                           
)
