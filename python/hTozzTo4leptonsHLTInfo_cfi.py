import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsHLTInfo = cms.EDProducer("HZZ4LeptonsHLTInfo",
  TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
  debug = cms.untracked.bool(False)                                
)

