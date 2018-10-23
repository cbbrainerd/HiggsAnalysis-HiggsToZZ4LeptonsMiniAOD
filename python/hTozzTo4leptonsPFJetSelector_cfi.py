import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsPFJetSelector = cms.EDProducer("HZZ4LeptonsPFJetSelector",
    #isLoosePFJetID            = cms.bool(True),
    isTightPFJetID            = cms.bool(True),
    #PFJetCollection = cms.InputTag("slimmedJets")
    PFJetCollection = cms.InputTag("updatedPatJetsUpdatedJEC")
)


