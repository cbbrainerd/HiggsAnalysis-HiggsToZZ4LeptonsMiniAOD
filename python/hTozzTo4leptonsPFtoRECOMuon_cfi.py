import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsPFtoRECOMuon = cms.EDProducer("HZZ4LeptonsPFtoRECOMuon",
    pfCollection = cms.InputTag("packedPFCandidates"),                               
)


