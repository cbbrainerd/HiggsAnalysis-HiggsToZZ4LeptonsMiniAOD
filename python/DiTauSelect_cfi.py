import FWCore.ParameterSet.Config as cms
ditauvispairs = cms.EDProducer("DiTauSelect",
                     MuonTag = cms.untracked.InputTag("hTozzTo4leptonsMuonSelector"),                
                     ElectronTag = cms.untracked.InputTag("hTozzTo4leptonsElectronSelector"),
                     pfMetTag = cms.untracked.InputTag("pfMet")
                     )
