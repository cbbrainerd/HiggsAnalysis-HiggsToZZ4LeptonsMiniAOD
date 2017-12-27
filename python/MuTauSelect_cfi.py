import FWCore.ParameterSet.Config as cms
mutauvispairs = cms.EDProducer("MuTauSelect",
                     MuonTag = cms.untracked.InputTag("hTozzTo4leptonsMuonSelector"),                
                     ElectronTag = cms.untracked.InputTag("hTozzTo4leptonsElectronSelector"),
                     #TrigTag = cms.untracked.InputTag("TriggerResults","","REDIGI311X"),
                     #MuonTrig = cms.untracked.string("HLT_Mu9"),
                     pfMetTag = cms.untracked.InputTag("pfMet")
                     )
