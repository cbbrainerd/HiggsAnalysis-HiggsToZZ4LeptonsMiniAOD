import FWCore.ParameterSet.Config as cms

edmLumi = cms.EDFilter("EdmLumi",
   lumi_by_LS_all_csv = cms.untracked.string("HiggsAnalysis/HiggsToZZ4Leptons/data/lumi_by_LS_all.csv"),
   method = cms.string("vertex"),
   #prescale_by_run    = cms.untracked.string("HiggsAnalysis/HiggsToZZ4Leptons/data/HLT_Mu3.txt"),
   scale = cms.double(4.29e+28)
)
edmLumiHF = edmLumi.clone(method = "hf")
edmLumiOnline = edmLumi.clone(
   method = "online",
   lumi_by_LS_all_csv ="HiggsAnalysis/HiggsToZZ4Leptons/data/lumi_by_LS_all.csv",
   scale = cms.double(1)
)
#Lumi_Path = cms.Path(edmLumi + edmLumiHF)
Lumi_Path_Online = cms.Path(edmLumiOnline)
