import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsHLTAnalysis = cms.EDProducer("HZZ4LeptonsHLTAnalysis",
    ElectronCollectionLabel = cms.InputTag("gsfElectrons"),
    # HLTPaths = cms.vstring('HLT_IsoMu11', 'HLT_Mu15', 'HLT_DoubleMu3', 'HLT_IsoEle15_L1I', 'HLT_IsoEle18_L1R','HLT_DoubleIsoEle10_L1I', 'HLT_DoubleIsoEle12_L1R', 'HLT_MinBiasBSC'),
    HLTPaths = cms.vstring('HLT_L2Mu0','HLT_L2Mu3','HLT_L1Mu20','HLT_L2Mu9','HLT_L2Mu11','HLT_L1Mu14_L1SingleEG10','HLT_L1Mu14_L1SingleJet6U','HLT_L1Mu14_L1ETM30','HLT_L2DoubleMu0','HLT_L1DoubleMuOpen','HLT_DoubleMu0','HLT_DoubleMu3','HLT_Mu3','HLT_Mu5','HLT_Mu9','HLT_IsoMu3','HLT_Mu0_L1MuOpen','HLT_Mu0_Track0_Jpsi','HLT_Mu3_L1MuOpen','HLT_Mu3_Track0_Jpsi','HLT_Mu5_L1MuOpen','HLT_Mu5_Track0_Jpsi','HLT_Mu0_L2Mu0','HLT_Mu3_L2Mu0','HLT_Mu5_L2Mu0','HLT_Photon10_L1R','HLT_Photon15_L1R','HLT_Photon15_LooseEcalIso_L1R','HLT_Photon20_L1R','HLT_Photon30_L1R_8E29','HLT_DoublePhoton4_Jpsi_L1R','HLT_DoublePhoton4_Upsilon_L1R','HLT_DoublePhoton4_eeRes_L1R','HLT_DoublePhoton5_eeRes_L1R','HLT_DoublePhoton5_Jpsi_L1R','HLT_DoublePhoton5_Upsilon_L1R','HLT_DoublePhoton5_L1R','HLT_DoublePhoton10_L1R','HLT_DoubleEle5_SW_L1R','HLT_Ele20_LW_L1R','HLT_Ele15_SiStrip_L1R','HLT_Ele15_SC10_LW_L1R','HLT_Ele15_LW_L1R','HLT_Ele10_LW_EleId_L1R','HLT_Ele10_LW_L1R','HLT_Photon15_TrackIso_L1R'),
    andOr = cms.bool(True),
    MuonCollectionLabel = cms.InputTag("muons"),
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
)


