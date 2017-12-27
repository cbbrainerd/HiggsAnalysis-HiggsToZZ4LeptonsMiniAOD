import FWCore.ParameterSet.Config as cms


## 
# Filter to select 2l2tau events
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi
hTozzTo4leptonsMCGenFilter2l2tau = HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi.hTozzTo4leptonsMCGenFilter.clone()
hTozzTo4leptonsMCGenFilter2l2tau.HZZ4LeptonsMCFilterLeptonFlavour=cms.int32(7)#(7 to 12)

# ParticleListDrawer
from  HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenParticleListDrawer_cfi import * 
hTozzTo4leptonsMCGenParticleListDrawer2l2tau = HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenParticleListDrawer_cfi.hTozzTo4leptonsMCGenParticleListDrawer.clone()

# Save MC truth: 2l2tau
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCDumper_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCDumper_cfi
hTozzTo4leptonsMCDumper2l2tau = HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCDumper_cfi.hTozzTo4leptonsMCDumper.clone()
#hTozzTo4leptonsMCDumper2l2tau.motherPdgId          = cms.int32(25)
#hTozzTo4leptonsMCDumper2l2tau.firstdaughtersPdgId  = cms.vint32( 23, 23)
#hTozzTo4leptonsMCDumper2l2tau.seconddaughtersPdgId = cms.vint32( 11, -11, 13, -13 )

# Build best candidate from MC truth: 2l2tau
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsBestCandidateProducer_cfi import *
hTozzTo4leptonsMCBestCandidateProducer=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsBestCandidateProducer_cfi.hTozzTo4leptonsBestCandidateProducer.clone()
hTozzTo4leptonsMCBestCandidateProducer.decaychannel = cms.string('2l2tau')
hTozzTo4leptonsMCBestCandidateProducer.RECOcollName = cms.VInputTag(cms.InputTag("hTozzTo4leptonsMCDumper2l2tau:hToZZTo4LeptonsMCDumperCompositeMother"))
hTozzTo4leptonsMCBestCandidateProducer.decayChain = cms.string('hToZZTo4LeptonsMCBestCandidate')

# CP producer: 2l2tau
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi import *
hTozzTo4leptonsMCCP2l2tau=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi.hTozzTo4leptonsCP.clone()
hTozzTo4leptonsMCCP2l2tau.decayChain = cms.string('hToZZTo4LeptonsMCBestCandidateCP')
hTozzTo4leptonsMCCP2l2tau.RECOcollName = cms.VInputTag(cms.InputTag("hTozzTo4leptonsMCBestCandidateProducer:hToZZTo4LeptonsMCBestCandidateMother"), cms.InputTag("hTozzTo4leptonsMCBestCandidateProducer:hToZZTo4LeptonsMCBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsMCBestCandidateProducer:hToZZTo4LeptonsMCBestCandidateBoson1"))

# Filter of low mass peak events
# from HiggsAnalysis.HiggsToZZ4Leptons.hTozzT4LeptonsLowMassPeakFilter_cfi import *

# module to apply skim
from HiggsAnalysis.Skimming.higgsToZZ4Leptons_Sequences_cff import *
higgsToZZ4LeptonsHLTAnalysisData=HiggsAnalysis.Skimming.higgsToZZ4LeptonsHLTAnalysis_cfi.higgsToZZ4LeptonsHLTAnalysis.clone()


## Data - Double Mu Double Electron
higgsToZZ4LeptonsHLTAnalysisData.HLTPaths= cms.vstring('HLT_DoubleMu7_v1','HLT_L1DoubleMu0_v1','HLT_DoubleMu4_Acoplanarity03_v1','HLT_L2DoubleMu23_NoVertex_v1','HLT_L2DoubleMu0_v2','HLT_DoubleMu3_v3','HLT_DoubleMu6_v1','HLT_Mu8_Jet40_v3','HLT_TripleMu5_v2','HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2','HLT_Ele17_CaloIdL_CaloIsoVL_v2','HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2','HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2','HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2','HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2','HLT_TripleEle10_CaloIdL_TrkIdVL_v2','HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v2','HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2','HLT_Ele8_CaloIdL_CaloIsoVL_v2','HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2','HLT_Ele8_CaloIdL_TrkIdVL_v2','HLT_Ele8_v2');

## MC samples 3_11_2
higgsToZZ4LeptonsHLTAnalysisData.HLTPaths= cms.vstring('HLT_L1DoubleMuOpen','HLT_L2DoubleMu0','HLT_L2DoubleMu20_NoVertex_v1','HLT_DoubleMu0','HLT_DoubleMu0_Quarkonium_v1','HLT_DoubleMu3_v2','HLT_DoubleMu5_v1','HLT_Ele10_SW_L1R_v2','HLT_Ele12_SW_TighterEleId_L1R_v2','HLT_Ele17_SW_L1R_v2','HLT_Ele17_SW_Isol_L1R_v2','HLT_Ele17_SW_TighterEleIdIsol_L1R_v3','HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2','HLT_Ele22_SW_L1R_v2','HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2','HLT_Ele22_SW_TighterEleId_L1R_v3','HLT_Ele32_SW_TighterEleId_L1R_v2','HLT_Photon20_Isol_Cleaned_L1R_v1','HLT_DoubleEle17_SW_L1R_v1');
                                                       



#higgsToZZ4LeptonsHLTAnalysisData.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
higgsToZZ4LeptonsHLTAnalysisData.TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI311X")
#higgsToZZ4LeptonsHLTAnalysisData.TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI41X")

higgsToZZ4LeptonsSkimTriLeptonProducerData=HiggsAnalysis.Skimming.higgsToZZ4LeptonsSkimTriLeptonProducer_cfi.higgsToZZ4LeptonsSkimTriLeptonProducer.clone()
higgsToZZ4LeptonsSkimTriLeptonProducerData.nLeptonMinimum=cms.int32(3)
higgsToZZ4LeptonsSkimTriLeptonProducerData.softMinimumPt= cms.double(5.0)
higgsToZZ4LeptonsSkimTriLeptonProducerData.stiffMinimumPt= cms.double(10.0)
higgsToZZ4LeptonsSkimFilterData=HiggsAnalysis.Skimming.higgsToZZ4LeptonsSkimFilter_cfi.higgsToZZ4LeptonsSkimFilter.clone()
higgsToZZ4LeptonsSkimFilterData.useHLT  = cms.untracked.bool(False) 
higgsToZZ4LeptonsSkimFilterData.HLTinst = cms.string('higgsToZZ4LeptonsHLTAnalysisData')
higgsToZZ4LeptonsSkimFilterData.useDiLeptonSkim = cms.untracked.bool(True)
higgsToZZ4LeptonsSkimFilterData.SkimDiLeptonflag = cms.string('SkimDiLeptonB')

higgsToZZ4LeptonsSkimFilterData.SkimTriLeptoninst  = cms.string('higgsToZZ4LeptonsSkimTriLeptonProducerData')

higgsToZZ4LeptonsSequenceData = cms.Sequence(higgsToZZ4LeptonsHLTAnalysisData           +
                                             higgsToZZ4LeptonsBuildLeptons              +
                                             higgsToZZ4LeptonsSkimDiLeptonProducer      +
                                             higgsToZZ4LeptonsSkimTriLeptonProducerData +
                                             higgsToZZ4LeptonsSkimFilterData )


# module to run HLT analysis
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHLTAnalysis_cfi import *
hTozzTo4leptonsHLTAnalysisData=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHLTAnalysis_cfi.hTozzTo4leptonsHLTAnalysis.clone()
#hTozzTo4leptonsHLTAnalysisData.HLTPaths= cms.vstring('HLT_L1MuOpen', 'HLT_L1MuOpen_DT','HLT_L1Mu', 'HLT_L1Mu20', 'HLT_L2Mu0_NoVertex','HLT_L2Mu0', 'HLT_L2Mu3', 'HLT_L2Mu9', 'HLT_L2Mu25', 'HLT_Mu3', 'HLT_Mu5', 'HLT_Mu7', 'HLT_Mu9', 'HLT_Mu11','HLT_IsoMu9','HLT_Mu20_NoVertex','HLT_L1DoubleMuOpen', 'HLT_L2DoubleMu0', 'HLT_DoubleMu0','HLT_DoubleMu3','HLT_Mu0_L1MuOpen', 'HLT_Mu3_L1MuOpen', 'HLT_Mu5_L1MuOpen', 'HLT_Mu0_L2Mu0', 'HLT_Mu5_L2Mu0', 'HLT_Mu0_Track0_Jpsi', 'HLT_Mu0_TkMu0_OST_Jpsi', 'HLT_Mu3_Track3_Jpsi','HLT_Mu3_TkMu0_OST_Jpsi','HLT_Mu5_Track0_Jpsi', 'HLT_Mu5_TkMu0_OST_Jpsi', 'HLT_L1SingleEG2', 'HLT_L1SingleEG5','HLT_L1SingleEG8', 'HLT_L1DoubleEG5', 'HLT_Ele10_SW_L1R', 'HLT_Ele12_SW_TightEleId_L1R', 'HLT_Ele12_SW_TightEleIdIsol_L1R',  'HLT_Ele17_SW_L1R', 'HLT_Ele17_SW_CaloEleId_L1R', 'HLT_Ele17_SW_LooseEleId_L1R','HLT_Ele17_SW_EleId_L1R', 'HLT_Ele22_SW_CaloEleId_L1R', 'HLT_Ele40_SW_L1R', 'HLT_DoubleEle4_SW_eeRes_L1R','HLT_DoubleEle10_SW_L1R','HLT_Photon10_Cleaned_L1R','HLT_Photon15_Cleaned_L1R','HLT_Photon20_NoHE_L1R', 'HLT_Photon20_Cleaned_L1R','HLT_Photon30_Cleaned_L1R','HLT_Photon50_NoHE_L1R','HLT_Photon50_NoHE_Cleaned_L1R','HLT_DoublePhoton5_CEP_L1R','HLT_DoublePhoton5_L1R','HLT_DoublePhoton10_L1R', 'HLT_DoublePhoton15_L1R','HLT_DoublePhoton17_L1R')

# Data Double Muon  Electron
#hTozzTo4leptonsHLTAnalysisData.HLTPaths= cms.vstring('HLT_DoubleMu7_v1','HLT_L1DoubleMu0_v1','HLT_DoubleMu4_Acoplanarity03_v1','HLT_L2DoubleMu23_NoVertex_v1','HLT_L2DoubleMu0_v2','HLT_DoubleMu3_v3','HLT_DoubleMu6_v1','HLT_Mu8_Jet40_v3','HLT_TripleMu5_v2','HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2','HLT_Ele17_CaloIdL_CaloIsoVL_v2','HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2','HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2','HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2','HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2','HLT_TripleEle10_CaloIdL_TrkIdVL_v2','HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v2','HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2','HLT_Ele8_CaloIdL_CaloIsoVL_v2','HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2','HLT_Ele8_CaloIdL_TrkIdVL_v2','HLT_Ele8_v2');

# MC 3_11_2
hTozzTo4leptonsHLTAnalysisData.HLTPaths= cms.vstring('HLT_L1DoubleMuOpen',
                                                     'HLT_L2DoubleMu0',
                                                     'HLT_L2DoubleMu20_NoVertex_v1',
                                                     'HLT_DoubleMu0',
                                                     'HLT_DoubleMu0_Quarkonium_v1',
                                                     'HLT_DoubleMu3_v2',
                                                     'HLT_DoubleMu5_v1',                                                     
                                                     'HLT_Ele10_SW_L1R_v2',                                                     
                                                     'HLT_Ele12_SW_TighterEleId_L1R_v2',
                                                     'HLT_Ele17_SW_L1R_v2',
                                                     'HLT_Ele17_SW_Isol_L1R_v2',
                                                     'HLT_Ele17_SW_TighterEleIdIsol_L1R_v3',
                                                     'HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2',
                                                     'HLT_Ele22_SW_L1R_v2',
                                                     'HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2',
                                                     'HLT_Ele22_SW_TighterEleId_L1R_v3',
                                                     'HLT_Ele32_SW_TighterEleId_L1R_v2',
                                                     'HLT_Photon20_Isol_Cleaned_L1R_v1',
                                                     'HLT_DoubleEle17_SW_L1R_v1',                                                     
                                                     'HLT_Mu9',
                                                     'HLT_Mu15',
                                                     'HLT_Ele10_LW_EleId_L1R',
                                                     'HLT_Ele10_SW_EleId_L1R',
                                                     'HLT_Ele12_SW_TightEleIdIso1_L1R',
                                                     'HLT_Ele12_SW_TighterEleIdIso1_L1R',
                                                     'HLT_IsoEle12_PFTau15',
                                                     'HLT_IsoMu12',
                                                     'HLT_IsoMu12_LooseIsoPFTau10',
                                                     'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15',
                                                     'HLT_Mu10_Ele10_CaloIdL1',
                                                     'L1SingleMu7',
                                                     'L1SingleEG5',
                                                     'L1SingleEG8',
                                                     'L1SingleEG12',
                                                     );

#hTozzTo4leptonsHLTAnalysisData.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
hTozzTo4leptonsHLTAnalysisData.TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI311X")
#hTozzTo4leptonsHLTAnalysisData.TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI41X")
#from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHLTAnalysisRootTree_cfi import *
#from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHLTAnalysisSequences_cff import *

hTozzTo4leptonsHLTInfo = cms.EDProducer("HZZ4LeptonsHLTInfo",
                                        TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI311X")
                                        #TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI41X")
                                        )

# TriggerResults:SkimEarlyData
hTozzTo4leptonsSkimEarlyDataAnalysis=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHLTAnalysis_cfi.hTozzTo4leptonsHLTAnalysis.clone()
hTozzTo4leptonsSkimEarlyDataAnalysis.TriggerResultsTag = cms.InputTag("TriggerResults","","Skim")
hTozzTo4leptonsSkimEarlyDataAnalysis.HLTPaths =cms.vstring('Skim_highEnergyMuons','Skim_highEnergyElectrons','Skim_recoWMNfromPf','Skim_recoWMNfromTc','Skim_recoWENfromPf','Skim_recoWENfromTc','Skim_diMuonsJPsi','Skim_diMuonsZ','Skim_diElectronsZ','Skim_triLeptonsMuMuMu','Skim_triLeptonsMuMuEl','Skim_triLeptonsMuElEl','Skim_triLeptonsElElEl','Skim_quadLeptons4Mu','Skim_quadLeptons2Mu2El','Skim_quadLeptons4El')


# Use Early Skim and change input collections
useSkimEarlyData='false'

if useSkimEarlyData == 'true':
    # electrons
    from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronSequences_cff import *
    ELECTRON_BASE_CUT=("(pt > 5 &&" +
                       " fbrem > 0 &&" +
                       " eSuperClusterOverP < 3 &&" +
                       " hcalOverEcal < 0.15 &&" +
                       " abs(deltaPhiSuperClusterTrackAtVtx) < 0.10 &&" +
                       " abs(deltaEtaSuperClusterTrackAtVtx) < 0.02 &&" +
                       " (( isEB && sigmaIetaIeta < 0.015) ||" +
                       "  (!isEB && sigmaIetaIeta < 0.035)) )");
    hTozzTo4leptonsElectronSelector = cms.EDFilter("GsfElectronRefSelector",
                                         src = cms.InputTag("gsfElectrons"),
                                         cut = cms.string(ELECTRON_BASE_CUT),
                                         )
    hTozzTo4leptonsElectronSequence=cms.Sequence(hTozzTo4leptonsElectronIdSequence + hTozzTo4leptonsElectronSelector)

    # muons
    from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonSelector_cfi import *
    TM_ARBITRATION = "numberOfMatches('SegmentAndTrackArbitration')>0";
    MUON_BASE_CUT="(isGlobalMuon || (isTrackerMuon && "+TM_ARBITRATION+"))"
    hTozzTo4leptonsMuonSelector = cms.EDFilter("MuonRefSelector",
                                     src = cms.InputTag("muons"),
                                     cut = cms.string(MUON_BASE_CUT),
                                     )    
else:
    
    # Electron relaxed selection
    from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronSequences_cff import *
    hTozzTo4leptonsElectronSelector=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronSelector_cfi.hTozzTo4leptonsElectronSelector.clone()
    hTozzTo4leptonsElectronSelector.electronEtaMax=cms.double(2.5)
    hTozzTo4leptonsElectronSelector.electronPtMin=cms.double(0.)
    hTozzTo4leptonsElectronSelector.useEleID=cms.bool(False)
    #hTozzTo4leptonsElectronSequence=cms.Sequence(hTozzTo4leptonsElectronIdSequence + hTozzTo4leptonsElectronSelector)
    hTozzTo4leptonsElectronSequence=cms.Sequence(hTozzTo4leptonsElectronSelector)
    
    # Muon relaxed selection
    from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonSelector_cfi import *
    hTozzTo4leptonsMuonSelector=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonSelector_cfi.hTozzTo4leptonsMuonSelector.clone()
    hTozzTo4leptonsMuonSelector.isGlobalMuon=cms.bool(True)
    hTozzTo4leptonsMuonSelector.isTrackerMuon=cms.bool(False)
    hTozzTo4leptonsMuonSelector.muonPtMin=cms.double(0.)
    hTozzTo4leptonsMuonSelector.muonEtaMax=cms.double(2.5)
    hTozzTo4leptonsMuonSequence=cms.Sequence(hTozzTo4leptonsMuonSelector)
    



#*******************
#2 Leptons No Presel 
#*******************

from HiggsAnalysis.HiggsToZZ4Leptons.DiTauSelect_cfi import *

# zToEE
from HiggsAnalysis.HiggsToZZ4Leptons.zToEE_cfi import *
zToEEclone= HiggsAnalysis.HiggsToZZ4Leptons.zToEE_cfi.zToEE.clone()
zToEEclone.decay = cms.string('hTozzTo4leptonsElectronSelector@+ hTozzTo4leptonsElectronSelector@-')
zToEEclone.cut = cms.string('mass > 0 && daughter(0).pt>10 && daughter(1).pt>10')
zToEEclone.checkCharge = cms.bool(True)

zToEEcleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                              srcObject = cms.InputTag("zToEEclone"),
                              srcObjectsToRemove = cms.InputTag("zToEEclone"),
                              #deltaRMin = cms.double(0.2),
                              module_label=cms.string("EEcleaned")                                              
                              )

# zToMuMu
from HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi import *                       
zToMuMuclone=HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi.zToMuMu.clone()
zToMuMuclone.decay=cms.string('hTozzTo4leptonsMuonSelector@+ hTozzTo4leptonsMuonSelector@-')
zToMuMuclone.cut = cms.string('mass > 0 && daughter(0).pt>10 && daughter(1).pt>10')
zToMuMuclone.checkCharge = cms.bool(True)

zToMuMucleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                srcObject = cms.InputTag("zToMuMuclone"),
                                srcObjectsToRemove = cms.InputTag("zToMuMuclone"),
                                #deltaRMin = cms.double(0.2),
                                module_label=cms.string("MuMucleaned")                                              
                                )

# zToCrossLeptons
# MuE (muon is possitive)
from HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi import *  
zToMuE=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
zToMuE.decay = cms.string("hTozzTo4leptonsMuonSelector@+ hTozzTo4leptonsElectronSelector@-")
zToMuE.cut = cms.string('mass > 0 && daughter(0).pt>5 && daughter(1).pt>5')
zToMuE.checkCharge = cms.bool(True)

zToMuEcleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                srcObject = cms.InputTag("zToMuE"),
                                srcObjectsToRemove = cms.InputTag("zToMuE"),
                                #deltaRMin = cms.double(0.2),
                                module_label=cms.string("MuEcleaned")                                              
                                )
# EMu (electron is possitive)
from HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi import *  
zToEMu=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
zToEMu.decay = cms.string("hTozzTo4leptonsMuonSelector@- hTozzTo4leptonsElectronSelector@+")
zToEMu.cut = cms.string('mass > 0 && daughter(0).pt>5 && daughter(1).pt>5')
zToEMu.checkCharge = cms.bool(True)

zToEMucleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                               srcObject = cms.InputTag("zToEMu"),
                               srcObjectsToRemove = cms.InputTag("zToEMu"),
                               #deltaRMin = cms.double(0.2),
                               module_label=cms.string("EMucleaned")                                              
                               )

zToCrossLeptonsMuEcleaned = cms.EDProducer("CandViewMerger",
                                           src = cms.VInputTag( "zToMuEcleaned", "zToEMucleaned")
                                           )
zToCrossLeptonsMuEPFMET=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
zToCrossLeptonsMuEPFMET.decay = cms.string('zToCrossLeptonsMuEcleaned pfMet')

# zToMuTau
zToMuTaucleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                 srcObject = cms.InputTag("ditauvispairs:visMuHPSTaucands"),
                                 srcObjectsToRemove = cms.InputTag("ditauvispairs:visMuHPSTaucands"),
                                 #deltaRMin = cms.double(0.2),
                                 module_label=cms.string("MuTaucleaned")                                              
                                 )
zToMuTauPFMET=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
zToMuTauPFMET.decay = cms.string('zToMuTaucleaned pfMet')

# zToETau
zToETaucleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                srcObject = cms.InputTag("ditauvispairs:visEleHPSTaucands"),
                                srcObjectsToRemove = cms.InputTag("ditauvispairs:visEleHPSTaucands"),
                                #deltaRMin = cms.double(0.2),
                                module_label=cms.string("ETaucleaned")                                              
                                )
zToETauPFMET=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
zToETauPFMET.decay = cms.string('zToETaucleaned pfMet')

# zToTauTau
zToTauTaucleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                srcObject = cms.InputTag("ditauvispairs:visHPSTauHPSTaucands"),
                                srcObjectsToRemove = cms.InputTag("ditauvispairs:visHPSTauHPSTaucands"),
                                #deltaRMin = cms.double(0.2),
                                module_label=cms.string("TauTaucleaned")                                              
                                )

#*******************
#4 Leptons No Presel 
#*******************

# EEMM
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi import *
hTozzToEEMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToEEMM.decay = cms.string('zToEEcleaned zToMuMucleaned')
hTozzToEEMM.cut = cms.string('mass > 100')
hTozzToEEMM.checkCharge = cms.bool(False)

hTozzToEEMMcleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                    srcObject = cms.InputTag("hTozzToEEMM"),
                                    srcObjectsToRemove = cms.InputTag("hTozzToEEMM"),
                                    module_label=cms.string("EEMMcleaned")                                              
                                    )
hTozzToEEMMPFMET=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToEEMMPFMET.decay = cms.string('hTozzToEEMM pfMet')

# MMEE
hTozzToMMEE=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToMMEE.decay = cms.string('zToMuMucleaned zToEEcleaned')
hTozzToMMEE.cut = cms.string('mass > 100')
hTozzToMMEE.checkCharge = cms.bool(False)

hTozzToMMEEcleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                    srcObject = cms.InputTag("hTozzToMMEE"),
                                    srcObjectsToRemove = cms.InputTag("hTozzToMMEE"),
                                    module_label=cms.string("MMEEcleaned")                                              
                                    )
hTozzToMMEEPFMET=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToMMEEPFMET.decay = cms.string('hTozzToMMEE pfMet')

# MMMM
hTozzToMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToMMMM.decay = cms.string('zToMuMucleaned zToMuMucleaned')
hTozzToMMMM.cut = cms.string('mass > 100')
hTozzToMMMM.checkCharge = cms.bool(False)

hTozzToMMMMcleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                    srcObject = cms.InputTag("hTozzToMMMM"),
                                    srcObjectsToRemove = cms.InputTag("hTozzToMMMM"),
                                    module_label=cms.string("MMMMcleaned")                                              
                                    )
hTozzToMMMMPFMET=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToMMMMPFMET.decay = cms.string('hTozzToMMMM pfMet')

# EEEE
hTozzToEEEE=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToEEEE.decay = cms.string('zToEEcleaned zToEEcleaned')
hTozzToEEEE.cut = cms.string('mass > 100')
hTozzToEEEE.checkCharge = cms.bool(False)

hTozzToEEEEcleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                    srcObject = cms.InputTag("hTozzToEEEE"),
                                    srcObjectsToRemove = cms.InputTag("hTozzToEEEE"),
                                    module_label=cms.string("EEEEcleaned")                                              
                                    )
hTozzToEEEEPFMET=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToEEEEPFMET.decay = cms.string('hTozzToEEEE pfMet')

# MMME
hTozzToMMME=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToMMME.decay = cms.string('zToMuMucleaned zToCrossLeptonsMuEcleaned')
hTozzToMMME.cut = cms.string('mass > 100')
hTozzToMMME.checkCharge = cms.bool(False)

hTozzToMMMEcleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                    srcObject = cms.InputTag("hTozzToMMME"),
                                    srcObjectsToRemove = cms.InputTag("hTozzToMMME"),
                                    module_label=cms.string("MMMEcleaned")                                              
                                    )
hTozzToMMMEPFMET=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToMMMEPFMET.decay = cms.string('hTozzToMMME pfMet')

# EEME
hTozzToEEEM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToEEEM.decay = cms.string('zToEEcleaned zToCrossLeptonsMuEcleaned')
hTozzToEEEM.cut = cms.string('mass > 100')
hTozzToEEEM.checkCharge = cms.bool(False)

hTozzToEEEMcleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                    srcObject = cms.InputTag("hTozzToEEEM"),
                                    srcObjectsToRemove = cms.InputTag("hTozzToEEEM"),
                                    module_label=cms.string("EEEMcleaned")                                              
                                    )
hTozzToEEEMPFMET=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToEEEMPFMET.decay = cms.string('hTozzToEEEM pfMet')

# MMMT
hTozzToMMMT=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToMMMT.decay = cms.string('zToMuMucleaned zToMuTaucleaned')
hTozzToMMMT.cut = cms.string('mass > 100')
hTozzToMMMT.checkCharge = cms.bool(False)

hTozzToMMMTcleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                    srcObject = cms.InputTag("hTozzToMMMT"),
                                    srcObjectsToRemove = cms.InputTag("hTozzToMMMT"),
                                    module_label=cms.string("MMMTcleaned")                                              
                                    )
hTozzToMMMTPFMET=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToMMMTPFMET.decay = cms.string('hTozzToMMMT pfMet')

# EEET
hTozzToEEET=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToEEET.decay = cms.string('zToEEcleaned zToETaucleaned')
#hTozzToEEET.decay = cms.string('zToEEcleaned ditauvispairs:visEleHPSTaucands')
hTozzToEEET.cut = cms.string('mass > 100')
hTozzToEEET.checkCharge = cms.bool(False)

hTozzToEEETcleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                    srcObject = cms.InputTag("hTozzToEEET"),
                                    srcObjectsToRemove = cms.InputTag("hTozzToEEET"),
                                    module_label=cms.string("EEETcleaned")                                              
                                    )
hTozzToEEETPFMET=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
#hTozzToEEETPFMET.decay = cms.string('hTozzToEEETcleaned pfMet')
hTozzToEEETPFMET.decay = cms.string('hTozzToEEET pfMet')

# MMTT
hTozzToMMTT=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToMMTT.decay = cms.string('zToMuMucleaned zToTauTaucleaned')
hTozzToMMTT.cut = cms.string('mass > 100')
hTozzToMMTT.checkCharge = cms.bool(False)

hTozzToMMTTcleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                    srcObject = cms.InputTag("hTozzToMMTT"),
                                    srcObjectsToRemove = cms.InputTag("hTozzToMMTT"),
                                    module_label=cms.string("MMTTcleaned")                                              
                                    )
hTozzToMMTTPFMET=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToMMTTPFMET.decay = cms.string('hTozzToMMTT pfMet')

# EETT
hTozzToEETT=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToEETT.decay = cms.string('zToEEcleaned zToTauTaucleaned')
hTozzToEETT.cut = cms.string('mass > 100')
hTozzToEETT.checkCharge = cms.bool(False)

hTozzToEETTcleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                    srcObject = cms.InputTag("hTozzToEETT"),
                                    srcObjectsToRemove = cms.InputTag("hTozzToEETT"),
                                    module_label=cms.string("EETTcleaned")                                              
                                    )
hTozzToEETTPFMET=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToEETTPFMET.decay = cms.string('hTozzToEETT pfMet')

# MMET
hTozzToMMET=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToMMET.decay = cms.string('zToMuMucleaned zToETaucleaned')
hTozzToMMET.cut = cms.string('mass > 100')
hTozzToMMET.checkCharge = cms.bool(False)

hTozzToMMETcleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                    srcObject = cms.InputTag("hTozzToMMET"),
                                    srcObjectsToRemove = cms.InputTag("hTozzToMMET"),
                                    module_label=cms.string("MMETcleaned")                                              
                                    )
hTozzToMMETPFMET=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToMMETPFMET.decay = cms.string('hTozzToMMET pfMet')

# EEMT
hTozzToEEMT=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToEEMT.decay = cms.string('zToEEcleaned zToMuTaucleaned')
hTozzToEEMT.cut = cms.string('mass > 100')
hTozzToEEMT.checkCharge = cms.bool(False)

hTozzToEEMTcleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                    srcObject = cms.InputTag("hTozzToEEMT"),
                                    srcObjectsToRemove = cms.InputTag("hTozzToEEMT"),
                                    module_label=cms.string("EEMTcleaned")                                              
                                    )
hTozzToEEMTPFMET=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo2l2tau_cfi.hTozzTo2l2tau.clone()
hTozzToEEMTPFMET.decay = cms.string('hTozzToEEMT pfMet')

# zToMuMu_SS, zToEE_SS and zToCrossLeptons
from HiggsAnalysis.HiggsToZZ4Leptons.zToMuMuss_cfi import *
zToMuMussplus=HiggsAnalysis.HiggsToZZ4Leptons.zToMuMuss_cfi.zToMuMussplus.clone()
zToMuMussplus.decay = cms.string('hTozzTo4leptonsMuonSelector@+ hTozzTo4leptonsMuonSelector@+')
zToMuMussplus.cut = cms.string('mass > 0 && (daughter(0).charge>0 && daughter(1).charge>0)')
zToMuMussminus=HiggsAnalysis.HiggsToZZ4Leptons.zToMuMuss_cfi.zToMuMussminus.clone()
zToMuMussminus.decay = cms.string('hTozzTo4leptonsMuonSelector@- hTozzTo4leptonsMuonSelector@-')
zToMuMussminus.cut = cms.string('mass > 0 && (daughter(0).charge<0 && daughter(1).charge<0)')
zToMuMussmerge=HiggsAnalysis.HiggsToZZ4Leptons.zToMuMuss_cfi.zToMuMussmerge.clone()
zToMuMussmerge.src = cms.VInputTag( "zToMuMussplus", "zToMuMussminus")
zToMuMuss=cms.Sequence(zToMuMussplus+zToMuMussminus+zToMuMussmerge)

from HiggsAnalysis.HiggsToZZ4Leptons.zToEEss_cfi import *
zToEEssplus=HiggsAnalysis.HiggsToZZ4Leptons.zToEEss_cfi.zToEEssplus.clone()
zToEEssplus.decay = cms.string('hTozzTo4leptonsElectronSelector@+ hTozzTo4leptonsElectronSelector@+')
zToEEssplus.cut = cms.string('mass > 0 && (daughter(0).charge>0 && daughter(1).charge>0)')
zToEEssminus=HiggsAnalysis.HiggsToZZ4Leptons.zToEEss_cfi.zToEEssminus.clone()
zToEEssminus.decay = cms.string('hTozzTo4leptonsElectronSelector@- hTozzTo4leptonsElectronSelector@-')
zToEEssminus.cut = cms.string('mass > 0 && (daughter(0).charge<0 && daughter(1).charge<0)')
zToEEssmerge=HiggsAnalysis.HiggsToZZ4Leptons.zToEEss_cfi.zToEEssmerge.clone()
zToEEssmerge.src = cms.VInputTag( "zToEEssplus", "zToEEssminus")
zToEEss=cms.Sequence(zToEEssplus+zToEEssminus+zToEEssmerge)


from HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi import *  
zToCrossLeptons=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
zToCrossLeptons.decay = cms.string("hTozzTo4leptonsMuonSelector hTozzTo4leptonsElectronSelector")
zToCrossLeptons.checkCharge = cms.bool(False)
zToCrossLeptons.cut = cms.string('mass > 0')

dileptons = cms.EDProducer("CandViewMerger",
       src = cms.VInputTag( "zToEE", "zToMuMu","zToEEssmerge","zToMuMussmerge","zToCrossLeptons")
)

#*******************
#4 Leptons No Presel 
#*******************


# hTozzToEEMuMu
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi import *
hTozzTo4leptons=HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi.zToMuMu.clone()
hTozzTo4leptons.decay = cms.string('zToEE zToMuMu')
hTozzTo4leptons.cut = cms.string('mass > 0')

# hTozzToMuMuMuMu
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi import *
hTozzTo4leptonsMMMM=HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi.zToMuMu.clone()
hTozzTo4leptonsMMMM.decay = cms.string('zToMuMu zToMuMu')
hTozzTo4leptonsMMMM.cut = cms.string('mass > 0')

# hTozzToEEEE
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi import *
hTozzTo4leptonsEEEE=HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi.zToMuMu.clone()
hTozzTo4leptonsEEEE.decay = cms.string('zToEE zToEE')
hTozzTo4leptonsEEEE.cut = cms.string('mass > 0')


# other 4 leptons combinations with SS Z and Opposite Flavour Z
quadLeptons4Mu=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
quadLeptons4Mu.decay = cms.string('zToMuMussmerge zToMuMussmerge')
quadLeptons4Mu.cut = cms.string('mass > 0') 

quadLeptons2Mu2E=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
quadLeptons2Mu2E.decay = cms.string('zToMuMussmerge zToEEssmerge')
quadLeptons2Mu2E.cut = cms.string('mass > 0') 

quadLeptons4E=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
quadLeptons4E.decay = cms.string('zToEEssmerge zToEEssmerge')
quadLeptons4E.cut = cms.string('mass > 0') 

#one Z SS and on Z OS
quadLeptonsSSOSele=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
quadLeptonsSSOSele.decay = cms.string('zToEE zToEEssmerge')
quadLeptonsSSOSele.cut = cms.string('mass > 0') 

quadLeptonsSSOSmu=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
quadLeptonsSSOSmu.decay = cms.string('zToMuMu zToMuMussmerge')
quadLeptonsSSOSmu.cut = cms.string('mass > 0') 

quadLeptonsSSOSmuele=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
quadLeptonsSSOSmuele.decay = cms.string('zToMuMu zToEEssmerge')
quadLeptonsSSOSmuele.cut = cms.string('mass > 0') 

quadLeptonsSSOSelemu=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
quadLeptonsSSOSelemu.decay = cms.string('zToEE zToMuMussmerge')
quadLeptonsSSOSelemu.cut = cms.string('mass > 0')


#3Mu+1E, 3E+1Mu
quadLeptons3Mu1E1Z =HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi.hTozzTo4leptons.clone()
quadLeptons3Mu1E1Z.decay = cms.string('zToMuMu zToCrossLeptons')
quadLeptons3Mu1E1Z.checkCharge = cms.bool(False)
quadLeptons3Mu1E1Z.cut = cms.string('mass > 0')

quadLeptons3Mu1E0Z =HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi.hTozzTo4leptons.clone()
quadLeptons3Mu1E0Z.decay = cms.string('zToMuMussmerge zToCrossLeptons')
quadLeptons3Mu1E0Z.checkCharge = cms.bool(False)
quadLeptons3Mu1E0Z.cut = cms.string('mass > 0')

quadLeptons3E1Mu1Z =HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi.hTozzTo4leptons.clone()
quadLeptons3E1Mu1Z.decay = cms.string('zToEE zToCrossLeptons')
quadLeptons3E1Mu1Z.checkCharge = cms.bool(False)
quadLeptons3E1Mu1Z.cut = cms.string('mass > 0')

quadLeptons3E1Mu0Z =HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi.hTozzTo4leptons.clone()
quadLeptons3E1Mu0Z.decay = cms.string('zToEEssmerge zToCrossLeptons')
quadLeptons3E1Mu0Z.checkCharge = cms.bool(False)
quadLeptons3E1Mu0Z.cut = cms.string('mass > 0')

#ME + ME
quadLeptonsCrossZ=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi.hTozzTo4leptons.clone()
quadLeptonsCrossZ.decay = cms.string('zToCrossLeptons zToCrossLeptons')
quadLeptonsCrossZ.checkCharge = cms.bool(False)
quadLeptonsCrossZ.cut = cms.string('mass > 0')



#*******************
#3 Leptons No Presel 
#*******************

triLeptonsMuMuMu=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
triLeptonsMuMuMu.decay = cms.string('hTozzTo4leptonsMuonSelector hTozzTo4leptonsMuonSelector hTozzTo4leptonsMuonSelector')
triLeptonsMuMuMu.cut = cms.string("mass > 0");

triLeptonsMuMuE=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
triLeptonsMuMuE.decay = cms.string('hTozzTo4leptonsMuonSelector hTozzTo4leptonsMuonSelector hTozzTo4leptonsElectronSelector')
triLeptonsMuMuE.cut = cms.string("mass > 0");

triLeptonsMuEE=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
triLeptonsMuEE.decay = cms.string('hTozzTo4leptonsMuonSelector hTozzTo4leptonsElectronSelector hTozzTo4leptonsElectronSelector')
triLeptonsMuEE.cut = cms.string("mass > 0");

triLeptonsEEE=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
triLeptonsEEE.decay = cms.string('hTozzTo4leptonsElectronSelector hTozzTo4leptonsElectronSelector hTozzTo4leptonsElectronSelector')
triLeptonsEEE.cut = cms.string("mass > 0");


# Veto electorn and muons for isolation
vetoMuons =  cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("isGlobalMuon && isTrackerMuon && pt>5")
)

vetoElectrons =  cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag("gsfElectrons"),
    cut = cms.string("pt>7 && gsfTrack().trackerExpectedHitsInner().numberOfHits<2")
)


# Electron loose isolation

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronIsolationProducerEgamma_cfi import *
hTozzTo4leptonsElectronIsolationProducerEgamma.threshold = cms.double(99999.)

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronIsolationEgammaSequences_cff import *
hTozzTo4leptonsElectronIsolationSequenceEgamma=cms.Sequence(hTozzTo4leptonsElectronIsolationDepositSequence + hTozzTo4leptonsElectronIsolationProducerEgamma)

# Electron Id
from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi import *
eidVeryLoose.src = "hTozzTo4leptonsElectronSelector"
eidLoose.src = "hTozzTo4leptonsElectronSelector"
eidMedium.src = "hTozzTo4leptonsElectronSelector"
eidTight.src  = "hTozzTo4leptonsElectronSelector"

from RecoEgamma.ElectronIdentification.cutsInCategoriesHZZElectronIdentificationV06_cfi import *
eidHZZVeryLoose.src = "hTozzTo4leptonsElectronSelector"
eidHZZLoose.src = "hTozzTo4leptonsElectronSelector"
eidHZZMedium.src = "hTozzTo4leptonsElectronSelector"
eidHZZHyperTight1.src  = "hTozzTo4leptonsElectronSelector"

# Muon loose isolation
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonIsolationSequences_cff import *
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonIsolationProducerMu_cfi import *
hTozzTo4leptonsMuonIsolationProducerMu.threshold=cms.double(99999.)

# Common preselection 
# from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCommonPreselectionSequences_cff import *

# zToEE loose isolated
from HiggsAnalysis.HiggsToZZ4Leptons.zToEE_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.zToEE_cfi
zToEELooseIsol = HiggsAnalysis.HiggsToZZ4Leptons.zToEE_cfi.zToEE.clone()
zToEELooseIsol.decay = ('hTozzTo4leptonsElectronSelector@+ hTozzTo4leptonsElectronSelector@-')
zToEELooseIsol.cut = cms.string('mass > 0')


# zToMuMu loose isolated
from HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi
zToMuMuLooseIsol = HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi.zToMuMu.clone()
zToMuMuLooseIsol.decay = ('hTozzTo4leptonsMuonSelector@+ hTozzTo4leptonsMuonSelector@-')
zToMuMuLooseIsol.cut = cms.string('mass > 0')

# hTozzToEEMuMu loose isolated
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi
hTozzTo4leptonsLooseIsol=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi.hTozzTo4leptons.clone()
hTozzTo4leptonsLooseIsol.decay = ('zToEELooseIsol zToMuMuLooseIsol')
hTozzTo4leptonsLooseIsol.cut = cms.string('mass > 0')


# hTozzToMuMuMuMu loose isolated
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi import *
hTozzTo4leptonsMMMMLooseIsol=HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi.zToMuMu.clone()
hTozzTo4leptonsMMMMLooseIsol.decay = cms.string('zToMuMuLooseIsol zToMuMuLooseIsol')
hTozzTo4leptonsMMMMLooseIsol.cut = cms.string('mass > 0')

# hTozzToEEEE loose isolated
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi import *
hTozzTo4leptonsEEEELooseIsol=HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi.zToMuMu.clone()
hTozzTo4leptonsEEEELooseIsol.decay = cms.string('zToEELooseIsol zToEELooseIsol')
hTozzTo4leptonsEEEELooseIsol.cut = cms.string('mass > 0')


# 2e2mu best candidate producer
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsBestCandidateProducer_cfi import *

# 4mu best candidate producer
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsBestCandidateProducer_cfi import *
hTozzTo4leptonsBestCandidateProducerMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsBestCandidateProducer_cfi.hTozzTo4leptonsBestCandidateProducer.clone()
hTozzTo4leptonsBestCandidateProducerMMMM.decaychannel = cms.string('4mu')
hTozzTo4leptonsBestCandidateProducerMMMM.RECOcollName = cms.VInputTag(cms.InputTag("hTozzTo4leptonsMMMMLooseIsol"))
hTozzTo4leptonsBestCandidateProducerMMMM.decayChain = cms.string('hToZZTo4LeptonsBestCandidate')

# 4e best candidate producer
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsBestCandidateProducer_cfi import *
hTozzTo4leptonsBestCandidateProducerEEEE=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsBestCandidateProducer_cfi.hTozzTo4leptonsBestCandidateProducer.clone()
hTozzTo4leptonsBestCandidateProducerEEEE.decaychannel = cms.string('4e')
hTozzTo4leptonsBestCandidateProducerEEEE.RECOcollName = cms.VInputTag(cms.InputTag("hTozzTo4leptonsEEEELooseIsol"))
hTozzTo4leptonsBestCandidateProducerEEEE.decayChain = cms.string('hToZZTo4LeptonsBestCandidate')



# higgs frame producer 2e2mu, 4mu, 4e
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHiggsFrame_cfi import *

hTozzTo4leptonsHiggsFrameMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHiggsFrame_cfi.hTozzTo4leptonsHiggsFrame.clone()
hTozzTo4leptonsHiggsFrameMMMM.prodinst = cms.string('hTozzTo4leptonsBestCandidateProducerMMMM')
hTozzTo4leptonsHiggsFrameMMMM.RECOcollName = cms.VInputTag(cms.InputTag("hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hToZZTo4LeptonsBestCandidateBoson1"))

hTozzTo4leptonsHiggsFrameEEEE=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHiggsFrame_cfi.hTozzTo4leptonsHiggsFrame.clone()
hTozzTo4leptonsHiggsFrameEEEE.prodinst = cms.string('hTozzTo4leptonsBestCandidateProducerEEEE')
hTozzTo4leptonsHiggsFrameEEEE.RECOcollName = cms.VInputTag(cms.InputTag("hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hToZZTo4LeptonsBestCandidateBoson1"))


# CP producer 2e2mu, 4mu, 4e
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi import *
hTozzTo4leptonsCP.RECOcollNameLabFrame = cms.VInputTag(cms.InputTag("hTozzTo4leptonsBestCandidateProducer:hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hTozzTo4leptonsBestCandidateProducer:hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsBestCandidateProducer:hToZZTo4LeptonsBestCandidateBoson1"))

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi import *
hTozzTo4leptonsCPMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi.hTozzTo4leptonsCP.clone()
hTozzTo4leptonsCPMMMM.RECOcollNameLabFrame = cms.VInputTag(cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateBoson1"))

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi import *
hTozzTo4leptonsCPEEEE=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi.hTozzTo4leptonsCP.clone()
hTozzTo4leptonsCPEEEE.RECOcollNameLabFrame = cms.VInputTag(cms.InputTag("hTozzTo4leptonsBestCandidateProducerEEEE:hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerEEEE:hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerEEEE:hToZZTo4LeptonsBestCandidateBoson1"))

# 3D IP KF
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi import *
hTozzTo4leptonsIpToVtxProducerKF=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi.hTozzTo4leptonsIpToVtxProducer.clone()
hTozzTo4leptonsIpToVtxProducerKF.VertexLabel = cms.InputTag("offlinePrimaryVertices")

# Deterministic annealing
# from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi import *

# 3D IP DA
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi import *
hTozzTo4leptonsIpToVtxProducer.VertexLabel = cms.InputTag("offlinePrimaryVertices")

# 2D IP DA
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsTipLipToVtxProducer_cfi import *
hTozzTo4leptonsTipLipToVtxProducer.VertexLabel = cms.InputTag("offlinePrimaryVertices")

# Geometrical Discriminator
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsGeomDiscrimProducer_cfi import *
#hTozzTo4leptonsGeomDiscrimProducer.RECOcollName=cms.InputTag("hTozzTo4leptonsBestCandidateProducer:hToZZTo4LeptonsBestCandidateMother")
hTozzTo4leptonsGeomDiscrimProducerMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsGeomDiscrimProducer_cfi.hTozzTo4leptonsGeomDiscrimProducer.clone()
hTozzTo4leptonsGeomDiscrimProducerMMMM.RECOcollName=cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateMother")
hTozzTo4leptonsGeomDiscrimProducerEEEE=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsGeomDiscrimProducer_cfi.hTozzTo4leptonsGeomDiscrimProducer.clone()
hTozzTo4leptonsGeomDiscrimProducerEEEE.RECOcollName=cms.InputTag("hTozzTo4leptonsBestCandidateProducerEEEE:hToZZTo4LeptonsBestCandidateMother")

# Constrained fit: input 4l:2e2mu
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsConstraintFitProducer_cfi import *
hTozzTo4leptonsConstraintFitProducer.VertexLabel = cms.InputTag("offlinePrimaryVertices")
# 4mu
hTozzTo4leptonsConstraintFitProducerMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsConstraintFitProducer_cfi.hTozzTo4leptonsConstraintFitProducer.clone()
hTozzTo4leptonsConstraintFitProducerMMMM.RECOcollName =cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateMother")
hTozzTo4leptonsConstraintFitProducerMMMM.VertexLabel = cms.InputTag("offlinePrimaryVertices")
# 4e
hTozzTo4leptonsConstraintFitProducerEEEE=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsConstraintFitProducer_cfi.hTozzTo4leptonsConstraintFitProducer.clone()
hTozzTo4leptonsConstraintFitProducerEEEE.RECOcollName =cms.InputTag("hTozzTo4leptonsBestCandidateProducerEEEE:hToZZTo4LeptonsBestCandidateMother")
hTozzTo4leptonsConstraintFitProducerEEEE.VertexLabel = cms.InputTag("offlinePrimaryVertices")

# 3D IP with GDvertexFitter
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi import *
hTozzTo4leptonsIpToVtxProducerGD=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi.hTozzTo4leptonsIpToVtxProducer.clone()
hTozzTo4leptonsIpToVtxProducerGD.VertexLabel = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:GDFitVertex")

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi import *
hTozzTo4leptonsIpToVtxProducerGDMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi.hTozzTo4leptonsIpToVtxProducer.clone()
hTozzTo4leptonsIpToVtxProducerGDMMMM.VertexLabel = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:GDFitVertex")

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi import *
hTozzTo4leptonsIpToVtxProducerGDEEEE=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi.hTozzTo4leptonsIpToVtxProducer.clone()
hTozzTo4leptonsIpToVtxProducerGDEEEE.VertexLabel = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:GDFitVertex")

# 3D IP with Standard vertexFitter
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi import *
hTozzTo4leptonsIpToVtxProducerStd=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi.hTozzTo4leptonsIpToVtxProducer.clone()
hTozzTo4leptonsIpToVtxProducerStd.VertexLabel = cms.InputTag("hTozzTo4leptonsConstraintFitProducer:StandardFitVertex")

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi import *
hTozzTo4leptonsIpToVtxProducerStdMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi.hTozzTo4leptonsIpToVtxProducer.clone()
hTozzTo4leptonsIpToVtxProducerStdMMMM.VertexLabel = cms.InputTag("hTozzTo4leptonsConstraintFitProducerMMMM:StandardFitVertex")

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi import *
hTozzTo4leptonsIpToVtxProducerStdEEEE=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi.hTozzTo4leptonsIpToVtxProducer.clone()
hTozzTo4leptonsIpToVtxProducerStdEEEE.VertexLabel = cms.InputTag("hTozzTo4leptonsConstraintFitProducerEEEE:StandardFitVertex")

# 3D IP with Kinematic vertexFitter
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi import *
hTozzTo4leptonsIpToVtxProducerKin=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi.hTozzTo4leptonsIpToVtxProducer.clone()
hTozzTo4leptonsIpToVtxProducerKin.VertexLabel = cms.InputTag("hTozzTo4leptonsConstraintFitProducer:KinematicFitVertex")

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi import *
hTozzTo4leptonsIpToVtxProducerKinMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi.hTozzTo4leptonsIpToVtxProducer.clone()
hTozzTo4leptonsIpToVtxProducerKinMMMM.VertexLabel = cms.InputTag("hTozzTo4leptonsConstraintFitProducerMMMM:KinematicFitVertex")

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi import *
hTozzTo4leptonsIpToVtxProducerKinEEEE=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi.hTozzTo4leptonsIpToVtxProducer.clone()
hTozzTo4leptonsIpToVtxProducerKinEEEE.VertexLabel = cms.InputTag("hTozzTo4leptonsConstraintFitProducerEEEE:KinematicFitVertex")

# COMMON ROOT TREE
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCommonRootTree_cfi  import *
hTozzTo4leptonsCommonRootTreePresel=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCommonRootTree_cfi.hTozzTo4leptonsCommonRootTree.clone()
hTozzTo4leptonsCommonRootTreePresel.decaychannel = cms.string('2l2tau')
hTozzTo4leptonsCommonRootTreePresel.rootFileName = cms.untracked.string('roottree_leptons.root')
# hlt
hTozzTo4leptonsCommonRootTreePresel.fillHLTinfo = cms.untracked.bool(True)                                           
hTozzTo4leptonsCommonRootTreePresel.HLTAnalysisinst = cms.string('hTozzTo4leptonsHLTAnalysisData')


# Data DoubleElectron and Muon
#hTozzTo4leptonsCommonRootTreePresel.flagHLTnames=cms.VInputTag(cms.InputTag("flagHLTDoubleMu7v1"),cms.InputTag("flagHLTL1DoubleMu0v1"),cms.InputTag("flagHLTDoubleMu4Acoplanarity03v1"),cms.InputTag("flagHLTL2DoubleMu23NoVertexv1"),cms.InputTag("flagHLTL2DoubleMu0v2"),cms.InputTag("flagHLTDoubleMu3v3"),cms.InputTag("flagHLTDoubleMu6v1"),cms.InputTag("flagHLTMu8Jet40v3"),cms.InputTag("flagHLTTripleMu5v2"),cms.InputTag("flagHLTEle17CaloIdLCaloIsoVLEle8CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTEle17CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTEle8CaloIdLCaloIsoVLJet40v2"),cms.InputTag("flagHLTEle17CaloIdTTrkIdVLCaloIsoVLTrkIsoVLEle8CaloIdTTrkIdVLCaloIsoVLTrkIsoVLv2"),cms.InputTag("flagHLTPhoton20CaloIdVTIsoTEle8CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTEle32CaloIdLCaloIsoVLSC17v2"),cms.InputTag("flagHLTTripleEle10CaloIdLTrkIdVLv2"),cms.InputTag("flagHLTEle17CaloIdLCaloIsoVLEle15HFLv2"),cms.InputTag("flagHLTEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8Mass30v2"),cms.InputTag("flagHLTEle8CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTDoubleEle10CaloIdLTrkIdVLEle10v2"),cms.InputTag("flagHLTEle8CaloIdLTrkIdVLv2"),cms.InputTag("flagHLTEle8v2"),cms.InputTag("flagHLTaccept"))

# MC 3_11_2
hTozzTo4leptonsCommonRootTreePresel.flagHLTnames=cms.VInputTag(cms.InputTag("flagHLTL1DoubleMuOpen"),cms.InputTag("flagHLTL2DoubleMu0"),cms.InputTag("flagHLTL2DoubleMu20NoVertexv1"),cms.InputTag("flagHLTDoubleMu0"),cms.InputTag("flagHLTDoubleMu0Quarkoniumv1"),cms.InputTag("flagHLTDoubleMu3v2"),cms.InputTag("flagHLTDoubleMu5v1"),cms.InputTag("flagHLTEle10SWL1Rv2"),cms.InputTag("flagHLTEle12SWTighterEleIdL1Rv2"),cms.InputTag("flagHLTEle17SWL1Rv2"),cms.InputTag("flagHLTEle17SWIsolL1Rv2"),cms.InputTag("flagHLTEle17SWTighterEleIdIsolL1Rv3"),cms.InputTag("flagHLTEle17SWTightCaloEleIdEle8HEL1Rv2"),cms.InputTag("flagHLTEle22SWL1Rv2"),cms.InputTag("flagHLTEle22SWTighterCaloIdIsolL1Rv2"),cms.InputTag("flagHLTEle22SWTighterEleIdL1Rv3"),cms.InputTag("flagHLTEle32SWTighterEleIdL1Rv2"),cms.InputTag("flagHLTPhoton20IsolCleanedL1Rv1"),cms.InputTag("flagHLTDoubleEle17SWL1Rv1"),cms.InputTag("flagHLTaccept"))


# skimEarlyData
if useSkimEarlyData == 'true':
    hTozzTo4leptonsCommonRootTreePresel.useSkimEarlyData = cms.untracked.bool(True)
else:
    hTozzTo4leptonsCommonRootTreePresel.useSkimEarlyData = cms.untracked.bool(False) 
hTozzTo4leptonsCommonRootTreePresel.SkimEarlyDataAnalysisinst = cms.string('hTozzTo4leptonsSkimEarlyDataAnalysis')
hTozzTo4leptonsCommonRootTreePresel.flagSkimEarlyDatanames=cms.VInputTag(cms.InputTag("flagSkimhighEnergyMuons"),cms.InputTag("flagSkimhighEnergyElectrons"),cms.InputTag("flagSkimrecoWMNfromPf"),cms.InputTag("flagSkimrecoWMNfromTc"),cms.InputTag("flagSkimrecoWENfromPf"),cms.InputTag("flagSkimrecoWENfromTc"),cms.InputTag("flagSkimdiMuonsJPsi"),cms.InputTag("flagSkimdiMuonsZ"),cms.InputTag("flagSkimdiElectronsZ"),cms.InputTag("flagSkimtriLeptonsMuMuMu"),cms.InputTag("flagSkimtriLeptonsMuMuEl"),cms.InputTag("flagSkimtriLeptonsMuElEl"),cms.InputTag("flagSkimtriLeptonsElElEl"),cms.InputTag("flagSkimquadLeptons4Mu"),cms.InputTag("flagSkimquadLeptons2Mu2El"),cms.InputTag("flagSkimquadLeptons4El"))
# presel
hTozzTo4leptonsCommonRootTreePresel.flaginst = cms.string('hTozzTo4leptonsCommonPreselection')
hTozzTo4leptonsCommonRootTreePresel.flagtags = cms.vstring('PreselAtleast2Ele','PreselAtleast2Mu','PreselAtleast1ZEE','PreselAtleast1ZMuMu','PreselAtleast1H','PreselLoose2IsolEle','PreselLoose2IsolMu')
# MC truth
hTozzTo4leptonsCommonRootTreePresel.fillMCTruth  = cms.untracked.bool(False)
#hTozzTo4leptonsCommonRootTreePresel.MCcollName = cms.VInputTag(cms.InputTag("hTozzTo4leptonsMCDumper2e2mu:hToZZTo4LeptonsMCDumperMother"), cms.InputTag("hTozzTo4leptonsMCDumper2e2mu:hToZZTo4LeptonsMCDumperBoson0"), cms.InputTag("hTozzTo4leptonsMCDumper2e2mu:hToZZTo4LeptonsMCDumperBoson1"), cms.InputTag("hTozzTo4leptonsMCDumper2e2mu:hToZZTo4LeptonsMCDumperLepton0"), cms.InputTag("hTozzTo4leptonsMCDumper2e2mu:hToZZTo4LeptonsMCDumperLepton1"), cms.InputTag("hTozzTo4leptonsMCDumper2e2mu:hToZZTo4LeptonsMCDumperLepton2"), cms.InputTag("hTozzTo4leptonsMCDumper2e2mu:hToZZTo4LeptonsMCDumperLepton3"))
hTozzTo4leptonsCommonRootTreePresel.MCcollName = cms.VInputTag(cms.InputTag("hTozzTo4leptonsMCBestCandidateProducer:hToZZTo4LeptonsMCBestCandidateMother"), cms.InputTag("hTozzTo4leptonsMCBestCandidateProducer:hToZZTo4LeptonsMCBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsMCBestCandidateProducer:hToZZTo4LeptonsMCBestCandidateBoson1"))
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameBest2e2mu= cms.VInputTag(cms.InputTag("hTozzTo4leptonsBestCandidateProducer:hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hTozzTo4leptonsBestCandidateProducer:hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsBestCandidateProducer:hToZZTo4LeptonsBestCandidateBoson1"))
#hTozzTo4leptonsCommonRootTreePresel.RECOcollName=cms.VInputTag(cms.InputTag("hTozzTo4leptonsLooseIsol"),cms.InputTag("zToMuMuLooseIsol"), cms.InputTag("zToEELooseIsol"))
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameBest4mu= cms.VInputTag(cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateBoson1"))
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameBest4e= cms.VInputTag(cms.InputTag("hTozzTo4leptonsBestCandidateProducerEEEE:hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerEEEE:hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerEEEE:hToZZTo4LeptonsBestCandidateBoson1"))
hTozzTo4leptonsCommonRootTreePresel.useAdditionalRECO  = cms.untracked.bool(False)
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameZ      =cms.VInputTag(cms.InputTag("zToMuMu"), cms.InputTag("zToEE"))
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameZss    = cms.VInputTag(cms.InputTag("zToMuMussmerge"),cms.InputTag("zToEEssmerge"),cms.InputTag("zToCrossLeptons"))
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameDiLep  = cms.InputTag("dileptons")
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameEEMM   = cms.VInputTag(cms.InputTag("hTozzTo4leptonsLooseIsol"))
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameMMMM   = cms.VInputTag(cms.InputTag("hTozzTo4leptonsMMMMLooseIsol"))
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameEEEE   = cms.VInputTag(cms.InputTag("hTozzTo4leptonsEEEELooseIsol"))
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameLLL    = cms.VInputTag(cms.InputTag("triLeptonsMuMuMu"),cms.InputTag("triLeptonsMuMuE"),cms.InputTag("triLeptonsMuEE"),cms.InputTag("triLeptonsEEE"))
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameLLLLss = cms.VInputTag(cms.InputTag("quadLeptons4Mu"),cms.InputTag("quadLeptons2Mu2E"),cms.InputTag("quadLeptons4E"))
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameLLLLssos = cms.VInputTag(cms.InputTag("quadLeptonsSSOSele"),cms.InputTag("quadLeptonsSSOSmu"),cms.InputTag("quadLeptonsSSOSelemu"),cms.InputTag("quadLeptonsSSOSmuele"))
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameLLLl   = cms.VInputTag(cms.InputTag("quadLeptons3Mu1E0Z"),cms.InputTag("quadLeptons3Mu1E1Z"),cms.InputTag(
"quadLeptons3E1Mu0Z"),cms.InputTag("quadLeptons3E1Mu1Z") )
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameLLLL   = cms.InputTag("allLLLL")



hTozzTo4leptonsCommonRootTreePresel.MuonsLabel     = cms.InputTag("hTozzTo4leptonsMuonSelector")
hTozzTo4leptonsCommonRootTreePresel.MuonsMapLabel  = cms.InputTag("hTozzTo4leptonsMuonSelector")
hTozzTo4leptonsCommonRootTreePresel.MuonsTkMapLabel = cms.InputTag("muIsoFromDepsTkOptimized")
hTozzTo4leptonsCommonRootTreePresel.MuonsEcalMapLabel = cms.InputTag("muIsoFromDepsEcalOptimized")
hTozzTo4leptonsCommonRootTreePresel.MuonsHcalMapLabel = cms.InputTag("muIsoFromDepsHcalOptimized")
hTozzTo4leptonsCommonRootTreePresel.MuonsLabelVert = cms.InputTag("hTozzTo4leptonsMuonSelector")
hTozzTo4leptonsCommonRootTreePresel.MuonsMapLabelVert = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexMuMap")
hTozzTo4leptonsCommonRootTreePresel.MuonsMapLabelVertValue    = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexValueMuMap")
hTozzTo4leptonsCommonRootTreePresel.MuonsMapLabelVertError    = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexErrorMuMap")
hTozzTo4leptonsCommonRootTreePresel.MuonsSTIPMapLabelVert     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipMuMap")
hTozzTo4leptonsCommonRootTreePresel.MuonsSLIPMapLabelVert     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipMuMap")
hTozzTo4leptonsCommonRootTreePresel.MuonsSTIPMapLabelVertValue = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipValueMuMap")
hTozzTo4leptonsCommonRootTreePresel.MuonsSLIPMapLabelVertValue = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipValueMuMap")
hTozzTo4leptonsCommonRootTreePresel.MuonsSTIPMapLabelVertError = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipErrorMuMap")
hTozzTo4leptonsCommonRootTreePresel.MuonsSLIPMapLabelVertError = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipErrorMuMap")

hTozzTo4leptonsCommonRootTreePresel.ElectronsEgmLabel        = cms.InputTag("hTozzTo4leptonsElectronSelector")
hTozzTo4leptonsCommonRootTreePresel.ElectronsEgmTkMapLabel   = cms.InputTag("eleIsoFromDepsTkOptimized")
hTozzTo4leptonsCommonRootTreePresel.ElectronsEgmEcalMapLabel = cms.InputTag("eleIsoFromDepsEcalFromHitsByCrystalOptimized")
hTozzTo4leptonsCommonRootTreePresel.ElectronsEgmHcalMapLabel = cms.InputTag("eleIsoFromDepsHcalFromTowersOptimized")        
hTozzTo4leptonsCommonRootTreePresel.ElectronsLabelVert         = cms.InputTag("hTozzTo4leptonsElectronSelector")
hTozzTo4leptonsCommonRootTreePresel.ElectronsMapLabelVert      = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexEleMap")
hTozzTo4leptonsCommonRootTreePresel.ElectronsMapLabelVertValue = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexValueEleMap")
hTozzTo4leptonsCommonRootTreePresel.ElectronsMapLabelVertError = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexErrorEleMap")


# CP variables
hTozzTo4leptonsCommonRootTreePresel.MCCP_2e2mu_PhiLabel          = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPPhi")
hTozzTo4leptonsCommonRootTreePresel.MCCP_2e2mu_Phi1Label         = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPPhi1")
hTozzTo4leptonsCommonRootTreePresel.MCCP_2e2mu_Phi2Label         = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPPhi2")
hTozzTo4leptonsCommonRootTreePresel.MCCP_2e2mu_phi1RFLabel       = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPphi1RF")
hTozzTo4leptonsCommonRootTreePresel.MCCP_2e2mu_phi2RFLabel       = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPphi2RF")
hTozzTo4leptonsCommonRootTreePresel.MCCP_2e2mu_cosThetaStarLabel = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPcosThetaStar")
hTozzTo4leptonsCommonRootTreePresel.MCCP_2e2mu_cosTheta1Label    = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPcosTheta1")
hTozzTo4leptonsCommonRootTreePresel.MCCP_2e2mu_cosTheta2Label    = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPcosTheta2")


hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_PhiLabel          = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPPhi")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_Phi1Label         = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPPhi1")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_Phi2Label         = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPPhi2")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_phi1RFLabel       = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPphi1RF")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_phi2RFLabel       = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPphi2RF")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_cosThetaStarLabel = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPcosThetaStar")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_cosTheta1Label    = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPcosTheta1")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_cosTheta2Label    = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPcosTheta2")

hTozzTo4leptonsCommonRootTreePresel.CP4mu_PhiLabel          = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPPhi")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_Phi1Label         = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPPhi1")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_Phi2Label         = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPPhi2")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_phi1RFLabel       = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPphi1RF")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_phi2RFLabel       = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPphi2RF")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_cosThetaStarLabel = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPcosThetaStar")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_cosTheta1Label    = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPcosTheta1")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_cosTheta2Label    = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPcosTheta2")

hTozzTo4leptonsCommonRootTreePresel.CP4e_PhiLabel          = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPPhi")
hTozzTo4leptonsCommonRootTreePresel.CP4e_Phi1Label         = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPPhi1")
hTozzTo4leptonsCommonRootTreePresel.CP4e_Phi2Label         = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPPhi2")
hTozzTo4leptonsCommonRootTreePresel.CP4e_phi1RFLabel       = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPphi1RF")
hTozzTo4leptonsCommonRootTreePresel.CP4e_phi2RFLabel       = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPphi2RF")
hTozzTo4leptonsCommonRootTreePresel.CP4e_cosThetaStarLabel = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPcosThetaStar")
hTozzTo4leptonsCommonRootTreePresel.CP4e_cosTheta1Label    = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPcosTheta1")
hTozzTo4leptonsCommonRootTreePresel.CP4e_cosTheta2Label    = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPcosTheta2")



######Conversion
from HiggsAnalysis.HiggsToZZ4Leptons.ConvValueMapProd_cfi  import *

#### FastJet corrections
from RecoJets.JetProducers.kt4PFJets_cfi import *
import RecoJets.JetProducers.kt4PFJets_cfi
kt6corPFJets=RecoJets.JetProducers.kt4PFJets_cfi.kt4PFJets.clone()
kt6corPFJets.rParam = cms.double(0.6)
kt6corPFJets.doRhoFastjet = cms.bool(True) 
kt6corPFJets.Rho_EtaMax = cms.double(2.5)


from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *
patTrigger.processName=cms.string( "REDIGI311X")

muonTriggerMatchHLT = cms.EDProducer( 'PATTriggerMatcherDRDPtLessByR',
    src     = cms.InputTag( 'hTozzTo4leptonsMuonSelector' ),
    matched = cms.InputTag( 'patTrigger' ),
    matchedCuts = cms.string( 'path( "HLT_DoubleMu*" )' ),                                 
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( '*' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring( '*' ),
    collectionTags = cms.vstring( 'hltL3MuonCandidates' ),
    maxDPtRel = cms.double( 1. ),
    maxDeltaR = cms.double( 0.2 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( False )
)

muonTriggerMatchHLTasym = cms.EDProducer( 'PATTriggerMatcherDRDPtLessByR',
    src     = cms.InputTag( 'hTozzTo4leptonsMuonSelector' ),
    matched = cms.InputTag( 'patTrigger' ),
    matchedCuts = cms.string( 'path( "HLT_Mu13_Mu8_v*" )' ),                                 
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( '*' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring( '*' ),
    collectionTags = cms.vstring( 'hltL3MuonCandidates' ),
    maxDPtRel = cms.double( 1. ),
    maxDeltaR = cms.double( 0.2 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( False )
)


electronTriggerMatchHLT = cms.EDProducer( 'PATTriggerMatcherDRDPtLessByR',
    src     = cms.InputTag( 'hTozzTo4leptonsElectronSelector' ),
    matched = cms.InputTag( 'patTrigger' ),
    matchedCuts = cms.string( 'path( "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*" ) || path( "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*" ) || path( "HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v*" )' ),                                 
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( '*' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring( '*' ),
    collectionTags = cms.vstring( '*' ),
    maxDPtRel = cms.double( 1. ),
    maxDeltaR = cms.double( 0.2 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( False )
)

#from PhysicsTools.PatAlgos.triggerLayer1.triggerMatcher_cfi import *
#cleanMuonTriggerMatchHLTDoubleIsoMu3.src= cms.InputTag( "hTozzTo4leptonsMuonSelector" )
#cleanMuonTriggerMatchHLTDoubleIsoMu3.matchedCuts = cms.string( 'path( "HLT_DoubleMu3_v*" )' )

#from PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi import *
#patTriggerEvent.patTriggerMatches = cms.VInputTag("muonTriggerMatchHLT")
#patTriggerEvent.processName= cms.string( 'REDIGI311X' )               # default; change only, if you know exactly, what you are doing!



hTozzTo4leptonsSelectionSequence2l2tau = cms.Sequence(
#        hTozzTo4leptonsMCGenFilter2l2tau             +
#        hTozzTo4leptonsMCGenParticleListDrawer2l2tau +
#        hTozzTo4leptonsMCDumper2l2tau                + hTozzTo4leptonsMCDumper4mu + hTozzTo4leptonsMCDumper4e +
#        hTozzTo4leptonsMCBestCandidateProducer      + hTozzTo4leptonsMCBestCandidateProducer4mu + hTozzTo4leptonsMCBestCandidateProducer4e +
#        hTozzTo4leptonsMCCP2l2tau                    + hTozzTo4leptonsMCCP4mu + hTozzTo4leptonsMCCP4e +
#        higgsToZZ4LeptonsSequenceData               +
        hTozzTo4leptonsHLTAnalysisData              +
        hTozzTo4leptonsHLTInfo                      +
#       hTozzTo4leptonsSkimEarlyDataAnalysis        +
	hTozzTo4leptonsElectronSelector             +
	hTozzTo4leptonsMuonSelector                 +
        #
       zToEE                                       +
       zToMuMu                                     +
       hTozzTo4leptons                             +
       hTozzTo4leptonsMMMM                         +
       hTozzTo4leptonsEEEE                         +
#       additional collection
        zToEEss                                     +
        zToMuMuss                                   +
        zToCrossLeptons                             +
        dileptons                                   +
 ##        triLeptonsMuMuMu                            +
##         triLeptonsMuMuE                             +
##         triLeptonsMuEE                              +
##         triLeptonsEEE                               +
##         quadLeptons4Mu                              +
##         quadLeptons2Mu2E                            +
##         quadLeptons4E                               +
##         quadLeptons3Mu1E0Z                          +
##         quadLeptons3Mu1E1Z                          +
##         quadLeptons3E1Mu0Z                          +
##         quadLeptons3E1Mu1Z                          +
##         quadLeptonsCrossZ                           +
##         quadLeptonsSSOSmu                           +
##         quadLeptonsSSOSele                          +
##         quadLeptonsSSOSmuele                        +
##         quadLeptonsSSOSelemu                        +
        #
        ditauvispairs                               +
        zToEEclone                                  +
        zToEEcleaned                                +
        zToMuMuclone                                +
        zToMuMucleaned                              +
        zToMuE                                      +
        zToMuEcleaned                               +
        zToEMu                                      +
        zToEMucleaned                               +
        zToCrossLeptonsMuEcleaned                   +
        zToCrossLeptonsMuEPFMET                     +
        zToMuTaucleaned                             +
        zToMuTauPFMET                               +
        zToETaucleaned                              +
        zToETauPFMET                                +
        zToTauTaucleaned                            +
        hTozzToEEMM                                 +
        hTozzToEEMMcleaned                          +
        hTozzToEEMMPFMET                            +
        hTozzToMMEE                                 +
        hTozzToMMEEcleaned                          +
        hTozzToMMEEPFMET                            +
        hTozzToMMMM                                 +
        hTozzToMMMMcleaned                          +
        hTozzToMMMMPFMET                            +
        hTozzToEEEE                                 +
        hTozzToEEEEcleaned                          +
        hTozzToEEEEPFMET                            +
        hTozzToEEEM                                 +
        hTozzToEEEMcleaned                          +
        hTozzToEEEMPFMET                            +
        hTozzToMMME                                 +
        hTozzToMMMEcleaned                          +
        hTozzToMMMEPFMET                            +
        hTozzToMMMT                                 +
        hTozzToMMMTcleaned                          +
        hTozzToMMMTPFMET                            +
        hTozzToEEET                                 +
        hTozzToEEETcleaned                          +
        hTozzToEEETPFMET                            +
        hTozzToMMTT                                 +
        hTozzToMMTTcleaned                          +
        hTozzToMMTTPFMET                            +
        hTozzToEETT                                 +
        hTozzToEETTcleaned                          +
        hTozzToEETTPFMET                            +
        hTozzToMMET                                 +
        hTozzToMMETcleaned                          +
        hTozzToMMETPFMET                            +
        hTozzToEEMT                                 +
        hTozzToEEMTcleaned                          +
        hTozzToEEMTPFMET                            +

#       additional collection
#        zToEEss                                     +
#        zToMuMuss                                   +
#        zToCrossLeptons                             +
#        dileptons                                   +
#        triLeptonsMuMuMu                            +
#        triLeptonsMuMuE                             +
#        triLeptonsMuEE                              +
#        triLeptonsEEE                               +
#        quadLeptons4Mu                              +
#        quadLeptons2Mu2E                            +
#        quadLeptons4E                               +
#        quadLeptons3Mu1E0Z                          +
#        quadLeptons3Mu1E1Z                          +
#        quadLeptons3E1Mu0Z                          +
#        quadLeptons3E1Mu1Z                          +
#        quadLeptonsCrossZ                           +
#        quadLeptonsSSOSmu                           +
#        quadLeptonsSSOSele                          +
#        quadLeptonsSSOSmuele                        +
#        quadLeptonsSSOSelemu                        +       
#        allLLLL                                     +
        vetoMuons                                   +
        vetoElectrons                               +
        hTozzTo4leptonsElectronIsolationSequenceEgamma   +
        eidVeryLoose                                + 
        eidLoose                                    + 
        eidMedium                                   +
        eidTight                                    + 
        eidHZZVeryLoose                             + 
        eidHZZLoose                                 + 
        eidHZZMedium                                +
        eidHZZHyperTight1                           + 
        hTozzTo4leptonsMuonIsolationSequence        +
        hTozzTo4leptonsMuonIsolationProducerMu      +   
       zToEELooseIsol                              +
       zToMuMuLooseIsol                            +
       hTozzTo4leptonsLooseIsol                    +
       hTozzTo4leptonsMMMMLooseIsol                +
       hTozzTo4leptonsEEEELooseIsol                +
        hTozzTo4leptonsBestCandidateProducer        +
       hTozzTo4leptonsBestCandidateProducerMMMM    +
       hTozzTo4leptonsBestCandidateProducerEEEE    +
#        hTozzTo4leptonsCP                           +
#        hTozzTo4leptonsCPMMMM                       +
#        hTozzTo4leptonsCPEEEE                       +
        hTozzTo4leptonsIpToVtxProducerKF            +
#        offlinePrimaryVerticesDA                    +
        hTozzTo4leptonsIpToVtxProducer              +
        hTozzTo4leptonsTipLipToVtxProducer          +
        # not needed
       hTozzTo4leptonsGeomDiscrimProducer          +
       hTozzTo4leptonsGeomDiscrimProducerMMMM      +
       hTozzTo4leptonsGeomDiscrimProducerEEEE      +
       hTozzTo4leptonsConstraintFitProducer        +
       hTozzTo4leptonsConstraintFitProducerMMMM    +
       hTozzTo4leptonsConstraintFitProducerEEEE    +
       hTozzTo4leptonsIpToVtxProducerGD            +
       hTozzTo4leptonsIpToVtxProducerGDMMMM        +
       hTozzTo4leptonsIpToVtxProducerGDEEEE        +
       hTozzTo4leptonsIpToVtxProducerStd           +
       hTozzTo4leptonsIpToVtxProducerStdMMMM       +
       hTozzTo4leptonsIpToVtxProducerStdEEEE       +
       hTozzTo4leptonsIpToVtxProducerKin           +
       hTozzTo4leptonsIpToVtxProducerKinMMMM       +
       hTozzTo4leptonsIpToVtxProducerKinEEEE       +
#        hTozzTo4leptonsConstraintFitProducerDiLeptons     +
#        hTozzTo4leptonsConstraintFitProducerTriLeptonsMMM +
#        hTozzTo4leptonsConstraintFitProducerTriLeptonsMME +
#        hTozzTo4leptonsConstraintFitProducerTriLeptonsEEE +
#        hTozzTo4leptonsConstraintFitProducerTriLeptonsMEE +
        ConvValueMapProd                            +
        kt6corPFJets                                +
        patTrigger                                  +
        muonTriggerMatchHLT                         +
        muonTriggerMatchHLTasym                    +
        electronTriggerMatchHLT                     +
        hTozzTo4leptonsCommonRootTreePresel        
	)
                                                 

