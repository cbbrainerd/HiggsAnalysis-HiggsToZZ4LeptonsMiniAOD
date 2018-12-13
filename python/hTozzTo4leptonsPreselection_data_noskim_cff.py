import FWCore.ParameterSet.Config as cms

# Generic MC Truth analysis
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsGenSequence_cff import *

## 
# Filter to select 2e2mu events
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi
hTozzTo4leptonsMCGenFilter2e2mu = HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi.hTozzTo4leptonsMCGenFilter.clone()
hTozzTo4leptonsMCGenFilter2e2mu.HZZ4LeptonsMCFilterLeptonFlavour=cms.int32(3)

# ParticleListDrawer
from  HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenParticleListDrawer_cfi import * 
hTozzTo4leptonsMCGenParticleListDrawer2e2mu = HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenParticleListDrawer_cfi.hTozzTo4leptonsMCGenParticleListDrawer.clone()

# Save MC truth: 
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCDumper_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCDumper_cfi
#hTozzTo4leptonsMCDumper.status=cms.vint32(62)

# CP producer: 
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi import *
hTozzTo4leptonsMCCP=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi.hTozzTo4leptonsCP.clone()
hTozzTo4leptonsMCCP.RECOcollName = cms.InputTag("hTozzTo4leptonsMCDumper")

# PF muons
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsPFtoRECOMuon_cfi import *

# PF photons
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsPFfsrPhoton_cfi import *



from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHLTInfo_cfi import *
hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
hTozzTo4leptonsHLTInfo.debug = cms.untracked.bool(False)

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHLTAnalysisFilter_cfi import *


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
    hTozzTo4leptonsElectronSelector = cms.EDFilter("PatElectronRefSelector",
                                         src = cms.InputTag("slimmedElectrons"),
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
    
    # Electron Preselector
    from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronSequences_cff import *
    hTozzTo4leptonsElectronPreSelector=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronSelector_cfi.hTozzTo4leptonsElectronSelector.clone()
    hTozzTo4leptonsElectronPreSelector.electronCollection=cms.InputTag("slimmedElectrons")
    hTozzTo4leptonsElectronPreSelector.electronEtaMax=cms.double(2.5)
    hTozzTo4leptonsElectronPreSelector.electronPtMin=cms.double(5.)
    hTozzTo4leptonsElectronPreSelector.useEleID=cms.bool(False)
    
    # Electron ordering in pT
    hTozzTo4leptonsElectronOrdering = cms.EDProducer("HZZ4LeptonsElectronOrdering",
                                                     electronCollection = cms.InputTag("hTozzTo4leptonsElectronPreSelector"),
                                                     )
    
#@#    # Electron scale calibration
#@#    # from EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi import *
    
#@#    from EgammaAnalysis.ElectronTools.calibrationTablesRun2 import correctionType
#@#    from EgammaAnalysis.ElectronTools.calibrationTablesRun2 import files
    
#@#    calibratedPatElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2",
#@#                                            
#@#                                            # input collections
#@#                                            electrons = cms.InputTag('slimmedElectrons'),
#@#                                            gbrForestName = cms.vstring('electron_eb_ECALTRK_lowpt', 'electron_eb_ECALTRK',
#@#                                                                        'electron_ee_ECALTRK_lowpt', 'electron_ee_ECALTRK',
#@#                                                                        'electron_eb_ECALTRK_lowpt_var', 'electron_eb_ECALTRK_var',
#@#                                                                        'electron_ee_ECALTRK_lowpt_var', 'electron_ee_ECALTRK_var'),
#@#                                            
#@#                                            # data or MC corrections
#@#                                            # if isMC is false, data corrections are applied
#@#                                            isMC = cms.bool(False),
#@#                                            autoDataType = cms.bool(True),
#@#                                            
#@#                                            # set to True to get special "fake" smearing for synchronization. Use JUST in case of synchronization
#@#                                            isSynchronization = cms.bool(False),
#@#                                            
#@#                                            correctionFile = cms.string(files[correctionType]),
#@#                                            recHitCollectionEB = cms.InputTag('reducedEgamma:reducedEBRecHits'),
#@#                                            recHitCollectionEE = cms.InputTag('reducedEgamma:reducedEERecHits')
#@#                                            
#@#                                            )
#@#
#@#    
#@#    
#@#    calibratedPatElectrons.electrons = cms.InputTag('hTozzTo4leptonsElectronOrdering')
#@#    calibratedPatElectrons.correctionFile = cms.string(files["Run2017_17Nov2017_v1"]) # reham run2 2017
#@#    calibratedPatElectrons.isMC = cms.bool(True)
    
    # Electron relaxed selection
    from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronSequences_cff import *
    hTozzTo4leptonsElectronSelector=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronSelector_cfi.hTozzTo4leptonsElectronSelector.clone()
   # hTozzTo4leptonsElectronSelector.electronCollection=cms.InputTag("calibratedPatElectrons")
   #@# hTozzTo4leptonsElectronSelector.electronCollection=cms.InputTag("slimmedElectrons::MonoHiggs")# this updated after new 2017 smearing instructions
    hTozzTo4leptonsElectronSelector.electronCollection=cms.InputTag("slimmedElectrons")
    hTozzTo4leptonsElectronSelector.electronEtaMax=cms.double(2.5)
    hTozzTo4leptonsElectronSelector.electronPtMin=cms.double(7.)
    hTozzTo4leptonsElectronSelector.useEleID=cms.bool(False)
    #hTozzTo4leptonsElectronSequence=cms.Sequence(hTozzTo4leptonsElectronIdSequence + hTozzTo4leptonsElectronSelector)
    hTozzTo4leptonsElectronSequence=cms.Sequence(hTozzTo4leptonsElectronSelector)
    
    # Muon Calibration
    #Kalman calibrator
    #from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonCalibrator_cfi import *
    #hTozzTo4leptonsMuonCalibrator = cms.EDProducer("KalmanPATMuonCorrector",
    #                                     src = cms.InputTag("slimmedMuons"),
    #                                     identifier = cms.string(""),
    #                                     isMC = cms.bool(False),
    #                                     isSynchronization = cms.bool(False),
    #                                     )


#    #Rochester calibrator

    myMuonsCalib = cms.EDFilter("CandViewShallowCloneProducer",
                               src = cms.InputTag("slimmedMuons"),
                                cut = cms.string('')
                                )

    goodMuonMCMatchCalib = cms.EDProducer ("MCMatcher",
                            #src     = cms.InputTag("myMuonsCalib"),      # RECO objects to match
                            src     = cms.InputTag("slimmedMuons"),      # RECO objects to match
                            matched = cms.InputTag("prunedGenParticles"), # mc-truth particle collection
                            mcPdgId     = cms.vint32(13),           # one or more PDG ID (13 = muon); absolute values (see below)
                            checkCharge = cms.bool(True),           # True = require RECO and MC objects to have the same charge
                            mcStatus = cms.vint32(1),               # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                            maxDeltaR = cms.double(0.15),           # MAX deltaR for the match
                            maxDPtRel = cms.double(0.5),            # MAX deltaPt/Pt for the match
                            resolveAmbiguities = cms.bool(True),    # Forbid two RECO objects to match to the same GEN object
                            resolveByMatchQuality = cms.bool(True)  # False = just match input in order; True = pick lowest deltaR pair first
                                  )

    from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonRochesterCalibrator_cfi import * #Reham to add Rochester calibrator 

    # Muon ghost cleaning
    from HiggsAnalysis.HiggsToZZ4Leptons.muonCleanerBySegments_cfi import *
    #cleanPatMuonsBySegments.src = cms.InputTag("hTozzTo4leptonsMuonCalibrator") # Reham to add new calibrator
    cleanPatMuonsBySegments.src = cms.InputTag("hTozzTo4leptonsMuonRochesterCalibrator") #Reham Rochester calibrator for 2017 data

    # Muon relaxed selection
    from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonSelector_cfi import *
    hTozzTo4leptonsMuonSelector=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonSelector_cfi.hTozzTo4leptonsMuonSelector.clone()
    hTozzTo4leptonsMuonSelector.muonCollection = cms.InputTag("cleanPatMuonsBySegments")
    #hTozzTo4leptonsMuonSelector.muonCollection = cms.InputTag("hTozzTo4leptonsMuonCalibrator") #Reham to add new calibrator
    #hTozzTo4leptonsMuonSelector.muonCollection = cms.InputTag("hTozzTo4leptonsMuonRochesterCalibrator")
    hTozzTo4leptonsMuonSelector.isGlobalMuon=cms.bool(True) 
    hTozzTo4leptonsMuonSelector.isTrackerMuon=cms.bool(True)  
    hTozzTo4leptonsMuonSelector.muonPtMin=cms.double(5. ) 
    hTozzTo4leptonsMuonSelector.muonEtaMax=cms.double(2.4) 
    hTozzTo4leptonsMuonSequence=cms.Sequence(hTozzTo4leptonsMuonSelector) 

##################################

# Veto electrons and muons for isolation
vetoMuons =  cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("(isGlobalMuon || isTrackerMuon) && pt>1.")
)

vetoElectrons =  cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag("slimmedElectrons"),
    #cut = cms.string("pt>7 && gsfTrack().trackerExpectedHitsInner().numberOfHits<2")
    cut = cms.string("pt>7")                             
)

##@### New MVA Electron ID
##@##from RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi import *
##@###electronMVAValueMapProducer.src = cms.InputTag("hTozzTo4leptonsElectronOrdering")


##@### MVA Electron ID 80X
##@##from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
##@##from RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cff import *
##@###from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_V1_cff import *
##@###from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff import *
##@##from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff import * # reham new ele ID for 2017


# MVA Photon ID 80X
#from RecoEgamma.PhotonIdentification.Identification.mvaTLEID_Fall15_V1_cff import * #reham commented 


## Electron Regression
# from EGamma.EGammaAnalysisTools.electronRegressionEnergyProducer_cfi import *
# eleRegressionEnergy.electronTag    = cms.InputTag("hTozzTo4leptonsElectronSelector")
## E/p combination with new E/P
# from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsRegressionElectronProducer_cfi import *

#need to change CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi

# Electron PF isolation
#@#from CommonTools.ParticleFlow.PFBRECO_cff import *
from CommonTools.ParticleFlow.Isolation.pfElectronIsolationPFBRECO_cff import *
elPFIsoDepositChargedPFBRECO.src    = cms.InputTag("hTozzTo4leptonsElectronSelector")
elPFIsoDepositChargedAllPFBRECO.src = cms.InputTag("hTozzTo4leptonsElectronSelector")
elPFIsoDepositNeutralPFBRECO.src    = cms.InputTag("hTozzTo4leptonsElectronSelector")
elPFIsoDepositGammaPFBRECO.src      = cms.InputTag("hTozzTo4leptonsElectronSelector")
elPFIsoDepositPUPFBRECO.src         = cms.InputTag("hTozzTo4leptonsElectronSelector")
#elPFIsoValueGamma03PFIdPFBRECO.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.)','EcalEndcaps:ConeVeto(0.)')

# Muon PF isolation
# from CommonTools.ParticleFlow.PFBRECO_cff import *
from CommonTools.ParticleFlow.Isolation.pfMuonIsolationPFBRECO_cff import *
muPFIsoDepositChargedPFBRECO.src    = cms.InputTag("hTozzTo4leptonsMuonSelector")
muPFIsoDepositChargedAllPFBRECO.src = cms.InputTag("hTozzTo4leptonsMuonSelector")
muPFIsoDepositNeutralPFBRECO.src    = cms.InputTag("hTozzTo4leptonsMuonSelector")
muPFIsoDepositGammaPFBRECO.src      = cms.InputTag("hTozzTo4leptonsMuonSelector")
muPFIsoDepositPUPFBRECO.src         = cms.InputTag("hTozzTo4leptonsMuonSelector")

# Photon PF
from CommonTools.ParticleFlow.Isolation.pfPhotonIsolationPFBRECO_cff import *
from CommonTools.ParticleFlow.Isolation.photonPFIsolationDepositsPFBRECO_cff import *

phPFIsoDepositChargedPFBRECO.src    = cms.InputTag("hTozzTo4leptonsPFfsrPhoton")
phPFIsoDepositChargedPFBRECO.ExtractorPSet.inputCandView = cms.InputTag("pfAllChargedHadronsPFBRECO")

phPFIsoDepositChargedAllPFBRECO.src = cms.InputTag("hTozzTo4leptonsPFfsrPhoton")
phPFIsoDepositChargedAllPFBRECO.ExtractorPSet.inputCandView = cms.InputTag("pfAllChargedParticlesPFBRECO")

phPFIsoDepositNeutralPFBRECO.src    = cms.InputTag("hTozzTo4leptonsPFfsrPhoton")
phPFIsoDepositNeutralPFBRECO.ExtractorPSet.inputCandView = cms.InputTag("pfAllNeutralHadronsPFBRECO")

phPFIsoDepositGammaPFBRECO.src      = cms.InputTag("hTozzTo4leptonsPFfsrPhoton")
phPFIsoDepositGammaPFBRECO.ExtractorPSet.inputCandView = cms.InputTag("pfAllPhotonsPFBRECO")

phPFIsoDepositPUPFBRECO.src         = cms.InputTag("hTozzTo4leptonsPFfsrPhoton")
phPFIsoDepositPUPFBRECO.ExtractorPSet.inputCandView = cms.InputTag("pfPileUpAllChargedParticlesPFBRECO")

from CommonTools.ParticleFlow.Isolation.photonPFIsolationValuesPFBRECO_cff import *
phPFIsoValueCharged03PFIdPFBRECO.deposits[0].src = cms.InputTag("phPFIsoDepositChargedPFBRECO")
phPFIsoValueCharged03PFIdPFBRECO.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.0001)','EcalEndcaps:ConeVeto(0.0001)','Threshold(0.2)')
phPFIsoValueCharged03PFIdPFBRECO.deposits[0].deltaR = cms.double(0.3)

phPFIsoValueChargedAll03PFIdPFBRECO.deposits[0].src = cms.InputTag("phPFIsoDepositChargedAllPFBRECO")
phPFIsoValueChargedAll03PFIdPFBRECO.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.0001)','EcalEndcaps:ConeVeto(0.0001)','Threshold(0.2)')
phPFIsoValueChargedAll03PFIdPFBRECO.deposits[0].vetos.deltaR = cms.double(0.3)

phPFIsoValueNeutral03PFIdPFBRECO.deposits[0].src = cms.InputTag("phPFIsoDepositNeutralPFBRECO")
phPFIsoValueNeutral03PFIdPFBRECO.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.01)','Threshold(0.5)')
phPFIsoValueNeutral03PFIdPFBRECO.deposits[0].deltaR = cms.double(0.3)

phPFIsoValueGamma03PFIdPFBRECO.deposits[0].src = cms.InputTag("phPFIsoDepositGammaPFBRECO")
phPFIsoValueGamma03PFIdPFBRECO.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.01)','Threshold(0.5)')
phPFIsoValueGamma03PFIdPFBRECO.deposits[0].deltaR = cms.double(0.3)

phPFIsoValuePU03PFIdPFBRECO.deposits[0].src = cms.InputTag("phPFIsoDepositPUPFBRECO")
phPFIsoValuePU03PFIdPFBRECO.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.0001)','EcalEndcaps:ConeVeto(0.0001)','Threshold(0.2)')
phPFIsoValuePU03PFIdPFBRECO.deposits[0].deltaR = cms.double(0.3)


# 3D IP DA

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi import *
hTozzTo4leptonsIpToVtxProducer.VertexLabel = cms.InputTag("offlineSlimmedPrimaryVertices")


# Matching sequence # for matching information Reham 
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMatchingSequence_cff  import *


# COMMON ROOT TREE
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCommonRootTree_cfi  import *
hTozzTo4leptonsCommonRootTreePresel=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCommonRootTree_cfi.hTozzTo4leptonsCommonRootTree.clone()
hTozzTo4leptonsCommonRootTreePresel.decaychannel = cms.string('2e2mu')
hTozzTo4leptonsCommonRootTreePresel.rootFileName = cms.untracked.string('roottree_leptons.root')
# hlt
hTozzTo4leptonsCommonRootTreePresel.fillHLTinfo = cms.untracked.bool(False)                                           
hTozzTo4leptonsCommonRootTreePresel.HLTAnalysisinst = cms.string('hTozzTo4leptonsHLTAnalysisData')


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
hTozzTo4leptonsCommonRootTreePresel.fillMCTruth  = cms.untracked.bool(True) # Reham change to true to add gen match informations
hTozzTo4leptonsCommonRootTreePresel.MCcollName = cms.InputTag("hTozzTo4leptonsMCDumper")
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameBest2e2mu= cms.VInputTag(cms.InputTag("hTozzTo4leptonsBestCandidateProducer:hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hTozzTo4leptonsBestCandidateProducer:hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsBestCandidateProducer:hToZZTo4LeptonsBestCandidateBoson1"))
#hTozzTo4leptonsCommonRootTreePresel.RECOcollName=cms.VInputTag(cms.InputTag("hTozzTo4leptonsLooseIsol"),cms.InputTag("zToMuMuLooseIsol"), cms.InputTag("zToEELooseIsol"))
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameBest4mu= cms.VInputTag(cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateBoson1"))
hTozzTo4leptonsCommonRootTreePresel.RECOcollNameBest4e= cms.VInputTag(cms.InputTag("hTozzTo4leptonsBestCandidateProducerEEEE:hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerEEEE:hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerEEEE:hToZZTo4LeptonsBestCandidateBoson1"))
#hTozzTo4leptonsCommonRootTreePresel.useAdditionalRECO  = cms.untracked.bool(False)


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

#@#hTozzTo4leptonsCommonRootTreePresel.ElectronsEgmLabel        = cms.InputTag("hTozzTo4leptonsElectronSelector")
hTozzTo4leptonsCommonRootTreePresel.ElectronsEgmTkMapLabel   = cms.InputTag("eleIsoFromDepsTkOptimized")
hTozzTo4leptonsCommonRootTreePresel.ElectronsEgmEcalMapLabel = cms.InputTag("eleIsoFromDepsEcalFromHitsByCrystalOptimized")
hTozzTo4leptonsCommonRootTreePresel.ElectronsEgmHcalMapLabel = cms.InputTag("eleIsoFromDepsHcalFromTowersOptimized")        
#@#hTozzTo4leptonsCommonRootTreePresel.ElectronsLabelVert         = cms.InputTag("hTozzTo4leptonsElectronSelector")
hTozzTo4leptonsCommonRootTreePresel.ElectronsMapLabelVert      = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexEleMap")
hTozzTo4leptonsCommonRootTreePresel.ElectronsMapLabelVertValue = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexValueEleMap")
hTozzTo4leptonsCommonRootTreePresel.ElectronsMapLabelVertError = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexErrorEleMap")


# CP variables
hTozzTo4leptonsCommonRootTreePresel.MCCP_PhiLabel          = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsCPPhi")
hTozzTo4leptonsCommonRootTreePresel.MCCP_Phi1Label         = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsCPPhi1")
hTozzTo4leptonsCommonRootTreePresel.MCCP_Phi2Label         = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsCPPhi2")
hTozzTo4leptonsCommonRootTreePresel.MCCP_phi1RFLabel       = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsCPphi1RF")
hTozzTo4leptonsCommonRootTreePresel.MCCP_phi2RFLabel       = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsCPphi2RF")
hTozzTo4leptonsCommonRootTreePresel.MCCP_cosThetaStarLabel = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsCPcosThetaStar")
hTozzTo4leptonsCommonRootTreePresel.MCCP_cosTheta1Label    = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsCPcosTheta1")
hTozzTo4leptonsCommonRootTreePresel.MCCP_cosTheta2Label    = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsCPcosTheta2")
hTozzTo4leptonsCommonRootTreePresel.MCCP_MELALabel         = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsCPMELA")


hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_PhiLabel          = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPPhi")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_Phi1Label         = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPPhi1")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_Phi2Label         = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPPhi2")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_phi1RFLabel       = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPphi1RF")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_phi2RFLabel       = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPphi2RF")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_cosThetaStarLabel = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPcosThetaStar")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_cosTheta1Label    = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPcosTheta1")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_cosTheta2Label    = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPcosTheta2")
hTozzTo4leptonsCommonRootTreePresel.CP2e2mu_MELALabel         = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsCPMELA")

hTozzTo4leptonsCommonRootTreePresel.CP4mu_PhiLabel          = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPPhi")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_Phi1Label         = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPPhi1")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_Phi2Label         = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPPhi2")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_phi1RFLabel       = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPphi1RF")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_phi2RFLabel       = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPphi2RF")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_cosThetaStarLabel = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPcosThetaStar")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_cosTheta1Label    = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPcosTheta1")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_cosTheta2Label    = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPcosTheta2")
hTozzTo4leptonsCommonRootTreePresel.CP4mu_MELALabel         = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPMELA")

hTozzTo4leptonsCommonRootTreePresel.CP4e_PhiLabel          = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPPhi")
hTozzTo4leptonsCommonRootTreePresel.CP4e_Phi1Label         = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPPhi1")
hTozzTo4leptonsCommonRootTreePresel.CP4e_Phi2Label         = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPPhi2")
hTozzTo4leptonsCommonRootTreePresel.CP4e_phi1RFLabel       = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPphi1RF")
hTozzTo4leptonsCommonRootTreePresel.CP4e_phi2RFLabel       = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPphi2RF")
hTozzTo4leptonsCommonRootTreePresel.CP4e_cosThetaStarLabel = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPcosThetaStar")
hTozzTo4leptonsCommonRootTreePresel.CP4e_cosTheta1Label    = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPcosTheta1")
hTozzTo4leptonsCommonRootTreePresel.CP4e_cosTheta2Label    = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPcosTheta2")
hTozzTo4leptonsCommonRootTreePresel.CP4e_MELALabel         = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsCPMELA")

######Conversion
from HiggsAnalysis.HiggsToZZ4Leptons.ConvValueMapProd_cfi  import *

#PFJet ID
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsPFJetSelector_cfi import *
#hTozzTo4leptonsPFJetSelector.PFJetCollection = cms.InputTag("ak4PFJetsCHS")


#//////////////////////////////////////////////////

from RecoJets.JetProducers.QGTagger_cfi import *
QGTaggerMC=QGTagger.clone()
QGTaggerMC.srcJets = cms.InputTag('ak4PFJetsCorrection')
QGTaggerMC.srcJets = cms.InputTag('hTozzTo4leptonsPFJetSelector')     
QGTaggerMC.jetsLabel = cms.string('QGL_AK4PFchs')


QGTaggerDATA=QGTagger.clone()
QGTaggerDATA.srcJets = cms.InputTag('ak4PFJetsCorrectionData')  
QGTaggerDATA.srcJets = cms.InputTag('hTozzTo4leptonsPFJetSelector')  
QGTaggerDATA.jetsLabel = cms.string('QGL_AK4PFchs')

###########################################


hTozzTo4leptonsSelectionSequenceData = cms.Sequence(
    ##test        hTozzTo4leptonsPFtoRECOMuon                 +
    hTozzTo4leptonsPFfsrPhoton                  +
    hTozzTo4leptonsHLTInfo                      +
    hTozzTo4leptonsHLTAnalysisFilter            +
    #@#hTozzTo4leptonsElectronPreSelector       +
    #@#hTozzTo4leptonsElectronOrdering          +
    #calibratedPatElectrons                     +
    #@#hTozzTo4leptonsElectronSelector          +
    #@#electronMVAValueMapProducer              + #Reham stop sequence here to get the output root file from MVA module and see the map name
    #@#hTozzTo4leptonsMuonCalibrator            + #coment kalman calibrator and add Rochester
   #@@@# myMuonsCalib                                + # for rochester
   #@@@# goodMuonMCMatchCalib                        + #for rochester
    hTozzTo4leptonsMuonRochesterCalibrator      +   #for rochester
    cleanPatMuonsBySegments                     +  
    hTozzTo4leptonsMuonSelector                 +  
    # PF isolation for electrons and muons
    ##test   pfParticleSelectionPFBRECOSequence  + 
    ##test   pfElectronIsolationPFBRECOSequence  +      
    ##test   muonPFIsolationPFBRECOSequence      +
    ##test   pfPhotonIsolationPFBRECOSequence    +
    hTozzTo4leptonsPFJetSelector                 +                             
    QGTaggerMC                                   + 
    QGTaggerDATA                                      
    )


hTozzTo4leptonsSelectionSequenceMC = cms.Sequence(
    ##test        hTozzTo4leptonsPFtoRECOMuon                 +
    hTozzTo4leptonsPFfsrPhoton                  +
    hTozzTo4leptonsHLTInfo                      +
    hTozzTo4leptonsHLTAnalysisFilter            +
    #@#hTozzTo4leptonsElectronPreSelector       +
    #@#hTozzTo4leptonsElectronOrdering          +
    #calibratedPatElectrons                     +
    #@#hTozzTo4leptonsElectronSelector          +
    #@#electronMVAValueMapProducer              + #Reham stop sequence here to get the output root file from MVA module and see the map name
    #@#hTozzTo4leptonsMuonCalibrator            + #coment kalman calibrator and add Rochester
    myMuonsCalib                                + # for rochester
    goodMuonMCMatchCalib                        + #for rochester
    hTozzTo4leptonsMuonRochesterCalibrator      +   #for rochester
    cleanPatMuonsBySegments                     +  
    hTozzTo4leptonsMuonSelector                 +  
    # PF isolation for electrons and muons
    ##test   pfParticleSelectionPFBRECOSequence  + 
    ##test   pfElectronIsolationPFBRECOSequence  +      
    ##test   muonPFIsolationPFBRECOSequence      +
    ##test   pfPhotonIsolationPFBRECOSequence    +
    hTozzTo4leptonsPFJetSelector                 +                             
    QGTaggerMC                                   + 
    QGTaggerDATA                                      
    )



