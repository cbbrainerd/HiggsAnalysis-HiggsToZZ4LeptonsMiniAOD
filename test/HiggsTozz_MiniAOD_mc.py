import FWCore.ParameterSet.Config as cms

process = cms.Process('MonoHiggs')

# Complete Preselection Sequence for 4l analysis

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/Geometry/GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')#reham
process.load('Configuration.StandardSequences.MagneticField_cff') #reham
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration/EventContent/EventContent_cff')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v7', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v13', '')#Reham Tag recommended for JEC 2017

# Random generator
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
    )
)

process.load('HiggsAnalysis.HiggsToZZ4Leptons.bunchSpacingProducer_cfi')
#process.load('HiggsAnalysis.HiggsToZZ4Leptons.metFiltersMiniAOD_cff')

process.load('RecoMET.METFilters.metFilters_cff')
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')

process.Path_BunchSpacingproducer=cms.Path(process.bunchSpacingProducer)

process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseIsoFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)                                                          
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.primaryVertexFilter.vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices')
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.BadPFMuonFilter.muons  = cms.InputTag("slimmedMuons")
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.muons  = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates") #reham
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter) # Reham added for 2017
#process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)## 

#///////////////////////////////
#new MET filter 2017 Reham & 2018  To be used (UNDER TEST)

#baddetEcallist2017 = cms.vuint32(
#    [872439604,872422825,872420274,872423218,
#     872423215,872416066,872435036,872439336,
#     872420273,872436907,872420147,872439731,
#     872436657,872420397,872439732,872439339,
#     872439603,872422436,872439861,872437051,
#     872437052,872420649])

#process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
#    "EcalBadCalibFilter",
#    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
#    ecalMinEt        = cms.double(50.),
#    baddetEcal    = baddetEcallist2017, #use baddetEcallist2018  for 2018 analysis
#    taggingMode = cms.bool(True),
#    debug = cms.bool(False)
#    )


#/////////////////////////////////////////            

process.goodOfflinePrimaryVerticestwo = cms.EDFilter("VertexSelector",
                                            src = cms.InputTag('offlineSlimmedPrimaryVertices'),
					    cut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
                                            filter = cms.bool(True)
                                        )
        

#@#Muon Rochester correction Reham

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsMuonRochesterCalibrator_cfi')
process.hTozzTo4leptonsMuonRochesterCalibrator.isData = cms.bool(False)
process.hTozzTo4leptonsMuonRochesterCalibrator.MCTruth = cms.bool(True)

#Kalman calibrator

#@#process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsMuonCalibrator_cfi')
#@#process.hTozzTo4leptonsMuonCalibrator.isData = cms.bool(False) 

#@#process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
#@#process.calibratedElectrons.isMC = cms.bool(True)

#/////////////////////////////////////////////////////////////
#Reham JEC

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(

process,
jetSource = cms.InputTag('slimmedJets'),
labelName = 'UpdatedJEC',
jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')

)

process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

#///////////////////////////////////////////////////////////

#Reham to update the MET after updating the JEC

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD(process,
                           isData=False, #(or False),
                           postfix = "TEST"
                           )

#/////////////////////////////////////////////////////

#Reham to add new instructiond for electron energy correction and smearing PLUS electron ID 

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=True,
                       eleIDModules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
                                       'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                      # runVID=True, #saves CPU time by not needlessly re-running VID
                      era='2017-Nov17ReReco')  


#/////////////////////////////////////////////////////

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim_cff')  

#@#process.calibratedPatElectrons.isMC = cms.bool(True)#Reham Run2 2017

process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.use2011EA = cms.untracked.bool(False)
process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = True
process.hTozzTo4leptonsCommonRootTreePresel.fillHLTinfo = cms.untracked.bool(False)                                           
process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f10QL3Filtered20Q')
process.hTozzTo4leptonsCommonRootTreePresel.triggerEleFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f10QL3Filtered20Q')
process.hTozzTo4leptonsCommonRootTreePresel.fillMCTruth  = cms.untracked.bool(True)    
process.hTozzTo4leptonsCommonRootTreePresel.isVBF  = cms.bool(False)
#//@
#This variable isData to apply muon calibrator inside commonRooTree.h and get the error on muon pT
process.hTozzTo4leptonsCommonRootTreePresel.isData = cms.bool(False)

process.hTozzTo4leptonsCommonRootTreePresel.noiseFilterTag = cms.InputTag("TriggerResults","","PAT")

#for LHE informations for Jets
#process.hTozzTo4leptonsCommonRootTreePresel.LHEProduct = cms.InputTag("source") #when commented, default=externalLHEProducer
process.hTozzTo4leptonsCommonRootTreePresel.LHEProduct = cms.InputTag("externalLHEProducer")# this inputTag depend on input mc sample 

process.genanalysis= cms.Sequence(
  process.hTozzTo4leptonsGenSequence                  *
  #       process.hTozzTo4leptonsMCGenFilter2e2mu             *
  #       process.hTozzTo4leptonsMCGenParticleListDrawer2e2mu *
  process.hTozzTo4leptonsMCDumper                     
 # process.hTozzTo4leptonsMCCP                         
  )


process.hTozzTo4leptonsSelectionPath = cms.Path(
 # process.ecalBadCalibReducedMINIAODFilter * #new met filter to be used
  process.goodOfflinePrimaryVerticestwo     *
  process.genanalysis *
  process.jecSequence* #Reham to add JEC
  process.fullPatMetSequenceTEST * #Reham To update MET after update JEC
  process.egammaPostRecoSeq * #Reham to include electron smearing due to kink at 50 Gev in electron pt spectrum from old electron scale and smearing 
  #process.hTozzTo4leptonsSelectionSequenceData * 
  process.hTozzTo4leptonsSelectionSequenceMC * 
  process.hTozzTo4leptonsMatchingSequence * # for MC matching Reham
  process.hTozzTo4leptonsCommonRootTreePresel   
  )

#//////////////////////////////////////////////////////////////////

#process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff') #reham comment to run in crab
#from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *   #reham need to comment in run in crab
#process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()  #reham need to comment in run in crab
#process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "MC_Synch_2017hTozzToLeptons.root"  #reham need to comment in run in crab

#process.o = cms.EndPath (process.hTozzTo4leptonsSelectionOutputModuleNew ) #reham comment in run in crab 
process.schedule = cms.Schedule( process.Path_BunchSpacingproducer,
                                 #@#process.Flag_HBHENoiseFilter,
                                 #@#process.Flag_HBHENoiseIsoFilter,
                                 #@#process.Flag_globalSuperTightHalo2016Filter,
                                 #@#process.Flag_EcalDeadCellTriggerPrimitiveFilter,
                                 #@#process.Flag_goodVertices,
                                 #process.Flag_eeBadScFilter,####Data only 
                                 #@#process.Flag_BadPFMuonFilter, ### reham uncomment
                                 #@#process.Flag_BadChargedCandidateFilter,#### reham uncomment
                                 #process.Flag_ecalBadCalibFilter, #new 2017 changed with Reduced
                                 #process.ecalBadCalibReducedMINIAODFilter,
                                 process.hTozzTo4leptonsSelectionPath
                                 #process.o #comment to run in crab
                                 )


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source ("PoolSource",
                             
  fileNames = cms.untracked.vstring(
'/store/mc/RunIIFall17MiniAOD/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/40000/205E2EB6-2600-E811-A8D9-A0369FC5E090.root', #2017 Synchronization
#'/store/mc/RunIIFall17MiniAOD/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v2/00000/E8505BB6-5F07-E811-B009-002590DE6E88.root',#2017 Synchronization
#'/store/mc/RunIIFall17MiniAOD/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/10000/80B92986-8501-E811-99BB-002590200900.root'#2017 Synchronization
#'file:RunIIFall17MiniAOD_VBF_HToZZTo4L_M125_13TeV_E8505BB6-5F07-E811-B009-002590DE6E88.root'
 )#,
#skipEvents= cms.untracked.uint32(18000)
)

## # Endpath
#process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew ) 
