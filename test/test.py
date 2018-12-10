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
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v6', '') # Reham Tag recommended for JEC 2017

############################################################

process.load('HiggsAnalysis.HiggsToZZ4Leptons.bunchSpacingProducer_cfi')
#process.load('HiggsAnalysis.HiggsToZZ4Leptons.metFiltersMiniAOD_cff')

process.load('RecoMET.METFilters.metFilters_cff')
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')

################################################################

process.Path_BunchSpacingproducer=cms.Path(process.bunchSpacingProducer)

process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseIsoFilter)
## process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)                                                               
## process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)                                           
## process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)      
## process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)   
## process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)                                                             
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
## process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter) 
process.primaryVertexFilter.vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices')
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)

################################################
## process.Flag_trackingFailureFilter = cms.Path(process.goodVertices + process.trackingFailureFilter)                                    
#process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.BadPFMuonFilter.muons  = cms.InputTag("slimmedMuons")
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.BadChargedCandidateFilter.muons  = cms.InputTag("slimmedMuons")
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)

#///////////////////////////////
#new MET filter 2017 Reham

process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter) #new 2017

################################################

process.goodOfflinePrimaryVerticestwo = cms.EDFilter("VertexSelector",
                                            src = cms.InputTag('offlineSlimmedPrimaryVertices'),
					    cut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
                                            filter = cms.bool(True)
                                        )
        


process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsMuonCalibrator_cfi')
process.hTozzTo4leptonsMuonCalibrator.isData = cms.bool(True) 
# process.hTozzTo4leptonsMuonCalibrator.iisMC = cms.bool(False)

####################################################

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(

process,
jetSource = cms.InputTag('slimmedJets'),
labelName = 'UpdatedJEC',
jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')

)

process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

################################################

##update MET after update jet 

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD(process,
                           isData=True, #(or False),
                           )

#/////////////////////////////////////////////////////

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim2_cff') 

process.calibratedPatElectrons.isMC = cms.bool(False) #reham run2 2017


process.hTozzTo4leptonsPFfsrPhoton.src = cms.InputTag("packedPFCandidates")
process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.use2011EA = cms.untracked.bool(False)
process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = False
process.hTozzTo4leptonsCommonRootTreePresel.fillHLTinfo = cms.untracked.bool(False)                                           
process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f10QL3Filtered20Q')
process.hTozzTo4leptonsCommonRootTreePresel.triggerEleFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f10QL3Filtered20Q')
  #process.hTozzTo4leptonsCommonRootTreePresel.triggerFilterAsym = cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8')
process.hTozzTo4leptonsCommonRootTreePresel.fillMCTruth  = cms.untracked.bool(False)    
process.hTozzTo4leptonsCommonRootTreePresel.isVBF  = cms.bool(False)
#//@
#This variable isData to apply muon calibrator inside commonRooTree.h and get the error on muon pT
process.hTozzTo4leptonsCommonRootTreePresel.isData = cms.bool(True)
#for MC only but put need to run and not crash
process.hTozzTo4leptonsCommonRootTreePresel.LHEProduct = cms.InputTag("externalLHEProducer")# this inputTag depend on input mc sample 

process.genanalysis= cms.Sequence(
  # process.hTozzTo4leptonsGenSequence                  *
  #       process.hTozzTo4leptonsMCGenFilter2e2mu             *
  #       process.hTozzTo4leptonsMCGenParticleListDrawer2e2mu *
  process.hTozzTo4leptonsMCDumper                     
 # process.hTozzTo4leptonsMCCP                         
  )

process.hTozzTo4leptonsSelectionPath = cms.Path(
  process.goodOfflinePrimaryVerticestwo     *
#  process.genanalysis * 
  process.jecSequence *#Reham to add JEC
 # process.fullPatMetSequenceTEST * #Reham To update MET after update JEC
 # process.testCands *
  process.fullPatMetSequence * #Reham To update MET after update JEC
  process.hTozzTo4leptonsSelectionSequenceData *# Reham to add again
  #process.pileupJetIdUpdated * #REham to add puileup jet ID
  #process.patJetCorrFactorsReapplyJEC * #Reham to add JEC for PileUP
  #process.updatedJets * #Reham to add JEC for PileUP
#  process.hTozzTo4leptonsMatchingSequence *
  process.hTozzTo4leptonsCommonRootTreePresel 
  )

#///////////////////////////////////////////////////////
#Reham comment to test
#quark/gluon tagging

#process.load("CondCore.CondDB.CondDB_cfi")
#qgDatabaseVersion = '80X'
#process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
#                                        DBParameters = cms.PSet(messageLevel = cms.untracked.int32(1)),
#                                        timetype = cms.string('runnumber'),
#                                        toGet = cms.VPSet(
#                                          cms.PSet(
#                                             record = cms.string('QGLikelihoodRcd'),
#                                             tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_AK4PFchs'),
#                                             label  = cms.untracked.string('QGL_AK4PFchs')
#                                             ),
#                                          ),
#                                          connect = cms.string('sqlite:QGL_'+qgDatabaseVersion+'.db')
#)
#process.es_prefer_qg = cms.ESPrefer('PoolDBESSource','QGPoolDBESSource')

#///////////////////////////////////////////

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *   #reham need to comment in run in crab
process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()  #reham need to comment in run in crab
process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "Data_2017_DoubleMuon_RunB_hTozzToLeptons.root"  #reham need to comment in run in crab

process.o = cms.EndPath (process.hTozzTo4leptonsSelectionOutputModuleNew ) #reham comment in run in crab
process.schedule = cms.Schedule( process.Path_BunchSpacingproducer,                            
                                 process.Flag_HBHENoiseFilter,
                                 process.Flag_HBHENoiseIsoFilter,
                                 process.Flag_globalTightHalo2016Filter,
                                 process.Flag_EcalDeadCellTriggerPrimitiveFilter,
                                 process.Flag_goodVertices,
#                                 process.Flag_eeBadScFilter,
#                                 process.Flag_BadPFMuonFilter,
#                                 process.Flag_BadChargedCandidateFilter,
                                  process.Flag_ecalBadCalibFilter,  #new 2017 
                                  process.hTozzTo4leptonsSelectionPath,
                                 process.o)
                                 


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source ("PoolSource",
                             
  fileNames = cms.untracked.vstring(
'file:data_DoubleMuon_2017_RunB_0852E0CB-E7D7-E711-B2DA-0025905C3DCE.root' 
#'file:Data_2017_DoubleMuon_RunB_hTozzToLeptons.root'
  )
)

## # Endpath
#process.o = cms.EndPath (process.hTozzTo4leptonsSelectionOutputModuleNew ) #reham comment in run in crab

 

####################################################


# process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
# process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# process.source = cms.Source ("PoolSource",
                             
#   fileNames = cms.untracked.vstring(
# 'file:data_DoubleMuon_2017_RunB_0852E0CB-E7D7-E711-B2DA-0025905C3DCE.root' 
#   )
# )

