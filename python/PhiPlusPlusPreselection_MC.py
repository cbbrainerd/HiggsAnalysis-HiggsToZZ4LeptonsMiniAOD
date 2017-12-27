import FWCore.ParameterSet.Config as cms

process = cms.Process('PhiPlusPlusPresel')

# Complete Preselection Sequence for Phi++ analysis

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

process.load("Configuration/StandardSequences/Reconstruction_cff")

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load("HiggsAnalysis.HiggsToZZ4Leptons.simpleEleIdSequence_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsPFIsolationProducer_cff")


process.GlobalTag.globaltag = "START42_V13::All"

# Preselection analysis sequence
process.load('HiggsAnalysis/HiggsToZZ4Leptons/phiplusplusPreselection_cff')

process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = True
process.higgsToZZ4LeptonsSkimFilterData.useHLT  = cms.untracked.bool(False)  
  
process.higgsToZZ4LeptonsHLTAnalysisData.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hTozzTo4leptonsHLTAnalysisData.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

process.patTrigger.processName=cms.string("HLT")
process.hTozzTo4leptonsHLTInfo.TriggerResultsTag=cms.InputTag("TriggerResults","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT")

process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltDiMuonL3PreFiltered7')
process.hTozzTo4leptonsCommonRootTreePresel.triggerFilterAsym = cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8')
  

process.hTozzTo4leptonsSelectionPath = cms.Path(process.hTozzTo4leptonsMuonSelector + process.simpleEleIdSequence + process.PFTau + process.phiplusplusSelectionSequence)

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *
process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()
process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "phiplusplus_MC.root"

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
'file:/lustre/cms/store/mc/Summer11/HPlusPlusHMinusMinusHTo4L_M-350_7TeV-pythia6/AODSIM/PU_S4_START42_V11-v1/0000/187C9AD7-BAA2-E011-B419-003048D460FC.root'
#'file:/lustre/cms/store/mc/Summer11/HPlusPlusHMinusHTo3L_M-250_7TeV-calchep-pythia6/AODSIM/PU_S4_START42_V11-v1/0000/368FF930-BDC5-E011-86D8-00215E221224.root'
                             )
                           )


# Endpath
# process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )

