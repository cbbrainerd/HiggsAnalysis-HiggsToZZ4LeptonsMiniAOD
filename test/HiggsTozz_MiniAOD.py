import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options=VarParsing('analysis')
options.register('dataset','',VarParsing.multiplicity.singleton,VarParsing.varType.string,'Dataset to process (e.g. DoubleMuon)')
options.register('year',-1,VarParsing.multiplicity.singleton,VarParsing.varType.int,'Year to run over (e.g. 2017)')
options.parseArguments()
dataset=options.dataset
isMC=(dataset=="MC")
AllowedDatasets=['DoubleMuon','SingleElectron','DoubleEG','MuonEG','SingleMuon','MC']
year=options.year

if dataset=='' or year==-1:
    print 'Must specify dataset and year to run over, e.g. cmsRun config.py dataset=DoubleMuon year=2018'
    raise SystemExit
elif dataset not in AllowedDatasets:
    print 'Unrecognized dataset "%s". Implemented datasets are:' % dataset
    for dataset in AllowedDatasets:
        print dataset
    raise SystemExit
else:
    print 'Running over dataset %s' % dataset

process = cms.Process('MonoHiggs')

# Complete Preselection Sequence for 4l analysis

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')#reham
process.load('Configuration.StandardSequences.MagneticField_cff') #reham
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.EventContent.EventContent_cff')


from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
GlobalTagByYearData={ 2017 : '94X_dataRun2_v11' , 2018 : '102X_dataRun2_v4' }
GlobalTagByYearMC={ 2017 : '94X_mc2017_realistic_v13' , 2018 : '102X_upgrade2018_realistic_v18' }
GlobalTagByYear=GlobalTagByYearMC if isMC else GlobalTagByYearData
process.GlobalTag = GlobalTag(process.GlobalTag, GlobalTagByYear[year], '') # Reham Tag recommended for JEC 2017

if isMC:
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

###ADDEDfor MET Walaa##
process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')
baddetEcallist = cms.vuint32(
                             [872439604,872422825,872420274,872423218,872423215,872416066,872435036,872439336, 872420273,872436907,872420147,872439731,872436657,872420397,872439732,872439339, 872439603,872422436,872439861,872437051,872437052,872420649,872421950,872437185, 872422564,872421566,872421695,872421955,872421567,872437184,872421951,872421694, 872437056,872437057,872437313,872438182,872438951,872439990,872439864,872439609, 872437181,872437182,872437053,872436794,872436667,872436536,872421541,872421413, 872421414,872421031,872423083,872421439])

process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
                                                        "EcalBadCalibFilter",
                                                        EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
                                                        ecalMinEt        = cms.double(50.),
                                                        baddetEcal    = baddetEcallist,
                                                        taggingMode = cms.bool(True),
                                                        debug = cms.bool(False)
                                                        )

process.Path_BunchSpacingproducer=cms.Path(process.bunchSpacingProducer)

process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseIsoFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)                                                         
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.primaryVertexFilter.vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices')
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)                                    
#process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter) ##for data only reham
process.BadPFMuonFilter.muons  = cms.InputTag("slimmedMuons")
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.muons  = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates") #reham
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter) # Reham added for 2017
#process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter) #new 2017

#///////////////////////////////
#new MET filter 2017 Reham new MET filter to be used (under test)

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
        
#@#Rochester correction

process.load('HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonRochesterCalibrator_cfi')
process.hTozzTo4leptonsMuonRochesterCalibrator.isData = cms.bool(not isMC)
process.hTozzTo4leptonsMuonRochesterCalibrator.MCTruth = cms.bool(isMC)

###############include Jet Walaa
import os
# Jet Energy Corrections
from CondCore.DBCommon.CondDBSetup_cfi import *
era = "Autumn18_V8_MC" if isMC else "Autumn18_RunABCD_V8_DATA"
#dBFile = os.environ.get('CMSSW_BASE')+"/src/HiggsAnalysis/HiggsToZZ4Leptons/test/"+era+".db"
dBFile = "Autumn18_V8_MC.db" if isMC else "Autumn18_RunABCD_V8_DATA.db"
process.jec = cms.ESSource("PoolDBESSource",
                           CondDBSetup,
                           connect = cms.string("sqlite_file:"+dBFile),
                           toGet =  cms.VPSet(
                                              
                                              cms.PSet(
                                                       record = cms.string("JetCorrectionsRecord"),
                                                       tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
                                                       label= cms.untracked.string("AK4PFchs")
                                                       ),
                                              )
                           )
process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')
process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")

process.jetCorrFactors = process.updatedPatJetCorrFactors.clone(
                                                                src = cms.InputTag("slimmedJets"),
                                                                levels = ['L1FastJet',
                                                                          'L2Relative',
                                                                          'L3Absolute'],
                                                                payload = 'AK4PFchs' )
process.slimmedJetsJEC = process.updatedPatJets.clone(
                                                      jetSource = cms.InputTag("slimmedJets"),
                                                      jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jetCorrFactors"))
                                                      )
### add pileup id and discriminant to patJetsReapplyJEC
process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdUpdated = process.pileupJetId.clone(
                                                       jets=cms.InputTag("slimmedJets"),
                                                       inputIsCorrected=False,
                                                       applyJec=True,
                                                       vertexes=cms.InputTag("offlineSlimmedPrimaryVertices")
                                                       )
process.slimmedJetsJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
process.slimmedJetsJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']
#######################################

#JER
if isMC:
    process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
    #dBJERFile = os.environ.get('CMSSW_BASE')+"/src/HiggsAnalysis/HiggsToZZ4Leptons/test/Autumn18_V1_MC.db"
    dBJERFile = "Autumn18_V1_MC.db"
    process.jer = cms.ESSource("PoolDBESSource",
                               CondDBSetup,
                               connect = cms.string("sqlite_file:"+dBJERFile),
                               toGet = cms.VPSet(
                                                 cms.PSet(
                                                          record = cms.string('JetResolutionRcd'),
                                                          tag    = cms.string('JR_Autumn18_V1_MC_PtResolution_AK4PFchs'),
                                                          label  = cms.untracked.string('AK4PFchs_pt')
                                                          ),
                                                 cms.PSet(
                                                          record = cms.string('JetResolutionRcd'),
                                                          tag    = cms.string('JR_Autumn18_V1_MC_PhiResolution_AK4PFchs'),
                                                          label  = cms.untracked.string('AK4PFchs_phi')
                                                          ),
                                                 cms.PSet(
                                                          record = cms.string('JetResolutionScaleFactorRcd'),
                                                          tag    = cms.string('JR_Autumn18_V1_MC_SF_AK4PFchs'),
                                                          label  = cms.untracked.string('AK4PFchs')
                                                          )
                                                 )
                               )
    process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD(process,
                           isData=True, #(or False),
                           postfix = "TEST"
                           )

#/////////////////////////////////////////////////////

#Reham to add new instructiond for electron energy correction and smearing PLUS electron ID 

postRecoSeqEras= { 2017 : '2017-Nov17ReReco' , 2018 : '2018-Prompt' }
postRecoSeqEleIDModules = { 
    2017 : [ 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff', 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff' ],
    2018 : [ 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Autumn18_ID_ISO_cff', 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff' ]
}

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       eleIDModules = postRecoSeqEleIDModules[year],
                       runVID=True, #saves CPU time by not needlessly re-running VID
                       runEnergyCorrections=True,
                      era=postRecoSeqEras[year])  


#/////////////////////////////////////////////////////


process.load('HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsPreselection_data_noskim_cff')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHLTAnalysisFilter_cfi import set_triggers
set_triggers(year=year,dataset=dataset)

#@#process.calibratedPatElectrons.isMC = cms.bool(False) #reham run2 2017
process.load('HiggsAnalysis.HiggsToZZ4Leptons.fsrPhotons_cff') #FSR Walaa

process.hTozzTo4leptonsPFfsrPhoton.src = cms.InputTag("packedPFCandidates")
process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.use2011EA = cms.untracked.bool(False)
process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT")
process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = isMC
process.hTozzTo4leptonsCommonRootTreePresel.fillHLTinfo = cms.untracked.bool(isMC)                                           
process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f10QL3Filtered20Q')
process.hTozzTo4leptonsCommonRootTreePresel.triggerEleFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f10QL3Filtered20Q')
  #process.hTozzTo4leptonsCommonRootTreePresel.triggerFilterAsym = cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8')
process.hTozzTo4leptonsCommonRootTreePresel.fillMCTruth  = cms.untracked.bool(False)    
process.hTozzTo4leptonsCommonRootTreePresel.year = cms.untracked.int32(year)
process.hTozzTo4leptonsCommonRootTreePresel.isVBF  = cms.bool(False)
#//@
#This variable isData to apply muon calibrator inside commonRooTree.h and get the error on muon pT
process.hTozzTo4leptonsCommonRootTreePresel.isData = cms.bool(True)

process.hTozzTo4leptonsCommonRootTreePresel.noiseFilterTag = cms.InputTag("TriggerResults","","PAT" if isMC else "RECO")

#for MC only but put need to run and not crash
process.hTozzTo4leptonsCommonRootTreePresel.LHEProduct = cms.InputTag("externalLHEProducer")# this inputTag depend on input mc sample 

if isMC:
    process.genanalysis= cms.Sequence(
      process.hTozzTo4leptonsGenSequence                  *
      #       process.hTozzTo4leptonsMCGenFilter2e2mu             *
      #       process.hTozzTo4leptonsMCGenParticleListDrawer2e2mu *
      process.hTozzTo4leptonsMCDumper
     # process.hTozzTo4leptonsMCCP                         
      )
else:
    process.genanalysis= cms.Sequence(
      process.hTozzTo4leptonsMCDumper
    )

process.hTozzTo4leptonsSelectionPath = cms.Path(
   # process.ecalBadCalibReducedMINIAODFilter  * New met filter to be used (under test)
    process.goodOfflinePrimaryVerticestwo *
    process.genanalysis *
    process.fsrPhotonSequence *
    process.jetCorrFactors *
    process.pileupJetIdUpdated *
    process.ecalBadCalibReducedMINIAODFilter *
    process.slimmedJetsJEC *
    process.fullPatMetSequenceTEST *
    process.egmGsfElectronIDSequence *
    process.egammaPostRecoSeq *
    (process.hTozzTo4leptonsSelectionSequenceMC*process.hTozzTo4leptonsMatchingSequence if isMC else process.hTozzTo4leptonsSelectionSequenceData)*
    process.hTozzTo4leptonsCommonRootTreePresel
    #process.jecSequence *#Reham to add JEC
    #process.fullPatMetSequenceTEST * #Reham To update MET after update JEC
    #process.egammaPostRecoSeq * #Reham to include electron smearing due to kink at 50 Gev in electron pt spectrum from old electron scale and smearing
    #process.hTozzTo4leptonsSelectionSequenceData *# Reham to add again
    #process.hTozzTo4leptonsCommonRootTreePresel 
    )

#///////////////////////////////////////////

process.load('HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *   #reham need to comment in run in crab
#process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()  #reham need to comment in run in crab
#process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "Data_2017_DoubleMuon_RunB_hTozzToLeptons.root"  #reham need to comment in run in crab

#process.o = cms.EndPath (process.hTozzTo4leptonsSelectionOutputModuleNew ) #reham comment in run in crab
process.schedule = cms.Schedule( process.Path_BunchSpacingproducer,                            
                                 #@#process.Flag_HBHENoiseFilter,
                                 #@#process.Flag_HBHENoiseIsoFilter,
                                 #@#process.Flag_globalSuperTightHalo2016Filter,
                                 #@#process.Flag_EcalDeadCellTriggerPrimitiveFilter,
                                 #@#process.Flag_goodVertices,
                                # process.Flag_eeBadScFilter,### data only
                                 #@#process.Flag_BadPFMuonFilter,#####
                                 #@#process.Flag_BadChargedCandidateFilter,###
                                 #process.Flag_ecalBadCalibFilter,  #new 2017 changed with reduced
                                 process.hTozzTo4leptonsSelectionPath,
                                 #process.o
)



#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck', #DEBUG
#    ignoreTotal=cms.untracked.int32(1), #DEBUG
#    showMallocInfo=cms.untracked.bool(True), #DEBUG
#    moduleMemorySummary=cms.untracked.bool(True) #DEBUG
#) #DEBUG

process.source = cms.Source ("PoolSource",
                             
  fileNames = cms.untracked.vstring(
#'file:data_DoubleMuon_2017_RunB_0852E0CB-E7D7-E711-B2DA-0025905C3DCE.root' 
#'root://cmsxrootd.fnal.gov//store/data/Run2017B/DoubleMuon/MINIAOD/17Nov2017-v1/30000/0852E0CB-E7D7-E711-B2DA-0025905C3DCE.root'
#'/store/data/Run2018C/DoubleMuon/MINIAOD/PromptReco-v3/000/320/065/00000/DC1C9ACB-2B90-E811-BD7F-FA163E635E53.root' #2018 data
#'file:DC1C9ACB-2B90-E811-BD7F-FA163E635E53.root',
'file:E8250D82-7AEE-A245-8BA2-DAC402BFF393.root', #2018MC
#'/store/data/Run2017C/DoubleMuon/MINIAOD/31Mar2018-v1/90000/047E2618-7738-E811-B77A-38EAA78D8AF4.root' #2017 data (for sync)
#'file:Data_2017_DoubleMuon_RunB_hTozzToLeptons.root'
  )
)
#from PhysicsTools.PythonAnalysis.LumiList import LumiList
#myLumis=LumiList(filename='goodList.json').getCMSSWString().split(',')
#process.source.lumisToProcess=cms.untracked.VLuminosityBlockRange()
#process.source.lumisToProcess.extend(myLumis)

#from IgTools.IgProf.IgProfTrigger import igprof
#process.load("IgTools.IgProf.IgProfTrigger")
#process.igprofPath = cms.Path(process.igprof)
#process.igprof.reportEventInterval     = cms.untracked.int32(25)
#process.igprof.reportToFileAtBeginJob  = cms.untracked.string("|gzip -c>igprof.begin-job.gz")
#process.igprof.reportToFileAtEvent = cms.untracked.string("|gzip -c>igprof.%I.%E.%L.%R.event.gz")
#process.schedule.append(process.igprofPath)
#
## # Endpath
#process.o = cms.EndPath (process.hTozzTo4leptonsSelectionOutputModuleNew ) #reham comment in run in crab
