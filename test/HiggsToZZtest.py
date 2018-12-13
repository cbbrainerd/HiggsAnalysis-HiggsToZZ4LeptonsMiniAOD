
import FWCore.ParameterSet.Config as cms

process = cms.Process('TestTest')


process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

###############
# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/Geometry/GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff') #reham
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration/EventContent/EventContent_cff')


from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v7', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v13', '')#Reham Tag recommended for JEC 2017
###########


process.printTree1 = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("prunedGenParticles"),
    maxEventsToPrint  = cms.untracked.int32(2)
)

###################################
### Configure electron IDs with VID
###################################

#from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

#switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

#my_id_modules = [ 
#        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
#                ]

#for idmod in my_id_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


###########################################
# Configure Energy correction for electron 
##########################################

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       eleIDModules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                       runVID=True, #saves CPU time by not needlessly re-running VID
                       era='2017-Nov17ReReco')

################

#process.hTozzTo4leptonsPath = cms.Path(process.egammaPostRecoSeq * process.egmGsfElectronIDSequence  * process.printTree1)
process.hTozzTo4leptonsPath = cms.Path( process.egammaPostRecoSeq *  process.printTree1)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('run.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('higgsToZZtest')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('hTozzTo4leptonsPath')
    )                               
 )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
'/store/mc/RunIIFall17MiniAOD/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/40000/205E2EB6-2600-E811-A8D9-A0369FC5E090.root' #2017 Synchronization
                             )
                           )


# Endpath
process.o = cms.EndPath ( process.output )


