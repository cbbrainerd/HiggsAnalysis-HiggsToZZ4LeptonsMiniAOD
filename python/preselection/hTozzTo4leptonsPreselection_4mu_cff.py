import FWCore.ParameterSet.Config as cms

## 
# Filter to select 4mu events
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi
hTozzTo4leptonsMCGenFilter4mu = HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi.hTozzTo4leptonsMCGenFilter.clone()
hTozzTo4leptonsMCGenFilter4mu.HZZ4LeptonsMCFilterLeptonFlavour=cms.int32(1)

# ParticleListDrawer
from  HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenParticleListDrawer_cfi import * 
hTozzTo4leptonsMCGenParticleListDrawer4mu = HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenParticleListDrawer_cfi.hTozzTo4leptonsMCGenParticleListDrawer.clone()

# Save MC truth
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCDumper_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCDumper_cfi
hTozzTo4leptonsMCDumper4mu = HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCDumper_cfi.hTozzTo4leptonsMCDumper.clone()
hTozzTo4leptonsMCDumper4mu.motherPdgId          = cms.int32(25)
hTozzTo4leptonsMCDumper4mu.firstdaughtersPdgId  = cms.vint32( 23, 23)
hTozzTo4leptonsMCDumper4mu.seconddaughtersPdgId = cms.vint32( 13, -13, 13, -13 )

# Build best candidate from MC truth
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsBestCandidateProducer_cfi import *
hTozzTo4leptonsMCBestCandidateProducer4mu=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsBestCandidateProducer_cfi.hTozzTo4leptonsBestCandidateProducer.clone()
hTozzTo4leptonsMCBestCandidateProducer4mu.decaychannel = cms.string('4mu')
hTozzTo4leptonsMCBestCandidateProducer4mu.RECOcollName = cms.VInputTag(cms.InputTag("hTozzTo4leptonsMCDumper4mu:hToZZTo4LeptonsMCDumperCompositeMother"))
hTozzTo4leptonsMCBestCandidateProducer4mu.decayChain = cms.string('hToZZTo4LeptonsMCBestCandidate')

# CP producer
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi import *
hTozzTo4leptonsMCCP4mu=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi.hTozzTo4leptonsCP.clone()
hTozzTo4leptonsMCCP4mu.decayChain = cms.string('hToZZTo4LeptonsMCBestCandidateCP')
hTozzTo4leptonsMCCP4mu.RECOcollNameLabFrame = cms.VInputTag(cms.InputTag("hTozzTo4leptonsMCBestCandidateProducer4mu:hToZZTo4LeptonsMCBestCandidateMother"), cms.InputTag("hTozzTo4leptonsMCBestCandidateProducer4mu:hToZZTo4LeptonsMCBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsMCBestCandidateProducer4mu:hToZZTo4LeptonsMCBestCandidateBoson1"))


# Filter of low mass peak events
# from HiggsAnalysis.HiggsToZZ4Leptons.hTozzT4LeptonsLowMassPeakFilter_cfi import *

# module to apply skim
from HiggsAnalysis.Skimming.higgsToZZ4Leptons_Sequences_cff import *
higgsToZZ4LeptonsHLTAnalysisData=HiggsAnalysis.Skimming.higgsToZZ4LeptonsHLTAnalysis_cfi.higgsToZZ4LeptonsHLTAnalysis.clone()
#higgsToZZ4LeptonsHLTAnalysisData.HLTPaths = cms.vstring('HLT_L1Mu14_L1ETM30','HLT_L1Mu14_L1SingleJet6U','HLT_L1Mu14_L1SingleEG10','HLT_L1Mu20','HLT_DoubleMu3','HLT_Mu3','HLT_Mu5','HLT_Mu9','HLT_L2Mu9','HLT_L2Mu11','HLT_L1Mu30','HLT_Mu7','HLT_L2Mu15')
higgsToZZ4LeptonsHLTAnalysisData.HLTPaths= cms.vstring('HLT_DoubleMu3','HLT_Mu9','HLT_DoubleEle10_SW_L1R','HLT_Ele10_SW_L1R')
higgsToZZ4LeptonsHLTAnalysisData.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
higgsToZZ4LeptonsSkimFilterData=HiggsAnalysis.Skimming.higgsToZZ4LeptonsSkimFilter_cfi.higgsToZZ4LeptonsSkimFilter.clone()
higgsToZZ4LeptonsSkimFilterData.useHLT  = cms.untracked.bool(False) 
higgsToZZ4LeptonsSkimFilterData.HLTinst = cms.string('higgsToZZ4LeptonsHLTAnalysisData')
higgsToZZ4LeptonsSkimFilterData.useDiLeptonSkim = cms.untracked.bool(False)


higgsToZZ4LeptonsSequenceData = cms.Sequence(higgsToZZ4LeptonsHLTAnalysisData           +
                                             higgsToZZ4LeptonsBuildLeptons              +
                                             higgsToZZ4LeptonsSkimDiLeptonProducer      +
                                             higgsToZZ4LeptonsSkimTriLeptonProducer     +
                                             higgsToZZ4LeptonsSkimFilterData )


# module to run HLT analysis
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHLTAnalysis_cfi import *
hTozzTo4leptonsHLTAnalysisData=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHLTAnalysis_cfi.hTozzTo4leptonsHLTAnalysis.clone()
hTozzTo4leptonsHLTAnalysisData.HLTPaths= cms.vstring('HLT_DoubleMu3','HLT_Mu9','HLT_DoubleEle10_SW_L1R','HLT_Ele17_SW_EleId_L1R')
hTozzTo4leptonsHLTAnalysisData.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

# Muon selection
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonSelector_cfi import *

# zToMuMu
from HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi import *

# hTozzToMMMM
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi
hTozzTo4leptonsMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi.hTozzTo4leptons.clone()
hTozzTo4leptonsMMMM.decay = cms.string('zToMuMu zToMuMu')

# Muon loose isolation
# from Configuration.StandardSequences.GeometryIdeal_cff import *
from Configuration.StandardSequences.GeometryDB_cff import *
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonIsolationSequences_cff import *


# Common preselection 
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCommonPreselectionSequences_cff import *
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCommonPreselectionSequences_cff
hTozzTo4leptonsCommonPreselectionMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCommonPreselectionSequences_cff.hTozzTo4leptonsCommonPreselection.clone()
hTozzTo4leptonsCommonPreselectionMMMM.decaychannel = cms.string('4mu')
hTozzTo4leptonsCommonPreselectionMMMM.MuonsLabel = cms.InputTag("hTozzTo4leptonsMuonSelector")
hTozzTo4leptonsCommonPreselectionMMMM.LeptonsLabel = cms.InputTag("allLeptons")
hTozzTo4leptonsCommonPreselectionMMMM.HLabel = cms.InputTag("hTozzTo4leptonsMMMM")
hTozzTo4leptonsCommonPreselectionMMMM.ZMuMuLabel = cms.InputTag("zToMuMu")
hTozzTo4leptonsCommonPreselectionMMMM.MuonsLooseIsolLabel = cms.InputTag("hTozzTo4leptonsMuonIsolationProducer")
hTozzTo4leptonsCommonPreselectionMMMM.cuts.nMu = cms.int32(4)
hTozzTo4leptonsCommonPreselectionMMMM.cuts.mumuMass = cms.double(12.0)
hTozzTo4leptonsCommonPreselectionMMMM.cuts.numberOfmumuCombs = cms.int32(2)
hTozzTo4leptonsCommonPreselectionMMMM.cuts.fourleptMass = cms.double(100.0)
hTozzTo4leptonsCommonPreselectionMMMM.cuts.numberOf4lCombs = cms.int32(1)
hTozzTo4leptonsCommonPreselectionMMMM.cuts.nlooseMu = cms.int32(4)

import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCommonPreselectionSequences_cff
hTozzTo4leptonsCommonPreselectionFilterMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCommonPreselectionSequences_cff.hTozzTo4leptonsCommonPreselectionFilter.clone()
hTozzTo4leptonsCommonPreselectionFilterMMMM.decaychannel = cms.string('4mu')
hTozzTo4leptonsCommonPreselectionFilterMMMM.preselinst = cms.string('hTozzTo4leptonsCommonPreselectionMMMM')
hTozzTo4leptonsCommonPreselectionFilterMMMM.preseltags = cms.vstring('PreselAtleast4Mu','PreselAtleast2ZMuMu','PreselAtleast1H','PreselLoose4IsolMu') 
hTozzTo4leptonsCommonPreselectionFilterMMMM.rootFileName = cms.string('preselect4mu.root')
hTozzTo4leptonsCommonPreselectionFilterMMMM.preSelectFileName = cms.string('preselect4mu.out')

hTozzTo4leptonsCommonPreselectionSequenceMMMM = cms.Sequence(hTozzTo4leptonsCommonPreselectionMMMM+hTozzTo4leptonsCommonPreselectionFilterMMMM)

# zToMuMu loose isolated
from HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi
zToMuMuLooseIsolMMMM = HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi.zToMuMu.clone()
zToMuMuLooseIsolMMMM.decay = ('hTozzTo4leptonsMuonIsolationProducer@+ hTozzTo4leptonsMuonIsolationProducer@-')
# hTozzToMuMuMuMu loose isolated
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi
hTozzTo4leptonsLooseIsolMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi.hTozzTo4leptons.clone()
hTozzTo4leptonsLooseIsolMMMM.decay = ('zToMuMuLooseIsolMMMM zToMuMuLooseIsolMMMM')

# best candidate producer
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsBestCandidateProducer_cfi import *
hTozzTo4leptonsBestCandidateProducerMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsBestCandidateProducer_cfi.hTozzTo4leptonsBestCandidateProducer.clone()
hTozzTo4leptonsBestCandidateProducerMMMM.decaychannel = cms.string('4mu')
hTozzTo4leptonsBestCandidateProducerMMMM.RECOcollName = cms.VInputTag(cms.InputTag("hTozzTo4leptonsLooseIsolMMMM"))
hTozzTo4leptonsBestCandidateProducerMMMM.decayChain = cms.string('hToZZTo4LeptonsBestCandidate')


# higgs frame producer
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHiggsFrame_cfi import *
hTozzTo4leptonsHiggsFrameMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHiggsFrame_cfi.hTozzTo4leptonsHiggsFrame.clone()
hTozzTo4leptonsHiggsFrameMMMM.prodinst = cms.string('hTozzTo4leptonsBestCandidateProducerMMMM')
hTozzTo4leptonsHiggsFrameMMMM.RECOcollName = cms.VInputTag(cms.InputTag("hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hToZZTo4LeptonsBestCandidateBoson1"))

# CP producer
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi import *
hTozzTo4leptonsCPMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi.hTozzTo4leptonsCP.clone()
hTozzTo4leptonsCPMMMM.RECOcollNameLabFrame = cms.VInputTag(cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateBoson1"))

# 3D IP
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi import *
hTozzTo4leptonsIpToVtxProducerMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsIpToVtxProducer_cfi.hTozzTo4leptonsIpToVtxProducer.clone()
hTozzTo4leptonsIpToVtxProducerMMMM.decaychannel = cms.string('4mu')
hTozzTo4leptonsIpToVtxProducerMMMM.BestCandidatesLeptons = cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateLeptons")
hTozzTo4leptonsIpToVtxProducerMMMM.MuonsLabel=cms.InputTag("hTozzTo4leptonsMuonSelector")
hTozzTo4leptonsIpToVtxProducerMMMM.useBeamSpot = cms.bool(False)

# 2D IP
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsTipLipToVtxProducer_cfi import *
hTozzTo4leptonsTipLipToVtxProducerMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsTipLipToVtxProducer_cfi.hTozzTo4leptonsTipLipToVtxProducer.clone()
hTozzTo4leptonsTipLipToVtxProducerMMMM.decaychannel = cms.string('4mu')
hTozzTo4leptonsTipLipToVtxProducerMMMM.BestCandidatesLeptons = cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateLeptons")
hTozzTo4leptonsTipLipToVtxProducerMMMM.MuonsLabel=cms.InputTag("hTozzTo4leptonsMuonSelector")
hTozzTo4leptonsTipLipToVtxProducerMMMM.useBeamSpot = cms.bool(False)

# Geometrical Discriminant
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsGeomDiscrimProducer_cfi import *
hTozzTo4leptonsGeomDiscrimProducerMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsGeomDiscrimProducer_cfi.hTozzTo4leptonsGeomDiscrimProducer.clone()
hTozzTo4leptonsGeomDiscrimProducerMMMM.decaychannel = cms.string('4mu')
hTozzTo4leptonsGeomDiscrimProducerMMMM.BeamSpotLabel= cms.InputTag("offlineBeamSpot")
hTozzTo4leptonsGeomDiscrimProducerMMMM.RECOcollName=cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateMother")

# Constrained fit: input 4mu
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsConstraintFitProducer_cfi import *
hTozzTo4leptonsConstraintFitProducerMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsConstraintFitProducer_cfi.hTozzTo4leptonsConstraintFitProducer.clone()
hTozzTo4leptonsConstraintFitProducerMMMM.RECOcollName =cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateMother")

        
# Commont Root tree
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCommonRootTree_cfi  import *
hTozzTo4leptonsCommonRootTreePreselMMMM=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCommonRootTree_cfi.hTozzTo4leptonsCommonRootTree.clone()
hTozzTo4leptonsCommonRootTreePreselMMMM.decaychannel = cms.string('4mu')
hTozzTo4leptonsCommonRootTreePreselMMMM.rootFileName = cms.untracked.string('roottree_4mu.root')
# hlt
#hTozzTo4leptonsCommonRootTreePreselMMMM.HLTAnalysisinst = cms.string('hTozzTo4leptonsHLTAnalysis')
#hTozzTo4leptonsCommonRootTreePreselMMMM.flagHLTnames=cms.VInputTag(cms.InputTag("flagHLTIsoMu11"),cms.InputTag("flagHLTMu15"),cms.InputTag("flagHLTDoubleMu3"),cms.InputTag("flagHLTIsoEle15L1I"),cms.InputTag("flagHLTIsoEle18L1R"),cms.InputTag("flagHLTDoubleIsoEle10L1I"), cms.InputTag("flagHLTDoubleIsoEle12L1R"), cms.InputTag("flagHLTaccept"))
hTozzTo4leptonsCommonRootTreePreselMMMM.HLTAnalysisinst = cms.string('hTozzTo4leptonsHLTAnalysisData')
#hTozzTo4leptonsCommonRootTreePresel.flagHLTnames=cms.VInputTag(cms.InputTag("flagHLTIsoMu11"),cms.InputTag("flagHLTMu15"),cms.InputTag("flagHLTDoubleMu3"),cms.InputTag("flagHLTIsoEle15L1I"),cms.InputTag("flagHLTIsoEle18L1R"),cms.InputTag("flagHLTDoubleIsoEle10L1I"), cms.InputTag("flagHLTDoubleIsoEle12L1R"), cms.InputTag("flagHLTaccept"))
hTozzTo4leptonsCommonRootTreePreselMMMM.flagHLTnames=cms.VInputTag(cms.InputTag("flagHLTDoubleMu3"),cms.InputTag("flagHLTMu9"),cms.InputTag("flagHLTDoubleEle5SWL1R"),cms.InputTag("flagHLTEle15LWL1R"),cms.InputTag("flagHLTaccept"))

# skim early data
hTozzTo4leptonsCommonRootTreePreselMMMM.useSkimEarlyData = cms.untracked.bool(False)
# presel
hTozzTo4leptonsCommonRootTreePreselMMMM.flaginst = cms.string('hTozzTo4leptonsCommonPreselectionMMMM')
hTozzTo4leptonsCommonRootTreePreselMMMM.flagtags = cms.vstring('PreselAtleast4Mu','PreselAtleast2ZMuMu','PreselAtleast1H','PreselLoose4IsolMu')
hTozzTo4leptonsCommonRootTreePreselMMMM.fillMCTruth  = cms.untracked.bool(False)
#hTozzTo4leptonsCommonRootTreePreselMMMM.MCcollName = cms.VInputTag(cms.InputTag("hTozzTo4leptonsMCDumper4mu:hToZZTo4LeptonsMCDumperMother"), cms.InputTag("hTozzTo4leptonsMCDumper4mu:hToZZTo4LeptonsMCDumperBoson0"), cms.InputTag("hTozzTo4leptonsMCDumper4mu:hToZZTo4LeptonsMCDumperBoson1"), cms.InputTag("hTozzTo4leptonsMCDumper4mu:hToZZTo4LeptonsMCDumperLepton0"), cms.InputTag("hTozzTo4leptonsMCDumper4mu:hToZZTo4LeptonsMCDumperLepton1"), cms.InputTag("hTozzTo4leptonsMCDumper4mu:hToZZTo4LeptonsMCDumperLepton2"), cms.InputTag("hTozzTo4leptonsMCDumper4mu:hToZZTo4LeptonsMCDumperLepton3"))
hTozzTo4leptonsCommonRootTreePreselMMMM.MCcollName = cms.VInputTag(cms.InputTag("hTozzTo4leptonsMCBestCandidateProducer4mu:hToZZTo4LeptonsMCBestCandidateMother"), cms.InputTag("hTozzTo4leptonsMCBestCandidateProducer4mu:hToZZTo4LeptonsMCBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsMCBestCandidateProducer4mu:hToZZTo4LeptonsMCBestCandidateBoson1"))
hTozzTo4leptonsCommonRootTreePreselMMMM.RECOcollNameBest4mu= cms.VInputTag(cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateMother"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateBoson0"), cms.InputTag("hTozzTo4leptonsBestCandidateProducerMMMM:hToZZTo4LeptonsBestCandidateBoson1"))
hTozzTo4leptonsCommonRootTreePreselMMMM.RECOcollNameBestRestFrame4mu = cms.VInputTag(cms.InputTag("hTozzTo4leptonsHiggsFrameMMMM:hToZZTo4LeptonsHiggsRestFrameMother"), cms.InputTag("hTozzTo4leptonsHiggsFrameMMMM:hToZZTo4LeptonsHiggsRestFrameBoson0"), cms.InputTag("hTozzTo4leptonsHiggsFrameMMMM:hToZZTo4LeptonsHiggsRestFrameBoson1"))
hTozzTo4leptonsCommonRootTreePreselMMMM.useAdditionalRECO  = cms.untracked.bool(False)
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsLabel     = cms.InputTag("hTozzTo4leptonsMuonSelector")
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsMapLabel  = cms.InputTag("hTozzTo4leptonsMuonIsolationProducer")
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsTkMapLabel = cms.InputTag("hTozzTo4leptonsMuonIsolationProducer:Tk")
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsEcalMapLabel = cms.InputTag("hTozzTo4leptonsMuonIsolationProducer:Ecal")
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsHcalMapLabel = cms.InputTag("hTozzTo4leptonsMuonIsolationProducer:Hcal")
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsLabelVert = cms.InputTag("hTozzTo4leptonsMuonSelector")
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsMapLabelVert = cms.InputTag("hTozzTo4leptonsIpToVtxProducerMMMM:VertexMuMap")
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsMapLabelVertValue    = cms.InputTag("hTozzTo4leptonsIpToVtxProducerMMMM:VertexValueMuMap")
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsMapLabelVertError    = cms.InputTag("hTozzTo4leptonsIpToVtxProducerMMMM:VertexErrorMuMap")
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsSTIPMapLabelVert     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducerMMMM:TipMuMap")
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsSLIPMapLabelVert     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducerMMMM:LipMuMap")
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsSTIPMapLabelVertValue = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducerMMMM:TipValueMuMap")
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsSLIPMapLabelVertValue = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducerMMMM:LipValueMuMap")
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsSTIPMapLabelVertError = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducerMMMM:TipErrorMuMap")
hTozzTo4leptonsCommonRootTreePreselMMMM.MuonsSLIPMapLabelVertError = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducerMMMM:LipErrorMuMap")
hTozzTo4leptonsCommonRootTreePreselMMMM.ftsigmaVertMMMM = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:ftsigma")
hTozzTo4leptonsCommonRootTreePreselMMMM.ftsigmalagVertMMMM = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:ftsigmalag")

#CP variables
hTozzTo4leptonsCommonRootTreePreselMMMM.MCCP_4mu_PhiLabel   = cms.InputTag("hTozzTo4leptonsMCCP4mu:hToZZTo4LeptonsMCBestCandidateCPPhi")
hTozzTo4leptonsCommonRootTreePreselMMMM.MCCP_4mu_Phi1Label   = cms.InputTag("hTozzTo4leptonsMCCP4mu:hToZZTo4LeptonsMCBestCandidateCPPhi1")
hTozzTo4leptonsCommonRootTreePreselMMMM.MCCP_4mu_Phi2Label   = cms.InputTag("hTozzTo4leptonsMCCP4mu:hToZZTo4LeptonsMCBestCandidateCPPhi2")
hTozzTo4leptonsCommonRootTreePreselMMMM.MCCP_4mu_phi1RFLabel   = cms.InputTag("hTozzTo4leptonsMCCP4mu:hToZZTo4LeptonsMCBestCandidateCPphi1RF")
hTozzTo4leptonsCommonRootTreePreselMMMM.MCCP_4mu_phi2RFLabel   = cms.InputTag("hTozzTo4leptonsMCCP4mu:hToZZTo4LeptonsMCBestCandidateCPphi2RF")
hTozzTo4leptonsCommonRootTreePreselMMMM.MCCP_4mu_cosThetaStarLabel = cms.InputTag("hTozzTo4leptonsMCCP4mu:hToZZTo4LeptonsMCBestCandidateCPcosThetaStar")
hTozzTo4leptonsCommonRootTreePreselMMMM.MCCP_4mu_cosTheta1Label = cms.InputTag("hTozzTo4leptonsMCCP4mu:hToZZTo4LeptonsMCBestCandidateCPcosTheta1")
hTozzTo4leptonsCommonRootTreePreselMMMM.MCCP_4mu_cosTheta2Label = cms.InputTag("hTozzTo4leptonsMCCP4mu:hToZZTo4LeptonsMCBestCandidateCPcosTheta2")

hTozzTo4leptonsCommonRootTreePreselMMMM.CP4mu_PhiLabel   = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPPhi")
hTozzTo4leptonsCommonRootTreePreselMMMM.CP4mu_Phi1Label   = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPPhi1")
hTozzTo4leptonsCommonRootTreePreselMMMM.CP4mu_Phi2Label   = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPPhi2")
hTozzTo4leptonsCommonRootTreePreselMMMM.CP4mu_phi1RFLabel   = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPphi1RF")
hTozzTo4leptonsCommonRootTreePreselMMMM.CP4mu_phi2RFLabel   = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPphi2RF")
hTozzTo4leptonsCommonRootTreePreselMMMM.CP4mu_cosThetaStarLabel = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPcosThetaStar")
hTozzTo4leptonsCommonRootTreePreselMMMM.CP4mu_cosTheta1Label = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPcosTheta1")
hTozzTo4leptonsCommonRootTreePreselMMMM.CP4mu_cosTheta2Label = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsCPcosTheta2")


hTozzTo4leptonsSelectionSequence4mu = cms.Sequence(
        hTozzTo4leptonsMCGenFilter4mu                  +
#        hTozzTo4leptonsMCGenParticleListDrawer4mu      +
#        hTozzTo4leptonsMCDumper4mu                     +
#        hTozzTo4leptonsMCBestCandidateProducer4mu      +
#        hTozzTo4leptonsMCCP4mu                         +
        higgsToZZ4LeptonsSequenceData                  +
        hTozzTo4leptonsHLTAnalysisData                 + 
	hTozzTo4leptonsMuonSelector                    + 
        zToMuMu                                        +
        hTozzTo4leptonsMMMM                            +
	hTozzTo4leptonsMuonIsolationSequence           + 
        hTozzTo4leptonsCommonPreselectionSequenceMMMM  +
        zToMuMuLooseIsolMMMM                           +
        hTozzTo4leptonsLooseIsolMMMM                   +
        hTozzTo4leptonsBestCandidateProducerMMMM       +
        hTozzTo4leptonsHiggsFrameMMMM                  +
        hTozzTo4leptonsCPMMMM                          +
        hTozzTo4leptonsIpToVtxProducerMMMM             +
        hTozzTo4leptonsTipLipToVtxProducerMMMM         +
 	hTozzTo4leptonsGeomDiscrimProducerMMMM         +
        hTozzTo4leptonsConstraintFitProducerMMMM       +
        hTozzTo4leptonsCommonRootTreePreselMMMM        
	)
                                                 

