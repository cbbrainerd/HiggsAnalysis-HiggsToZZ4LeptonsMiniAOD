import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.pfParticleSelection_cff import *
from CommonTools.ParticleFlow.Isolation.tools_cfi import * 

#Put all charged particles in charged hadron collection(electrons and muons)
pfAllChargedHadrons.pdgId = cms.vint32(211,-211,321,-321,999211,2212,-2212,11,-11,13,-13)

#For DB correction#######
pfPileUpCandidates = cms.EDProducer("TPPFCandidatesOnPFCandidates",
                                    enable =  cms.bool(True),
                                    verbose = cms.untracked.bool(False),
                                    name = cms.untracked.string("pileUpCandidates"),
                                    topCollection = cms.InputTag("pfNoPileUp"),
                                    bottomCollection = cms.InputTag("particleFlow")
                                    )

#enable PF no Pile Up
pfPileUp.Enable = cms.bool(True)

pileUpHadrons = cms.EDFilter("PdgIdPFCandidateSelector",
                             src = cms.InputTag("pfPileUpCandidates"),
                             pdgId = cms.vint32(211,-211,321,-321,999211,2212,-2212,11,-11,13,-13)
                             )
#########################

pfPostSequence = cms.Sequence(
    pfParticleSelectionSequence+
    pfPileUpCandidates+
    pileUpHadrons
    )


###Electron Isolation
electronPFIsoDepositAll  = isoDeposits.clone()
electronPFIsoDepositAll.src = cms.InputTag("hTozzTo4leptonsElectronSelector")
electronPFIsoDepositAll.ExtractorPSet.inputCandView = cms.InputTag("pfNoPileUp")

electronPFIsoDepositCharged  = isoDeposits.clone()
electronPFIsoDepositCharged.src = cms.InputTag("hTozzTo4leptonsElectronSelector")
electronPFIsoDepositCharged.ExtractorPSet.inputCandView = cms.InputTag("pfAllChargedHadrons")

electronPFIsoDepositNeutral  = isoDeposits.clone()
electronPFIsoDepositNeutral.src = cms.InputTag("hTozzTo4leptonsElectronSelector")
electronPFIsoDepositNeutral.ExtractorPSet.inputCandView = cms.InputTag("pfAllNeutralHadrons")

electronPFIsoDepositGamma  = isoDeposits.clone()
electronPFIsoDepositGamma.src = cms.InputTag("hTozzTo4leptonsElectronSelector")
electronPFIsoDepositGamma.ExtractorPSet.inputCandView = cms.InputTag("pfAllPhotons")

#Isodeposit from PileUp- For Vertex subtraction!!!!
electronPFIsoDepositPU  = isoDeposits.clone()
electronPFIsoDepositPU.src = cms.InputTag("hTozzTo4leptonsElectronSelector")
electronPFIsoDepositPU.ExtractorPSet.inputCandView = cms.InputTag("pileUpHadrons")



electronPFIsoDeposits = cms.Sequence(
    electronPFIsoDepositAll*
    electronPFIsoDepositCharged*
    electronPFIsoDepositPU*
    electronPFIsoDepositNeutral*
    electronPFIsoDepositGamma
    )


#And Values
electronPFIsoValueAll = cms.EDProducer("CandIsolatorFromDeposits",
                                 deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("electronPFIsoDepositAll"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring('0.001','Threshold(0.5)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
                                 )

electronPFIsoValueCharged = cms.EDProducer("CandIsolatorFromDeposits",
                                     deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("electronPFIsoDepositCharged"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring('0.0001','Threshold(0.0)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
                                     )

electronPFIsoValueNeutral = cms.EDProducer("CandIsolatorFromDeposits",
                                     deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("electronPFIsoDepositNeutral"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring('0.01','Threshold(0.5)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
                                     )

electronPFIsoValueGamma = cms.EDProducer("CandIsolatorFromDeposits",
                                   deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("electronPFIsoDepositGamma"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring('0.01','Threshold(0.5)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
                                   )

electronPFIsoValuePU = cms.EDProducer("CandIsolatorFromDeposits",
                                deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("electronPFIsoDepositPU"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring('0.0001','Threshold(0.5)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
                                )

electronPFIsoValuePULow = cms.EDProducer("CandIsolatorFromDeposits",
                                   deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("electronPFIsoDepositPU"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring('0.0001','Threshold(0.0)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
   )
    )
                                   )



electronPFIsoValues =  cms.Sequence( electronPFIsoValueAll
                               * electronPFIsoValueCharged
                               * electronPFIsoValueNeutral
                               * electronPFIsoValueGamma
                               * electronPFIsoValuePU
                               * electronPFIsoValuePULow
                               )





###Muon Isolation
muPFIsoDepositAll  = isoDeposits.clone()
muPFIsoDepositAll.src = cms.InputTag("hTozzTo4leptonsMuonSelector")
muPFIsoDepositAll.ExtractorPSet.inputCandView = cms.InputTag("pfNoPileUp")

muPFIsoDepositCharged  = isoDeposits.clone()
muPFIsoDepositCharged.src = cms.InputTag("hTozzTo4leptonsMuonSelector")
muPFIsoDepositCharged.ExtractorPSet.inputCandView = cms.InputTag("pfAllChargedHadrons")

muPFIsoDepositNeutral  = isoDeposits.clone()
muPFIsoDepositNeutral.src = cms.InputTag("hTozzTo4leptonsMuonSelector")
muPFIsoDepositNeutral.ExtractorPSet.inputCandView = cms.InputTag("pfAllNeutralHadrons")

muPFIsoDepositGamma  = isoDeposits.clone()
muPFIsoDepositGamma.src = cms.InputTag("hTozzTo4leptonsMuonSelector")
muPFIsoDepositGamma.ExtractorPSet.inputCandView = cms.InputTag("pfAllPhotons")

#Isodeposit from PileUp- For Vertex subtraction!!!!
muPFIsoDepositPU  = isoDeposits.clone()
muPFIsoDepositPU.src = cms.InputTag("hTozzTo4leptonsMuonSelector")
muPFIsoDepositPU.ExtractorPSet.inputCandView = cms.InputTag("pileUpHadrons")



muPFIsoDeposits = cms.Sequence(
    muPFIsoDepositAll*
    muPFIsoDepositCharged*
    muPFIsoDepositPU*
    muPFIsoDepositNeutral*
    muPFIsoDepositGamma
    )


#And Values
muPFIsoValueAll = cms.EDProducer("CandIsolatorFromDeposits",
                                 deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("muPFIsoDepositAll"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring('0.001','Threshold(0.5)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
                                 )

muPFIsoValueCharged = cms.EDProducer("CandIsolatorFromDeposits",
                                     deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("muPFIsoDepositCharged"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring('0.0001','Threshold(0.0)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
                                     )

muPFIsoValueNeutral = cms.EDProducer("CandIsolatorFromDeposits",
                                     deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("muPFIsoDepositNeutral"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring('0.01','Threshold(0.5)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
                                     )

muPFIsoValueGamma = cms.EDProducer("CandIsolatorFromDeposits",
                                   deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("muPFIsoDepositGamma"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring('0.01','Threshold(0.5)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
                                   )

muPFIsoValuePU = cms.EDProducer("CandIsolatorFromDeposits",
                                deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("muPFIsoDepositPU"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring('0.0001','Threshold(0.5)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
    )
    )
                                )

muPFIsoValuePULow = cms.EDProducer("CandIsolatorFromDeposits",
                                   deposits = cms.VPSet(
    cms.PSet(
    src = cms.InputTag("muPFIsoDepositPU"),
    deltaR = cms.double(0.4),
    weight = cms.string('1'),
    vetos = cms.vstring('0.0001','Threshold(0.0)'),
    skipDefaultVeto = cms.bool(True),
    mode = cms.string('sum')
   )
    )
                                   )



muPFIsoValues =  cms.Sequence( muPFIsoValueAll
                               * muPFIsoValueCharged
                               * muPFIsoValueNeutral
                               * muPFIsoValueGamma
                               * muPFIsoValuePU
                               * muPFIsoValuePULow
                               )
