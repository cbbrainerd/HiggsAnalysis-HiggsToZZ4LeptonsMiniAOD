import FWCore.ParameterSet.Config as cms

from RecoMuon.MuonIsolationProducers.muIsoDepositTk_cfi import *
import RecoMuon.MuonIsolationProducers.muIsoDepositTk_cfi

# Old way to compute tracker isolation with veto
muIsoDepositTkNew=RecoMuon.MuonIsolationProducers.muIsoDepositTk_cfi.muIsoDepositTk.clone()
muIsoDepositTkNew.IOPSet.inputMuonCollection = cms.InputTag("hTozzTo4leptonsMuonSelector")


from RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi  import *
import RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi 
muIsoDepositCalByAssociatorTowersNew=RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi.muIsoDepositCalByAssociatorTowers.clone()
muIsoDepositCalByAssociatorTowersNew.IOPSet.inputMuonCollection = cms.InputTag("hTozzTo4leptonsMuonSelector")

from RecoMuon.MuonIsolationProducers.muIsoDepositJets_cfi  import *
import RecoMuon.MuonIsolationProducers.muIsoDepositJets_cfi 
muIsoDepositJetsNew=RecoMuon.MuonIsolationProducers.muIsoDepositJets_cfi.muIsoDepositJets.clone()
muIsoDepositJetsNew.IOPSet.inputMuonCollection = cms.InputTag("hTozzTo4leptonsMuonSelector")

muIsoDeposits_muonsNew = cms.Sequence(muIsoDepositTkNew+muIsoDepositCalByAssociatorTowersNew+muIsoDepositJetsNew)
muIsolation_muonsNew = cms.Sequence(muIsoDeposits_muonsNew)
muIsolationNew = cms.Sequence(muIsolation_muonsNew)


# from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonIsolationProducer_cfi import *

# hTozzTo4leptonsMuonIsolationSequence = cms.Sequence( muIsolationNew + hTozzTo4leptonsMuonIsolationProducer )



muIsoFromDepsTkOptimized = cms.EDProducer("CandIsolatorFromDeposits",
   deposits = cms.VPSet(cms.PSet(
       mode = cms.string('sum'),
       src = cms.InputTag("muIsoDepositTkNew"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring(
              'vetoMuons:0.015',
              'vetoElectrons:0.015',
              'Threshold(1.0)'),
       skipDefaultVeto = cms.bool(True))
   )
)


muIsoFromDepsEcalOptimized = cms.EDProducer("CandIsolatorFromDeposits",
   deposits = cms.VPSet(cms.PSet(
       mode = cms.string('sum'),
       src = cms.InputTag("muIsoDepositCalByAssociatorTowersNew","ecal"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('vetoMuons:0.015',
                           'vetoElectrons:0.015',
                           'Threshold(1.0)'),
       skipDefaultVeto = cms.bool(True)
   ))
)

muIsoFromDepsHcalOptimized = cms.EDProducer("CandIsolatorFromDeposits",
   deposits = cms.VPSet(cms.PSet(
       mode = cms.string('sum'),
       src = cms.InputTag("muIsoDepositCalByAssociatorTowersNew","hcal"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('vetoMuons:0.015',
                           'vetoElectrons:0.015',
                           'Threshold(1.0)'),
       skipDefaultVeto = cms.bool(True)
   ))
)

hTozzTo4leptonsMuonIsolationSequence = cms.Sequence(
    muIsoDepositTkNew+muIsoDepositCalByAssociatorTowersNew+muIsoDepositJetsNew +
    muIsoFromDepsTkOptimized   +
    muIsoFromDepsEcalOptimized +
    muIsoFromDepsHcalOptimized
)

