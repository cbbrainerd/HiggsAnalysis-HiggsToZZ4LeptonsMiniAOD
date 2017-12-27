import FWCore.ParameterSet.Config as cms

from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

hTozzTo4leptonsTipLipToVtxProducer = cms.EDProducer("HZZ4LeptonsTipLipToVtxProducer",
    # BestCandidatesLeptons = cms.InputTag("hTozzTo4leptonsBestCandidateProducer:hToZZTo4LeptonsBestCandidateLeptons"),
#    ElectronsLabel = cms.InputTag("hTozzTo4leptonsElectronIsolationProducerEgamma"),
    ElectronsLabel = cms.InputTag("hTozzTo4leptonsElectronSelector"),
#    MuonsLabel = cms.InputTag("hTozzTo4leptonsMuonIsolationProducer"),
    MuonsLabel = cms.InputTag("hTozzTo4leptonsMuonSelector"),
    VertexLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
    BeamSpotLabel  = cms.InputTag("offlineBeamSpot"),
    useBeamSpot = cms.bool(False),
    decaychannel = cms.string('2e2mu')
)


