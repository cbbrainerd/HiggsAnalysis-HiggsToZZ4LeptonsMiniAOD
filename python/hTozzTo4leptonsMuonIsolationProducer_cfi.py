import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsMuonIsolationProducer = cms.EDProducer("HZZ4LeptonsMuonIsolationProducer",

    # primary vertex
    PVLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),

    # deposits
    ECALIsoDepositLabel = cms.InputTag("muIsoDepositCalByAssociatorTowersNew","ecal"),
    HCALIsoDepositLabel = cms.InputTag("muIsoDepositCalByAssociatorTowersNew","hcal"),
    HOCALIsoDepositLabel = cms.InputTag("muIsoDepositCalByAssociatorTowersNew","ho"),
    TrackerIsoDepositLabel = cms.InputTag("muIsoDepositTkNew"),

    # objects
    MuonsLabel = cms.InputTag("hTozzTo4leptonsMuonSelector"),
    TracksLabel = cms.InputTag("generalTracks"),
    ElectronsLabel = cms.InputTag("gsfElectrons"),

    # algo parameters
    isolationCone = cms.double(0.3),
    isolationConeVeto = cms.double(0.015),
    trkIsoWeight = cms.double(1.),
    ecalWeight = cms.double(1.),
    hcalWeight = cms.double(1.),

    # isolation cut
    isolationcut = cms.double(1.0),
   
    # use relative iso
    useRelativeIso = cms.bool(True)
)


