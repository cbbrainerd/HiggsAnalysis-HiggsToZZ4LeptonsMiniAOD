import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsCommonPreselection = cms.EDProducer("HZZ4LeptonsCommonPreselection",
    HLabel = cms.InputTag("hTozzTo4leptons"),
    ZEELabel = cms.InputTag("zToEE"),
    decaychannel = cms.string('2e2mu'),
    MuonsLooseIsolLabel = cms.InputTag("hTozzTo4leptonsMuonIsolationProducerMu"),
    MuonsLabel = cms.InputTag("hTozzTo4leptonsMuonSelector"),
    LeptonsLabel = cms.InputTag("allLeptons"),
    ZMMLabel = cms.InputTag("zToMuMu"),
    cuts = cms.PSet(
        fourleptMass = cms.double(100.0),
        nEle = cms.int32(2),
        nMu = cms.int32(2),
        eeMass = cms.double(12.0),
        numberOfeeCombs = cms.int32(1),
        numberOf4lCombs = cms.int32(1),
        nlooseMu = cms.int32(2),
        nlooseEle = cms.int32(2),
        mumuMass = cms.double(12.0),
        numberOfmumuCombs = cms.int32(1)
    ),
    ElectronsLooseIsolLabel = cms.InputTag("hTozzTo4leptonsElectronIsolationProducerEgamma"),
    ElectronsLabel = cms.InputTag("hTozzTo4leptonsElectronSelector")
)

hTozzTo4leptonsCommonPreselectionFilter = cms.EDFilter("HZZ4LeptonsCommonPreselectionFilter",
    preselinst = cms.string('hTozzTo4leptonsCommonPreselection'),
    preseltags = cms.vstring('PreselAtleast2Ele', 
        'PreselAtleast2Mu', 
        'PreselAtleast1ZEE', 
        'PreselAtleast1ZMuMu', 
        'PreselAtleast1H', 
        'PreselLoose2IsolEle', 
        'PreselLoose2IsolMu'),
    preSelectFileName = cms.string('preselect2e2mu.out'),
    rootFileName = cms.string('preselect2e2mu.root'),
    decaychannel = cms.string('2e2mu')
)

hTozzTo4leptonsCommonPreselectionSequence = cms.Sequence(hTozzTo4leptonsCommonPreselection+hTozzTo4leptonsCommonPreselectionFilter)

