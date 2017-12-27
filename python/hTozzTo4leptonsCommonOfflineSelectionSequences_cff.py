import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsCommonOfflineSelection = cms.EDProducer("HZZ4LeptonsCommonOfflineSelection",
    decaychannel = cms.string('2e2mu'),
    BestCandidatesLeptons = cms.InputTag("hTozzTo4leptonsBestCandidateProducer:hToZZTo4LeptonsBestCandidateLeptons"),
    useBestCandidate = cms.bool(True),                                                     
    # vertexing
    MuonsLabelVert = cms.InputTag("hTozzTo4leptonsMuonIsolationProducer"),
    MuonsMapLabelVert = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexMuMap"),
    ElectronsLabelVert = cms.InputTag("hTozzTo4leptonsHadIsolationProducer"),
    ElectronsMapLabelVert = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexEleMap"),
    vertVarCut = cms.vdouble(12., 8.),
    # tight isolation
    ElectronsLabel = cms.InputTag("hTozzTo4leptonsElectronIsolationProducer"),
    ElectronsMapLabel = cms.InputTag("hTozzTo4leptonsHadIsolationProducer"),
    MuonsLabel = cms.InputTag("hTozzTo4leptonsMuonIsolationProducer"),
    MuonsMapLabel = cms.InputTag("hTozzTo4leptonsMuonIsolationProducerOffsel"),
    isoVarTagElectrons = cms.InputTag("HadIsolationXele"),
    isoVarCutElectrons = cms.vdouble(0.35),                                                   
    isoVarTagMuons = cms.InputTag("MuonIsolationX"),
    isoVarCutMuons = cms.vdouble(30.0)                                                       
)

hTozzTo4leptonsCommonOfflineSelectionFilter = cms.EDFilter("HZZ4LeptonsCommonOfflineSelectionFilter",
    offseltags = cms.vstring('OffselTightCombIsolEle', 
        'OffselTightCombIsolMu', 
        'OffselVertComb'),
    offSelectFileName = cms.string('offselect2e2mu.out'),
    offselinst = cms.string('hTozzTo4leptonsCommonOfflineSelection'),
    rootFileName = cms.string('offselect2e2mu.root'),
    decaychannel = cms.string('2e2mu')
)

hTozzTo4leptonsCommonOffSelSequence = cms.Sequence(hTozzTo4leptonsCommonOfflineSelection+hTozzTo4leptonsCommonOfflineSelectionFilter)
