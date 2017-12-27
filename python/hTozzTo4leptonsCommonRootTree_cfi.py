import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsCommonRootTree = cms.EDAnalyzer("HZZ4LeptonsCommonRootTree",
    decaychannel = cms.string('2e2mu'),
    rootFileName = cms.untracked.string('roottree_2e2mu.root'),
    useRECOformat = cms.untracked.bool(False),                                         

    module_to_search =  cms.untracked.vstring('a','b'),                                     
    par_to_search = cms.untracked.string('filename'),
    # PU
    fillPUinfo = cms.untracked.bool(True),
    PileupSrc  = cms.InputTag("slimmedAddPileupInfo"),
      
    # Generator
    Generator  = cms.InputTag("generator"),                                           

    # HLT
    fillHLTinfo  = cms.untracked.bool(False),
    HLTInfoFired = cms.InputTag("hTozzTo4leptonsHLTInfo"),                                           
    HLTAnalysisinst = cms.string('hTozzTo4leptonsHLTAnalysis'),
    flagHLTnames=cms.VInputTag(cms.InputTag("flagHLTIsoMu11"), cms.InputTag("flagHLTMu15"),cms.InputTag("flagHLTDoubleMu3"),cms.InputTag("flagHLTIsoEle15L1I"),cms.InputTag("flagHLTIsoEle18L1R"),cms.InputTag("flagHLTDoubleIsoEle10L1I"), cms.InputTag("flagHLTDoubleIsoEle12L1R"), cms.InputTag("flagHLTaccept")),

    # Trigger matching                                           
    triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    triggerobjects = cms.InputTag("selectedPatTrigger"),

    triggerbits = cms.InputTag("TriggerResults","","HLT"),

    triggerFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q'),
    triggerMatchObject   =  cms.InputTag("muonTriggerMatchHLT"),
    triggerEleFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q'),                                               
    triggerHLTcollection =  cms.string("hltL3MuonCandidates"),                                              
    triggerMatchObjectEle =  cms.InputTag("electronTriggerMatchHLT"),                                           

    # Skim Early Data Higgs WG
    useSkimEarlyData=cms.untracked.bool(False),
    SkimEarlyDataAnalysisinst = cms.string('hTozzTo4leptonsSkimEarlyDataAnalysis'),
    flagSkimEarlyDatanames=cms.VInputTag(cms.InputTag("Flag_spikes"), cms.InputTag("Skim_highEnergyMuons"),cms.InputTag("Skim_highEnergyElectrons"),cms.InputTag("Skim_recoWMNfromPf"),cms.InputTag("Skim_recoWMNfromTc"),cms.InputTag("Skim_recoWENfromPf"),cms.InputTag("Skim_recoWENfromTc"),cms.InputTag("Skim_diMuonsJPsi"),cms.InputTag("Skim_diMuonsZ"),cms.InputTag("Skim_diElectronsZ"),cms.InputTag("Skim_triLeptonsMuMuMu"),cms.InputTag("Skim_triLeptonsMuMuEl"),cms.InputTag("Skim_triLeptonsMuElEl"),cms.InputTag("Skim_triLeptonsElElEl"),cms.InputTag("Skim_quadLeptons4Mu"),cms.InputTag("Skim_quadLeptons2Mu2El"),cms.InputTag("Skim_quadLeptons4El")),                                               

    # MC truth
    fillMCTruth    = cms.untracked.bool(False),
    MCcollName     = cms.InputTag("hTozzTo4leptonsMCDumper"),
    genParticles   = cms.InputTag("prunedGenParticles"),
    fourgenleptons = cms.InputTag("fourgenleptons"),
    digenZ         = cms.InputTag("digenZ"),
                                             
    # RECO                                               
    RECOcollNameBest2e2mu = cms.VInputTag(cms.InputTag("hToZZTo4LeptonsBestCandidateOffselMother"), cms.InputTag("hToZZTo4LeptonsBestCandidateOffselBoson0"), cms.InputTag("hToZZTo4LeptonsBestCandidateOffselBoson1")),
    RECOcollNameBest4mu = cms.VInputTag(cms.InputTag("hToZZTo4LeptonsBestCandidateOffselMother"), cms.InputTag("hToZZTo4LeptonsBestCandidateOffselBoson0"), cms.InputTag("hToZZTo4LeptonsBestCandidateOffselBoson1")),
    RECOcollNameBest4e = cms.VInputTag(cms.InputTag("hToZZTo4LeptonsBestCandidateOffselMother"), cms.InputTag("hToZZTo4LeptonsBestCandidateOffselBoson0"), cms.InputTag("hToZZTo4LeptonsBestCandidateOffselBoson1")),

    # Additional RECO    
    useAdditionalRECO  = cms.untracked.bool(True),
    use2011EA          = cms.untracked.bool(True),                                     
    RECOcollNameZ      = cms.VInputTag(cms.InputTag("zToMuMu"), cms.InputTag("zToEE")),
    zToMuMu            = cms.InputTag("zToMuMu"),
    zToEE              = cms.InputTag("zToEE"),
    RECOcollNameZss    = cms.VInputTag(cms.InputTag("zToMuMussmerge"),cms.InputTag("zToEEssmerge"),cms.InputTag("zToCrossLeptons")),
    zToMuMussmerge     = cms.InputTag("zToMuMussmerge"),
    zToEEssmerge       = cms.InputTag("zToEEssmerge"),
    zToCrossLeptons    = cms.InputTag("zToCrossLeptons"),
    RECOcollNameDiLep  = cms.InputTag("dileptons"),                                                                                     
#    RECOcollNameEEMM   = cms.VInputTag(cms.InputTag("hTozzTo4leptonsLooseIsol")),
#    RECOcollNameMMMM   = cms.VInputTag(cms.InputTag("hTozzTo4leptonsMMMMLooseIsol")),
#    RECOcollNameEEEE   = cms.VInputTag(cms.InputTag("hTozzTo4leptonsEEEELooseIsol")),
    RECOcollNameLLLLss = cms.VInputTag(cms.InputTag("quadLeptons4Mu"),cms.InputTag("quadLeptons2Mu2E"),cms.InputTag("quadLeptons4E")),
    RECOcollNameLLL    = cms.VInputTag(cms.InputTag("triLeptonsMuMuMu"),cms.InputTag("triLeptonsMuMuE"),cms.InputTag("triLeptonsMuEE"),cms.InputTag("triLeptonsEEE")),
    RECOcollNameLLLl   = cms.VInputTag(cms.InputTag("quadLeptons3Mu1E"),cms.InputTag("quadLeptons3E1Mu")),
    RECOcollNameLLLLssos = cms.VInputTag(cms.InputTag("quadLeptonsSSOSele"),cms.InputTag("quadLeptonsSSOSmu"),cms.InputTag("quadLeptonsSSOSelemu"),cms.InputTag("quadLeptonsSSOSmuele")),
#    RECOcollNameLLLL   = cms.InputTag("allLLLL"),

    RECOcollNameEEMM   = cms.InputTag("hTozzTo4leptonsLooseIsol"),
    RECOcollNameMMMM   = cms.InputTag("hTozzTo4leptonsMMMMLooseIsol"),
    RECOcollNameEEEE   = cms.InputTag("hTozzTo4leptonsEEEELooseIsol"),
    quadLeptons4Mu     = cms.InputTag("quadLeptons4Mu"),
    quadLeptons2Mu2E   = cms.InputTag("quadLeptons2Mu2E"),
    quadLeptons4E      = cms.InputTag("quadLeptons4E"),
    triLeptonsMuMuMu   = cms.InputTag("triLeptonsMuMuMu"),
    triLeptonsMuMuE    = cms.InputTag("triLeptonsMuMuE"),
    triLeptonsMuEE     = cms.InputTag("triLeptonsMuEE"),
    triLeptonsEEE      = cms.InputTag("triLeptonsEEE"), 
    quadLeptons3Mu1E   = cms.InputTag("quadLeptons3Mu1E"),
    quadLeptons3E1Mu   = cms.InputTag("quadLeptons3E1Mu"),
    quadLeptonsSSOSele = cms.InputTag("quadLeptonsSSOSele"),
    quadLeptonsSSOSmu  = cms.InputTag("quadLeptonsSSOSmu"),
    quadLeptonsSSOSelemu = cms.InputTag("quadLeptonsSSOSelemu"),
    quadLeptonsSSOSmuele = cms.InputTag("quadLeptonsSSOSmuele"),
    RECOcollNameLLLL   = cms.InputTag("allLLLL"),


                                               
    # isolation Tk, Ecal and Hcal
    SuperClustersLabel       = cms.InputTag("hTozzTo4leptonsMergedSuperClusters"),
    GsfTracksElectronsLabel  = cms.InputTag("electronGsfTracks"),
    ElectronsEgmLabel        = cms.InputTag("hTozzTo4leptonsElectronIsolationProducerEgamma"),                                               
    ElectronsEgmTkMapLabel   = cms.InputTag("eleIsoFromDepsTkOptimized"),
    ElectronsEgmEcalMapLabel = cms.InputTag("eleIsoFromDepsEcalFromHitsByCrystalOptimized"),
    ElectronsEgmHcalMapLabel = cms.InputTag("eleIsoFromDepsHcalFromTowersOptimized"),

#    MuonsLabel               = cms.InputTag("hTozzTo4leptonsMuonIsolationProducer"),
    MuonsCorrPtErrorMapLabel = cms.InputTag("hTozzTo4leptonsMuonCalibrator:CorrPtError"),

    # PF muons
   PFMuonsLabel             = cms.InputTag("slimmedMuons"),                                           

    # Particle Flow Isolation
    MuonPFIsoValueChargedAll    = cms.InputTag("muPFIsoValueChargedAll03PFBRECO"),
    MuonPFIsoValueCharged       = cms.InputTag("muPFIsoValueCharged03PFBRECO"),
    MuonPFIsoValueNeutral       = cms.InputTag("muPFIsoValueNeutral03PFBRECO"),
    MuonPFIsoValueGamma         = cms.InputTag("muPFIsoValueGamma03PFBRECO"),
    MuonPFIsoValuePU            = cms.InputTag("muPFIsoValuePU03PFBRECO"),

    ElectronPFIsoValueChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPFBRECO"),
    ElectronPFIsoValueCharged    = cms.InputTag("elPFIsoValueCharged03PFIdPFBRECO"),
    ElectronPFIsoValueNeutral    = cms.InputTag("elPFIsoValueNeutral03PFIdPFBRECO"),
    ElectronPFIsoValueGamma      = cms.InputTag("elPFIsoValueGamma03PFIdPFBRECO"),
    ElectronPFIsoValuePU         = cms.InputTag("elPFIsoValuePU03PFIdPFBRECO"),
    
    PFPhotonsLabel             = cms.InputTag("hTozzTo4leptonsPFfsrPhoton"),                                           
    PFpterrorLabel             = cms.InputTag("hTozzTo4leptonsPFfsrPhoton:ErrorMap"),

    PhotonPFIsoValueChargedAll = cms.InputTag("phPFIsoValueChargedAll03PFIdPFBRECO"),
    PhotonPFIsoValueCharged    = cms.InputTag("phPFIsoValueCharged03PFIdPFBRECO"),
    PhotonPFIsoValueNeutral    = cms.InputTag("phPFIsoValueNeutral03PFIdPFBRECO"),
    PhotonPFIsoValueGamma      = cms.InputTag("phPFIsoValueGamma03PFIdPFBRECO"),
    PhotonPFIsoValuePU         = cms.InputTag("phPFIsoValuePU03PFIdPFBRECO"),
                                           
    
    # vertexing w.r.t primary vertex DA
    MuonsLabelVert           = cms.InputTag("hTozzTo4leptonsMuonSelector"),    
    MuonsMapLabelVert        = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexMuMap"),
    MuonsMapLabelVertValue   = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexValueMuMap"),
    MuonsMapLabelVertError   = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexErrorMuMap"),

    # vertexing w.r.t primary vertex KF
    MuonsMapLabelVertKF       = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKF:VertexMuMap"),
    MuonsMapLabelVertValueKF  = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKF:VertexValueMuMap"),
    MuonsMapLabelVertErrorKF  = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKF:VertexErrorMuMap"),                                           

    # vertexing w.r.t GD, standard kalman and kinematic fit
    MuonsMapLabelVertGD       = cms.InputTag("hTozzTo4leptonsIpToVtxProducerGD:VertexMuMap"),
    MuonsMapLabelVertGDMMMM   = cms.InputTag("hTozzTo4leptonsIpToVtxProducerGDMMMM:VertexMuMap"),
    MuonsMapLabelVertStd      = cms.InputTag("hTozzTo4leptonsIpToVtxProducerStd:VertexMuMap"),
    MuonsMapLabelVertStdMMMM  = cms.InputTag("hTozzTo4leptonsIpToVtxProducerStdMMMM:VertexMuMap"),
    MuonsMapLabelVertKin      = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKin:VertexMuMap"),
    MuonsMapLabelVertKinMMMM  = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKinMMMM:VertexMuMap"),

    # vertexing w.r.t primary vertex DA                                     
    ElectronsLabelVert         = cms.InputTag("hTozzTo4leptonsElectronIsolationProducerEgamma"),
    ElectronsMapLabelVert      = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexEleMap"),
    ElectronsMapLabelVertValue = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexValueEleMap"), 
    ElectronsMapLabelVertError = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexErrorEleMap"),

    # vertexing w.r.t primary vertex KF
    ElectronsMapLabelVertKF      = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKF:VertexEleMap"),
    ElectronsMapLabelVertValueKF = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKF:VertexValueEleMap"), 
    ElectronsMapLabelVertErrorKF = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKF:VertexErrorEleMap"),                                          
                                               

    # vertexing w.r.t GD, standard kalman and kinematic fit
    ElectronsMapLabelVertGD      = cms.InputTag("hTozzTo4leptonsIpToVtxProducerGD:VertexEleMap"),
    ElectronsMapLabelVertGDEEEE  = cms.InputTag("hTozzTo4leptonsIpToVtxProducerGDEEEE:VertexEleMap"),
    ElectronsMapLabelVertStd     = cms.InputTag("hTozzTo4leptonsIpToVtxProducerStd:VertexEleMap"),
    ElectronsMapLabelVertStdEEEE = cms.InputTag("hTozzTo4leptonsIpToVtxProducerStdEEEE:VertexEleMap"),
    ElectronsMapLabelVertKin     = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKin:VertexEleMap"),
    ElectronsMapLabelVertKinEEEE = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKinEEEE:VertexEleMap"),
                                               

    # vertexing with respect to primary vertex                                           
    MuonsSTIPMapLabelVert          = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipMuMap"),
    MuonsSLIPMapLabelVert          = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipMuMap"),
    MuonsSTIPMapLabelVertValue     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipValueMuMap"),
    MuonsSLIPMapLabelVertValue     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipValueMuMap"),
    MuonsSTIPMapLabelVertError     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipErrorMuMap"),
    MuonsSLIPMapLabelVertError     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipErrorMuMap"),
    	
    ElectronsSTIPMapLabelVert      = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipEleMap"),
    ElectronsSLIPMapLabelVert      = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipEleMap"),   	
    ElectronsSTIPMapLabelVertValue = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipValueEleMap"),
    ElectronsSLIPMapLabelVertValue = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipValueEleMap"),
    ElectronsSTIPMapLabelVertError = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipErrorEleMap"),
    ElectronsSLIPMapLabelVertError = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipErrorEleMap"),

    # electron ID
    eleIDLabel                = cms.VInputTag(cms.InputTag("IsoleidClassLoose"),cms.InputTag("IsoleidClassMedium")),

    # electron regression
    eleRegressionEnergyErrorLabel  = cms.InputTag("eleRegressionEnergy:eneErrorRegForGsfEle"),
    eleRegressionEnergyLabel       = cms.InputTag("eleRegressionEnergy:eneRegForGsfEle"),          

    # MVA ele ID BDT
#    mvaElectronTag            = cms.InputTag("slimmedElectrons"),
    mvaElectronTag            = cms.InputTag("hTozzTo4leptonsElectronOrdering"),
    mvaTrigV0MapTag           = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
#    mvaNonTrigV0MapTag        = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
    mvaNonTrigV0MapTag        = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16V1Values"),
                                        
    # GD                                          
    ftsigmaVert               = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:ftsigma"),
    ftsigmalagVert            = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:ftsigmalag"),
    gdX_Vert                  = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdX"),
    gdY_Vert                  = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdY"),                                           
    gdZ_Vert                  = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdZ"),
    gdlagX_Vert               = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdlagX"),
    gdlagY_Vert               = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdlagY"),                                           
    gdlagZ_Vert               = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdlagZ"),
    gdlagProb_Vert            = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdlagProb"),
    gdlagNdof_Vert            = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdlagNdof"),                                           
    ftsigmaVertMMMM           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:ftsigma"),
    ftsigmalagVertMMMM        = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:ftsigmalag"),
    gdX_VertMMMM              = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdX"),
    gdY_VertMMMM              = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdY"),                                           
    gdZ_VertMMMM              = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdZ"),
    gdlagX_VertMMMM           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdlagX"),
    gdlagY_VertMMMM           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdlagY"),                                           
    gdlagZ_VertMMMM           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdlagZ"),
    gdlagProb_VertMMMM        = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdlagProb"),
    gdlagNdof_VertMMMM        = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdlagNdof"),                                             
    ftsigmaVertEEEE           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:ftsigma"),
    ftsigmalagVertEEEE        = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:ftsigmalag"),
    gdX_VertEEEE              = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdX"),
    gdY_VertEEEE              = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdY"),                                           
    gdZ_VertEEEE              = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdZ"),
    gdlagX_VertEEEE           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdlagX"),
    gdlagY_VertEEEE           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdlagY"),                                           
    gdlagZ_VertEEEE           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdlagZ"),
    gdlagProb_VertEEEE        = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdlagProb"),
    gdlagNdof_VertEEEE        = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdlagNdof"),                                             

    # ConstraintFit for 4l vertex                                        
    StandardFitVertex          = cms.InputTag("hTozzTo4leptonsConstraintFitProducer:StandardFitVertex"),
    StandardFitVertexMMMM      = cms.InputTag("hTozzTo4leptonsConstraintFitProducerMMMM:StandardFitVertex"),                                     
    StandardFitVertexEEEE      = cms.InputTag("hTozzTo4leptonsConstraintFitProducerEEEE:StandardFitVertex"),
    KinematicFitVertex         = cms.InputTag("hTozzTo4leptonsConstraintFitProducer:KinematicFitVertex"),
    KinematicFitVertexMMMM     = cms.InputTag("hTozzTo4leptonsConstraintFitProducerMMMM:KinematicFitVertex"),                                     
    KinematicFitVertexEEEE     = cms.InputTag("hTozzTo4leptonsConstraintFitProducerEEEE:KinematicFitVertex"),                                           
    RefittedMass               = cms.InputTag("hTozzTo4leptonsConstraintFitProducer:RefittedMass"),
    RefittedMassMMMM           = cms.InputTag("hTozzTo4leptonsConstraintFitProducerMMMM:RefittedMass"),
    RefittedMassEEEE           = cms.InputTag("hTozzTo4leptonsConstraintFitProducerEEEE:RefittedMass"),                                           

    # ConstraintFit for dilepton vertex
    StandardFitVertexDiLep     = cms.InputTag("hTozzTo4leptonsConstraintFitProducerDiLeptons:StandardFitVertex"),

    # ConstraintFit for 3l vertex                                               
    StandardFitVertexMMM       = cms.InputTag("hTozzTo4leptonsConstraintFitProducerTriLeptonsMMM:StandardFitVertex"),
    StandardFitVertexMME       = cms.InputTag("hTozzTo4leptonsConstraintFitProducerTriLeptonsMME:StandardFitVertex"),                                     
    StandardFitVertexEEE       = cms.InputTag("hTozzTo4leptonsConstraintFitProducerTriLeptonsEEE:StandardFitVertex"),
    StandardFitVertexMEE       = cms.InputTag("hTozzTo4leptonsConstraintFitProducerTriLeptonsMEE:StandardFitVertex"),                                              
                                                                                                                              

    # Other Objetcs
    PhotonsLabel       = cms.InputTag("slimmedPhotons"),
#    TracksLabel        = cms.InputTag("generalTracks"),
    JetsLabel          = cms.InputTag("hTozzTo4leptonsPFJetSelector"),
#    JetsLabel          = cms.InputTag("ak4PFJetsCorrection"),
#    JetsDataLabel      = cms.InputTag("ak4PFJetsCorrectionData"),
    JetsMVALabel       = cms.InputTag("hTozzTo4leptonsPFJetSelector"),
#    PuJetMvaMCfullDiscrLabel   = cms.InputTag("recoPuJetIdMvaMC:fullDiscriminant"),
#    PuJetMvaMCfullIdLabel      = cms.InputTag("recoPuJetIdMvaMC:full53xId"),                                              
#    PuJetMvaDatafullDiscrLabel = cms.InputTag("recoPuJetIdMvaData:fullDiscriminant"),
#    PuJetMvaDatafullIdLabel    = cms.InputTag("recoPuJetIdMvaData:full53xId"),                                                             
    # RhoJetsLabel       = cms.InputTag("kt6corPFJets:rho"),
    # RhoJetsLabel       = cms.InputTag("kt6PFJetsCentral:rho"),                                           
    # RhoJetsLabel       = cms.InputTag("kt4PFJetsNew:rho"),
    RhoJetsLabel       = cms.InputTag("fixedGridRhoFastjetAll"),
    VerticesLabel      = cms.InputTag("offlineSlimmedPrimaryVertices"),

    # GenJet
    GenJetLabel        = cms.InputTag("slimmedGenJets"),
    # Gen MET
#    GenMETLabel        = cms.InputTag("genMetTrue"),
    # Tracker MET                                           
#    TrackerMETLabel    = cms.InputTag("tcMet"),                                           
    # Calo MET
#    CaloMETLabel           = cms.InputTag("met"),                                           
#    CaloMET_NoHFLabel      = cms.InputTag("metNoHF"),	

 #   useAdditionalMET       = cms.untracked.bool(False),
#    CaloMET_HOLabel        = cms.InputTag("metHO"),
#    CaloMET_OptLabel       = cms.InputTag("metOpt"), 
#    CaloMET_OptNoHFLabel   = cms.InputTag("metOptNoHF"), 
#    CaloMET_OptNoHFHOLabel = cms.InputTag("metOptNoHFHO"), 
#    CaloMET_OptHOLabel     = cms.InputTag("metOptHO"), 
#    CaloMET_NoHFHOLabel    = cms.InputTag("metNoHFHO"),
 
    # PF MET
    PfMETLabel             = cms.InputTag("slimmedMETs"), 
    # HT MET                                          
#    HtMET_IC5Label         = cms.InputTag("htMetIC5"), 
#    HtMET_KT4Label         = cms.InputTag("htMetKT4"),
#    HtMET_KT6Label         = cms.InputTag("htMetKT6"),
    #HtMET_SC5Label         = cms.InputTag("htMetSC5"),
    #HtMET_SC7Label         = cms.InputTag("htMetSC7"),
#    HtMET_SC5Label         = cms.InputTag("htMetAK5"),
#    HtMET_SC7Label         = cms.InputTag("htMetAK7"),                                           
    # JES correction                                           
#    MET_JESCorIC5CaloJetLabel = cms.InputTag("metJESCorIC5CaloJet"),
#    MET_JESCorKT4CaloJetLabel = cms.InputTag("metJESCorKT4CaloJet"),
#    MET_JESCorKT6CaloJetLabel = cms.InputTag("metJESCorKT6CaloJet"),
#    MET_JESCorSC5CaloJetLabel = cms.InputTag("metJESCorSC5CaloJet"),
#    MET_JESCorSC7CaloJetLabel = cms.InputTag("metJESCorSC7CaloJet"),
    # MET correction for muons
    CorMETGlobalMuLabel       = cms.InputTag("corMetGlobalMuons"),

    # btagging                                           
    tCHighEff_bTagLabel  = cms.string("pfTrackCountingHighEffBJetTags"),
    tCHighPur_bTagLabel  = cms.string("pfTrackCountingHighPurBJetTags"),
    jPHighEff_bTagLabel  = cms.string("pfjetProbabilityBJetTags"),
    jBP_bTagLabel        = cms.string("pfjetBProbabilityBJetTags"),
    sSVHighEff_bTagLabel = cms.string("pfSimpleSecondaryVertexHighEffBJetTags"),
    sSVHighPur_bTagLabel = cms.string("pfSimpleSecondaryVertexHighPurBJetTags"),
    cSV_bTagLabel        = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    cSVMVA_bTagLabel     = cms.string("pfCombinedSecondaryVertexMVABJetTags"),
    sEByIP3d_bTagLabel   = cms.InputTag("softElectronByIP3dBJetTags"),
    sEByPt_bTagLabel     = cms.InputTag("softPFElectronByPtBJetTags"),
    sM_bTagLabel         = cms.InputTag("softPFMuonBJetTags"),
    sMByIP3d_bTagLabel   = cms.InputTag("softMuonByIP3dBJetTags"),
    sMByPt_bTagLabel     = cms.InputTag("softMuonByPtBJetTags"),
                                                
    # Conversion finder
    ConvMapDist          = cms.InputTag("ConvValueMapProd:dist"),                          
    ConvMapDcot          = cms.InputTag("ConvValueMapProd:dcot"),

    # Matching
    goodElectronMCMatch  = cms.InputTag("goodElectronMCMatch"),
    myElectrons          = cms.InputTag("myElectrons"),
    goodMuonMCMatch      = cms.InputTag("goodMuonMCMatch"),
    myMuons              = cms.InputTag("myMuons"),
    goodGammaMCMatch      = cms.InputTag("goodGammaMCMatch"),
    myGammas              = cms.InputTag("myGammas"),

    goodZtoMuMuMCMatch   = cms.InputTag("goodZtoMuMuMCMatch"),
    goodZtoEEMCMatch     = cms.InputTag("goodZtoEEMCMatch"),

    goodHiggsTozzToEEMMMCMatch = cms.InputTag("goodHiggsTozzToEEMMMCMatch"),
    goodHiggsTozzToMMMMMCMatch = cms.InputTag("goodHiggsTozzToMMMMMCMatch"),
    goodHiggsTozzToEEEEMCMatch = cms.InputTag("goodHiggsTozzToEEEEMCMatch"),
    
    # Beam Spot
    offlineBeamSpot       = cms.InputTag("offlineBeamSpot"),

    pfCands = cms.InputTag("packedPFCandidates"),
)
