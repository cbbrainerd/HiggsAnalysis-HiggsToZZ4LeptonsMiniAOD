import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *_*_*_*', 
        'keep edmHepMCProduct_*_*_*', 
        'keep edmGenInfoProduct_*_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep *_allElectrons_*_*', 
        'keep *_allMuons_*_*', 
        'keep *_allTracks_*_*', 
        'keep *_allLeptons_*_*', 
        'keep *_hTozzToEEMuMuElectronCand_*_*', 
        'keep *_hTozzToEEMuMuMuonCand_*_*', 
        'keep *_zToEE_*_*', 
        'keep *_zToMuMu_*_*', 
        'keep *_hTozzToEEMuMu_*_*', 
        'keep *_hTozzToEEMuMuHiggsFrame_*_*', 
        'keep *_hTozzToEEMuMuCP_*_*', 
        'keep *_hTozzToEEMuMuBestCandidateProducer_*_*', 
        'keep *_genParticleCandidates_*_*', 
        'keep *_zToMuMuGenParticlesMatch_*_*', 
        'keep *_zToEEGenParticlesMatch_*_*', 
        'keep *_hTozzToEEMuMuIsolation_*_*', 
        'keep *_hTozzToEEMuMuInnerBrem_*_*', 
        'keep *_globalMuons_*_*', 
        'keep recoMuons_muons_*_*', 
        'keep recoMuons_trackerMuons_*_*', 
        'keep recoMuons_hTozzToEEMuMuMuonSelector_*_*', 
        'keep *_islandBasicClusters_*_*', 
        'keep *_islandSuperClusters_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep hybridSuperClusters_hybridShapeAssoc_*_*', 
        'keep islandBasicClusters_islandEndcapShapeAssoc_*_*', 
        'keep *_correctedIslandEndcapSuperClusters_*_*', 
        'keep *_correctedIslandBarrelSuperClusters_*_*', 
        'keep *_correctedHybridSuperClusters_*_*', 
        'keep *_correctedEndcapSuperClustersWithPreshower_*_*', 
        'keep *_photons_*_*', 
        'keep *_correctedPhotons_*_*', 
        'keep recoConvertedPhotons_convertedPhotons_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_overlapElectronResolver_*_*', 
        'keep recoGsfElectrons_hTozzToEEMuMuElectronSelector_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep *_offlinePrimaryVerticesFromCTFTracks_*_*', 
        'keep *_genMet_*_*', 
        'keep *_met_*_*', 
        'keep *_iterativeCone5CaloJets_*_*', 
        'keep *_hTozzToEEMuMuMC*_*_*', 
        'keep *_MCToRecoMatched*_*_*', 
        'keep *_*_vtxInfo*_*', 
        'keep float_*_*_*', 
        'keep double_*_*_*', 
        'keep ints_*_*_*', 
        'keep bool_*_*_*', 
        'keep TrackingRecHitsOwned_*_*_*', 
        'keep TrajectorySeeds_*_*_*', 
        'keep recoGsfTracks_*_*_*', 
        'keep recoGsfTrackExtras_*_*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep *_recoPuJet*_*_*',
        'keep *_*_*_*')
)
hTozzTo4leptonsEventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('hTozzTo4leptonsPreselectionPath')
    )
)
hTozzTo4leptonsEventOffSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('hTozzTo4leptonsOfflineselectionPath')
    )
)
hTozzTo4leptonsEventCompleteSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('hTozzTo4leptonsSelectionPath')
    )
)

hTozzTo4leptonsEventCompleteSelectionTwoeTwomu = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('hTozzTo4leptonsSelectionPath2e2mu')
    )
)

hTozzTo4leptonsEventCompleteSelection4e = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('hTozzTo4leptonsSelectionPath4e')
    )
)

hTozzTo4leptonsEventCompleteSelection4mu = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('hTozzTo4leptonsSelectionPath4mu')
    )
)




