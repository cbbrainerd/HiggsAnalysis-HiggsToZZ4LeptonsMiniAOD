import FWCore.ParameterSet.Config as cms


# Matching sequence

myMuons = cms.EDFilter("CandViewShallowCloneProducer",
                               src = cms.InputTag("hTozzTo4leptonsMuonSelector"),                            
                               cut = cms.string('')
                                )

myElectrons = cms.EDFilter("CandViewShallowCloneProducer",
                                   src = cms.InputTag("slimmedElectrons"),
                                   cut = cms.string('')
                                )


myGammas = cms.EDFilter("CandViewShallowCloneProducer",
                                   src = cms.InputTag("hTozzTo4leptonsPFfsrPhoton"),
                                   cut = cms.string('')
                                )


goodMuonMCMatch = cms.EDProducer ("MCMatcher",
                                  src     = cms.InputTag("slimmedMuons"),      # RECO objects to match
                                  matched = cms.InputTag("prunedGenParticles"), # mc-truth particle collection
                                  mcPdgId     = cms.vint32(13),           # one or more PDG ID (13 = muon); absolute values (see below)
                                  checkCharge = cms.bool(True),           # True = require RECO and MC objects to have the same charge
                                  mcStatus = cms.vint32(1),               # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                  maxDeltaR = cms.double(0.15),           # MAX deltaR for the match
                                  maxDPtRel = cms.double(0.5),            # MAX deltaPt/Pt for the match
                                  resolveAmbiguities = cms.bool(True),    # Forbid two RECO objects to match to the same GEN object
                                  resolveByMatchQuality = cms.bool(True)  # False = just match input in order; True = pick lowest deltaR pair first
                                  )

goodElectronMCMatch = cms.EDProducer ("MCMatcher",
                                      #src     = cms.InputTag("myElectrons"),  # RECO objects to match  
                                      src     = cms.InputTag("slimmedElectrons"),  # RECO objects to match                                  
                                      matched = cms.InputTag("prunedGenParticles"), # mc-truth particle collection
                                      mcPdgId     = cms.vint32(11),           # one or more PDG ID (13 = muon); absolute values (see below)
                                      checkCharge = cms.bool(True),           # True = require RECO and MC objects to have the same charge
                                      mcStatus = cms.vint32(1),               # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                      maxDeltaR = cms.double(0.15),           # MAX deltaR for the match
                                      maxDPtRel = cms.double(0.5),            # MAX deltaPt/Pt for the match
                                      resolveAmbiguities = cms.bool(True),    # Forbid two RECO objects to match to the same GEN object
                                      resolveByMatchQuality = cms.bool(True)  # False = just match input in order; True = pick lowest deltaR pair first
                                      )

goodGammaMCMatch = cms.EDProducer ("MCMatcher",
                                   src     = cms.InputTag("myGammas"),     # RECO objects to match
                                   matched = cms.InputTag("prunedGenParticles"), # mc-truth particle collection
                                   mcPdgId     = cms.vint32(22),           # one or more PDG ID (13 = muon); absolute values (see below)
                                   checkCharge = cms.bool(True),           # True = require RECO and MC objects to have the same charge
                                   mcStatus = cms.vint32(1),               # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                   maxDeltaR = cms.double(0.5),            # MAX deltaR for the match
                                   maxDPtRel = cms.double(1.),             # MAX deltaPt/Pt for the match
                                   resolveAmbiguities = cms.bool(True),    # Forbid two RECO objects to match to the same GEN object
                                   resolveByMatchQuality = cms.bool(True)  # False = just match input in order; True = pick lowest deltaR pair first
                                   )

goodZtoMuMuMCMatch = cms.EDProducer ("MCMatcher",
                                     src     = cms.InputTag("zToMuMu"),                                     
                                     matched = cms.InputTag("prunedGenParticles"),                                     
                                     mcPdgId     = cms.vint32(23),                                     
                                     checkCharge = cms.bool(True),                                     
                                     mcStatus = cms.vint32(2),                                     
                                     maxDeltaR = cms.double(0.15),                                     
                                     maxDPtRel = cms.double(0.5),                                     
                                     resolveAmbiguities = cms.bool(True),                                     
                                     resolveByMatchQuality = cms.bool(True)                                     
                                     )


goodZtoEEMCMatch = cms.EDProducer ("MCMatcher",
                                   src     = cms.InputTag("zToEE"),                                   
                                   matched = cms.InputTag("prunedGenParticles"),                                   
                                   mcPdgId     = cms.vint32(23),                                   
                                   checkCharge = cms.bool(True),                                   
                                   mcStatus = cms.vint32(2),                                   
                                   maxDeltaR = cms.double(0.15),                                   
                                   maxDPtRel = cms.double(0.5),                                   
                                   resolveAmbiguities = cms.bool(True),                                   
                                   resolveByMatchQuality = cms.bool(True)                                   
                                   )


goodHiggsTozzToEEMMMCMatch = cms.EDProducer ("MCMatcher",
                                        src     = cms.InputTag("hTozzTo4leptonsLooseIsol"),                                             
                                        matched = cms.InputTag("prunedGenParticles"),                                             
                                        mcPdgId     = cms.vint32(25),                                             
                                        checkCharge = cms.bool(True),                                             
                                        mcStatus = cms.vint32(2),                                             
                                        maxDeltaR = cms.double(0.15),                                             
                                        maxDPtRel = cms.double(0.5),                                             
                                        resolveAmbiguities = cms.bool(False),                                             
                                        resolveByMatchQuality = cms.bool(True)                                             
                                       )



goodHiggsTozzToMMMMMCMatch = cms.EDProducer ("MCMatcher",
                                             src     = cms.InputTag("hTozzTo4leptonsMMMMLooseIsol"),                                             
                                             matched = cms.InputTag("prunedGenParticles"),                                             
                                             mcPdgId     = cms.vint32(25),                                             
                                             checkCharge = cms.bool(True),                                             
                                             mcStatus = cms.vint32(2),                                             
                                             maxDeltaR = cms.double(0.15),                                             
                                             maxDPtRel = cms.double(0.5),                                             
                                             resolveAmbiguities = cms.bool(False),                                             
                                             resolveByMatchQuality = cms.bool(True)                                             
                                             )




goodHiggsTozzToEEEEMCMatch = cms.EDProducer ("MCMatcher",
                                             src     = cms.InputTag("hTozzTo4leptonsEEEELooseIsol"),       
                                             matched = cms.InputTag("prunedGenParticles"), 
                                             mcPdgId     = cms.vint32(25),          
                                             checkCharge = cms.bool(True),           
                                             mcStatus = cms.vint32(2),               
                                             maxDeltaR = cms.double(0.15),                                             
                                             maxDPtRel = cms.double(0.5),                                             
                                             resolveAmbiguities = cms.bool(False),                                             
                                             resolveByMatchQuality = cms.bool(True)                                             
                                             )


hTozzTo4leptonsMatchingSequence = cms.Sequence(
    myMuons             +
    goodMuonMCMatch     +
    myElectrons         +
    goodElectronMCMatch +
    myGammas            +
    goodGammaMCMatch    
   # goodZtoMuMuMCMatch  +
   # goodZtoEEMCMatch    +
   # goodHiggsTozzToEEMMMCMatch +
   # goodHiggsTozzToMMMMMCMatch +
   # goodHiggsTozzToEEEEMCMatch 
    )
