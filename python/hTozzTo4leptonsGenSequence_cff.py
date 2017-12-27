import FWCore.ParameterSet.Config as cms


# select only Z, and save clones
genZ = cms.EDFilter("PdgIdAndStatusCandViewSelector",
                          src = cms.InputTag("prunedGenParticles"),
                           pdgId = cms.vint32( 23 ),
#                           status = cms.vint32(3)  
#                           status = cms.vint32(62)                                         
                          status = cms.vint32(22)                    
                           )

# select Gamma, and save clones
Gamma = cms.EDFilter("PdgIdAndStatusCandViewSelector",
                          src = cms.InputTag("prunedGenParticles"),
                          pdgId = cms.vint32( 22 ),
                          #status = cms.vint32(3)
                          status = cms.vint32(1)
                           )


# di-Z
digenZ = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string("genZ genZ"),
				   checkCharge= cms.bool(False),
                                   cut = cms.string("0< mass")
                                   )


# leptons
genleptons = cms.EDFilter("PdgIdAndStatusCandViewSelector",
                                 src = cms.InputTag("prunedGenParticles"),
                                 pdgId = cms.vint32(11,13,15),
                                 status = cms.vint32(1)                     
                                 )

# di-leptons
digenleptons = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string("genleptons@+ genleptons@-"),
                                   cut = cms.string("0< mass")
                                   )
## trileptons
trigenleptons = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string("genleptons@+ genleptons@- genleptons@+"),
                                    cut = cms.string("0< mass")
                                    )

#fourleptons
fourgenleptons = cms.EDProducer("CandViewShallowCloneCombiner",
                                     decay = cms.string("trigenleptons@+ genleptons@-"),
                                     cut = cms.string("0< mass")
                                     )


hTozzTo4leptonsGenSequence = cms.Sequence(
    genZ           +
    digenZ         +
    genleptons     +
    digenleptons   +
    trigenleptons  +
    fourgenleptons )
