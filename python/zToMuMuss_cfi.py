import FWCore.ParameterSet.Config as cms

zToMuMussplus = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('hTozzTo4leptonsMuonSelector@+ hTozzTo4leptonsMuonSelector@+'),
    cut = cms.string('0.0 < mass < 20000.0 && (daughter(0).charge>0 && daughter(1).charge>0)'),
    checkCharge = cms.bool(True)
)

zToMuMussminus = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('hTozzTo4leptonsMuonSelector@- hTozzTo4leptonsMuonSelector@-'),
    cut = cms.string('0.0 < mass < 20000.0 && (daughter(0).charge<0 && daughter(1).charge<0)'),
    checkCharge = cms.bool(True)
)

zToMuMussmerge = cms.EDProducer("CandViewMerger",
       src = cms.VInputTag( "zToMuMussplus", "zToMuMussminus")
) 

zToMuMuss=cms.Sequence(zToMuMussplus+zToMuMussminus+zToMuMussmerge)


