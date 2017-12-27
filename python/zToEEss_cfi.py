import FWCore.ParameterSet.Config as cms

zToEEssplus = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('hTozzTo4leptonsElectronSelector@+ hTozzTo4leptonsElectronSelector@+'),
    cut = cms.string('0.0 < mass < 20000.0 && (daughter(0).charge>0 && daughter(1).charge>0)')
)

zToEEssminus = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('hTozzTo4leptonsElectronSelector@- hTozzTo4leptonsElectronSelector@-'),
    cut = cms.string('0.0 < mass < 20000.0 && (daughter(0).charge<0 && daughter(1).charge<0)')
)

zToEEssmerge = cms.EDProducer("CandViewMerger",
       src = cms.VInputTag( "zToEEssplus", "zToEEssminus")
) 

zToEEss=cms.Sequence(zToEEssplus+zToEEssminus+zToEEssmerge)
