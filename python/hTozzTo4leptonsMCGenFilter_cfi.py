import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsMCGenFilter = cms.EDFilter("HZZ4LeptonsMCGenFilter",

    genParticles                     = cms.InputTag("genParticles"),
    DebugHZZ4LeptonsMCGenFilter      = cms.bool(False),
    # Four leptons e, mu, tau in the acceptance
    # LeptonFlavour
    # 0 = 4l including tau
    # 1 = 4 mu
    # 2 = 4 e
    # 3 = 2e 2mu
    # 4 = 2e 2tau
    # 5 = 2mu2tau
    # 6 = 4tau
    # 7 = 2e 2tau(LH or HL)
    # 8 = 2e 2tau(HH)
    # 9 = 2e 2tau(LL)
    # 10= 2mu2tau(LH or HL)
    # 11= 2mu2tau(HH)
    # 12= 2mu2tau(LL)
    # 13= z+light jets
    # 14= zbbbar 
    # 15= zccbar                 
    HZZ4LeptonsMCFilterLeptonFlavour = cms.int32(3),
    acceptance                       = cms.double(2.5),
    # zgen                             = cms.InputTag("Z"),
    # zbbgen                           = cms.InputTag("llbBcands"),
    # zccgen                           = cms.InputTag("llcccands")                                     
)


