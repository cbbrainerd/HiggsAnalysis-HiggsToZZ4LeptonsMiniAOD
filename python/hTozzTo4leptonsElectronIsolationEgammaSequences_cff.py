import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositTk_cff import *
eleIsoDepositTk.src = "hTozzTo4leptonsElectronSelector"

from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositEcalFromHits_cff import *
eleIsoDepositEcalFromHits.ExtractorPSet.barrelEcalHits = cms.InputTag("reducedEgamma","reducedEBRecHits")
eleIsoDepositEcalFromHits.ExtractorPSet.endcapEcalHits = cms.InputTag("reducedEgamma","reducedEERecHits")
eleIsoDepositEcalFromHits.src = "hTozzTo4leptonsElectronSelector"

from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositHcalFromTowers_cff import *
eleIsoDepositHcalFromTowers.src = "hTozzTo4leptonsElectronSelector"

## Old prescription for veto
# eleIsoFromDepsTkOptimized = cms.EDProducer("CandIsolatorFromDeposits",
#   deposits = cms.VPSet(cms.PSet(
#       mode = cms.string('sum'),
#       src = cms.InputTag("eleIsoDepositTk"),
#       weight = cms.string('1'),
#       deltaR = cms.double(0.3),
#       vetos = cms.vstring('muons:0.01', 
#                           'hTozzTo4leptonsElectronSelector:0.015', 
#                           'hTozzTo4leptonsElectronSelector:RectangularEtaPhiVeto(-0.005,0.005,-0.3,0.3)', 
#                           'EcalEndcaps:RectangularEtaPhiVeto(-0.005,0.005,-0.5,0.5)',
#                           'Threshold(0.7)'),
#       skipDefaultVeto = cms.bool(True)
#   ))
# )


eleIsoFromDepsTkOptimized = cms.EDProducer("CandIsolatorFromDeposits",
     deposits = cms.VPSet(cms.PSet(
         mode = cms.string('sum'),
         src = cms.InputTag("eleIsoDepositTk"),
         weight = cms.string('1'),
         deltaR = cms.double(0.3),
         vetos = cms.vstring(
             'vetoMuons:0.015',
             'vetoElectrons:RectangularEtaPhiVeto(-0.015,0.015,-0.5,0.5)',
             'RectangularEtaPhiVeto(-0.015,0.015,-0.5,0.5)',
             'Threshold(0.7)'),
         skipDefaultVeto = cms.bool(True) 
   ))
)


eleIsoFromDepsEcalFromHitsByCrystalOptimized = cms.EDProducer("CandIsolatorFromDeposits",
   deposits = cms.VPSet(cms.PSet(
       mode = cms.string('sum'),
       src = cms.InputTag("eleIsoDepositEcalFromHits"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('NumCrystalVeto(3.0)',
                           'EcalBarrel:NumCrystalEtaPhiVeto(1.0,9999.0)',
                           'EcalEndcaps:NumCrystalEtaPhiVeto(1.5,9999.0)',
                           'EcalBarrel:AbsThresholdFromTransverse(0.08)',
                           'EcalEndcaps:AbsThreshold(0.20)',
                           'muons:0.05',
                           'hTozzTo4leptonsElectronSelector:NumCrystalVeto(3.0)',
                           'hTozzTo4leptonsElectronSelector:NumCrystalEtaPhiVeto(1.5,15.0)'),
       skipDefaultVeto = cms.bool(True)
   ))
)

eleIsoFromDepsHcalFromTowersOptimized = cms.EDProducer("CandIsolatorFromDeposits",
   deposits = cms.VPSet(cms.PSet(
       src = cms.InputTag("eleIsoDepositHcalFromTowers"),
       deltaR = cms.double(0.4),
       weight = cms.string('1'),
       vetos = cms.vstring('0.15', 'muons:0.05'),
       skipDefaultVeto = cms.bool(True),
       mode = cms.string('sum')
   ))
) 

# from Configuration.Geometry.GeometryDB_cff import *

hTozzTo4leptonsElectronIsolationMakeIsoDeposits = cms.Sequence(
    eleIsoDepositTk+
    eleIsoDepositEcalFromHits+
    eleIsoDepositHcalFromTowers
)

hTozzTo4leptonsElectronIsolationMakeIsoValueMaps = cms.Sequence(
    eleIsoFromDepsTkOptimized+
    eleIsoFromDepsEcalFromHitsByCrystalOptimized+
    eleIsoFromDepsHcalFromTowersOptimized
)

hTozzTo4leptonsElectronIsolationDepositSequence = cms.Sequence(
    hTozzTo4leptonsElectronIsolationMakeIsoDeposits*
    hTozzTo4leptonsElectronIsolationMakeIsoValueMaps
)
