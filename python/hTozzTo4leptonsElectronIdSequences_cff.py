import FWCore.ParameterSet.Config as cms

from RecoEgamma.ElectronIdentification.electronIdCutBasedClassBasedExt_cfi import *
import RecoEgamma.ElectronIdentification.electronIdCutBasedClassBasedExt_cfi
eidClassLoose = RecoEgamma.ElectronIdentification.electronIdCutBasedClassBasedExt_cfi.eidCutBasedClassBasedExt.clone()
eidClassLoose.electronQuality.string = 'Eff95Cuts'
eidClassLoose.src = "gsfElectrons"
eidClassMedium = RecoEgamma.ElectronIdentification.electronIdCutBasedClassBasedExt_cfi.eidCutBasedClassBasedExt.clone()
eidClassMedium.electronQuality.string = 'Eff90Cuts'
eidClassLoose.src = "gsfElectrons"
hTozzTo4leptonsElectronIdSequence = cms.Sequence( eidClassLoose + eidClassMedium )
