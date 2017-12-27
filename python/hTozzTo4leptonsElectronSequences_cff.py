import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronIdSequences_cff import *
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronSelector_cfi import *

hTozzTo4leptonsElectronSequence = cms.Sequence(hTozzTo4leptonsElectronIdSequence + hTozzTo4leptonsElectronSelector)

