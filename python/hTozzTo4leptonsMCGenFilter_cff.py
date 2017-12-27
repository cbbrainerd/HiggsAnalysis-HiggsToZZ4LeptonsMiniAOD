import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi import *

## 
# Filter to select 2e2mu events
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi
hTozzTo4leptonsMCGenFilter2e2mu = HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi.hTozzTo4leptonsMCGenFilter.clone()
hTozzTo4leptonsMCGenFilter2e2mu.HZZ4LeptonsMCFilterLeptonFlavour=cms.int32(3)

hTozzTo4leptonsSelectionSequence2e2mu = cms.Sequence(
        hTozzTo4leptonsMCGenFilter2e2mu            
	)
                                                 
##
# Filter to select 4mu events
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi
hTozzTo4leptonsMCGenFilter4mu = HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi.hTozzTo4leptonsMCGenFilter.clone()
hTozzTo4leptonsMCGenFilter4mu.HZZ4LeptonsMCFilterLeptonFlavour=cms.int32(1)

hTozzTo4leptonsSelectionSequence4mu = cms.Sequence(
        hTozzTo4leptonsMCGenFilter4mu
        )

##
# Filter to select 4e events
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi
hTozzTo4leptonsMCGenFilter4e = HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi.hTozzTo4leptonsMCGenFilter.clone()
hTozzTo4leptonsMCGenFilter4e.HZZ4LeptonsMCFilterLeptonFlavour=cms.int32(2)

hTozzTo4leptonsSelectionSequence4e = cms.Sequence(
        hTozzTo4leptonsMCGenFilter4e
        )


