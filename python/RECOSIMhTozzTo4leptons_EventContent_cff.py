import FWCore.ParameterSet.Config as cms

from Configuration.EventContent.EventContent_cff import *
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_EventContent_cff import *
RECOSIMhTozzTo4leptonsEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring()
)
RECOSIMhTozzTo4leptonsEventContent.outputCommands.extend(RECOSIMEventContent.outputCommands)

