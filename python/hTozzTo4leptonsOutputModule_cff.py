import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_EventContent_cff import *
from HiggsAnalysis.HiggsToZZ4Leptons.RECOSIMhTozzTo4leptons_EventContent_cff import *
hTozzTo4leptonsOutputModule = cms.OutputModule("PoolOutputModule",
    hTozzTo4leptonsEventContent,
    hTozzTo4leptonsEventSelection,
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('hTozzTo4leptons'),
        dataTier = cms.untracked.string('USER')
    ),
    fileName = cms.untracked.string('hTozzTo4leptons.root')
)

hTozzTo4leptonsOffselOutputModule = cms.OutputModule("PoolOutputModule",
    hTozzTo4leptonsEventContent,
    hTozzTo4leptonsEventOffSelection,
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('hTozzTo4leptonsOffsel'),
        dataTier = cms.untracked.string('USER')
    ),
    fileName = cms.untracked.string('hTozzTo4leptonsOffsel.root')
)

hTozzTo4leptonsSelectionOutputModule = cms.OutputModule("PoolOutputModule",
    hTozzTo4leptonsEventCompleteSelection,
    hTozzTo4leptonsEventContent,
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('hTozzTo4leptonsSel'),
        dataTier = cms.untracked.string('USER')
    ),
    fileName = cms.untracked.string('hTozzTo4leptonsSel.root')
)


hTozzTo4leptonsSelectiontwoetwomuOutputModule = cms.OutputModule("PoolOutputModule",
    hTozzTo4leptonsEventCompleteSelectionTwoeTwomu,
    hTozzTo4leptonsEventContent,
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('hTozzTo4leptonsSel2e2mu'),
        dataTier = cms.untracked.string('USER')
    ),
    fileName = cms.untracked.string('hTozzTo4leptonsSel2e2mu.root')
)


hTozzTo4leptonsSelectionfoureOutputModule = cms.OutputModule("PoolOutputModule",
    hTozzTo4leptonsEventCompleteSelection4e,
    hTozzTo4leptonsEventContent,
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('hTozzTo4leptonsSel4e'),
        dataTier = cms.untracked.string('USER')
    ),
    fileName = cms.untracked.string('hTozzTo4leptonsSel4e.root')
)


hTozzTo4leptonsSelectionfourmuOutputModule = cms.OutputModule("PoolOutputModule",
    hTozzTo4leptonsEventCompleteSelection4mu,
    hTozzTo4leptonsEventContent,
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('hTozzTo4leptonsSel4mu'),
        dataTier = cms.untracked.string('USER')
    ),
    fileName = cms.untracked.string('hTozzTo4leptonsSel4mu.root')
)


hTozzTo4leptonsOutputModuleAll = cms.OutputModule("PoolOutputModule",
    hTozzTo4leptonsEventContent,
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('USER')
    ),
    fileName = cms.untracked.string('hTozzTo4leptons.root')
)

