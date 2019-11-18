import FWCore.ParameterSet.Config as cms

from triggerMenus_cfi import get_triggers

hTozzTo4leptonsHLTAnalysisFilter = cms.EDFilter("HZZ4LeptonsHLTAnalysisFilter",
    HLTInfoFired = cms.InputTag("hTozzTo4leptonsHLTInfo"),                                           
)

def set_triggers(year,dataset):
    print "Setting up trigger for year %i on dataset %s" % (year,dataset)
    triggers=get_triggers(year,dataset)
    print "Triggering on events that pass the triggers:"
    for trigger in triggers['pass_triggers']:
        print trigger
    veto=triggers['veto_triggers']
    if veto:
        print "Except those that pass the triggers:"
        for trigger in veto:
            print trigger
    else:
        print "No veto triggers."
    hTozzTo4leptonsHLTAnalysisFilter.pass_triggers=cms.vstring(triggers['pass_triggers'])
    hTozzTo4leptonsHLTAnalysisFilter.veto_triggers=cms.vstring(triggers['veto_triggers'])
    print "Trigger configuration done."
