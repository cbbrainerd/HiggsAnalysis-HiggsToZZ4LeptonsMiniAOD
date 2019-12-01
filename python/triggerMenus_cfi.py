#This class handles trigger vetos automatically: the first dataset to have triggers added has no vetoes, and further datasets veto all the previously added triggers
#Calling get_dataset for MC returns all the added triggers, with no vetoes
class TriggersByDataset:
    def __init__(self):
        self.triggerMenu={}
        self.all_triggers=[]
    def add_dataset(self,dataset,trigger_list):
        #list(x) copies the list x
        self.triggerMenu[dataset]={'veto_triggers' : list(self.all_triggers)}
        self.triggerMenu[dataset]['pass_triggers']=[]
        for trigger in trigger_list:
            self.triggerMenu[dataset]['pass_triggers'].append(trigger)
            self.all_triggers.append(trigger)
    def get_dataset(self,dataset):
        if(dataset=='MC'):
            return { 'pass_triggers' : list(self.all_triggers) , 'veto_triggers':[] }
        else:
            try:
                return self.triggerMenu[dataset]
            except KeyError:
                print 'Error: dataset "%s" has no defined trigger menu.' % dataset
                raise

Triggers={}
#2016 Trigger Menu
Triggers[2016]=TriggersByDataset()
Triggers[2016].add_dataset('DoubleEG',
    [
    #dielectron
    "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",
    #trielectron
    "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v",
    ]
)

Triggers[2016].add_dataset('DoubleMuon',
    [
    #dimuon
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
    "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
    "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
    #trimuon
    "HLT_TripleMu_12_10_5_v",
    ]
)

Triggers[2016].add_dataset('MuonEG',
    [
    #muele
    "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v",
    "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v",
    ]
)

Triggers[2016].add_dataset('SingleElectron',
    [
    #single electron
    "HLT_Ele25_eta2p1_WPTight_Gsf_v",
    "HLT_Ele27_WPTight_Gsf_v",
    "HLT_Ele27_eta2p1_WPLoose_Gsf_v",
    "HLT_Ele32_eta2p1_WPTight_Gsf_v",
    ]
)

Triggers[2016].add_dataset('SingleElectron',
    [
    #single muon
    "HLT_IsoMu20_v",
    "HLT_IsoTkMu20_v",
    "HLT_IsoMu22_v",
    "HLT_IsoTkMu22_v",
    "HLT_IsoMu24_v",
    "HLT_IsoTkMu24_v",
    ]
)

#2017 Trigger Menu
Triggers[2017]=TriggersByDataset()
Triggers[2017].add_dataset('DoubleEG',
    [
    #dielectron
     "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
     "HLT_DoubleEle33_CaloIdL_MW_v",
    #trielectron
     "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v",
    ]
)

Triggers[2017].add_dataset('DoubleMuon',
    [
    #dimuon
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",
    #trimuon
    "HLT_TripleMu_10_5_5_DZ_v",
    "HLT_TripleMu_12_10_5_v",
    ],
)

Triggers[2017].add_dataset('MuonEG',
    [
    #muele
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v",
    "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v",
    "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v",
    ]
)

Triggers[2017].add_dataset('SingleElectron',
    [
    #single electron
    "HLT_Ele35_WPTight_Gsf_v",
    "HLT_Ele38_WPTight_Gsf_v",
    "HLT_Ele40_WPTight_Gsf_v",
    ]
)

Triggers[2017].add_dataset('SingleMuon',
    [
    #single muon
    "HLT_IsoMu27_v",
    ]
)

Triggers[2018]=TriggersByDataset()
Triggers[2018].add_dataset('DoubleEG',
    [
    #dielectron
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_DoubleEle25_CaloIdL_MW_v",
    #all trielectron triggers are prescaled
    ]
)

Triggers[2018].add_dataset('DoubleMuon',
    [
    #dimuon
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",
    #trimuon
    "HLT_TripleMu_10_5_5_DZ_v",
    "HLT_TripleMu_12_10_5_v",
    ]
)

Triggers[2018].add_dataset('MuonEG',
    [
    #muele
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v",
    "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v",
    ]
)

Triggers[2018].add_dataset('SingleElectron',
    [
    #single electron
    "HLT_Ele32_WPTight_Gsf_v",
    ]
)
Triggers[2018].add_dataset('SingleMuon',
    [
    #single muon
    "HLT_IsoMu24_v",
    ]
)

def get_triggers(year,dataset):
    try:
        return Triggers[year].get_dataset(dataset)
    except KeyError:
        raise
        #raise KeyError("Trigger menu for dataset %s in year %i not defined." % (dataset,year))
