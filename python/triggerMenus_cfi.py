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
#2017 Trigger Menu
Triggers[2017]=TriggersByDataset()
Triggers[2017].add_dataset('DoubleEG',
    [
     "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_",
     "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL",
     "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
    ]
)

Triggers[2017].add_dataset('DoubleMuon',
    [
     "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
     "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
     "HLT_TripleMu_12_10_5",
     "HLT_TripleMu_10_5_5_D",
    ],
)

Triggers[2017].add_dataset('MuonEG',
    [
     "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
     "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
     "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
     "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
     "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
     "HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
     "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ",
    ]
)

Triggers[2017].add_dataset('SingleElectron',
    [
     'HLT_Ele35_WPTight_Gsf_v',
     'HLT_Ele38_WPTight_Gsf_v',
     'HLT_Ele40_WPTight_Gsf_v',
    ]
)

Triggers[2017].add_dataset('SingleMuon',
    [
     'HLT_IsoMu27',
    ]
)

Triggers[2018]=TriggersByDataset()
Triggers[2018].add_dataset('DoubleEG',
    [
     'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v',
     'HLT_DoubleEle25_CaloIdL_MW_v',
    ]
)

Triggers[2018].add_dataset('DoubleMuon',
    [
     'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v',
     'HLT_TripleMu_10_5_5_DZ_v',
     'HLT_TripleMu_12_10_5_v',
    ]
)

Triggers[2018].add_dataset('MuonEG',
    [
     "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
     "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
     "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
     "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
     "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v",
     "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v",
     "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v",
    ]
)

Triggers[2018].add_dataset('SingleElectron',
    [
     'HLT_Ele32_WPTight_Gsf_v',
    ]
)
Triggers[2018].add_dataset('SingleMuon',
    [
     'HLT_IsoMu24_v',
    ]
)

def get_triggers(year,dataset):
    try:
        return Triggers[year].get_dataset(dataset)
    except KeyError:
        raise
        #raise KeyError("Trigger menu for dataset %s in year %i not defined." % (dataset,year))
