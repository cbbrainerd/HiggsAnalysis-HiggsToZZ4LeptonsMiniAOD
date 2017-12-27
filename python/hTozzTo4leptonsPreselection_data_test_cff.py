import FWCore.ParameterSet.Config as cms


# Generic MC Truth analysis
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsGenSequence_cff import *

## 
# Filter to select 2e2mu events
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi
hTozzTo4leptonsMCGenFilter2e2mu = HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi.hTozzTo4leptonsMCGenFilter.clone()
hTozzTo4leptonsMCGenFilter2e2mu.HZZ4LeptonsMCFilterLeptonFlavour=cms.int32(3)

# ParticleListDrawer
from  HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenParticleListDrawer_cfi import *
hTozzTo4leptonsMCGenParticleListDrawer2e2mu = HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenParticleListDrawer_cfi.hTozzTo4leptonsMCGenParticleListDrawer.clone()

# Save MC truth: 
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCDumper_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCDumper_cfi
#hTozzTo4leptonsMCDumper.status=cms.vint32(62)

# CP producer: 
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi import *
hTozzTo4leptonsMCCP=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsCP_cfi.hTozzTo4leptonsCP.clone()
hTozzTo4leptonsMCCP.RECOcollName = cms.InputTag("hTozzTo4leptonsMCDumper")

# PF muons
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsPFtoRECOMuon_cfi import *

# PF photons
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsPFfsrPhoton_cfi import *

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHLTInfo_cfi import *
hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
hTozzTo4leptonsHLTInfo.debug = cms.untracked.bool(False)

from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHLTAnalysisFilter_cfi import *


useSkimEarlyData='false'

if useSkimEarlyData == 'true':
    # electrons
    from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronSequences_cff import *
    ELECTRON_BASE_CUT=("(pt > 5 &&" +
                       " fbrem > 0 &&" +
                       " eSuperClusterOverP < 3 &&" +
                       " hcalOverEcal < 0.15 &&" +
                       " abs(deltaPhiSuperClusterTrackAtVtx) < 0.10 &&" +
                       " abs(deltaEtaSuperClusterTrackAtVtx) < 0.02 &&" +
                       " (( isEB && sigmaIetaIeta < 0.015) ||" +
                       "  (!isEB && sigmaIetaIeta < 0.035)) )");
    hTozzTo4leptonsElectronSelector = cms.EDFilter("GsfElectronRefSelector",
                                         src = cms.InputTag("gedGsfElectrons"),
                                         cut = cms.string(ELECTRON_BASE_CUT),
                                         )
    hTozzTo4leptonsElectronSequence=cms.Sequence(hTozzTo4leptonsElectronIdSequence + hTozzTo4leptonsElectronSelector)

    # muons
    from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonSelector_cfi import *
    TM_ARBITRATION = "numberOfMatches('SegmentAndTrackArbitration')>0";
    MUON_BASE_CUT="(isGlobalMuon || (isTrackerMuon && "+TM_ARBITRATION+"))"
    hTozzTo4leptonsMuonSelector = cms.EDFilter("MuonRefSelector",
                                     src = cms.InputTag("muons"),
                                     cut = cms.string(MUON_BASE_CUT),
                                     )
else:

    # Electron Preselector
    from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronSequences_cff import *
    hTozzTo4leptonsElectronPreSelector=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronSelector_cfi.hTozzTo4leptonsElectronSelector.clone()
    hTozzTo4leptonsElectronPreSelector.electronCollection=cms.InputTag("gedGsfElectrons")
    hTozzTo4leptonsElectronPreSelector.electronEtaMax=cms.double(2.5)
    hTozzTo4leptonsElectronPreSelector.electronPtMin=cms.double(5.)
    hTozzTo4leptonsElectronPreSelector.useEleID=cms.bool(False)

    # Electron ordering in pT
    hTozzTo4leptonsElectronOrdering = cms.EDProducer("HZZ4LeptonsElectronOrdering",
     electronCollection = cms.InputTag("hTozzTo4leptonsElectronPreSelector"),
    )

    # Electron scale calibration
    from EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi import *
    calibratedElectrons.electrons = cms.InputTag('hTozzTo4leptonsElectronOrdering')
    calibratedElectrons.correctionFile = cms.string(files["76XReReco"])
    calibratedElectrons.isMC = cms.bool(True)

    # Electron relaxed selection
    from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronSequences_cff import *
    hTozzTo4leptonsElectronSelector=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsElectronSelector_cfi.hTozzTo4leptonsElectronSelector.clone()
    hTozzTo4leptonsElectronSelector.electronCollection=cms.InputTag("calibratedElectrons")
    hTozzTo4leptonsElectronSelector.electronEtaMax=cms.double(2.5)
    hTozzTo4leptonsElectronSelector.electronPtMin=cms.double(3.)
    hTozzTo4leptonsElectronSelector.useEleID=cms.bool(False)
    #hTozzTo4leptonsElectronSequence=cms.Sequence(hTozzTo4leptonsElectronIdSequence + hTozzTo4leptonsElectronSelector)
    hTozzTo4leptonsElectronSequence=cms.Sequence(hTozzTo4leptonsElectronSelector)

    # Muon Calibration
    from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonCalibrator_cfi import *

    # Muon ghost cleaning
    from HiggsAnalysis.HiggsToZZ4Leptons.muonCleanerBySegments_cfi import *
    cleanMuonsBySegments.src = cms.InputTag("hTozzTo4leptonsMuonCalibrator")

    # Muon relaxed selection
    from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonSelector_cfi import *
    hTozzTo4leptonsMuonSelector=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonSelector_cfi.hTozzTo4leptonsMuonSelector.clone()
    hTozzTo4leptonsMuonSelector.muonCollection = cms.InputTag("cleanMuonsBySegments")
    hTozzTo4leptonsMuonSelector.isGlobalMuon=cms.bool(False)
    hTozzTo4leptonsMuonSelector.isTrackerMuon=cms.bool(True)
    hTozzTo4leptonsMuonSelector.muonPtMin=cms.double(3.)
    hTozzTo4leptonsMuonSelector.muonEtaMax=cms.double(2.5)
    hTozzTo4leptonsMuonSequence=cms.Sequence(hTozzTo4leptonsMuonSelector)



#*******************
#2 Leptons No Presel 
#*******************

# zToEE
from HiggsAnalysis.HiggsToZZ4Leptons.zToEE_cfi import *
zToEE.decay = cms.string('hTozzTo4leptonsElectronSelector@+ hTozzTo4leptonsElectronSelector@-')
zToEE.cut = cms.string('mass > 0')

# zToMuMu
from HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi import *
zToMuMu.decay=cms.string('hTozzTo4leptonsMuonSelector@+ hTozzTo4leptonsMuonSelector@-')
zToMuMu.cut = cms.string('mass > 0')

# zToMuMu_SS, zToEE_SS and zToCrossLeptons
from HiggsAnalysis.HiggsToZZ4Leptons.zToMuMuss_cfi import *
zToMuMussplus.decay = cms.string('hTozzTo4leptonsMuonSelector@+ hTozzTo4leptonsMuonSelector@+')
zToMuMussplus.cut = cms.string('mass > 0 && (daughter(0).charge>0 && daughter(1).charge>0)')
zToMuMussminus.decay = cms.string('hTozzTo4leptonsMuonSelector@- hTozzTo4leptonsMuonSelector@-')
zToMuMussminus.cut = cms.string('mass > 0 && (daughter(0).charge<0 && daughter(1).charge<0)')
zToMuMussmerge.src = cms.VInputTag( "zToMuMussplus", "zToMuMussminus")
zToMuMuss=cms.Sequence(zToMuMussplus+zToMuMussminus+zToMuMussmerge)

from HiggsAnalysis.HiggsToZZ4Leptons.zToEEss_cfi import *
zToEEssplus.decay = cms.string('hTozzTo4leptonsElectronSelector@+ hTozzTo4leptonsElectronSelector@+')
zToEEssplus.cut = cms.string('mass > 0 && (daughter(0).charge>0 && daughter(1).charge>0)')
zToEEssminus.decay = cms.string('hTozzTo4leptonsElectronSelector@- hTozzTo4leptonsElectronSelector@-')
zToEEssminus.cut = cms.string('mass > 0 && (daughter(0).charge<0 && daughter(1).charge<0)')
zToEEssmerge.src = cms.VInputTag( "zToEEssplus", "zToEEssminus")
zToEEss=cms.Sequence(zToEEssplus+zToEEssminus+zToEEssmerge)


from HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi import *
zToCrossLeptons.decay = cms.string("hTozzTo4leptonsMuonSelector hTozzTo4leptonsElectronSelector")
zToCrossLeptons.checkCharge = cms.bool(False)
zToCrossLeptons.cut = cms.string('mass > 0')

dileptons = cms.EDProducer("CandViewMerger",
       src = cms.VInputTag( "zToEE", "zToMuMu","zToEEssmerge","zToMuMussmerge","zToCrossLeptons")
)

#*******************
#4 Leptons No Presel 
#*******************


# hTozzToEEMuMu
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi import *
hTozzTo4leptons=HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi.zToMuMu.clone()
hTozzTo4leptons.decay = cms.string('zToEE zToMuMu')
hTozzTo4leptons.cut = cms.string('mass > 0')

# hTozzToMuMuMuMu
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi import *
hTozzTo4leptonsMMMM=HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi.zToMuMu.clone()
hTozzTo4leptonsMMMM.decay = cms.string('zToMuMu zToMuMu')
hTozzTo4leptonsMMMM.cut = cms.string('mass > 0')

# hTozzToEEEE
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi import *
hTozzTo4leptonsEEEE=HiggsAnalysis.HiggsToZZ4Leptons.zToMuMu_cfi.zToMuMu.clone()
hTozzTo4leptonsEEEE.decay = cms.string('zToEE zToEE')
hTozzTo4leptonsEEEE.cut = cms.string('mass > 0')


# other 4 leptons combinations with SS Z and Opposite Flavour Z
quadLeptons4Mu=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
quadLeptons4Mu.decay = cms.string('zToMuMussmerge zToMuMussmerge')
quadLeptons4Mu.cut = cms.string('mass > 0')

quadLeptons2Mu2E=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
quadLeptons2Mu2E.decay = cms.string('zToMuMussmerge zToEEssmerge')
quadLeptons2Mu2E.cut = cms.string('mass > 0')

quadLeptons4E=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
quadLeptons4E.decay = cms.string('zToEEssmerge zToEEssmerge')
quadLeptons4E.cut = cms.string('mass > 0')

#one Z SS and on Z OS
quadLeptonsSSOSele=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
quadLeptonsSSOSele.decay = cms.string('zToEE zToEEssmerge')
quadLeptonsSSOSele.cut = cms.string('mass > 0')

quadLeptonsSSOSmu=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
quadLeptonsSSOSmu.decay = cms.string('zToMuMu zToMuMussmerge')
quadLeptonsSSOSmu.cut = cms.string('mass > 0')

quadLeptonsSSOSmuele=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
quadLeptonsSSOSmuele.decay = cms.string('zToMuMu zToEEssmerge')
quadLeptonsSSOSmuele.cut = cms.string('mass > 0')

quadLeptonsSSOSelemu=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
quadLeptonsSSOSelemu.decay = cms.string('zToEE zToMuMussmerge')
quadLeptonsSSOSelemu.cut = cms.string('mass > 0')

#3Mu+1E, 3E+1Mu
quadLeptons3Mu1E1Z =HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi.hTozzTo4leptons.clone()
quadLeptons3Mu1E1Z.decay = cms.string('zToMuMu zToCrossLeptons')
quadLeptons3Mu1E1Z.checkCharge = cms.bool(False)
quadLeptons3Mu1E1Z.cut = cms.string('mass > 0')

quadLeptons3Mu1E0Z =HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi.hTozzTo4leptons.clone()
quadLeptons3Mu1E0Z.decay = cms.string('zToMuMussmerge zToCrossLeptons')
quadLeptons3Mu1E0Z.checkCharge = cms.bool(False)
quadLeptons3Mu1E0Z.cut = cms.string('mass > 0')

quadLeptons3E1Mu1Z =HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi.hTozzTo4leptons.clone()
quadLeptons3E1Mu1Z.decay = cms.string('zToEE zToCrossLeptons')
quadLeptons3E1Mu1Z.checkCharge = cms.bool(False)
quadLeptons3E1Mu1Z.cut = cms.string('mass > 0')

quadLeptons3E1Mu0Z =HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi.hTozzTo4leptons.clone()
quadLeptons3E1Mu0Z.decay = cms.string('zToEEssmerge zToCrossLeptons')
quadLeptons3E1Mu0Z.checkCharge = cms.bool(False)
quadLeptons3E1Mu0Z.cut = cms.string('mass > 0')

#ME + ME
quadLeptonsCrossZ=HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptons_cfi.hTozzTo4leptons.clone()
quadLeptonsCrossZ.decay = cms.string('zToCrossLeptons zToCrossLeptons')
quadLeptonsCrossZ.checkCharge = cms.bool(False)
quadLeptonsCrossZ.cut = cms.string('mass > 0')

#*******************
#3 Leptons No Presel 
#*******************

triLeptonsMuMuMu=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
triLeptonsMuMuMu.decay = cms.string('hTozzTo4leptonsMuonSelector hTozzTo4leptonsMuonSelector hTozzTo4leptonsMuonSelector')
triLeptonsMuMuMu.cut = cms.string("mass > 0");

triLeptonsMuMuE=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
triLeptonsMuMuE.decay = cms.string('hTozzTo4leptonsMuonSelector hTozzTo4leptonsMuonSelector hTozzTo4leptonsElectronSelector')
triLeptonsMuMuE.cut = cms.string("mass > 0");

triLeptonsMuEE=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
triLeptonsMuEE.decay = cms.string('hTozzTo4leptonsMuonSelector hTozzTo4leptonsElectronSelector hTozzTo4leptonsElectronSelector')
triLeptonsMuEE.cut = cms.string("mass > 0");

triLeptonsEEE=HiggsAnalysis.HiggsToZZ4Leptons.zToCrossLeptons_cfi.zToCrossLeptons.clone()
triLeptonsEEE.decay = cms.string('hTozzTo4leptonsElectronSelector hTozzTo4leptonsElectronSelector hTozzTo4leptonsElectronSelector')
triLeptonsEEE.cut = cms.string("mass > 0");

#*************************
#All Combination No Presel 
#*************************


allLLLL = cms.EDProducer("CandViewMerger",
       src = cms.VInputTag("quadLeptons4Mu","quadLeptons2Mu2E","quadLeptons4E","quadLeptonsSSOSele", "quadLeptonsSSOSmu", "quadLeptonsSSOSmuele", "quadLeptonsSSOSelemu", "quadLeptons3Mu1E1Z","quadLeptons3Mu1E0Z","quadLeptons3E1Mu1Z","quadLeptons3E1Mu0Z", "hTozzTo4leptonsEEEE", "hTozzTo4leptonsMMMM", "hTozzTo4leptons", "quadLeptonsCrossZ" )
)


# Veto electrons and muons for isolation
vetoMuons =  cms.EDFilter("MuonRefSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("(isGlobalMuon || isTrackerMuon) && pt>1.")
)

vetoElectrons =  cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag("gedGsfElectrons"),
    #cut = cms.string("pt>7 && gsfTrack().trackerExpectedHitsInner().numberOfHits<2")
    cut = cms.string("pt>7")
)

# MVA Electron ID
from HiggsAnalysis.HiggsToZZ4Leptons.electronIdMVAProducer_CSA14_cfi import *
# mvaTrigV025nsPHYS14.electronTag    = cms.InputTag("gedGsfElectrons")
# mvaNonTrigV025nsPHYS14.electronTag = cms.InputTag("gedGsfElectrons")
mvaTrigV025nsPHYS14.electronTag    = cms.InputTag("hTozzTo4leptonsElectronOrdering")
mvaNonTrigV025nsPHYS14.electronTag = cms.InputTag("hTozzTo4leptonsElectronOrdering")

# New MVA Electron ID
from RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi import *

# Electron PF isolation
from CommonTools.ParticleFlow.PFBRECO_cff import *
from CommonTools.ParticleFlow.Isolation.pfElectronIsolationPFBRECO_cff import *
elPFIsoDepositChargedPFBRECO.src    = cms.InputTag("hTozzTo4leptonsElectronSelector")
elPFIsoDepositChargedAllPFBRECO.src = cms.InputTag("hTozzTo4leptonsElectronSelector")
elPFIsoDepositNeutralPFBRECO.src    = cms.InputTag("hTozzTo4leptonsElectronSelector")
elPFIsoDepositGammaPFBRECO.src      = cms.InputTag("hTozzTo4leptonsElectronSelector")
elPFIsoDepositPUPFBRECO.src         = cms.InputTag("hTozzTo4leptonsElectronSelector")

# Muon PF isolation
# from CommonTools.ParticleFlow.PFBRECO_cff import *
#from CommonTools.ParticleFlow.Isolation.pfMuonIsolationPFBRECO_cff import *
#muPFIsoDepositChargedPFBRECO.src    = cms.InputTag("hTozzTo4leptonsMuonSelector")
#muPFIsoDepositChargedAllPFBRECO.src = cms.InputTag("hTozzTo4leptonsMuonSelector")
#muPFIsoDepositNeutralPFBRECO.src    = cms.InputTag("hTozzTo4leptonsMuonSelector")
#muPFIsoDepositGammaPFBRECO.src      = cms.InputTag("hTozzTo4leptonsMuonSelector")
#muPFIsoDepositPUPFBRECO.src         = cms.InputTag("hTozzTo4leptonsMuonSelector")

# Photon PF
