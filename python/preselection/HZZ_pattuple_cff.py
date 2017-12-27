#  HZZ-PAT configuration fragment
#
#  PAT configuration for the hzz group - 42X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/hzzPatLayer1DefV10


import FWCore.ParameterSet.Config as cms

def addDefaultHZZPAT(process,mcInfo=True,HLTMenu='HLT',jetMetCorrections=['L1Offset', 'L2Relative', 'L3Absolute'],mcVersion='',theJetNames = ['AK5PF'],doValidation=False):
	#    loadPF2PAT(process,mcInfo,jetMetCorrections,extMatch,dohzzTopProjection,'PF')
	addTagInfos(process,jetMetCorrections)
	if not mcInfo:
		removeMCDependence(process)
	#    loadMCVersion(process,mcVersion,mcInfo)
	loadPAT(process,jetMetCorrections, mcInfo)
	addJetMET(process,theJetNames,jetMetCorrections,mcVersion)
	useDAVertices(process)
	loadPATTriggers(process,HLTMenu)
	
	#-- Counter for the number of processed events --------------------------------
	process.eventCountProducer = cms.EDProducer("EventCountProducer")
	
	process.load("PhysicsTools.PatAlgos.producersLayer1.photonProducer_cff")
	process.load("PhysicsTools.PatAlgos.selectionLayer1.photonSelector_cfi")

	# Full path
	process.hzzPatDefaultSequence = cms.Sequence( 
		process.cicEleIdSequence+
		process.hzzEleIdSequence+
		process.simpleEleIdSequence+
		process.eventCountProducer 
		* process.daVertices * process.patDefaultSequence 
	)

def addTagInfos(process,jetMetCorrections):
    from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection
    switchJetCollection( process,
                     jetCollection=cms.InputTag('ak5CaloJets'),
                     jetCorrLabel=('AK5Calo', jetMetCorrections))

def loadPAT(process,jetMetCorrections, mcInfo):
    #-- Changes for electron and photon ID ----------------------------------------
    # Turn off photon-electron cleaning (i.e., flag only)
	process.cleanPatPhotons.checkOverlaps.electrons.requireNoOverlaps = False

	# Add CiC electron ID's and WP eleID's
	process.load("TripHto4lep.Skim3lep.cicEleIdSequence_42_cff")
	process.load("TripHto4lep.Skim3lep.hzzEleIdSequence_cff")
	process.load("TripHto4lep.Skim3lep.simpleEleIdSequence_cff")

	process.patElectrons.electronIDSources = cms.PSet(
	eidVeryLoose= cms.InputTag("eidVeryLooseMC"),
    eidLoose= cms.InputTag("eidLooseMC"),
    eidMedium= cms.InputTag("eidMediumMC"),
    eidTight= cms.InputTag("eidTightMC"),
    eidSuperTight= cms.InputTag("eidSuperTightMC"),
    eidHyperTight1= cms.InputTag("eidHyperTight1MC"),
    eidHZZVeryLoose= cms.InputTag("eidHZZVeryLoose"),
    eidHZZLoose= cms.InputTag("eidHZZLoose"),
    eidHZZMedium= cms.InputTag("eidHZZMedium"),
    eidHZZTight= cms.InputTag("eidHZZTight"),
    eidHZZSuperTight= cms.InputTag("eidHZZSuperTight"),
#    eidHZZHyperTight1= cms.InputTag("eidHZZHyperTight1"),
    simpleEleId95relIso= cms.InputTag("simpleEleId95relIso"),
    simpleEleId90relIso= cms.InputTag("simpleEleId90relIso"),
    simpleEleId85relIso= cms.InputTag("simpleEleId85relIso"),
    simpleEleId80relIso= cms.InputTag("simpleEleId80relIso"),
	simpleEleId70relIso= cms.InputTag("simpleEleId70relIso"),
	simpleEleId60relIso= cms.InputTag("simpleEleId60relIso"))
	process.patDefaultSequence.replace(process.patElectrons,process.simpleEleIdSequence+process.cicEleIdSequence+process.hzzEleIdSequence+process.patElectrons)

	# Switch to PF Jets
	from PhysicsTools.PatAlgos.tools.jetTools import *
	#switchJetCollection(process, jetCollection = cms.InputTag("ak5PFJets"), outputModule = 'out')
	switchJetCollection(process,cms.InputTag('ak5PFJets'),
		doJTA        = True,
		doBTagging   = True,
		jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute'])),
		doType1MET   = True,
		genJetCollection=cms.InputTag("ak5GenJets"),
		doJetID      = True
		)
	process.patJetCorrFactors.rho = cms.InputTag('')

def loadPATTriggers(process,HLTMenu):
	#-- Trigger matching ----------------------------------------------------------
	from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger, switchOnTriggerMatchEmbedding
#	from PhysicsTools.PatAlgos.triggerLayer1.triggerMatcher_cfi import cleanPhotonTriggerMatchHLTPhoton20CleanedL1R
	#    process.patPhotonMatch = cleanPhotonTriggerMatchHLTPhoton20CleanedL1R.clone(matchedCuts = cms.string( photonMatches ))
	# firing trigger objects used in succeeding HLT path 'HLT_Photon20_Cleaned_L1R'
	process.patPhotonMatch = cms.EDProducer(
	"PATTriggerMatcherDRDPtLessByR"                 # match by DeltaR only, best match by DeltaR
	, src     = cms.InputTag( "selectedPatPhotons" )
	, matched = cms.InputTag( "patTrigger" )          # default producer label as defined in PhysicsTools/PatAlgos/python/triggerLayer1/triggerProducer_cfi.py
	, matchedCuts = cms.string( 'path( "HLT_Photon26_CaloIdL_IsoVL_Photon18_v*" ) || path( "HLT_Photon20_CaloIdVL*" ) || path( "HLT_Photon30_CaloIdVL_v*" )|| path( "HLT_Photon50_CaloIdVL*" )|| path( "HLT_Photon75_CaloIdVL*" )|| path( "HLT_Photon90_CaloIdVL*" )' )
	#, andOr                      = cms.bool( False )  # AND
	#, filterIdsEnum              = cms.vstring( '*' ) # wildcard, overlaps with 'filterIds'
	#, filterIds                  = cms.vint32( 0 )    # wildcard, overlaps with 'filterIdsEnum'
	#, filterLabels               = cms.vstring( '*' ) # wildcard
	#, pathNames                  = cms.vstring(
	#'HLT_Photon20_Cleaned_L1R'
	#)
	#, pathLastFilterAcceptedOnly = cms.bool( True )   # select only trigger objects used in last filters of succeeding paths
	#, collectionTags             = cms.vstring( '*' ) # wildcard
	, maxDPtRel = cms.double( 0.5 )
	, maxDeltaR = cms.double( 0.5 )
	, resolveAmbiguities    = cms.bool( True )        # only one match per trigger object
	, resolveByMatchQuality = cms.bool( True )        # take best match found per reco object: by DeltaR here (s. above)
	)
	
	switchOnTrigger(process)
	process.patTrigger.addL1Algos = cms.bool( True )
	switchOnTrigger( process ) # to fix event content
	switchOnTriggerMatchEmbedding( process, ['patPhotonMatch'])
	#	removeCleaningFromTriggerMatching( process, outputModule = '')

def addJetMET(process,theJetNames,jetMetCorrections,mcVersion):
	# use PFMET, not the other stuff
	from PhysicsTools.PatAlgos.tools.pfTools import *
	switchToPFMET(process,input=cms.InputTag('pfMet'))

def removeMCDependence( process ):
    #-- Remove MC dependence ------------------------------------------------------
    from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
    removeMCMatching(process, ['All'])

def useDAVertices(process):
    #-- Deterministic Annealing Vertices ---------------------------------------------------
    from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi import offlinePrimaryVerticesDA
    from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDAWithBS_cfi import offlinePrimaryVerticesDA as offlinePrimaryVerticesDAWithBS
    process.offlinePrimaryVertices = offlinePrimaryVerticesDA.clone()
    process.offlinePrimaryVerticesWithBS = offlinePrimaryVerticesDAWithBS.clone()

    #-- Keep Gap Method Vertices for Comparison --------------------------------------------
    try:
        from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesGap_cfi import offlinePrimaryVerticesGap
        from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesGapWithBS_cfi import offlinePrimaryVerticesGapWithBS
    except ImportError:
        from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import offlinePrimaryVertices as offlinePrimaryVerticesGap
        from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesWithBS_cfi import offlinePrimaryVerticesWithBS as offlinePrimaryVerticesGapWithBS
        print "Could not find OfflinePrimaryVerticesGap. (Please ignore this warning if using CMSSW_4_1_X or lower)"
        
    process.offlinePrimaryVerticesGap = offlinePrimaryVerticesGap.clone()
    process.offlinePrimaryVerticesGapWithBS = offlinePrimaryVerticesGapWithBS.clone()

    process.daVertices = cms.Sequence(
        process.offlinePrimaryVertices
        + process.offlinePrimaryVerticesWithBS
        + process.offlinePrimaryVerticesGap
        + process.offlinePrimaryVerticesGapWithBS
        )

def getHZZ_pattuple_outputCommands( process ):
	from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent, patExtraAodEventContent, patTriggerEventContent, patTriggerStandAloneEventContent, patEventContentTriggerMatch
	keepList = []
    	hzzAddEventContent = [ # PAT Objects
    	
		'keep *_selectedPatJets*_*_*',           ## keep refactorized pat jet elements
		'drop patJets_selectedPatJets*_*_*',     ## drop the actual selected pat jets, they're redundant
		'drop *_selectedPatJets_pfCandidates_*', ## drop for default patJets which are CaloJets
		'drop *_*PF_caloTowers_*',               ## drop collections not needed for the corresponding jet types
		'drop *_*JPT_pfCandidates_*',            ## drop collections not needed for the corresponding jet types
		'drop *_*Calo_pfCandidates_*',           ## drop collections not needed for the corresponding jet types
		'keep *_cleanPatPhotons*_*_*',
		'keep *_cleanPatElectrons*_*_*',
		'keep *_cleanPatMuons*_*_*',
		'keep *_cleanPatTaus*_*_*',
		'keep *_cleanPatJets*_*_*',
		'keep *_patMET*_*_*',
		'keep *_patJets*_*_*',
		'keep *_cleanPatHemispheres*_*_*',
		'keep *_cleanPatPFParticles*_*_*',
		'keep *_cleanPatTrackCands*_*_*',

		'keep patTriggerAlgorithms_patTrigger_*_*',
		'keep patTriggerConditions_patTrigger_*_*',
		'keep patTriggerObjects_patTrigger_*_*',
		'keep patTriggerFilters_patTrigger_*_*',
		'keep patTriggerPaths_patTrigger_*_*',
		'keep *_patTriggerEvent_*_*'

		'keep *_eid*_*_*'
        #'keep *_triggerMatched*_*_*',         
	# Keep PF2PAT output
        'keep *_selectedPatMuonsPF_*_*',         
        'keep *_selectedPatElectronsPF_*_*',         
        'keep *_selectedPatTausPF_*_*',         
        'keep *_selectedPatJetsPF_*_*',
	#L1 trigger info         
		'keep L1GlobalTriggerObjectMapRecord_*_*_*',
        'keep L1GlobalTriggerReadoutRecord_*_*_*',
        # Generator information
        'keep recoGenJets_*GenJets*_*_*',
        'keep recoGenMETs_*_*_*',
	#Number of processed events
		'keep recoPFCandidates_particleFlow_*_*',
        #'keep *_gsfElectrons_*_*',    #Keep electron core
        'keep *_photonCore_*_*',        #Keep electron core
        'keep recoConversions_conversions_*_*',
        'keep recoTracks_*onversions_*_*',
        'keep HcalNoiseSummary_*_*_*', #Keep the one in RECO
       ] 
	keepList.extend(patEventContent)
	keepList.extend(patExtraAodEventContent)
	keepList.extend(patTriggerEventContent)
	keepList.extend(patEventContentTriggerMatch)
	keepList.extend(hzzAddEventContent)
	return keepList

