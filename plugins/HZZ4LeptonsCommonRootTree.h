#ifndef HZZ4LeptonsCommonRootTree_h
#define HZZ4LeptonsCommonRootTree_h

/** \class  HZZ4LeptonsCommonRootTree
 *
 *  Root Tree for H->ZZ->4l analysis.
 *
 *  Author: N. De Filippis - Politecnico and INFN Bari
 *          
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/DataKeyTags.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/FileBlock.h"

// MCTruth
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// Data format
#include "DataFormats/Common/interface/Handle.h" 
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"    
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
//#include "DataFormats/EgammaCandidates/interface/Photon.h"
//#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
//#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
//#include "DataFormats/METReco/interface/PFMET.h"
//#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "PhysicsTools/PatAlgos/interface/PATUserDataHelper.h"


// Trigger
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"

// Geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

// Maps
#include "DataFormats/Common/interface/ValueMap.h"

// Pileup
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

// user include files
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionBaseClass.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

// Transient tracks
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

//PU Jet ID
#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"

//Full Error
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToCartesian.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"

#include <TMatrixD.h>

//chisquare
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

class MultiTrajectoryStateMode ;
class EgammaTowerIsolation ;


// Class to create TTree variables
#include <TFile.h> 
#include <TTree.h> 

// Namespaces
using namespace reco;
using namespace std;
using namespace pat;




class HZZ4LeptonsCommonRootTree : public edm::EDAnalyzer {
  
 public:
  
  HZZ4LeptonsCommonRootTree(const edm::ParameterSet& pset);
  
  ~HZZ4LeptonsCommonRootTree();
  
  void analyze(const edm::Event& e, const edm::EventSetup& c);
  void beginJob();
  void endJob();
  
  void respondToOpenInputFile(edm::FileBlock const& fb) {
    inputfileName = fb.fileName();
    cout << "Input Filename is=" << inputfileName.c_str() << endl;
  }
	
  void ReadParameters(const edm::ParameterSet& pset){ 
    std::cout << "Reading parameters from cfg" << std::endl;
    
    typedef std::vector<edm::InputTag> vtag;
    // Get the various input parameters
    decaychannel              = pset.getParameter<std::string>("decaychannel");
    rootFileName              = pset.getUntrackedParameter<std::string>("rootFileName");
    useRECOformat             = pset.getUntrackedParameter<bool>("useRECOformat");

    // search param in cfg
    module_to_search=pset.getUntrackedParameter<vector<std::string> >("module_to_search");
    par_to_search= pset.getUntrackedParameter<std::string>("par_to_search");
    

    // Get PU simulation info
    fillPUinfo                = pset.getUntrackedParameter<bool>("fillPUinfo");
    PileupSrc_                = consumes<std::vector<PileupSummaryInfo> >(pset.getParameter<edm::InputTag>("PileupSrc"));

    // Generator
    generator_                = consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("Generator"));

    // Get HLT flags
    fillHLTinfo               = pset.getUntrackedParameter<bool>("fillHLTinfo");
    HLTInfoFired              = pset.getParameter<edm::InputTag>("HLTInfoFired");
    HLTAnalysisinst           = pset.getParameter<string>("HLTAnalysisinst");
    flagHLTnames              = pset.getParameter<vtag>("flagHLTnames");
    // Get HLT matching
//AOD    triggerEvent              = consumes<trigger::TriggerEvent >(pset.getParameter<edm::InputTag>("triggerEvent"));
    triggerObjects_              = consumes<pat::TriggerObjectStandAloneCollection>(pset.getParameter<edm::InputTag>("triggerobjects"));

    triggerBits_              = consumes<edm::TriggerResults>(pset.getParameter<edm::InputTag>("triggerbits"));
   
    triggerFilter             = pset.getParameter<std::string>("triggerFilter");
    triggerMatchObject        = consumes<edm::Association<std::vector<pat::TriggerObjectStandAlone> > >(pset.getParameter<edm::InputTag>("triggerMatchObject"));
    triggerMatchObjectEle     = pset.getParameter<edm::InputTag>("triggerMatchObjectEle");
    triggerHLTcollection      = pset.getParameter<std::string>("triggerHLTcollection");

    // Get flags
    flaginst                  = pset.getParameter<std::string>("flaginst");
    flagtags                  = pset.getParameter<std::vector<std::string> >("flagtags");
    // MCtruth tags
    fillMCTruth               = pset.getUntrackedParameter<bool>("fillMCTruth");
    MCcollName                = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("MCcollName"));
    genParticles_             = consumes<vector<reco::GenParticle> >(pset.getParameter<edm::InputTag>("genParticles"));
    fourgenleptons_           = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("fourgenleptons"));
    digenZ_                   = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("digenZ"));

    
    // RECO tags
    RECOcollNameBest2e2mu     = pset.getParameter<vtag>("RECOcollNameBest2e2mu");
    RECOcollNameBest4mu       = pset.getParameter<vtag>("RECOcollNameBest4mu");
    RECOcollNameBest4e        = pset.getParameter<vtag>("RECOcollNameBest4e");
 
    // RECO additional tags
    useAdditionalRECO         = pset.getUntrackedParameter<bool>("useAdditionalRECO");
    RECOcollNameZ             = pset.getParameter<vtag>("RECOcollNameZ");
    RECOcollNameZss           = pset.getParameter<vtag>("RECOcollNameZss");
    RECOcollNameDiLep_         = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameDiLep"));

    RECOcollNameDiLep         = pset.getParameter<edm::InputTag>("RECOcollNameDiLep");
    RECOcollNameMMMM_         = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameMMMM"));
    RECOcollNameEEEE          = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameEEEE")); 
    RECOcollNameEEMM          = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameEEMM"));
    RECOcollNameLLLLss        = pset.getParameter<vtag>("RECOcollNameLLLLss");
    RECOcollNameLLLLssos      = pset.getParameter<vtag>("RECOcollNameLLLLssos");
    RECOcollNameLLL           = pset.getParameter<vtag>("RECOcollNameLLL");
    RECOcollNameLLLl          = pset.getParameter<vtag>("RECOcollNameLLLl");
    RECOcollNameLLLL         = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollNameLLLL"));   

//    RECOcollNameMMMM          = pset.getParameter<vtag>("RECOcollNameMMMM");
//    RECOcollNameEEEE          = pset.getParameter<vtag>("RECOcollNameEEEE");
//    RECOcollNameEEMM          = pset.getParameter<vtag>("RECOcollNameEEMM");
    RECOcollNameLLLLss        = pset.getParameter<vtag>("RECOcollNameLLLLss");
    RECOcollNameLLLLssos      = pset.getParameter<vtag>("RECOcollNameLLLLssos");
    RECOcollNameLLL           = pset.getParameter<vtag>("RECOcollNameLLL");
    RECOcollNameLLLl          = pset.getParameter<vtag>("RECOcollNameLLLl");
//    RECOcollNameLLLL          = pset.getParameter<edm::InputTag>("RECOcollNameLLLL");

    // electrons and muons tags
    use2011EA                 = pset.getUntrackedParameter<bool>("use2011EA");

 
    muonTag_                  = consumes<edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("MuonsLabel"));
    muonCorrPtErrorMapTag_    = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsCorrPtErrorMapLabel"));
    
    muonPFTag_                  = consumes<edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("PFMuonsLabel"));
/*
    muonPFIsoValueChargedAllTag_= consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValueChargedAll"));
    muonPFIsoValueChargedTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValueCharged"));
    muonPFIsoValueNeutralTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValueNeutral"));
    muonPFIsoValueGammaTag_     = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValueGamma"));
    muonPFIsoValuePUTag_        = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValuePU"));
 */   
    
    clusterCollectionTag_    = pset.getParameter<edm::InputTag>("SuperClustersLabel");
    gsftrackCollection_      = pset.getParameter<edm::InputTag>("GsfTracksElectronsLabel");

    electronEgmTag_          = consumes<edm::View<pat::Electron> >(pset.getParameter<edm::InputTag>("ElectronsEgmLabel"));
    electronEgmTkMapTag_     = pset.getParameter<edm::InputTag>("ElectronsEgmTkMapLabel");
    electronEgmEcalMapTag_   = pset.getParameter<edm::InputTag>("ElectronsEgmEcalMapLabel");
    electronEgmHcalMapTag_   = pset.getParameter<edm::InputTag>("ElectronsEgmHcalMapLabel");
    
    mvaElectronTag_          = consumes<edm::View<pat::Electron> >(pset.getParameter<edm::InputTag>("mvaElectronTag"));
    mvaTrigV0MapTag_         = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("mvaTrigV0MapTag"));
    mvaNonTrigV0MapTag_      = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("mvaNonTrigV0MapTag"));
/*    
    electronPFIsoValueChargedAllTag_= consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueChargedAll"));
    electronPFIsoValueChargedTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueCharged"));
    electronPFIsoValueNeutralTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueNeutral"));
    electronPFIsoValueGammaTag_     = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueGamma"));
    electronPFIsoValuePUTag_        = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValuePU"));
  */  
//    eleRegressionEnergyErrorTag_= pset.getParameter<edm::InputTag>("eleRegressionEnergyErrorLabel");
//    eleRegressionEnergyTag_     = pset.getParameter<edm::InputTag>("eleRegressionEnergyLabel");
    
    
    // PF photons
//    pfphotonsTag_                 = consumes<edm::View<reco::PFCandidate>>(pset.getParameter<edm::InputTag>("PFPhotonsLabel"));
    pfTag_                 = consumes<edm::View<pat::PackedCandidate>>(pset.getParameter<edm::InputTag>("pfCands"));
/*    photonPFIsoValueChargedAllTag_= consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValueChargedAll"));
    photonPFIsoValueChargedTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValueCharged"));
    photonPFIsoValueNeutralTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValueNeutral"));
    photonPFIsoValueGammaTag_     = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValueGamma"));
    photonPFIsoValuePUTag_        = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValuePU"));
*/   

    // vertexing 
    // 3D w.r.t primary vertex DA
    muonTag_Vert           = pset.getParameter<edm::InputTag>("MuonsLabelVert");
    muonMapTag_Vert        = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVert"));
    muonMapTag_VertValue   = consumes<edm::ValueMap<float> >( pset.getParameter<edm::InputTag>("MuonsMapLabelVertValue"));
    muonMapTag_VertError   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVertError"));
    // KF
    muonMapTag_VertKF        = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVertKF"));
    muonMapTag_VertValueKF   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVertValueKF"));
    muonMapTag_VertErrorKF   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsMapLabelVertErrorKF"));
    
    
    // GD, Std and Kin vertex
    muonMapTag_VertGD      = pset.getParameter<edm::InputTag>("MuonsMapLabelVertGD");
    muonMapTag_VertGDMMMM  = pset.getParameter<edm::InputTag>("MuonsMapLabelVertGDMMMM");
    muonMapTag_VertStd     = pset.getParameter<edm::InputTag>("MuonsMapLabelVertStd");
    muonMapTag_VertStdMMMM = pset.getParameter<edm::InputTag>("MuonsMapLabelVertStdMMMM");
    muonMapTag_VertKin     = pset.getParameter<edm::InputTag>("MuonsMapLabelVertKin");
    muonMapTag_VertKinMMMM = pset.getParameter<edm::InputTag>("MuonsMapLabelVertKinMMMM");
    
    // STIP SLIP
    muonSTIPMapTag_Vert   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSTIPMapLabelVert"));
    muonSLIPMapTag_Vert   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSLIPMapLabelVert"));
    
    muonSTIPMapTag_VertValue   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSTIPMapLabelVertValue"));
    muonSLIPMapTag_VertValue   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSLIPMapLabelVertValue"));
    muonSTIPMapTag_VertError   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsSTIPMapLabelVertError"));
    muonSLIPMapTag_VertError   =consumes<edm::ValueMap<float> >( pset.getParameter<edm::InputTag>("MuonsSLIPMapLabelVertError"));
    
    // 3D w.r.t primary vertex DA
    electronTag_Vert          = pset.getParameter<edm::InputTag>("ElectronsLabelVert");
    electronMapTag_Vert       = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ElectronsMapLabelVert"));
    electronMapTag_VertValue  = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ElectronsMapLabelVertValue"));
    electronMapTag_VertError  = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ElectronsMapLabelVertError"));
    // KF
    electronMapTag_VertKF       = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ElectronsMapLabelVertKF"));
    electronMapTag_VertValueKF  = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ElectronsMapLabelVertValueKF"));
    electronMapTag_VertErrorKF  = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ElectronsMapLabelVertErrorKF"));
    
    
    // GD, Std and Kin vertex
    electronMapTag_VertGD      = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertGD");
    electronMapTag_VertGDEEEE  = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertGDEEEE");
    electronMapTag_VertStd     = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertStd");
    electronMapTag_VertStdEEEE = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertStdEEEE");
    electronMapTag_VertKin     = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertKin");
    electronMapTag_VertKinEEEE = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertKinEEEE");
    
    // STIP SLIP
    electronSTIPMapTag_Vert   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ElectronsSTIPMapLabelVert"));
    electronSLIPMapTag_Vert   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ElectronsSLIPMapLabelVert"));
    
    electronSTIPMapTag_VertValue   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ElectronsSTIPMapLabelVertValue"));
    electronSLIPMapTag_VertValue   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ElectronsSLIPMapLabelVertValue"));
    electronSTIPMapTag_VertError   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ElectronsSTIPMapLabelVertError"));
    electronSLIPMapTag_VertError   = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ElectronsSLIPMapLabelVertError"));
   

    
    // geom. discr.
    ftsigma_Vert           = pset.getParameter<edm::InputTag>("ftsigmaVert");
    ftsigmalag_Vert        = pset.getParameter<edm::InputTag>("ftsigmalagVert");
    gdX_Vert               = pset.getParameter<edm::InputTag>("gdX_Vert");
    gdY_Vert               = pset.getParameter<edm::InputTag>("gdY_Vert");
    gdZ_Vert               = pset.getParameter<edm::InputTag>("gdZ_Vert");
    gdlagX_Vert            = pset.getParameter<edm::InputTag>("gdlagX_Vert");
    gdlagY_Vert            = pset.getParameter<edm::InputTag>("gdlagY_Vert");
    gdlagZ_Vert            = pset.getParameter<edm::InputTag>("gdlagZ_Vert");
    gdlagProb_Vert         = pset.getParameter<edm::InputTag>("gdlagProb_Vert");
    gdlagNdof_Vert         = pset.getParameter<edm::InputTag>("gdlagNdof_Vert");

    ftsigma_VertMMMM       = pset.getParameter<edm::InputTag>("ftsigmaVertMMMM");
    ftsigmalag_VertMMMM    = pset.getParameter<edm::InputTag>("ftsigmalagVertMMMM");
    gdX_VertMMMM           = pset.getParameter<edm::InputTag>("gdX_VertMMMM");
    gdY_VertMMMM           = pset.getParameter<edm::InputTag>("gdY_VertMMMM");
    gdZ_VertMMMM           = pset.getParameter<edm::InputTag>("gdZ_VertMMMM");
    gdlagX_VertMMMM        = pset.getParameter<edm::InputTag>("gdlagX_VertMMMM");
    gdlagY_VertMMMM        = pset.getParameter<edm::InputTag>("gdlagY_VertMMMM");
    gdlagZ_VertMMMM        = pset.getParameter<edm::InputTag>("gdlagZ_VertMMMM");
    gdlagProb_VertMMMM     = pset.getParameter<edm::InputTag>("gdlagProb_VertMMMM");
    gdlagNdof_VertMMMM     = pset.getParameter<edm::InputTag>("gdlagNdof_VertMMMM");

    ftsigma_VertEEEE       = pset.getParameter<edm::InputTag>("ftsigmaVertEEEE");
    ftsigmalag_VertEEEE    = pset.getParameter<edm::InputTag>("ftsigmalagVertEEEE");
    gdX_VertEEEE           = pset.getParameter<edm::InputTag>("gdX_VertEEEE");
    gdY_VertEEEE           = pset.getParameter<edm::InputTag>("gdY_VertEEEE");
    gdZ_VertEEEE           = pset.getParameter<edm::InputTag>("gdZ_VertEEEE");
    gdlagX_VertEEEE        = pset.getParameter<edm::InputTag>("gdlagX_VertEEEE");
    gdlagY_VertEEEE        = pset.getParameter<edm::InputTag>("gdlagY_VertEEEE");
    gdlagZ_VertEEEE        = pset.getParameter<edm::InputTag>("gdlagZ_VertEEEE");
    gdlagProb_VertEEEE     = pset.getParameter<edm::InputTag>("gdlagProb_VertEEEE");
    gdlagNdof_VertEEEE     = pset.getParameter<edm::InputTag>("gdlagNdof_VertEEEE");

    //ConstraintVertex 4l
    StandardFitVertex        = pset.getParameter<edm::InputTag>("StandardFitVertex");
    StandardFitVertexMMMM    = pset.getParameter<edm::InputTag>("StandardFitVertexMMMM");
    StandardFitVertexEEEE    = pset.getParameter<edm::InputTag>("StandardFitVertexEEEE");
    KinematicFitVertex       = pset.getParameter<edm::InputTag>("KinematicFitVertex");
    KinematicFitVertexMMMM   = pset.getParameter<edm::InputTag>("KinematicFitVertexMMMM");
    KinematicFitVertexEEEE   = pset.getParameter<edm::InputTag>("KinematicFitVertexEEEE");
    RefittedMass             = pset.getParameter<edm::InputTag>("RefittedMass");
    RefittedMassMMMM         = pset.getParameter<edm::InputTag>("RefittedMassMMMM");
    RefittedMassEEEE         = pset.getParameter<edm::InputTag>("RefittedMassEEEE");

    //ConstraintVertex 3l
    StandardFitVertexMMM     = pset.getParameter<edm::InputTag>("StandardFitVertexMMM");
    StandardFitVertexMME     = pset.getParameter<edm::InputTag>("StandardFitVertexMME");
    StandardFitVertexEEE     = pset.getParameter<edm::InputTag>("StandardFitVertexEEE");
    StandardFitVertexMEE     = pset.getParameter<edm::InputTag>("StandardFitVertexMEE");

    // ConstraintVertex Dileptons
    StandardFitVertexDiLep   = pset.getParameter<edm::InputTag>("StandardFitVertexDiLep");

    //electronID
    eleIDTag_                 = pset.getParameter<vtag>("eleIDLabel");

    // Other objets
    photonsTag_              = consumes<edm::View<pat::Photon> >(pset.getParameter<edm::InputTag>("PhotonsLabel"));
//    tracksTag_               = consumes<vector<reco::Track> >(pset.getParameter<edm::InputTag>("TracksLabel"));
    jetsTag_                 = consumes<vector<pat::Jet> >(pset.getParameter<edm::InputTag>("JetsLabel"));
 //   jetsDataTag_             = consumes<vector<pat::Jet> >(pset.getParameter<edm::InputTag>("JetsDataLabel"));
    jetsMVATag_              = consumes<vector<pat::Jet> >(pset.getParameter<edm::InputTag>("JetsMVALabel"));
//    PuJetMvaMCfullDiscr_     = consumes<edm::ValueMap<float> > (pset.getParameter<edm::InputTag>("PuJetMvaMCfullDiscrLabel"));
//    PuJetMvaMCfullId_        = consumes<edm::ValueMap<float> > (pset.getParameter<edm::InputTag>("PuJetMvaMCfullIdLabel"));
//    PuJetMvaDatafullDiscr_   = consumes<edm::ValueMap<float> > (pset.getParameter<edm::InputTag>("PuJetMvaDatafullDiscrLabel"));
//    PuJetMvaDatafullId_      = consumes<edm::ValueMap<float> > (pset.getParameter<edm::InputTag>("PuJetMvaDatafullIdLabel"));
 
    rhojetsTag_              = consumes<double>(pset.getParameter<edm::InputTag>("RhoJetsLabel"));
    verticesTag_             = consumes<vector<reco::Vertex> >(pset.getParameter<edm::InputTag>("VerticesLabel"));

    // GenJet 
    genjetTag_              = consumes<vector<reco::GenJet> >(pset.getParameter<edm::InputTag>("GenJetLabel")); // GenJet

/*    
    // MET reco
    genmetTag_              = consumes<vector<reco::GenMET> >(pset.getParameter<edm::InputTag>("GenMETLabel")); // GenMET
    trackermetTag_          = pset.getParameter<edm::InputTag>("TrackerMETLabel"); // TrackerMET
    // Calo MET
    calometTag_             = pset.getParameter<edm::InputTag>("CaloMETLabel");
    calometnohfTag_         = pset.getParameter<edm::InputTag>("CaloMET_NoHFLabel");
    useAdditionalMET_       = pset.getUntrackedParameter<bool>("useAdditionalMET");
    calomethoTag_           = pset.getParameter<edm::InputTag>("CaloMET_HOLabel");
    calometoptTag_          = pset.getParameter<edm::InputTag>("CaloMET_OptLabel");
    calometoptnohfTag_      = pset.getParameter<edm::InputTag>("CaloMET_OptNoHFLabel");
    calometoptnohfhoTag_    = pset.getParameter<edm::InputTag>("CaloMET_OptNoHFHOLabel");
    calometopthoTag_        = pset.getParameter<edm::InputTag>("CaloMET_OptHOLabel");
    calometnohfhoTag_       = pset.getParameter<edm::InputTag>("CaloMET_NoHFHOLabel");
*/    
    // PF MET
    pfmetTag_               = consumes<vector<pat::MET> >(pset.getParameter<edm::InputTag>("PfMETLabel"));
/*
    // htMET
    htmetic5Tag_            = pset.getParameter<edm::InputTag>("HtMET_IC5Label");
    htmetkt4Tag_            = pset.getParameter<edm::InputTag>("HtMET_KT4Label");
    htmetkt6Tag_            = pset.getParameter<edm::InputTag>("HtMET_KT6Label");
    htmetsc5Tag_            = pset.getParameter<edm::InputTag>("HtMET_SC5Label");
    htmetsc7Tag_            = pset.getParameter<edm::InputTag>("HtMET_SC7Label");
    // JES correction MET
    jescormetic5Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorIC5CaloJetLabel");
    jescormetkt4Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorKT4CaloJetLabel");
    jescormetkt6Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorKT6CaloJetLabel");
    jescormetsc5Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorSC5CaloJetLabel");
    jescormetsc7Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorSC7CaloJetLabel");
    // Type I Muon Correction MET
    cormetMuTag_            = pset.getParameter<edm::InputTag>("CorMETGlobalMuLabel");
*/    

    //quark/gluon tagging
    //QGTaggerMC_LHTag_        = consumes<edm::ValueMap<float> > (pset.getParameter<edm::InputTag>("QGTaggerMC_LHLabel"));
    //QGTaggerData_LHTag_      = consumes<edm::ValueMap<float> > (pset.getParameter<edm::InputTag>("QGTaggerData_LHLabel"));
    
    // btagging
    //   bDiscriminators_         =pset.getParameter<std::vector<std::string> >("bDiscriminators");

   tCHighEff_bTag_         =pset.getParameter<std::string>("tCHighEff_bTagLabel");
   tCHighPur_bTag_         =pset.getParameter<std::string>("tCHighPur_bTagLabel");
   cSV_bTag_         =pset.getParameter<std::string>("cSV_bTagLabel");
//    tCHighPur_bTag_         = consumes<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,std::vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference> >(pset.getParameter<edm::InputTag>("tCHighPur_bTagLabel"));
//    cSV_bTag_               = consumes<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,std::vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference> >(pset.getParameter<edm::InputTag>("cSV_bTagLabel"));

    // Conversion finder
//    ConvMapDistTag_       = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ConvMapDist"));
//    ConvMapDcotTag_       = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ConvMapDcot"));
    
    // Matching
/*
    goodElectronMCMatch_  = consumes<edm::Association<std::vector<reco::GenParticle> >>(pset.getParameter<edm::InputTag>("goodElectronMCMatch"));
    myElectrons_          = consumes<reco::CandidateCollection>(pset.getParameter<edm::InputTag>("myElectrons"));
    goodMuonMCMatch_      = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodMuonMCMatch"));
    myMuons_              = consumes<reco::CandidateCollection>(pset.getParameter<edm::InputTag>("myMuons"));
    goodGammaMCMatch_     = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodGammaMCMatch"));
    myGammas_             = consumes<reco::CandidateCollection>(pset.getParameter<edm::InputTag>("myGammas"));
*/

    // Beam Spot
    offlineBeamSpot_      = consumes<reco::BeamSpot>(pset.getParameter<edm::InputTag>("offlineBeamSpot"));
  }
  
  
  void DefineBranches(TTree *Tree_){
    // Run event
    Tree_->Branch("Run",&irun,"irun/i");
    Tree_->Branch("Event",&ievt,"ievt/i");
    Tree_->Branch("LumiSection",&ils,"ils/i");
    Tree_->Branch("Avginstlumi",&Avginstlumi,"Avginstlumi/F");
    
    // MC
    Tree_->Branch("num_PU_vertices",&num_PU_vertices,"num_PU_vertices/I");
    Tree_->Branch("PU_BunchCrossing",&PU_BunchCrossing,"PU_BunchCrossing/I");
    Tree_->Branch("MC_weighting",&MC_weighting,"MC_weighting/F");

    // HLT 
    Tree_->Branch("RECO_nMuHLTMatch",&RECO_nMuHLTMatch,"RECO_nMuHLTMatch/I");
    Tree_->Branch("RECOMU_PT_MuHLTMatch",RECOMU_PT_MuHLTMatch,"RECOMU_PT_MuHLTMatch[100]/F");
    Tree_->Branch("RECOMU_ETA_MuHLTMatch",RECOMU_ETA_MuHLTMatch,"RECOMU_ETA_MuHLTMatch[100]/F");
    Tree_->Branch("RECOMU_PHI_MuHLTMatch",RECOMU_PHI_MuHLTMatch,"RECOMU_PHI_MuHLTMatch[100]/F");

    Tree_->Branch("RECO_nEleHLTMatch",&RECO_nEleHLTMatch,"RECO_nEleHLTMatch/I");
    Tree_->Branch("RECOELE_PT_EleHLTMatch",RECOELE_PT_EleHLTMatch,"RECOELE_PT_EleHLTMatch[100]/F");
    Tree_->Branch("RECOELE_ETA_EleHLTMatch",RECOELE_ETA_EleHLTMatch,"RECOELE_ETA_EleHLTMatch[100]/F");
    Tree_->Branch("RECOELE_PHI_EleHLTMatch",RECOELE_PHI_EleHLTMatch,"RECOELE_PHI_EleHLTMatch[100]/F");

    Tree_->Branch("HLTPathsFired",HLTPathsFired,"HLTPathsFired/C");

   
    // MC block 
    Tree_->Branch("MC_E", MC_E, "MC_E[7]/F"); 
    Tree_->Branch("MC_PT", MC_PT, "MC_PT[7]/F"); 
    Tree_->Branch("MC_ETA", MC_ETA, "MC_ETA[7]/F"); 
    Tree_->Branch("MC_THETA", MC_THETA, "MC_THETA[7]/F");
    Tree_->Branch("MC_PHI", MC_PHI, "MC_PHI[7]/F");
    Tree_->Branch("MC_MASS", MC_MASS, "MC_MASS[7]/F");
    Tree_->Branch("MC_PDGID", &MC_PDGID, "MC_PDGID[7]/F");
   

    //  Leptons stable ordered in PT
    Tree_->Branch("MC_LEPT_PT", MC_LEPT_PT, "MC_LEPT_PT[4]/F"); 
    Tree_->Branch("MC_LEPT_ETA", MC_LEPT_ETA, "MC_LEPT_ETA[4]/F"); 
    Tree_->Branch("MC_LEPT_PHI", MC_LEPT_PHI, "MC_LEPT_PHI[4]/F"); 
    Tree_->Branch("MC_LEPT_THETA", MC_LEPT_THETA, "MC_LEPT_THETA[4]/F"); 
    Tree_->Branch("MC_LEPT_PDGID", MC_LEPT_PDGID, "MC_LEPT_PDGID[4]/F");

    // MC Z1, Z2 (first index) and daughter leptons and photon (second index, order: L1, L2, P1, P2 with PT ordering of the leptons and P1 associated to L1)
    Tree_->Branch("MC_Z_PT", MC_Z_PT, "MC_Z_PT[2][5]/F");
    Tree_->Branch("MC_Z_ETA", MC_Z_ETA, "MC_Z_ETA[2][5]/F");
    Tree_->Branch("MC_Z_PHI", MC_Z_PHI, "MC_Z_PHI[2][5]/F");
    Tree_->Branch("MC_Z_THETA", MC_Z_THETA, "MC_Z_THETA[2][5]/F");
    Tree_->Branch("MC_Z_MASS", MC_Z_MASS, "MC_Z_MASS[2][5]/F");
    Tree_->Branch("MC_Z_PDGID", MC_Z_PDGID, "MC_Z_PDGID[2][5]/F");

    // 4l from stable particles
    Tree_->Branch("MC_fourl_MASS",MC_fourl_MASS,"MC_fourl_MASS[50][5]/F");
    Tree_->Branch("MC_fourl_PT",MC_fourl_PT,"MC_fourl_PT[50][5]/F");
    Tree_->Branch("MC_fourl_PDGID",MC_fourl_PDGID,"MC_fourl_PDGID[50][5]/F");

    // diZ
    Tree_->Branch("MC_ZZ_MASS",MC_ZZ_MASS,"MC_ZZ_MASS[4][7]/F");
    Tree_->Branch("MC_ZZ_PT",MC_ZZ_PT,"MC_ZZ_PT[4][7]/F");
    Tree_->Branch("MC_ZZ_ETA",MC_ZZ_ETA,"MC_ZZ_ETA[4][7]/F");
    Tree_->Branch("MC_ZZ_PHI",MC_ZZ_PHI,"MC_ZZ_PHI[4][7]/F");
    Tree_->Branch("MC_ZZ_THETA",MC_ZZ_THETA,"MC_ZZ_THETA[4][7]/F");
    Tree_->Branch("MC_ZZ_PDGID",MC_ZZ_PDGID,"MC_ZZ_PDGID[4][7]/F");


    // GenJet
    Tree_->Branch( "MC_GENJET_PT",  MC_GENJET_PT,  "MC_GENJET_PT[100]/F");
    Tree_->Branch( "MC_GENJET_ETA", MC_GENJET_ETA, "MC_GENJET_ETA[100]/F");
    Tree_->Branch( "MC_GENJET_PHI", MC_GENJET_PHI, "MC_GENJET_PHI[100]/F");

    // GenMET  
    Tree_->Branch("MC_GENMET", &genmet, "MC_GENMET/F");
    
  
     // RECORF block 2e2mu
    
    Tree_->Branch("RECORF_2e2mu_cosTheta1_spin",RECORF_2e2mu_cosTheta1_spin,"RECORF_2e2mu_cosTheta1_spin[100]/D");
    Tree_->Branch("RECORF_2e2mu_cosTheta2_spin",RECORF_2e2mu_cosTheta2_spin,"RECORF_2e2mu_cosTheta2_spin[100]/D");
    Tree_->Branch("RECORF_2e2mu_cosThetaStar_spin",RECORF_2e2mu_cosThetaStar_spin,"RECORF_2e2mu_cosThetaStar_spin[100]/D");
    Tree_->Branch("RECORF_2e2mu_Phi_spin",RECORF_2e2mu_Phi_spin,"RECORF_2e2mu_Phi_spin[100]/D");
    Tree_->Branch("RECORF_2e2mu_Phi1_spin",RECORF_2e2mu_Phi1_spin,"RECORF_2e2mu_Phi1_spin[100]/D");
    Tree_->Branch("RECORF_2e2mu_Phi2_spin",RECORF_2e2mu_Phi2_spin,"RECORF_2e2mu_Phi2_spin[100]/D");
    Tree_->Branch("RECORF_2e2mu_phi1RF_spin",RECORF_2e2mu_phi1RF_spin,"RECORF_2e2mu_phi1RF_spin[100]/D");
    Tree_->Branch("RECORF_2e2mu_phi2RF_spin",RECORF_2e2mu_phi2RF_spin,"RECORF_2e2mu_phi2RF_spin[100]/D");
    Tree_->Branch("RECORF_2e2mu_MELA",RECORF_2e2mu_MELA,"RECORF_2e2mu_MELA[100]/D");

   
    Tree_->Branch("RECORF_4e_cosTheta1_spin", RECORF_4e_cosTheta1_spin,"RECORF_4e_cosTheta1_spin[100]/D");
    Tree_->Branch("RECORF_4e_cosTheta2_spin", RECORF_4e_cosTheta2_spin,"RECORF_4e_cosTheta2_spin[100]/D");
    Tree_->Branch("RECORF_4e_cosThetaStar_spin",RECORF_4e_cosThetaStar_spin,"RECORF_4e_cosThetaStar_spin[100]/D");
    Tree_->Branch("RECORF_4e_Phi_spin", RECORF_4e_Phi_spin,"RECORF_4e_Phi_spin[100]/D");
    Tree_->Branch("RECORF_4e_Phi1_spin", RECORF_4e_Phi1_spin,"RECORF_4e_Phi1_spin[100]/D");
    Tree_->Branch("RECORF_4e_Phi2_spin", RECORF_4e_Phi2_spin,"RECORF_4e_Phi2_spin[100]/D");
    Tree_->Branch("RECORF_4e_phi1RF_spin", RECORF_4e_phi1RF_spin,"RECORF_4e_phi1RF_spin[100]/D");
    Tree_->Branch("RECORF_4e_phi2RF_spin", RECORF_4e_phi2RF_spin,"RECORF_4e_phi2RF_spin[100]/D");
    Tree_->Branch("RECORF_4e_MELA", RECORF_4e_MELA,"RECORF_4e_MELA[100]/D");
     
    Tree_->Branch("RECORF_4mu_cosTheta1_spin", RECORF_4mu_cosTheta1_spin,"RECORF_4mu_cosTheta1_spin[100]/D");
    Tree_->Branch("RECORF_4mu_cosTheta2_spin", RECORF_4mu_cosTheta2_spin,"RECORF_4mu_cosTheta2_spin[100]/D");
    Tree_->Branch("RECORF_4mu_cosThetaStar_spin", RECORF_4mu_cosThetaStar_spin,"RECORF_4mu_cosThetaStar_spin[100]/D");
    Tree_->Branch("RECORF_4mu_Phi_spin", RECORF_4mu_Phi_spin,"RECORF_4mu_Phi_spin[100]/D");
    Tree_->Branch("RECORF_4mu_Phi1_spin", RECORF_4mu_Phi1_spin,"RECORF_4mu_Phi1_spin[100]/D");
    Tree_->Branch("RECORF_4mu_Phi2_spin", RECORF_4mu_Phi2_spin,"RECORF_4mu_Phi2_spin[100]/D");
    Tree_->Branch("RECORF_4mu_phi1RF_spin", RECORF_4mu_phi1RF_spin,"RECORF_4mu_phi1RF_spin[100]/D");
    Tree_->Branch("RECORF_4mu_phi2RF_spin", RECORF_4mu_phi2RF_spin,"RECORF_4mu_phi2RF_spin[100]/D");
    Tree_->Branch("RECORF_4mu_MELA", RECORF_4mu_MELA,"RECORF_4mu_MELA[100]/D");



    // RECO additional block for reconstructed higgs, Z and their daughters
    Tree_->Branch("RECO_ZMM_MASS", RECO_ZMM_MASS, "RECO_ZMM_MASS[50]/F");
    Tree_->Branch("RECO_ZEE_MASS", RECO_ZEE_MASS, "RECO_ZEE_MASS[50]/F");
    Tree_->Branch("RECO_DiLep_MASS", RECO_DiLep_MASS, "RECO_DiLep_MASS[50]/F");
    Tree_->Branch("RECO_ZMM_PT", RECO_ZMM_PT, "RECO_ZMM_PT[3][50]/F");
    Tree_->Branch("RECO_ZEE_PT", RECO_ZEE_PT, "RECO_ZEE_PT[3][50]/F");  
    Tree_->Branch("RECO_DiLep_PT", RECO_DiLep_PT, "RECO_DiLep_PT[3][50]/F");  
    Tree_->Branch("RECO_ZMM_ETA", RECO_ZMM_ETA, "RECO_ZMM_ETA[3][50]/F");
    Tree_->Branch("RECO_ZEE_ETA", RECO_ZEE_ETA, "RECO_ZEE_ETA[3][50]/F");
    Tree_->Branch("RECO_DiLep_ETA", RECO_DiLep_ETA, "RECO_DiLep_ETA[3][50]/F");  
    Tree_->Branch("RECO_ZMM_PHI", RECO_ZMM_PHI, "RECO_ZMM_PHI[3][50]/F");
    Tree_->Branch("RECO_ZEE_PHI", RECO_ZEE_PHI, "RECO_ZEE_PHI[3][50]/F");
    Tree_->Branch("RECO_DiLep_PHI", RECO_DiLep_PHI, "RECO_DiLep_PHI[3][50]/F");  
                              
    Tree_->Branch("RECO_ZMMss_MASS", RECO_ZMMss_MASS, "RECO_ZMMss_MASS[50]/F");
    Tree_->Branch("RECO_ZEEss_MASS", RECO_ZEEss_MASS, "RECO_ZEEss_MASS[50]/F");
    Tree_->Branch("RECO_ZEM_MASS", RECO_ZEM_MASS, "RECO_ZEM_MASS[50]/F");
    Tree_->Branch("RECO_ZMMss_PT", RECO_ZMMss_PT, "RECO_ZMMss_PT[3][50]/F");
    Tree_->Branch("RECO_ZEEss_PT", RECO_ZEEss_PT, "RECO_ZEEss_PT[3][50]/F");
    Tree_->Branch("RECO_ZEM_PT", RECO_ZEM_PT, "RECO_ZEM_PT[3][50]/F");
    Tree_->Branch("RECO_ZMMss_ETA", RECO_ZMMss_ETA, "RECO_ZMMss_ETA[3][50]/F");
    Tree_->Branch("RECO_ZEEss_ETA", RECO_ZEEss_ETA, "RECO_ZEEss_ETA[3][50]/F");
    Tree_->Branch("RECO_ZEM_ETA", RECO_ZEM_ETA, "RECO_ZEM_ETA[3][50]/F");
    Tree_->Branch("RECO_ZMMss_PHI", RECO_ZMMss_PHI, "RECO_ZMMss_PHI[3][50]/F");
    Tree_->Branch("RECO_ZEEss_PHI", RECO_ZEEss_PHI, "RECO_ZEEss_PHI[3][50]/F");
    Tree_->Branch("RECO_ZEM_PHI", RECO_ZEM_PHI, "RECO_ZEM_PHI[3][50]/F");

    
    Tree_->Branch("RECO_MMMM_MASS", RECO_MMMM_MASS, "RECO_MMMM_MASS[7][100]/F");
    Tree_->Branch("RECO_MMMM_PT", RECO_MMMM_PT, "RECO_MMMM_PT[7][100]/F");
    Tree_->Branch("RECO_MMMM_ETA", RECO_MMMM_ETA, "RECO_MMMM_ETA[7][100]/F");
    Tree_->Branch("RECO_MMMM_PHI", RECO_MMMM_PHI, "RECO_MMMM_PHI[7][100]/F");
    Tree_->Branch("RECO_MMMM_MASS_REFIT", RECO_MMMM_MASS_REFIT, "RECO_MMMM_MASS_REFIT[100]/F");

    Tree_->Branch("RECO_EEEE_MASS", RECO_EEEE_MASS, "RECO_EEEE_MASS[7][100]/F");
    Tree_->Branch("RECO_EEEE_PT", RECO_EEEE_PT, "RECO_EEEE_PT[7][100]/F"); 
    Tree_->Branch("RECO_EEEE_ETA", RECO_EEEE_ETA, "RECO_EEEE_ETA[7][100]/F");
    Tree_->Branch("RECO_EEEE_PHI", RECO_EEEE_PHI, "RECO_EEEE_PHI[7][100]/F");
    Tree_->Branch("RECO_EEEE_MASS_REFIT", RECO_EEEE_MASS_REFIT, "RECO_EEEE_MASS_REFIT[100]/F");

    Tree_->Branch("RECO_EEMM_MASS", RECO_EEMM_MASS, "RECO_EEMM_MASS[7][100]/F");
    Tree_->Branch("RECO_EEMM_PT", RECO_EEMM_PT, "RECO_EEMM_PT[7][100]/F");
    Tree_->Branch("RECO_EEMM_ETA", RECO_EEMM_ETA, "RECO_EEMM_ETA[7][100]/F");
    Tree_->Branch("RECO_EEMM_PHI", RECO_EEMM_PHI, "RECO_EEMM_PHI[7][100]/F");
    Tree_->Branch("RECO_EEMM_MASS_REFIT", RECO_EEMM_MASS_REFIT, "RECO_EEMM_MASS_REFIT[100]/F");

    Tree_->Branch("RECO_LLL0_MASS", RECO_LLL0_MASS, "RECO_LLL0_MASS[50]/F"); 
    Tree_->Branch("RECO_LLL1_MASS", RECO_LLL1_MASS, "RECO_LLL1_MASS[50]/F"); 
    Tree_->Branch("RECO_LLL2_MASS", RECO_LLL2_MASS, "RECO_LLL2_MASS[50]/F"); 
    Tree_->Branch("RECO_LLL3_MASS", RECO_LLL3_MASS, "RECO_LLL3_MASS[50]/F"); 
    Tree_->Branch("RECO_LLL0_PT", RECO_LLL0_PT, "RECO_LLL0_PT[4][50]/F"); 
    Tree_->Branch("RECO_LLL1_PT", RECO_LLL1_PT, "RECO_LLL1_PT[4][50]/F"); 
    Tree_->Branch("RECO_LLL2_PT", RECO_LLL2_PT, "RECO_LLL2_PT[4][50]/F"); 
    Tree_->Branch("RECO_LLL3_PT", RECO_LLL3_PT, "RECO_LLL3_PT[4][50]/F"); 

    Tree_->Branch("RECO_LLLl0_MASS", RECO_LLLl0_MASS, "RECO_LLLl0_MASS[20]/F"); 
    Tree_->Branch("RECO_LLLl1_MASS", RECO_LLLl1_MASS, "RECO_LLLl1_MASS[20]/F"); 
    Tree_->Branch("RECO_LLLl0_PT", RECO_LLLl0_PT, "RECO_LLLl0_PT[5][20]/F"); 
    Tree_->Branch("RECO_LLLl1_PT", RECO_LLLl1_PT, "RECO_LLLl1_PT[5][20]/F"); 

    Tree_->Branch("RECO_LLLL0ss_MASS", RECO_LLLL0ss_MASS, "RECO_LLLL0ss_MASS[20]/F"); 
    Tree_->Branch("RECO_LLLL0ss_PT", RECO_LLLL0ss_PT, "RECO_LLLL0ss_PT[5][20]/F"); 
    Tree_->Branch("RECO_LLLL1ss_MASS", RECO_LLLL1ss_MASS, "RECO_LLLL1ss_MASS[20]/F"); 
    Tree_->Branch("RECO_LLLL1ss_PT", RECO_LLLL1ss_PT, "RECO_LLLL1ss_PT[5][20]/F"); 
    Tree_->Branch("RECO_LLLL2ss_MASS", RECO_LLLL2ss_MASS, "RECO_LLLL2ss_MASS[20]/F"); 
    Tree_->Branch("RECO_LLLL2ss_PT", RECO_LLLL2ss_PT, "RECO_LLLL2ss_PT[5][20]/F"); 
       
    //Tree_->Branch("RECOcollNameLLLLssos_MASS",RECOcollNameLLLLssos_MASS,"RECOcollNameLLLLssos_MASS[20]/F");
    //Tree_->Branch("RECOcollNameLLLLssos_PT",RECOcollNameLLLLssos_PT,"RECOcollNameLLLLssos_PT[5][20]/F");

    Tree_->Branch("RECO_LLLL_MASS", RECO_LLLL_MASS, "RECO_LLLL_MASS[7][50]/F");
    Tree_->Branch("RECO_LLLL_PT", RECO_LLLL_PT, "RECO_LLLL_PT[7][50]/F");
    Tree_->Branch("RECO_LLLL_ETA", RECO_LLLL_ETA, "RECO_LLLL_ETA[7][50]/F");
    Tree_->Branch("RECO_LLLL_PHI", RECO_LLLL_PHI, "RECO_LLLL_PHI[7][50]/F");
 
   
    // Electron block
    Tree_->Branch("RECOELE_E", RECOELE_E, "RECOELE_E[100]/F"); 
    Tree_->Branch("RECOELE_PT",RECOELE_PT,"RECOELE_PT[100]/F");
    Tree_->Branch("RECOELE_PTError",RECOELE_PTError,"RECOELE_PTError[100]/F");
    Tree_->Branch("RECOELE_P",RECOELE_P,"RECOELE_P[100]/F");
    Tree_->Branch("RECOELE_ETA",RECOELE_ETA,"RECOELE_ETA[100]/F"); 
    Tree_->Branch("RECOELE_THETA",RECOELE_THETA,"RECOELE_THETA[100]/F"); 
    Tree_->Branch("RECOELE_PHI",RECOELE_PHI,"RECOELE_PHI[100]/F"); 
    Tree_->Branch("RECOELE_MASS",RECOELE_MASS,"RECOELE_MASS[100]/F"); 
    Tree_->Branch("RECOELE_CHARGE",RECOELE_CHARGE,"RECOELE_CHARGE[100]/F");   
    // Core attributes
    Tree_->Branch("RECOELE_isEcalDriven",   RECOELE_isEcalDriven,   "RECOELE_isEcalDriven[100]/b");   
    Tree_->Branch("RECOELE_isTrackerDriven",RECOELE_isTrackerDriven,"RECOELE_isTrackerDriven[100]/b");  
    Tree_->Branch("RECOELE_gsftrack_NPixHits",RECOELE_gsftrack_NPixHits,"RECOELE_gsftrack_NPixHits[100]/f");
    Tree_->Branch("RECOELE_gsftrack_NStripHits",RECOELE_gsftrack_NStripHits,"RECOELE_gsftrack_NStripHits[100]/f");
    Tree_->Branch("RECOELE_gsftrack_chi2",  RECOELE_gsftrack_chi2, "RECOELE_gsftrack_chi2[100]/F");
    Tree_->Branch("RECOELE_gsftrack_dxyB",  RECOELE_gsftrack_dxyB,"RECOELE_gsftrack_dxyB[100]/F");
    Tree_->Branch("RECOELE_gsftrack_dxy",   RECOELE_gsftrack_dxy,"RECOELE_gsftrack_dxy[100]/F");
    Tree_->Branch("RECOELE_gsftrack_dxyError", RECOELE_gsftrack_dxyError,"RECOELE_gsftrack_dxyError[100]/F");
    Tree_->Branch("RECOELE_gsftrack_dzB",   RECOELE_gsftrack_dzB,"RECOELE_gsftrack_dzB[100]/F");
    Tree_->Branch("RECOELE_gsftrack_dz",    RECOELE_gsftrack_dz,"RECOELE_gsftrack_dz[100]/F");
    Tree_->Branch("RECOELE_gsftrack_dzError", RECOELE_gsftrack_dzError,"RECOELE_gsftrack_dzError[100]/F");
    Tree_->Branch("RECOELE_gsftrack_losthits", RECOELE_gsftrack_losthits,"RECOELE_gsftrack_losthits[100]/I");
    Tree_->Branch("RECOELE_gsftrack_validhits",RECOELE_gsftrack_validhits,"RECOELE_gsftrack_validhits[100]/I");
    Tree_->Branch("RECOELE_gsftrack_expected_inner_hits",RECOELE_gsftrack_expected_inner_hits,"RECOELE_gsftrack_expected_inner_hits[100]/I") ; 
    Tree_->Branch("RECOELE_scl_E",  RECOELE_scl_E,"RECOELE_scl_E[100]/F");
    Tree_->Branch("RECOELE_scl_Et", RECOELE_scl_Et,"RECOELE_scl_Et[100]/F");
    Tree_->Branch("RECOELE_scl_Eta",RECOELE_scl_Eta,"RECOELE_scl_Eta[100]/F");
    Tree_->Branch("RECOELE_scl_Phi",RECOELE_scl_Phi,"RECOELE_scl_Phi[100]/F");    
    Tree_->Branch("RECOELE_ecalEnergy",RECOELE_ecalEnergy,"RECOELE_ecalEnergy[100]/F");

    // Track-Cluster Matching    
    Tree_->Branch("RECOELE_ep",          RECOELE_ep,"RECOELE_ep[100]/F");
    Tree_->Branch("RECOELE_eSeedp",      RECOELE_eSeedp,"RECOELE_eSeedp[100]/F");
    Tree_->Branch("RECOELE_eSeedpout",   RECOELE_eSeedpout,"RECOELE_eSeedpout[100]/F");
    Tree_->Branch("RECOELE_eElepout",    RECOELE_eElepout,"RECOELE_eElepout[100]/F");
    Tree_->Branch("RECOELE_deltaEtaIn",  RECOELE_deltaEtaIn,"RECOELE_deltaEtaIn[100]/F");
    Tree_->Branch("RECOELE_deltaEtaSeed",RECOELE_deltaEtaSeed,"RECOELE_deltaEtaSeed[100]/F");
    Tree_->Branch("RECOELE_deltaEtaEle", RECOELE_deltaEtaEle,"RECOELE_deltaEtaEle[100]/F");
    Tree_->Branch("RECOELE_deltaPhiIn",  RECOELE_deltaPhiIn,"RECOELE_deltaPhiIn[100]/F");
    Tree_->Branch("RECOELE_deltaPhiSeed",RECOELE_deltaPhiSeed,"RECOELE_deltaPhiSeed[100]/F");
    Tree_->Branch("RECOELE_deltaPhiEle", RECOELE_deltaPhiEle,"RECOELE_deltaPhiEle[100]/F");
    // Fiducial flags
    Tree_->Branch("RECOELE_isbarrel",    RECOELE_isbarrel,   "RECOELE_isbarrel[100]/I");   
    Tree_->Branch("RECOELE_isendcap",    RECOELE_isendcap,   "RECOELE_isendcap[100]/I");   
    Tree_->Branch("RECOELE_isGap",       RECOELE_isGap,      "RECOELE_isGap[100]/I");
    Tree_->Branch("RECOELE_isEBetaGap",  RECOELE_isEBetaGap, "RECOELE_isEBetaGap[100]/I");   
    Tree_->Branch("RECOELE_isEBphiGap",  RECOELE_isEBphiGap, "RECOELE_isEBphiGap[100]/I");   
    Tree_->Branch("RECOELE_isEEdeeGap",  RECOELE_isEEdeeGap, "RECOELE_isEEdeeGap[100]/I");   
    Tree_->Branch("RECOELE_isEEringGap", RECOELE_isEEringGap,"RECOELE_isEEringGap[100]/I");   
    // Shower shape
    Tree_->Branch("RECOELE_sigmaIetaIeta", RECOELE_sigmaIetaIeta, "RECOELE_sigmaIetaIeta[100]/F");   
    Tree_->Branch("RECOELE_sigmaEtaEta",   RECOELE_sigmaEtaEta,   "RECOELE_sigmaEtaEta[100]/F");   
    Tree_->Branch("RECOELE_e15",           RECOELE_e15,           "RECOELE_e15[100]/F");   
    Tree_->Branch("RECOELE_e25max",        RECOELE_e25max,        "RECOELE_e25max[100]/F");   
    Tree_->Branch("RECOELE_e55",           RECOELE_e55,           "RECOELE_e55[100]/F");   
    Tree_->Branch("RECOELE_he",            RECOELE_he,            "RECOELE_he[100]/F");   
    Tree_->Branch("RECOELE_r9",            RECOELE_r9,            "RECOELE_r9[100]/F");   
    // Particle flow
    Tree_->Branch("RECOELE_mva", RECOELE_mva,"RECOELE_mva[100]/F");   
    // Brem & Classifaction
    Tree_->Branch("RECOELE_fbrem",  RECOELE_fbrem,  "RECOELE_fbrem[100]/F");   
    Tree_->Branch("RECOELE_nbrems", RECOELE_nbrems, "RECOELE_nbrems[100]/I");   
    //  golden/bigbrem/(narrow)/showering/crack
    Tree_->Branch("RECOELE_Class",  RECOELE_Class,  "RECOELE_Class[100]/I");  
    //fBrem addition
    Tree_->Branch("RECOELE_fbrem_mode", &RECOELE_fbrem_mode,"RECOELE_fbrem_mode[100]/D");
    Tree_->Branch("RECOELE_fbrem_mean", &RECOELE_fbrem_mean,"RECOELE_fbrem_mean[100]/D");
    
    // Isolation
    Tree_->Branch("RECOELE_EGMTRACKISO",RECOELE_EGMTRACKISO,"RECOELE_EGMTRACKISO[100]/F");  
    Tree_->Branch("RECOELE_EGMHCALISO",RECOELE_EGMHCALISO,"RECOELE_EGMHCALISO[100]/F");  
    Tree_->Branch("RECOELE_EGMECALISO",RECOELE_EGMECALISO,"RECOELE_EGMECALISO[100]/F"); 
    Tree_->Branch("RECOELE_EGMX",RECOELE_EGMX,"RECOELE_EGMX[100]/F"); 

    // PF isolation
    Tree_->Branch("RECOELE_PFchAllPart",  RECOELE_PFchAllPart,  "RECOELE_PFchAllPart[100]/D");
    Tree_->Branch("RECOELE_PFchHad",         RECOELE_PFchHad,  "RECOELE_PFchHad[100]/D");
    Tree_->Branch("RECOELE_PFneuHad",        RECOELE_PFneuHad, "RECOELE_PFneuHad[100]/D");
    Tree_->Branch("RECOELE_PFphoton",        RECOELE_PFphoton, "RECOELE_PFphoton[100]/D");
    Tree_->Branch("RECOELE_PFPUchAllPart",   RECOELE_PFPUchAllPart,"RECOELE_PFPUchAllPart[100]/D");
    Tree_->Branch("RECOELE_PFX_dB",          RECOELE_PFX_dB,   "RECOELE_PFX_dB[100]/D");
    Tree_->Branch("RECOELE_PFX_rho",         RECOELE_PFX_rho,  "RECOELE_PFX_rho[100]/D");

    // Electron Regression
    Tree_->Branch("RECOELE_regEnergy",RECOELE_regEnergy,"RECOELE_regEnergy[100]/D");
    Tree_->Branch("RECOELE_regEnergyError",RECOELE_regEnergyError,"RECOELE_regEnergyError[100]/D");
    
    // Vertexing DA and KF
    Tree_->Branch("RECOELE_SIP",RECOELE_SIP,"RECOELE_SIP[100]/F"); 
    Tree_->Branch("RECOELE_IP",RECOELE_IP,"RECOELE_IP[100]/F"); 
    Tree_->Branch("RECOELE_IPERROR",RECOELE_IPERROR,"RECOELE_IPERROR[100]/F"); 
    Tree_->Branch("RECOELE_SIP_KF",RECOELE_SIP_KF,"RECOELE_SIP_KF[100]/F"); 
    Tree_->Branch("RECOELE_IP_KF",RECOELE_IP_KF,"RECOELE_IP_KF[100]/F"); 
    Tree_->Branch("RECOELE_IPERROR_KF",RECOELE_IPERROR_KF,"RECOELE_IPERROR_KF[100]/F"); 

     // GD vertex
    Tree_->Branch("RECOELE_SIP_GD",RECOELE_SIP_GD,"RECOELE_SIP_GD[100]/F"); //2e2mu
    Tree_->Branch("RECOELE_SIP_GDEEEE",RECOELE_SIP_GDEEEE,"RECOELE_SIP_GDEEEE[100]/F");  //4e
    // Std vertex
    Tree_->Branch("RECOELE_SIP_Std",RECOELE_SIP_Std,"RECOELE_SIP_Std[100]/F"); //2e2mu
    Tree_->Branch("RECOELE_SIP_StdEEEE",RECOELE_SIP_StdEEEE,"RECOELE_SIP_StdEEEE[100]/F");  //4e
    // Kin vertex
    Tree_->Branch("RECOELE_SIP_Kin",RECOELE_SIP_Kin,"RECOELE_SIP_Kin[100]/F"); //2e2mu
    Tree_->Branch("RECOELE_SIP_KinEEEE",RECOELE_SIP_KinEEEE,"RECOELE_SIP_KinEEEE[100]/F");  //4e


    Tree_->Branch("RECOELE_STIP",RECOELE_STIP,"RECOELE_STIP[100]/F"); 
    Tree_->Branch("RECOELE_SLIP",RECOELE_SLIP,"RECOELE_SLIP[100]/F"); 
    Tree_->Branch("RECOELE_TIP",RECOELE_TIP,"RECOELE_TIP[100]/F"); 
    Tree_->Branch("RECOELE_LIP",RECOELE_LIP,"RECOELE_LIP[100]/F"); 
    Tree_->Branch("RECOELE_TIPERROR",RECOELE_TIPERROR,"RECOELE_TIPERROR[100]/F"); 
    Tree_->Branch("RECOELE_LIPERROR",RECOELE_LIPERROR,"RECOELE_LIPERROR[100]/F"); 

    Tree_->Branch("RECOELE_sclRawE",ele_sclRawE,"RECOELE_sclRawE[100]/D") ;
    Tree_->Branch("RECOELE_sclX",ele_sclX,"RECOELE_sclX[100]/D") ; 
    Tree_->Branch("RECOELE_sclY",ele_sclY,"RECOELE_sclY[100]/D") ; 
    Tree_->Branch("RECOELE_sclZ",ele_sclZ,"RECOELE_sclZ[100]/D") ;
    Tree_->Branch("RECOELE_seedSubdet1",ele_seedSubdet1,"RECOELE_seedSubdet1[100]/D") ;
    Tree_->Branch("RECOELE_seedDphi1",ele_seedDphi1,"RECOELE_seedDphi[100]/D") ; 
    Tree_->Branch("RECOELE_seedDrz1",ele_seedDrz1,"RECOELE_seedDrz1[100]/D") ;
    Tree_->Branch("RECOELE_seedSubdet2",ele_seedSubdet2,"RECOELE_seedSubdet2[100]/D") ;
    Tree_->Branch("RECOELE_seedDphi2",ele_seedDphi2,"RECOELE_seedDphi2[100]/D") ; 
    Tree_->Branch("RECOELE_seedDrz2",ele_seedDrz2,"RECOELE_seedDrz2[100]/D") ;
    Tree_->Branch("RECOELE_eidVeryLoose",ele_eidVeryLoose,"RECOELE_eidVeryLoose[100]/D") ; 
    Tree_->Branch("RECOELE_eidLoose",ele_eidLoose,"RECOELE_eidLoose[100]/D") ; 
    Tree_->Branch("RECOELE_eidMedium",ele_eidMedium,"RECOELE_eidMedium[100]/D") ; 
    Tree_->Branch("RECOELE_eidTight",ele_eidTight,"RECOELE_eidTight[100]/D") ; 
    Tree_->Branch("RECOELE_eidHZZVeryLoose",ele_eidHZZVeryLoose,"RECOELE_eidHZZVeryLoose[100]/D") ; 
    Tree_->Branch("RECOELE_eidHZZLoose",ele_eidHZZLoose,"RECOELE_eidHZZLoose[100]/D") ; 
    Tree_->Branch("RECOELE_eidHZZMedium",ele_eidHZZMedium,"RECOELE_eidHZZMedium[100]/D") ; 
    Tree_->Branch("RECOELE_eidHZZTight",ele_eidHZZTight,"RECOELE_eidHZZTight[100]/D") ; 
    Tree_->Branch("RECOELE_mvaTrigV0",RECOELE_mvaTrigV0,"RECOELE_mvaTrigV0[100]/D") ;     
    Tree_->Branch("RECOELE_mvaNonTrigV0",RECOELE_mvaNonTrigV0,"RECOELE_mvaNonTrigV0[100]/D") ; 
    Tree_->Branch("RECOELE_COV",RECOELE_COV,"RECOELE_COV[100][3][3]/D"); 

    Tree_->Branch("RECOELE_TLE_ParentSC_X",RECOELE_TLE_ParentSC_X,"RECOELE_TLE_ParentSC_X[20]/F");
    Tree_->Branch("RECOELE_TLE_ParentSC_Y",RECOELE_TLE_ParentSC_Y,"RECOELE_TLE_ParentSC_Y[20]/F");
    Tree_->Branch("RECOELE_TLE_ParentSC_Z",RECOELE_TLE_ParentSC_Z,"RECOELE_TLE_ParentSC_X[20]/F");

    // Muon block
    Tree_->Branch("RECOMU_isPFMu", RECOMU_isPFMu, "RECOMU_isPFMu[100]/b");
    Tree_->Branch("RECOMU_isGlobalMu", RECOMU_isGlobalMu, "RECOMU_isGlobalMu[100]/b");
    Tree_->Branch("RECOMU_isStandAloneMu", RECOMU_isStandAloneMu, "RECOMU_isStandAloneMu[100]/b");
    Tree_->Branch("RECOMU_isTrackerMu", RECOMU_isTrackerMu, "RECOMU_isTrackerMu[100]/b");
    Tree_->Branch("RECOMU_isCaloMu", RECOMU_isCaloMu, "RECOMU_isCaloMu[100]/b");
    Tree_->Branch("RECOMU_isTrackerHighPtMu", RECOMU_isTrackerHighPtMu, "RECOMU_isTrackerHighPtMu[100]/b");
    Tree_->Branch("RECOMU_isME0Muon", RECOMU_isME0Muon, "RECOMU_isME0Muon[100]/b");
    Tree_->Branch("RECOMU_E",RECOMU_E,"RECOMU_E[100]/F"); 
    Tree_->Branch("RECOMU_PT",RECOMU_PT,"RECOMU_PT[100]/F"); 
    Tree_->Branch("RECOMU_P",RECOMU_P,"RECOMU_P[100]/F"); 
    Tree_->Branch("RECOMU_ETA",RECOMU_ETA,"RECOMU_ETA[100]/F"); 
    Tree_->Branch("RECOMU_THETA",RECOMU_THETA,"RECOMU_THETA[100]/F"); 
    Tree_->Branch("RECOMU_PHI",RECOMU_PHI,"RECOMU_PHI[100]/F"); 
    Tree_->Branch("RECOMU_MASS",RECOMU_MASS,"RECOMU_MASS[100]/F"); 
    Tree_->Branch("RECOMU_CHARGE",RECOMU_CHARGE,"RECOMU_CHARGE[100]/F");

    Tree_->Branch("RECOMU_COV",RECOMU_COV,"RECOMU_COV[100][3][3]/D"); 

    Tree_->Branch("RECOMU_TRACKISO",RECOMU_TRACKISO,"RECOMU_TRACKISO[100]/F");  
    Tree_->Branch("RECOMU_TRACKISO_SUMPT",RECOMU_TRACKISO_SUMPT,"RECOMU_TRACKISO_SUMPT[100]/F");  
    Tree_->Branch("RECOMU_HCALISO",RECOMU_HCALISO,"RECOMU_HCALISO[100]/F");  
    Tree_->Branch("RECOMU_ECALISO",RECOMU_ECALISO,"RECOMU_ECALISO[100]/F"); 
    Tree_->Branch("RECOMU_X",RECOMU_X,"RECOMU_X[100]/F");

    Tree_->Branch("RECOMU_PFchHad",  RECOMU_PFchHad,  "RECOMU_PFchHad[100]/D");
    Tree_->Branch("RECOMU_PFneuHad", RECOMU_PFneuHad, "RECOMU_PFneuHad[100]/D");
    Tree_->Branch("RECOMU_PFphoton", RECOMU_PFphoton, "RECOMU_PFphoton[100]/D");
    Tree_->Branch("RECOMU_PFPUchAllPart",RECOMU_PFPUchAllPart,"RECOMU_PFPUchAllPart[100]/D");
    Tree_->Branch("RECOMU_PFX_dB",   RECOMU_PFX_dB,   "RECOMU_PFX_dB[100]/D");
    Tree_->Branch("RECOMU_PFX_rho",  RECOMU_PFX_rho,  "RECOMU_PFX_rho[100]/D");


    // photon
    Tree_->Branch("RECOPFPHOT_PFchHad",  RECOPFPHOT_PFchHad,  "RECOPFPHOT_PFchHad[20]/D");
    Tree_->Branch("RECOPFPHOT_PFneuHad", RECOPFPHOT_PFneuHad, "RECOPFPHOT_PFneuHad[20]/D");
    Tree_->Branch("RECOPFPHOT_PFphoton", RECOPFPHOT_PFphoton, "RECOPFPHOT_PFphoton[20]/D");
    Tree_->Branch("RECOPFPHOT_PFPUchAllPart",RECOPFPHOT_PFPUchAllPart,"RECOPFPHOT_PFPUchAllPart[20]/D");
    Tree_->Branch("RECOPFPHOT_PFX_rho",  RECOPFPHOT_PFX_rho,  "RECOPFPHOT_PFX_rho[20]/D");

    // vertexing DA and KF
    Tree_->Branch("RECOMU_SIP",RECOMU_SIP,"RECOMU_SIP[100]/F"); 
    Tree_->Branch("RECOMU_IP",RECOMU_IP,"RECOMU_IP[100]/F"); 
    Tree_->Branch("RECOMU_IPERROR",RECOMU_IPERROR,"RECOMU_IPERROR[100]/F"); 
    Tree_->Branch("RECOMU_SIP_KF",RECOMU_SIP_KF,"RECOMU_SIP_KF[100]/F"); 
    Tree_->Branch("RECOMU_IP_KF",RECOMU_IP_KF,"RECOMU_IP_KF[100]/F"); 
    Tree_->Branch("RECOMU_IPERROR_KF",RECOMU_IPERROR_KF,"RECOMU_IPERROR_KF[100]/F"); 

    // GD vertex
    Tree_->Branch("RECOMU_SIP_GD",RECOMU_SIP_GD,"RECOMU_SIP_GD[100]/F"); //2e2mu
    Tree_->Branch("RECOMU_SIP_GDMMMM",RECOMU_SIP_GDMMMM,"RECOMU_SIP_GDMMMM[100]/F");  //4mu
    // Std vertex
    Tree_->Branch("RECOMU_SIP_Std",RECOMU_SIP_Std,"RECOMU_SIP_Std[100]/F"); //2e2mu
    Tree_->Branch("RECOMU_SIP_StdMMMM",RECOMU_SIP_StdMMMM,"RECOMU_SIP_StdMMMM[100]/F");  //4mu
    // Kin vertex
    Tree_->Branch("RECOMU_SIP_Kin",RECOMU_SIP_Kin,"RECOMU_SIP_Kin[100]/F"); //2e2mu
    Tree_->Branch("RECOMU_SIP_KinMMMM",RECOMU_SIP_KinMMMM,"RECOMU_SIP_KinMMMM[100]/F");  //4mu



    Tree_->Branch("RECOMU_STIP",RECOMU_STIP,"RECOMU_STIP[100]/F"); 
    Tree_->Branch("RECOMU_SLIP",RECOMU_SLIP,"RECOMU_SLIP[100]/F"); 
    Tree_->Branch("RECOMU_TIP",RECOMU_TIP,"RECOMU_TIP[100]/F"); 
    Tree_->Branch("RECOMU_LIP",RECOMU_LIP,"RECOMU_LIP[100]/F"); 
    Tree_->Branch("RECOMU_TIPERROR",RECOMU_TIPERROR,"RECOMU_TIPERROR[100]/F"); 
    Tree_->Branch("RECOMU_LIPERROR",RECOMU_LIPERROR,"RECOMU_LIPERROR[100]/F"); 
    
 

    Tree_->Branch("RECOMU_caloCompatibility",RECOMU_caloCompatibility,"RECOMU_caloCompatibility[100]/F");
    Tree_->Branch("RECOMU_segmentCompatibility",RECOMU_segmentCompatibility,"RECOMU_segmentCompatibility[100]/F"); 
    Tree_->Branch("RECOMU_numberOfMatches",RECOMU_numberOfMatches,"RECOMU_numberOfMatches[100]/i");
    Tree_->Branch("RECOMU_numberOfMatchedStations",RECOMU_numberOfMatchedStations,"RECOMU_numberOfMatchedStations[100]/i");
    Tree_->Branch("RECOMU_glbmuPromptTight",RECOMU_glbmuPromptTight,"RECOMU_glbmuPromptTight[100]/b");
 
    // track variables from muons:
    Tree_->Branch( "RECOMU_trkmuArbitration", RECOMU_trkmuArbitration, "RECOMU_trkmuArbitration[100]/b");
    Tree_->Branch( "RECOMU_trkmu2DCompatibilityLoose", RECOMU_trkmu2DCompatibilityLoose, "RECOMU_trkmu2DCompatibilityLoose[100]/b");
    Tree_->Branch( "RECOMU_trkmu2DCompatibilityTight", RECOMU_trkmu2DCompatibilityTight, "RECOMU_trkmu2DCompatibilityTight[100]/b");
    Tree_->Branch( "RECOMU_trkmuOneStationLoose", RECOMU_trkmuOneStationLoose, "RECOMU_trkmuOneStationLoose[100]/b");
    Tree_->Branch( "RECOMU_trkmuOneStationTight", RECOMU_trkmuOneStationTight, "RECOMU_trkmuOneStationTight[100]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationLoose", RECOMU_trkmuLastStationLoose, "RECOMU_trkmuLastStationLoose[100]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationTight", RECOMU_trkmuLastStationTight, "RECOMU_trkmuLastStationTight[100]/b");
    Tree_->Branch( "RECOMU_trkmuOneStationAngLoose", RECOMU_trkmuOneStationAngLoose, "RECOMU_trkmuOneStationAngLoose[100]/b");
    Tree_->Branch( "RECOMU_trkmuOneStationAngTight", RECOMU_trkmuOneStationAngTight, "RECOMU_trkmuOneStationAngTight[100]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationAngLoose", RECOMU_trkmuLastStationAngLoose, "RECOMU_trkmuLastStationAngLoose[100]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationAngTight", RECOMU_trkmuLastStationAngTight, "RECOMU_trkmuLastStationAngTight[100]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationOptimizedLowPtLoose",RECOMU_trkmuLastStationOptimizedLowPtLoose , "RECOMU_trkmuLastStationOptimizedLowPtLoose[100]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationOptimizedLowPtTight",RECOMU_trkmuLastStationOptimizedLowPtTight , "RECOMU_trkmuLastStationOptimizedLowPtTight[100]/b");

    Tree_->Branch( "RECOMU_mutrkPT", RECOMU_mutrkPT, "RECOMU_mutrkPT[100]/F");
    Tree_->Branch( "RECOMU_mutrkPTError", RECOMU_mutrkPTError, "RECOMU_mutrkPTError[100]/F");
    Tree_->Branch( "RECOMU_mutrkDxy", RECOMU_mutrkDxy, "RECOMU_mutrkDxy[100]/F");
    Tree_->Branch( "RECOMU_mutrkDxyError", RECOMU_mutrkDxyError, "RECOMU_mutrkDxyError[100]/F");
    Tree_->Branch( "RECOMU_mutrkDxyB", RECOMU_mutrkDxyB, "RECOMU_mutrkDxyB[100]/F");
    Tree_->Branch( "RECOMU_mutrkDz", RECOMU_mutrkDz, "RECOMU_mutrkDz[100]/F");
    Tree_->Branch( "RECOMU_mutrkDzError", RECOMU_mutrkDzError, "RECOMU_mutrkDzError[100]/F");
    Tree_->Branch( "RECOMU_mutrkDzB", RECOMU_mutrkDzB, "RECOMU_mutrkDzB[100]/F");
    Tree_->Branch( "RECOMU_mutrkChi2PerNdof", RECOMU_mutrkChi2PerNdof, "RECOMU_mutrkChi2PerNdof[100]/F");
    Tree_->Branch( "RECOMU_mutrkCharge", RECOMU_mutrkCharge, "RECOMU_mutrkCharge[100]/F");
    Tree_->Branch( "RECOMU_mutrkNHits", RECOMU_mutrkNHits, "RECOMU_mutrkNHits[100]/F");
    Tree_->Branch( "RECOMU_mutrkNStripHits", RECOMU_mutrkNStripHits, "RECOMU_mutrkNStripHits[100]/F");
    Tree_->Branch( "RECOMU_mutrkNPixHits", RECOMU_mutrkNPixHits, "RECOMU_mutrkNPixHits[100]/F");
    Tree_->Branch( "RECOMU_mutrkNMuonHits", RECOMU_mutrkNMuonHits, "RECOMU_mutrkNMuonHits[100]/F");
    Tree_->Branch( "RECOMU_mutrktrackerLayersWithMeasurement",RECOMU_mutrktrackerLayersWithMeasurement,"RECOMU_mutrktrackerLayersWithMeasurement[100]/F");
    
    Tree_->Branch( "RECOMU_muInnertrkDxy", RECOMU_muInnertrkDxy, "RECOMU_muInnertrkDxy[100]/F");
    Tree_->Branch( "RECOMU_muInnertrkDxyError", RECOMU_muInnertrkDxyError, "RECOMU_muInnertrkDxyError[100]/F");
    Tree_->Branch( "RECOMU_muInnertrkDxyB", RECOMU_muInnertrkDxyB, "RECOMU_muInnertrkDxyB[100]/F");
    Tree_->Branch( "RECOMU_muInnertrkDz", RECOMU_muInnertrkDz, "RECOMU_muInnertrkDz[100]/F");
    Tree_->Branch( "RECOMU_muInnertrkDzError", RECOMU_muInnertrkDzError, "RECOMU_muInnertrkDzError[100]/F");
    Tree_->Branch( "RECOMU_muInnertrkDzB", RECOMU_muInnertrkDzB, "RECOMU_muInnertrkDzB[100]/F");
    Tree_->Branch( "RECOMU_muInnertrkChi2PerNdof", RECOMU_muInnertrkChi2PerNdof, "RECOMU_muInnertrkChi2PerNdof[100]/F");
    Tree_->Branch( "RECOMU_muInnertrktrackerLayersWithMeasurement",RECOMU_muInnertrktrackerLayersWithMeasurement,"RECOMU_muInnertrktrackerLayersWithMeasurement[100]/F");
    Tree_->Branch( "RECOMU_muInnertrkPT", RECOMU_muInnertrkPT, "RECOMU_muInnertrkPT[100]/F");
    Tree_->Branch( "RECOMU_muInnertrkPTError", RECOMU_muInnertrkPTError, "RECOMU_muInnertrkPTError[100]/F");
    Tree_->Branch( "RECOMU_muInnertrkCharge", RECOMU_muInnertrkCharge, "RECOMU_muInnertrkCharge[100]/F");
    Tree_->Branch( "RECOMU_muInnertrkNHits", RECOMU_muInnertrkNHits, "RECOMU_muInnertrkNHits[100]/F");
    Tree_->Branch( "RECOMU_muInnertrkNStripHits", RECOMU_muInnertrkNStripHits, "RECOMU_muInnertrkNStripHits[100]/F");
    Tree_->Branch( "RECOMU_muInnertrkNPixHits", RECOMU_muInnertrkNPixHits, "RECOMU_muInnertrkNPixHits[100]/F");
    // best tracks for 13 TeV analysis
    Tree_->Branch( "RECOMU_mubesttrkType", RECOMU_mubesttrkType, "RECOMU_mubesttrkType[100]/I");
    Tree_->Branch( "RECOMU_mubesttrkDxy", RECOMU_mubesttrkDxy, "RECOMU_mubesttrkDxy[100]/F");
    Tree_->Branch( "RECOMU_mubesttrkDxyError", RECOMU_mubesttrkDxyError, "RECOMU_mubesttrkDxyError[100]/F");
    Tree_->Branch( "RECOMU_mubesttrkDz", RECOMU_mubesttrkDz, "RECOMU_mubesttrkDz[100]/F");
    Tree_->Branch( "RECOMU_mubesttrkDzError", RECOMU_mubesttrkDzError, "RECOMU_mubesttrkDzError[100]/F");
    Tree_->Branch( "RECOMU_mubesttrkPTError", RECOMU_mubesttrkPTError, "RECOMU_mubesttrkPTError[100]/F");

    // Geom. Discri.
    Tree_->Branch("ftsigma",        &ftsigma,        "ftsigma[100]/D");
    Tree_->Branch("gdX",            &gdX,            "gdX[100]/D");
    Tree_->Branch("gdY",            &gdY,            "gdY[100]/D");
    Tree_->Branch("gdZ",            &gdZ,            "gdZ[100]/D");
    Tree_->Branch("ftsigmalag",     &ftsigmalag,     "ftsigmalag[100]/D");
    Tree_->Branch("gdlagX",         &gdlagX,         "gdlagX[100]/D");
    Tree_->Branch("gdlagY",         &gdlagY,         "gdlagY[100]/D");
    Tree_->Branch("gdlagZ",         &gdlagZ,         "gdlagZ[100]/D");
    Tree_->Branch("gdlagProb",      &gdlagProb,      "gdlagProb[100]/D");
    Tree_->Branch("gdlagNdof",      &gdlagNdof,      "gdlagNdof[100]/D");
    Tree_->Branch("ftsigmaMMMM",    &ftsigmaMMMM,    "ftsigmaMMMM[100]/D");
    Tree_->Branch("gdXMMMM",        &gdXMMMM,        "gdXMMMM[100]/D");
    Tree_->Branch("gdYMMMM",        &gdYMMMM,        "gdYMMMM[100]/D");
    Tree_->Branch("gdZMMMM",        &gdZMMMM,        "gdZMMMM[100]/D");
    Tree_->Branch("ftsigmalagMMMM", &ftsigmalagMMMM, "ftsigmalagMMMM[100]/D");
    Tree_->Branch("gdlagXMMMM",     &gdlagXMMMM,     "gdlagXMMMM[100]/D");
    Tree_->Branch("gdlagYMMMM",     &gdlagYMMMM,     "gdlagYMMMM[100]/D");
    Tree_->Branch("gdlagZMMMM",     &gdlagZMMMM,     "gdlagZMMMM[100]/D");
    Tree_->Branch("gdlagProbMMMM",  &gdlagProbMMMM,  "gdlagProbMMMM[100]/D");
    Tree_->Branch("gdlagNdofMMMM",  &gdlagNdofMMMM,  "gdlagNdofMMMM[100]/D");
    Tree_->Branch("ftsigmaEEEE",    &ftsigmaEEEE,    "ftsigmaEEEE[100]/D");
    Tree_->Branch("gdXEEEE",        &gdXEEEE,        "gdXEEEE[100]/D");
    Tree_->Branch("gdYEEEE",        &gdYEEEE,        "gdYEEEE[100]/D");
    Tree_->Branch("gdZEEEE",        &gdZEEEE,        "gdZEEEE[100]/D");
    Tree_->Branch("ftsigmalagEEEE", &ftsigmalagEEEE, "ftsigmalagEEEE[100]/D");
    Tree_->Branch("gdlagXEEEE",     &gdlagXEEEE,     "gdlagXEEEE[100]/D");
    Tree_->Branch("gdlagYEEEE",     &gdlagYEEEE,     "gdlagYEEEE[100]/D");
    Tree_->Branch("gdlagZEEEE",     &gdlagZEEEE,     "gdlagZEEEE[100]/D");
    Tree_->Branch("gdlagProbEEEE",  &gdlagProbEEEE,  "gdlagProbEEEE[100]/D");
    Tree_->Branch("gdlagNdofEEEE",  &gdlagNdofEEEE,  "gdlagNdofEEEE[100]/D");
    
    // ConstraintFit 4l
    Tree_->Branch("StdFitVertexX",        StdFitVertexX,        "StdFitVertexX[100]/D");
    Tree_->Branch("StdFitVertexY",        StdFitVertexY,        "StdFitVertexY[100]/D");
    Tree_->Branch("StdFitVertexZ",        StdFitVertexZ,        "StdFitVertexZ[100]/D");
    Tree_->Branch("StdFitVertexChi2r",    StdFitVertexChi2r,    "StdFitVertexChi2r[100]/D");
    Tree_->Branch("StdFitVertexProb",     StdFitVertexProb,     "StdFitVertexProb[100]/D");
    Tree_->Branch("StdFitVertexTrack_PT", StdFitVertexTrack_PT, "StdFitVertexTrack_PT[4][100]/F");
    Tree_->Branch("StdFitVertexTrack_ETA",StdFitVertexTrack_ETA,"StdFitVertexTrack_ETA[4][100]/F");
    Tree_->Branch("StdFitVertexTrack_PHI",StdFitVertexTrack_PHI,"StdFitVertexTrack_PHI[4][100]/F");
    Tree_->Branch("KinFitVertexX",        KinFitVertexX,        "KinFitVertexX[100]/D");
    Tree_->Branch("KinFitVertexY",        KinFitVertexY,        "KinFitVertexY[100]/D");
    Tree_->Branch("KinFitVertexZ",        KinFitVertexZ,        "KinFitVertexZ[100]/D");
    Tree_->Branch("KinFitVertexChi2r",    KinFitVertexChi2r,    "KinFitVertexChi2r[100]/D");
    Tree_->Branch("KinFitVertexProb",     KinFitVertexProb,     "KinFitVertexProb[100]/D");

    Tree_->Branch("StdFitVertexXMMMM",        StdFitVertexXMMMM,        "StdFitVertexXMMMM[100]/D");
    Tree_->Branch("StdFitVertexYMMMM",        StdFitVertexYMMMM,        "StdFitVertexYMMMM[100]/D");
    Tree_->Branch("StdFitVertexZMMMM",        StdFitVertexZMMMM,        "StdFitVertexZMMMM[100]/D");
    Tree_->Branch("StdFitVertexChi2rMMMM",    StdFitVertexChi2rMMMM,    "StdFitVertexChi2rMMMM[100]/D");
    Tree_->Branch("StdFitVertexProbMMMM",     StdFitVertexProbMMMM,     "StdFitVertexProbMMMM[100]/D");
    Tree_->Branch("StdFitVertexTrackMMMM_PT", StdFitVertexTrackMMMM_PT, "StdFitVertexTrackMMMM_PT[4][100]/F");
    Tree_->Branch("StdFitVertexTrackMMMM_ETA",StdFitVertexTrackMMMM_ETA,"StdFitVertexTrackMMMM_ETA[4][100]/F");
    Tree_->Branch("StdFitVertexTrackMMMM_PHI",StdFitVertexTrackMMMM_PHI,"StdFitVertexTrackMMMM_PHI[4][100]/F");
    Tree_->Branch("KinFitVertexXMMMM",        KinFitVertexXMMMM,        "KinFitVertexXMMMM[100]/D");
    Tree_->Branch("KinFitVertexYMMMM",        KinFitVertexYMMMM,        "KinFitVertexYMMMM[100]/D");
    Tree_->Branch("KinFitVertexZMMMM",        KinFitVertexZMMMM,        "KinFitVertexZMMMM[100]/D");
    Tree_->Branch("KinFitVertexChi2rMMMM",    KinFitVertexChi2rMMMM,    "KinFitVertexChi2rMMMM[100]/D");
    Tree_->Branch("KinFitVertexProbMMMM",     KinFitVertexProbMMMM,     "KinFitVertexProbMMMM[100]/D");
    

    Tree_->Branch("StdFitVertexXEEEE",        StdFitVertexXEEEE,        "StdFitVertexXEEEE[100]/D");
    Tree_->Branch("StdFitVertexYEEEE",        StdFitVertexYEEEE,        "StdFitVertexYEEEE[100]/D");
    Tree_->Branch("StdFitVertexZEEEE",        StdFitVertexZEEEE,        "StdFitVertexZEEEE[100]/D");
    Tree_->Branch("StdFitVertexChi2rEEEE",    StdFitVertexChi2rEEEE,    "StdFitVertexChi2rEEEE[100]/D");
    Tree_->Branch("StdFitVertexProbEEEE",     StdFitVertexProbEEEE,     "StdFitVertexProbEEEE[100]/D");
    Tree_->Branch("StdFitVertexTrackEEEE_PT", StdFitVertexTrackEEEE_PT, "StdFitVertexTrackEEEE_PT[4][100]/F");
    Tree_->Branch("StdFitVertexTrackEEEE_ETA",StdFitVertexTrackEEEE_ETA,"StdFitVertexTrackEEEE_ETA[4][100]/F");
    Tree_->Branch("StdFitVertexTrackEEEE_PHI",StdFitVertexTrackEEEE_PHI,"StdFitVertexTrackEEEE_PHI[4][100]/F");
    Tree_->Branch("KinFitVertexXEEEE",        KinFitVertexXEEEE,        "KinFitVertexXEEEE[100]/D");
    Tree_->Branch("KinFitVertexYEEEE",        KinFitVertexYEEEE,        "KinFitVertexYEEEE[100]/D");
    Tree_->Branch("KinFitVertexZEEEE",        KinFitVertexZEEEE,        "KinFitVertexZEEEE[100]/D");
    Tree_->Branch("KinFitVertexChi2rEEEE",    KinFitVertexChi2rEEEE,    "KinFitVertexChi2rEEEE[100]/D");
    Tree_->Branch("KinFitVertexProbEEEE",     KinFitVertexProbEEEE,     "KinFitVertexProbEEEE[100]/D");

    // constrintFit 3l
    Tree_->Branch("StdFitVertexChi2rMMM",     StdFitVertexChi2rMMM,    "StdFitVertexChi2rMMM[50]/D");
    Tree_->Branch("StdFitVertexProbMMM",      StdFitVertexProbMMM,     "StdFitVertexProbMMM[50]/D");
    Tree_->Branch("StdFitVertexChi2rMME",     StdFitVertexChi2rMME,    "StdFitVertexChi2rMME[50]/D");
    Tree_->Branch("StdFitVertexProbMME",      StdFitVertexProbMME,     "StdFitVertexProbMME[50]/D");
    Tree_->Branch("StdFitVertexChi2rEEE",     StdFitVertexChi2rEEE,    "StdFitVertexChi2rEEE[50]/D");
    Tree_->Branch("StdFitVertexProbEEE",      StdFitVertexProbEEE,     "StdFitVertexProbEEE[50]/D");
    Tree_->Branch("StdFitVertexChi2rMEE",     StdFitVertexChi2rMEE,    "StdFitVertexChi2rMEE[50]/D");
    Tree_->Branch("StdFitVertexProbMEE",      StdFitVertexProbMEE,     "StdFitVertexProbMEE[50]/D");


     // constrintFit Dileptons
    Tree_->Branch("StdFitVertexChi2rDiLep",   StdFitVertexChi2rDiLep,    "StdFitVertexChi2rDiLep[40]/D");
    Tree_->Branch("StdFitVertexProbDiLep",    StdFitVertexProbDiLep,     "StdFitVertexProbDiLep[40]/D");

    // Conversions
    Tree_->Branch("ConvMapDist",              ConvMapDist,              "ConvMapDist[100]/F");
    Tree_->Branch("ConvMapDcot",              ConvMapDcot,              "ConvMapDcot[100]/F");



    //MatchingMC:
    //Muons:
    Tree_->Branch("RECOMU_MatchingMCTruth", RECOMU_MatchingMCTruth, "RECOMU_MatchingMCTruth[100]/b");
    Tree_->Branch("RECOMU_MatchingMCpT", RECOMU_MatchingMCpT, "RECOMU_MatchingMCpT[100]/F");
    Tree_->Branch("RECOMU_MatchingMCEta", RECOMU_MatchingMCEta, "RECOMU_MatchingMCEta[100]/F");
    Tree_->Branch("RECOMU_MatchingMCPhi", RECOMU_MatchingMCPhi, "RECOMU_MatchingMCPhi[100]/F");

    //Electrons:
    Tree_->Branch("RECOELE_MatchingMCTruth", RECOELE_MatchingMCTruth, "RECOELE_MatchingMCTruth[100]/b");
    Tree_->Branch("RECOELE_MatchingMCpT", RECOELE_MatchingMCpT, "RECOELE_MatchingMCpT[100]/F");
    Tree_->Branch("RECOELE_MatchingMCEta", RECOELE_MatchingMCEta, "RECOELE_MatchingMCEta[100]/F");
    Tree_->Branch("RECOELE_MatchingMCPhi", RECOELE_MatchingMCPhi, "RECOELE_MatchingMCPhi[100]/F");

    //Gamma:
    Tree_->Branch("RECOPHOT_MatchingMCTruth", RECOPHOT_MatchingMCTruth, "RECOPHOT_MatchingMCTruth[50]/b");
    Tree_->Branch("RECOPHOT_MatchingMCpT", RECOPHOT_MatchingMCpT, "RECOPHOT_MatchingMCpT[50]/F");
    Tree_->Branch("RECOPHOT_MatchingMCEta", RECOPHOT_MatchingMCEta, "RECOPHOT_MatchingMCEta[50]/F");
    Tree_->Branch("RECOPHOT_MatchingMCPhi", RECOPHOT_MatchingMCPhi, "RECOPHOT_MatchingMCPhi[50]/F");

    //ZtoMuMu:
    Tree_->Branch("RECOzMuMu_MatchingMCTruth", RECOzMuMu_MatchingMCTruth, "RECOzMuMu_MatchingMCTruth[50]/b");
    Tree_->Branch("RECOzMuMu_MatchingMCpT", RECOzMuMu_MatchingMCpT, "RECOzMuMu_MatchingMCpT[50]/F");
    Tree_->Branch("RECOzMuMu_MatchingMCmass", RECOzMuMu_MatchingMCmass, "RECOzMuMu_MatchingMCmass[50]/F");
    Tree_->Branch("RECOzMuMu_MatchingMCEta", RECOzMuMu_MatchingMCEta, "RECOzMuMu_MatchingMCEta[50]/F");
    Tree_->Branch("RECOzMuMu_MatchingMCPhi", RECOzMuMu_MatchingMCPhi, "RECOzMuMu_MatchingMCPhi[50]/F");

    //ZtoEE:
    Tree_->Branch("RECOzEE_MatchingMCTruth", RECOzEE_MatchingMCTruth, "RECOzEE_MatchingMCTruth[50]/b");
    Tree_->Branch("RECOzEE_MatchingMCpT", RECOzEE_MatchingMCpT, "RECOzEE_MatchingMCpT[50]/F");
    Tree_->Branch("RECOzEE_MatchingMCmass", RECOzEE_MatchingMCmass, "RECOzEE_MatchingMCmass[50]/F");
    Tree_->Branch("RECOzEE_MatchingMCEta", RECOzEE_MatchingMCEta, "RECOzEE_MatchingMCEta[50]/F");
    Tree_->Branch("RECOzEE_MatchingMCPhi", RECOzEE_MatchingMCPhi, "RECOzEE_MatchingMCPhi[50]/F");

    //HtoZtoEEEE:
    Tree_->Branch("RECOHzzEEEE_MatchingMCTruth", RECOHzzEEEE_MatchingMCTruth, "RECOHzzEEEE_MatchingMCTruth[100]/b");
    Tree_->Branch("RECOHzzEEEE_MatchingMCpT", RECOHzzEEEE_MatchingMCpT, "RECOHzzEEEE_MatchingMCpT[100]/F");
    Tree_->Branch("RECOHzzEEEE_MatchingMCmass", RECOHzzEEEE_MatchingMCmass, "RECOHzzEEEE_MatchingMCmass[100]/F");
    Tree_->Branch("RECOHzzEEEE_MatchingMCEta", RECOHzzEEEE_MatchingMCEta, "RECOHzzEEEE_MatchingMCEta[100]/F");
    Tree_->Branch("RECOHzzEEEE_MatchingMCPhi", RECOHzzEEEE_MatchingMCPhi, "RECOHzzEEEE_MatchingMCPhi[100]/F");

    //HtoZtoEEMM:
    Tree_->Branch("RECOHzzEEMM_MatchingMCTruth", RECOHzzEEMM_MatchingMCTruth, "RECOHzzEEMM_MatchingMCTruth[100]/b");
    Tree_->Branch("RECOHzzEEMM_MatchingMCpT", RECOHzzEEMM_MatchingMCpT, "RECOHzzEEMM_MatchingMCpT[100]/F");
    Tree_->Branch("RECOHzzEEMM_MatchingMCmass", RECOHzzEEMM_MatchingMCmass, "RECOHzzEEMM_MatchingMCmass[100]/F");
    Tree_->Branch("RECOHzzEEMM_MatchingMCEta", RECOHzzEEMM_MatchingMCEta, "RECOHzzEEMM_MatchingMCEta[100]/F");
    Tree_->Branch("RECOHzzEEMM_MatchingMCPhi", RECOHzzEEMM_MatchingMCPhi, "RECOHzzEEMM_MatchingMCPhi[100]/F");

    //HtoZtoMMMM:
    Tree_->Branch("RECOHzzMMMM_MatchingMCTruth", RECOHzzMMMM_MatchingMCTruth, "RECOHzzMMMM_MatchingMCTruth[100]/b");
    Tree_->Branch("RECOHzzMMMM_MatchingMCpT", RECOHzzMMMM_MatchingMCpT, "RECOHzzMMMM_MatchingMCpT[100]/F");
    Tree_->Branch("RECOHzzMMMM_MatchingMCmass", RECOHzzMMMM_MatchingMCmass, "RECOHzzMMMM_MatchingMCmass[100]/F");
    Tree_->Branch("RECOHzzMMMM_MatchingMCEta", RECOHzzMMMM_MatchingMCEta, "RECOHzzMMMM_MatchingMCEta[100]/F");
    Tree_->Branch("RECOHzzMMMM_MatchingMCPhi", RECOHzzMMMM_MatchingMCPhi, "RECOHzzMMMM_MatchingMCPhi[100]/F");



    //Global Event 
    Tree_->Branch( "RECO_NMU", &RECO_NMU, "RECO_NMU/I"); 
    Tree_->Branch( "RECO_NELE", &RECO_NELE, "RECO_NELE/I"); 
    
    // Tracks
    Tree_->Branch( "RECO_NTRACK", &RECO_NTRACK, "RECO_NTRACK/I");
    Tree_->Branch( "RECO_TRACK_PT", &RECO_TRACK_PT, "RECO_TRACK_PT[200]/F");
    Tree_->Branch( "RECO_TRACK_ETA", &RECO_TRACK_ETA, "RECO_TRACK_ETA[200]/F");
    Tree_->Branch( "RECO_TRACK_PHI", &RECO_TRACK_PHI, "RECO_TRACK_PHI[200]/F");
    Tree_->Branch( "RECO_TRACK_CHI2", &RECO_TRACK_CHI2, "RECO_TRACK_CHI2[200]/F");
    Tree_->Branch( "RECO_TRACK_CHI2RED", &RECO_TRACK_CHI2RED, "RECO_TRACK_CHI2RED[200]/F");
    Tree_->Branch( "RECO_TRACK_CHI2PROB", &RECO_TRACK_CHI2PROB, "RECO_TRACK_CHI2PROB[200]/F");
    Tree_->Branch( "RECO_TRACK_NHITS", &RECO_TRACK_NHITS, "RECO_TRACK_NHITS[200]/I");
    Tree_->Branch( "RECO_TRACK_DXY", &RECO_TRACK_DXY, "RECO_TRACK_DXY[200]/F");
    Tree_->Branch( "RECO_TRACK_DXYERR", &RECO_TRACK_DXYERR, "RECO_TRACK_DXYERR[200]/F");
    Tree_->Branch( "RECO_TRACK_DZ", &RECO_TRACK_DZ, "RECO_TRACK_DZ[200]/F");
    Tree_->Branch( "RECO_TRACK_DZERR", &RECO_TRACK_DZERR, "RECO_TRACK_DZERR[200]/F");
    
    // Photons
    Tree_->Branch("RECO_NPHOT", &RECO_NPHOT, "RECO_NPHOT/I");
    Tree_->Branch("RECOPHOT_PT",RECOPHOT_PT,"RECOPHOT_PT[20]/F"); 
    Tree_->Branch("RECOPHOT_ETA",RECOPHOT_ETA,"RECOPHOT_ETA[20]/F"); 
    Tree_->Branch("RECOPHOT_PHI",RECOPHOT_PHI,"RECOPHOT_PHI[20]/F"); 
    Tree_->Branch("RECOPHOT_THETA",RECOPHOT_THETA,"RECOPHOT_THETA[20]/F"); 
    Tree_->Branch("RECOPHOT_TLE_ParentSC_X",RECOPHOT_TLE_ParentSC_X,"RECOPHOT_TLE_ParentSC_X[20]/F");
    Tree_->Branch("RECOPHOT_TLE_ParentSC_Y",RECOPHOT_TLE_ParentSC_Y,"RECOPHOT_TLE_ParentSC_Y[20]/F");
    Tree_->Branch("RECOPHOT_TLE_ParentSC_Z",RECOPHOT_TLE_ParentSC_Z,"RECOPHOT_TLE_ParentSC_X[20]/F");

    Tree_->Branch("RECO_NPFPHOT", &RECO_NPFPHOT, "RECO_NPFPHOT/I");
    Tree_->Branch("RECOPFPHOT_PT",RECOPFPHOT_PT,"RECOPFPHOT_PT[20]/F"); 
    Tree_->Branch("RECOPFPHOT_PTError",RECOPFPHOT_PTError,"RECOPFPHOT_PTError[20]/F");  
    Tree_->Branch("RECOPFPHOT_ETA",RECOPFPHOT_ETA,"RECOPFPHOT_ETA[20]/F"); 
    Tree_->Branch("RECOPFPHOT_PHI",RECOPFPHOT_PHI,"RECOPFPHOT_PHI[20]/F"); 
    Tree_->Branch("RECOPFPHOT_THETA",RECOPFPHOT_THETA,"RECOPFPHOT_THETA[20]/F"); 
    
    
    //Beam Spot position
    Tree_->Branch("BeamSpot_X",&BeamSpot_X,"BeamSpot_X/D");
    Tree_->Branch("BeamSpot_Y",&BeamSpot_Y,"BeamSpot_Y/D");
    Tree_->Branch("BeamSpot_Z",&BeamSpot_Z,"BeamSpot_Z/D");
    // Vertices
    Tree_->Branch( "RECO_NVTX", &RECO_NVTX, "RECO_NVTX/I");
    Tree_->Branch( "RECO_VERTEX_x", RECO_VERTEX_x, "RECO_VERTEX_x[15]/F");
    Tree_->Branch( "RECO_VERTEX_y", RECO_VERTEX_y, "RECO_VERTEX_y[15]/F");
    Tree_->Branch( "RECO_VERTEX_z", RECO_VERTEX_z, "RECO_VERTEX_z[15]/F");
    Tree_->Branch( "RECO_VERTEX_ndof", RECO_VERTEX_ndof, "RECO_VERTEX_ndof[15]/F");
    Tree_->Branch( "RECO_VERTEX_chi2", RECO_VERTEX_chi2, "RECO_VERTEX_chi2[15]/F");
    Tree_->Branch( "RECO_VERTEX_ntracks", RECO_VERTEX_ntracks, "RECO_VERTEX_ntracks[15]/I");
    Tree_->Branch( "RECO_VERTEXPROB", RECO_VERTEXPROB, "RECO_VERTEXPROB[15]/F");
    Tree_->Branch( "RECO_VERTEX_isValid", RECO_VERTEX_isValid, "RECO_VERTEX_isValid[15]/b");
    Tree_->Branch( "RECO_VERTEX_TRACK_PT",RECO_VERTEX_TRACK_PT,"RECO_VERTEX_TRACK_PT[15][100]/F");
    
    // PFJets
    Tree_->Branch( "RECO_PFJET_N",   &RECO_PFJET_N,   "RECO_PFJET_N/I");
    Tree_->Branch( "RECO_PFJET_CHARGE",  RECO_PFJET_CHARGE,  "RECO_PFJET_CHARGE[200]/I");
    Tree_->Branch( "RECO_PFJET_ET",  RECO_PFJET_ET,  "RECO_PFJET_ET[200]/F");
    Tree_->Branch( "RECO_PFJET_PT",  RECO_PFJET_PT,  "RECO_PFJET_PT[200]/F");
    Tree_->Branch( "RECO_PFJET_ETA", RECO_PFJET_ETA, "RECO_PFJET_ETA[200]/F");
    Tree_->Branch( "RECO_PFJET_PHI", RECO_PFJET_PHI, "RECO_PFJET_PHI[200]/F");
    Tree_->Branch( "RECO_PFJET_PUID", RECO_PFJET_PUID, "RECO_PFJET_PUID[200]/I");
    Tree_->Branch( "RECO_PFJET_PUID_MVA", RECO_PFJET_PUID_MVA, "RECO_PFJET_PUID_MVA[200]/F");
    Tree_->Branch( "RHO_ele", &RHO_ele, "RHO_ele/D");
    Tree_->Branch( "RHO_mu", &RHO_mu, "RHO_mu/D");
    
    //CaloMET
    Tree_->Branch( "RECO_CALOMET",          &calomet,          "RECO_CALOMET/F");
 /*    Tree_->Branch( "RECO_CALOMETHO",        &calometho,        "RECO_CALOMETHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETNOHFHO",    &calometnohfho,    "RECO_CALOMETNOHFHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETNOHF",      &calometnohf,      "RECO_CALOMETNOHF/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPTHO",     &calometoptho,     "RECO_CALOMETOPTHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPTNOHFHO", &calometoptnohfho, "RECO_CALOMETOPTNOHFHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPTNOHF",   &calometoptnohf,   "RECO_CALOMETOPTNOHF/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPT",       &calometopt,       "RECO_CALOMETOPT/F"); */
    //Particle Flow MET
    Tree_->Branch( "RECO_PFMET", &pfmet, "RECO_PFMET/F");
    Tree_->Branch( "RECO_PFMET_X", &pfmet_x, "RECO_PFMET_X/F");
    Tree_->Branch( "RECO_PFMET_Y", &pfmet_y, "RECO_PFMET_Y/F");
    Tree_->Branch( "RECO_PFMET_PHI", &pfmet_phi, "RECO_PFMET_PHI/F");
    Tree_->Branch( "RECO_PFMET_THETA", &pfmet_theta, "RECO_PFMET_THETA/F");
    //Track Corrected MET
    Tree_->Branch( "RECO_TCMET", &tcmet, "RECO_TCMET/F");
    //Type I correction MET
    Tree_->Branch( "RECO_CORMETMUONS",  &cormetmuons,  "RECO_CORMETMUONS/F");
   

    // Btagging jets and discriminators
    Tree_->Branch("tCHighEff_BTagJet_PT",tCHighEff_BTagJet_PT,"tCHighEff_BTagJet_PT[50]/F");
    Tree_->Branch("tCHighPur_BTagJet_PT", tCHighPur_BTagJet_PT,"tCHighPur_BTagJet_PT[50]/F");
    Tree_->Branch("cSV_BTagJet_PT",cSV_BTagJet_PT,"cSV_BTagJet_PT[50]/F");
    Tree_->Branch("tCHighEff_BTagJet_ETA",tCHighEff_BTagJet_ETA,"tCHighEff_BTagJet_ETA[50]/F");
    Tree_->Branch("tCHighPur_BTagJet_ETA", tCHighPur_BTagJet_ETA,"tCHighPur_BTagJet_ETA[50]/F");
    Tree_->Branch("cSV_BTagJet_ETA",cSV_BTagJet_ETA,"cSV_BTagJet_ETA[50]/F");
    Tree_->Branch("tCHighEff_BTagJet_PHI",tCHighEff_BTagJet_PHI,"tCHighEff_BTagJet_PHI[50]/F");
    Tree_->Branch("tCHighPur_BTagJet_PHI", tCHighPur_BTagJet_PHI,"tCHighPur_BTagJet_PHI[50]/F");
    Tree_->Branch("cSV_BTagJet_PHI",cSV_BTagJet_PHI,"cSV_BTagJet_PHI[50]/F");
    Tree_->Branch("cSV_BTagJet_ET",cSV_BTagJet_ET,"cSV_BTagJet_ET[50]/F");
    Tree_->Branch("tCHighEff_BTagJet_DISCR",tCHighEff_BTagJet_DISCR,"tCHighEff_BTagJet_DISCR[50]/F");
    Tree_->Branch("tCHighPur_BTagJet_DISCR", tCHighPur_BTagJet_DISCR,"tCHighPur_BTagJet_DISCR[50]/F");
    Tree_->Branch("cSV_BTagJet_DISCR",cSV_BTagJet_DISCR,"cSV_BTagJet_DISCR[50]/F");     
  }
  
  
  void Initialize(){
    
    irun=-999,ievt=-999,ils=-999;
    Avginstlumi=-999.;
    RHO=-999.,RHO_ele=-999.,RHO_mu=-999.;

    // PU
    num_PU_vertices=-999;
    PU_BunchCrossing=-999;

    
    // HLT flags   

    RECO_nMuHLTMatch=0;
    RECO_nEleHLTMatch=0;
    
    //for (int ii=0;ii<200;ii++){
    //  HLTPathsFired[ii]="";
    //}

    leptonscands2e2mu_= new (CandidateCollection);
    leptonscands2e2murf_= new (CandidateCollection);
    leptonscands4mu_= new (CandidateCollection);
    leptonscands4murf_= new (CandidateCollection);
    leptonscands4e_= new (CandidateCollection);
    leptonscands4erf_= new (CandidateCollection);

    leptonscands_Z0= new (CandidateCollection);
    leptonscands_Z1= new (CandidateCollection);
    leptonscands_Zss0= new (CandidateCollection);
    leptonscands_Zss1= new (CandidateCollection);
    leptonscands_Zcross= new (CandidateCollection);
    leptonscands_DiLep= new (CandidateCollection);
    leptonscands_MMMM= new (CandidateCollection);
    leptonscands_EEEE= new (CandidateCollection);
    leptonscands_EEMM= new (CandidateCollection);
    leptonscands_LLL0= new (CandidateCollection);
    leptonscands_LLL1= new (CandidateCollection);
    leptonscands_LLL2= new (CandidateCollection);
    leptonscands_LLL3= new (CandidateCollection);
    leptonscands_LLLLss0= new (CandidateCollection);
    leptonscands_LLLLss1= new (CandidateCollection);
    leptonscands_LLLLss2= new (CandidateCollection);
    leptonscands_LLLl0= new (CandidateCollection);
    leptonscands_LLLl1= new (CandidateCollection);
    leptonscands_LLLL= new (CandidateCollection);

    for  (int i=0; i<4;i++){ 
      MC_LEPT_PT[i]=-999.;
      MC_LEPT_ETA[i]=-999.;
      MC_LEPT_PHI[i]=-999.;
      MC_LEPT_THETA[i]=-999.;
      MC_LEPT_PDGID[i]=-999.;
    }
    for  (int i=0; i<4;i++){ 
      for  (int j=0; j<7;j++){ 
	MC_ZZ_MASS[i][j]=-999.;
	MC_ZZ_PT[i][j]=-999.;
	MC_ZZ_ETA[i][j]=-999.;
	MC_ZZ_PHI[i][j]=-999.;
	MC_ZZ_THETA[i][j]=-999.;
	MC_ZZ_PDGID[i][j]=-999.;
      }
    }

    for(int i=0; i<2; ++i) {
      for(int j=0; j<5; ++j) {
	MC_Z_PT[i][j]=-999;
	MC_Z_ETA[i][j]=-999;
	MC_Z_PHI[i][j]=-999;
	MC_Z_THETA[i][j]=-999;
	MC_Z_MASS[i][j]=-999;
	MC_Z_PDGID[i][j]=-999;
      }
    }

    for  (int i=0; i<50;i++){ 
      for  (int j=0; j<5;j++){ 
	MC_fourl_MASS[i][j]=-999.;
	MC_fourl_PT[i][j]=-999.;
	MC_fourl_PDGID[i][j]=-999.;
      }
    }

    for (int i=0; i<7;i++){
      MC_E[i]=-999.;
      MC_PT[i]=-999.;
      MC_ETA[i]=-999.;
      MC_THETA[i]=-999.;
      MC_PHI[i]=-999.;
      MC_MASS[i]=-999.;
      MC_PDGID[i]=-999.;
      


  

    }

   
    for (int i=0; i<100;i++){
      RECORF_2e2mu_cosTheta1_spin[i]=-999.;
      RECORF_2e2mu_cosTheta2_spin[i]=-999.;
      RECORF_2e2mu_cosThetaStar_spin[i]=-999.;
      RECORF_2e2mu_Phi_spin[i]=-999.;
      RECORF_2e2mu_Phi1_spin[i]=-999.;
      RECORF_2e2mu_Phi2_spin[i]=-999.;
      RECORF_2e2mu_phi1RF_spin[i]=-999.;
      RECORF_2e2mu_phi2RF_spin[i]=-999.;
      RECORF_2e2mu_MELA[i]=-999.;
      
      RECORF_4mu_cosTheta1_spin[i]=-999.;
      RECORF_4mu_cosTheta2_spin[i]=-999.;
      RECORF_4mu_cosThetaStar_spin[i]=-999.;
      RECORF_4mu_Phi_spin[i]=-999.;
      RECORF_4mu_Phi1_spin[i]=-999.;
      RECORF_4mu_Phi2_spin[i]=-999.;
      RECORF_4mu_phi1RF_spin[i]=-999.;
      RECORF_4mu_phi2RF_spin[i]=-999.;
      RECORF_4mu_MELA[i]=-999.;

      RECORF_4e_cosTheta1_spin[i]=-999.;
      RECORF_4e_cosTheta2_spin[i]=-999.;
      RECORF_4e_cosThetaStar_spin[i]=-999.;
      RECORF_4e_Phi_spin[i]=-999.;
      RECORF_4e_Phi1_spin[i]=-999.;
      RECORF_4e_Phi2_spin[i]=-999.;
      RECORF_4e_phi1RF_spin[i]=-999.;
      RECORF_4e_phi2RF_spin[i]=-999.;
      RECORF_4e_MELA[i]=-999.;
   
    }

    for (int i=0; i<10;i++){
      MCRF_cosTheta1_spin[i]=-999.;
      MCRF_cosTheta2_spin[i]=-999.;
      MCRF_cosThetaStar_spin[i]=-999.;
      MCRF_Phi_spin[i]=-999.;
      MCRF_Phi1_spin[i]=-999.;
      MCRF_Phi2_spin[i]=-999.;
      MCRF_phi1RF_spin[i]=-999.;
      MCRF_phi2RF_spin[i]=-999.;          	       
      MCRF_MELA[i]=-999.; 
   
    }

    genmet=-999.,calomet=-999.;  
    RECO_NMU=0,RECO_NELE=0;
    RECO_NTRACK=0;
    
    
    RECO_NPHOT=0,RECO_NJET=0,RECO_NVTX=0;
    RECO_NPFPHOT=0;

    for (int i=0; i<20;i++){
      RECOPHOT_PT[i]=-999.;
      RECOPHOT_ETA[i]=-999.;
      RECOPHOT_PHI[i]=-999.;
      RECOPHOT_THETA[i]=-999.;
      RECOPHOT_TLE_ParentSC_X[i]=-999.;
      RECOPHOT_TLE_ParentSC_Y[i]=-999.;
      RECOPHOT_TLE_ParentSC_Z[i]=-999.;
      RECOPFPHOT_PT[i]=-999.;
      RECOPFPHOT_PTError[i]=-999.;
      RECOPFPHOT_ETA[i]=-999.;
      RECOPFPHOT_PHI[i]=-999.;
      RECOPFPHOT_THETA[i]=-999.;

      RECOPFPHOT_PFX_rho[i]=-999.;

      RECOPFPHOT_PFchAllPart[i]=-999.;
      RECOPFPHOT_PFchHad[i]=-999.;
      RECOPFPHOT_PFneuHad[i]=999.;
      RECOPFPHOT_PFphoton[i]=-999.;
      RECOPFPHOT_PFPUchAllPart[i]=-999.;
      RECOPFPHOT_PFX_rho[i]=-999.;


    }
		
    for(int ivtx=0;ivtx<15;ivtx++) {
      RECO_VERTEX_x[ivtx] = -999;
      RECO_VERTEX_y[ivtx] = -999;
      RECO_VERTEX_z[ivtx] = -999;
      RECO_VERTEXPROB[ivtx]=-999.;
      RECO_VERTEX_ndof[ivtx] = -999;
      RECO_VERTEX_chi2[ivtx] = -999;
      RECO_VERTEX_ntracks[ivtx]=-999;
      for (int i=0;i<100;i++){
	RECO_VERTEX_TRACK_PT[ivtx][i]=-999;
      }
    }
    
    RECO_PFJET_N = 0;
    for (int ijets=0;ijets<200;ijets++) {
      RECO_PFJET_CHARGE[ijets]   = -999.; 
      RECO_PFJET_ET[ijets]   = -999.; 
      RECO_PFJET_PT[ijets]  = -999.; 
      RECO_PFJET_ETA[ijets] = -999.; 
      RECO_PFJET_PHI[ijets] = -999.;
      RECO_PFJET_PUID[ijets] = -999;
      RECO_PFJET_PUID_MVA[ijets] = -999.;
    }
    
/*     calomet=-999.; */
/*     calometopt=-999.; */
/*     calometoptnohf=-999.; */
/*     calometoptnohfho=-999.; */
/*     calometoptho=-999.; */
/*     calometnohf=-999.; */
/*     calometnohfho=-999.; */
/*     calometho=-999.; */
    pfmet=-999.;
    pfmet_x=-999.;
    pfmet_y=-999.;
    pfmet_phi=-999.;
    pfmet_theta=-999.;
    tcmet=-999.;
    cormetmuons=-999.;
    
    BeamSpot_X=-999.;
    BeamSpot_Y=-999.;
    BeamSpot_Z=-999.;
    
    
    for (int i=0; i<100;i++){
      MC_GENJET_PT[i]=-999;
      MC_GENJET_ETA[i]=-999;
      MC_GENJET_PHI[i]=-999;

      ftsigma[i]=-999.;
      ftsigmalag[i]=-999.;
      gdX[i]=-999.;
      gdY[i]=-999.;
      gdZ[i]=-999.;
      gdlagX[i]=-999.;
      gdlagY[i]=-999.;
      gdlagZ[i]=-999.;
      gdlagProb[i]=-999.;
      gdlagNdof[i]=-999.;
      ftsigmaMMMM[i]=-999.;
      ftsigmalagMMMM[i]=-999.;
      gdXMMMM[i]=-999.;
      gdYMMMM[i]=-999.;
      gdZMMMM[i]=-999.;
      gdlagXMMMM[i]=-999.;
      gdlagYMMMM[i]=-999.;
      gdlagZMMMM[i]=-999.;
      gdlagProbMMMM[i]=-999.;
      gdlagNdofMMMM[i]=-999.;
      ftsigmaEEEE[i]=-999.;
      ftsigmalagEEEE[i]=-999.;
      gdXEEEE[i]=-999.;
      gdYEEEE[i]=-999.;
      gdZEEEE[i]=-999.;
      gdlagXEEEE[i]=-999.;
      gdlagYEEEE[i]=-999.;
      gdlagZEEEE[i]=-999.;
      gdlagProbEEEE[i]=-999.;
      gdlagNdofEEEE[i]=-999.;

      StdFitVertexX[i]=-999.;
      StdFitVertexY[i]=-999.;
      StdFitVertexZ[i]=-999.;
      StdFitVertexChi2r[i]=-999.;
      StdFitVertexProb[i]=-999.;
      KinFitVertexX[i]=-999.;
      KinFitVertexY[i]=-999.;
      KinFitVertexZ[i]=-999.;
      KinFitVertexChi2r[i]=-999.;
      KinFitVertexProb[i]=-999.;

      StdFitVertexXMMMM[i]=-999.;
      StdFitVertexYMMMM[i]=-999.;
      StdFitVertexZMMMM[i]=-999.;
      StdFitVertexChi2rMMMM[i]=-999.;
      StdFitVertexProbMMMM[i]=-999.;
      KinFitVertexXMMMM[i]=-999.;
      KinFitVertexYMMMM[i]=-999.;
      KinFitVertexZMMMM[i]=-999.;
      KinFitVertexChi2rMMMM[i]=-999.;
      KinFitVertexProbMMMM[i]=-999.;

      StdFitVertexXEEEE[i]=-999.;
      StdFitVertexYEEEE[i]=-999.;
      StdFitVertexZEEEE[i]=-999.;
      StdFitVertexChi2rEEEE[i]=-999.;
      StdFitVertexProbEEEE[i]=-999.;
      KinFitVertexXEEEE[i]=-999.;
      KinFitVertexYEEEE[i]=-999.;
      KinFitVertexZEEEE[i]=-999.;
      KinFitVertexChi2rEEEE[i]=-999.;
      KinFitVertexProbEEEE[i]=-999.;

      StdFitVertexChi2rMMM[i]=-999.;
      StdFitVertexProbMMM[i]=-999.;
      StdFitVertexChi2rMME[i]=-999.;
      StdFitVertexProbMME[i]=-999.;
      StdFitVertexChi2rEEE[i]=-999.;
      StdFitVertexProbEEE[i]=-999.;
      StdFitVertexChi2rMEE[i]=-999.;
      StdFitVertexProbMEE[i]=-999.;

      RECO_MMMM_MASS_REFIT[i]=-999.;
      RECO_EEMM_MASS_REFIT[i]=-999.;
      RECO_EEEE_MASS_REFIT[i]=-999.;

      for (int k=0;k<4;k++){
	StdFitVertexTrack_PT[k][i]=-999.;
	StdFitVertexTrack_ETA[k][i]=-999.;
	StdFitVertexTrack_PHI[k][i]=-999.;
	
	StdFitVertexTrackMMMM_PT[k][i]=-999.;
	StdFitVertexTrackMMMM_ETA[k][i]=-999.;
	StdFitVertexTrackMMMM_PHI[k][i]=-999.;
	
	StdFitVertexTrackEEEE_PT[k][i]=-999.;
	StdFitVertexTrackEEEE_ETA[k][i]=-999.;
	StdFitVertexTrackEEEE_PHI[k][i]=-999.;
      }
    }
    

    for (int i=0; i<40;i++){
      StdFitVertexChi2rDiLep[i]=-999.;
      StdFitVertexProbDiLep[i]=-999.;
    }
    
    for (int j=0; j<7;j++){
      for (int i=0; i<100;i++){
	RECO_MMMM_MASS[j][i]=-999.;
	RECO_MMMM_PT[j][i]=-999.;	
	RECO_EEMM_MASS[j][i]=-999.;
	RECO_EEMM_PT[j][i]=-999.;
	RECO_EEEE_MASS[j][i]=-999.;
	RECO_EEEE_PT[j][i]=-999.;

	RECO_MMMM_ETA[j][i]=-999.;
	RECO_MMMM_PHI[j][i]=-999.;
	RECO_EEMM_ETA[j][i]=-999.;
	RECO_EEMM_PHI[j][i]=-999.;
	RECO_EEEE_ETA[j][i]=-999.;
	RECO_EEEE_PHI[j][i]=-999.;

      }
    }

    for (int i=0; i<50;i++){

      RECO_ZMM_MASS[i]=-999.;
      RECO_ZMM_PT[0][i]=-999.;
      RECO_ZMM_PT[1][i]=-999.;
      RECO_ZMM_PT[2][i]=-999.;
      RECO_ZMM_ETA[0][i]=-999.;
      RECO_ZMM_ETA[1][i]=-999.;
      RECO_ZMM_ETA[2][i]=-999.;
      RECO_ZMM_PHI[0][i]=-999.;
      RECO_ZMM_PHI[1][i]=-999.;
      RECO_ZMM_PHI[2][i]=-999.;
      
      RECO_ZEE_MASS[i]=-999.;
      RECO_ZEE_PT[0][i]=-999.;
      RECO_ZEE_PT[1][i]=-999.;
      RECO_ZEE_PT[2][i]=-999.;
      RECO_ZEE_ETA[0][i]=-999.;
      RECO_ZEE_ETA[1][i]=-999.;
      RECO_ZEE_ETA[2][i]=-999.;
      RECO_ZEE_PHI[0][i]=-999.;
      RECO_ZEE_PHI[1][i]=-999.;
      RECO_ZEE_PHI[2][i]=-999.;
      
      RECO_ZMMss_MASS[i]=-999.;
      RECO_ZMMss_PT[0][i]=-999.;
      RECO_ZMMss_PT[1][i]=-999.;
      RECO_ZMMss_PT[2][i]=-999.;
      RECO_ZMMss_ETA[0][i]=-999.;
      RECO_ZMMss_ETA[1][i]=-999.;
      RECO_ZMMss_ETA[2][i]=-999.;
      RECO_ZMMss_PHI[0][i]=-999.;
      RECO_ZMMss_PHI[1][i]=-999.;
      RECO_ZMMss_PHI[2][i]=-999.;
      

      RECO_ZEEss_MASS[i]=-999.;
      RECO_ZEEss_PT[0][i]=-999.;
      RECO_ZEEss_PT[1][i]=-999.;
      RECO_ZEEss_PT[2][i]=-999.;
      RECO_ZEEss_ETA[0][i]=-999.;
      RECO_ZEEss_ETA[1][i]=-999.;
      RECO_ZEEss_ETA[2][i]=-999.;
      RECO_ZEEss_PHI[0][i]=-999.;
      RECO_ZEEss_PHI[1][i]=-999.;
      RECO_ZEEss_PHI[2][i]=-999.;


      RECO_ZEM_MASS[i]=-999.;
      RECO_ZEM_PT[0][i]=-999.;
      RECO_ZEM_PT[1][i]=-999.;
      RECO_ZEM_PT[2][i]=-999.;
      RECO_ZEM_ETA[0][i]=-999.;
      RECO_ZEM_ETA[1][i]=-999.;
      RECO_ZEM_ETA[2][i]=-999.;
      RECO_ZEM_PHI[0][i]=-999.;
      RECO_ZEM_PHI[1][i]=-999.;
      RECO_ZEM_PHI[2][i]=-999.;

      RECO_DiLep_MASS[i]=-999.;
      RECO_DiLep_PT[0][i]=-999.;
      RECO_DiLep_PT[1][i]=-999.;
      RECO_DiLep_PT[2][i]=-999.;
      RECO_DiLep_ETA[0][i]=-999.;
      RECO_DiLep_ETA[1][i]=-999.;
      RECO_DiLep_ETA[2][i]=-999.;
      RECO_DiLep_PHI[0][i]=-999.;
      RECO_DiLep_PHI[1][i]=-999.;
      RECO_DiLep_PHI[2][i]=-999.;
     
      RECO_LLL0_MASS[i]=-999;
      RECO_LLL1_MASS[i]=-999;
      RECO_LLL2_MASS[i]=-999;
      RECO_LLL3_MASS[i]=-999;
      for (int j=0; j<4;j++){
	RECO_LLL0_PT[j][i]=-999;
	RECO_LLL1_PT[j][i]=-999;
	RECO_LLL2_PT[j][i]=-999;
	RECO_LLL3_PT[j][i]=-999;
      }
      
    }

    for (int i=0; i<20;i++){
      RECO_LLLL0ss_MASS[i]=-999;
      RECO_LLLL1ss_MASS[i]=-999;
      RECO_LLLL2ss_MASS[i]=-999;

      RECO_LLLl0_MASS[i]=-999;
      RECO_LLLl1_MASS[i]=-999;
      
      for (int j=0; j<4;j++){
	RECO_LLLL1ss_PT[j][i]=-999;
	RECO_LLLL1ss_PT[j][i]=-999;
	RECO_LLLL2ss_PT[j][i]=-999;
	RECO_LLLl0_PT[j][i]=-999;
	RECO_LLLl1_PT[j][i]=-999;
      }
    }
    
    
    for (int j=0; j<7;j++){
      for (int i=0; i<50;i++){
	RECO_LLLL_MASS[j][i]=-999.;
	RECO_LLLL_PT[j][i]=-999.;
	RECO_LLLL_ETA[j][i]=-999.;
	RECO_LLLL_PHI[j][i]=-999.;
      }
    }

    for (int i=0; i<100;i++){
      RECOELE_E[i]     = -999.;
      RECOELE_PT[i]=-999.;
      RECOELE_PTError[i]=-999.;
      RECOELE_P[i]=-999.;
      RECOELE_PHI[i]=-999.;
      RECOELE_ETA[i]=-999.;
      RECOELE_THETA[i]=-999.;
      RECOELE_MASS[i]=-999.;
      RECOELE_CHARGE[i]=-999.;
      
      // Core attributes
      RECOELE_isEcalDriven[i]    = false;
      RECOELE_isTrackerDriven[i] = false;
      RECOELE_gsftrack_NPixHits[i]=-999.;
      RECOELE_gsftrack_NStripHits[i]=-999.;
      RECOELE_gsftrack_chi2[i]   = -999;
      RECOELE_gsftrack_dxyB[i]   = -999;
      RECOELE_gsftrack_dxy[i]    = -999;
      RECOELE_gsftrack_dxyError[i]    = -999;
      RECOELE_gsftrack_dzB[i]    = -999;
      RECOELE_gsftrack_dz[i]     = -999;
      RECOELE_gsftrack_dzError[i]   = -999;
      RECOELE_gsftrack_losthits[i]  = -999;
      RECOELE_gsftrack_validhits[i] = -999;
      RECOELE_gsftrack_expected_inner_hits[i]=-999;
      RECOELE_scl_E[i]   = -999;
      RECOELE_scl_Et[i]  = -999;
      RECOELE_scl_Eta[i] = -999;
      RECOELE_scl_Phi[i] = -999;
      RECOELE_ecalEnergy[i] = -999;

      // Track-Cluster matching attributes
      RECOELE_ep[i]             = -999.;
      RECOELE_eSeedp[i]         = -999.;
      RECOELE_eSeedpout[i]      = -999.;
      RECOELE_eElepout[i]       = -999.;
      //
      RECOELE_deltaEtaIn[i]     = -999.;
      RECOELE_deltaEtaSeed[i]   = -999.;
      RECOELE_deltaEtaEle[i]    = -999.;
      RECOELE_deltaPhiIn[i]     = -999.;
      RECOELE_deltaPhiSeed[i]   = -999.;
      RECOELE_deltaPhiEle[i]    = -999.;
      // Fiducial flags 
      RECOELE_isbarrel[i]    = 0;
      RECOELE_isendcap[i]    = 0;
      RECOELE_isGap[i]  = 0;
      RECOELE_isEBetaGap[i]  = 0;
      RECOELE_isEBphiGap[i]  = 0;
      RECOELE_isEEdeeGap[i]  = 0;
      RECOELE_isEEringGap[i] = 0;
      // Shower shape
      RECOELE_sigmaIetaIeta[i] = -999.;
      RECOELE_sigmaEtaEta[i]   = -999.;
      RECOELE_e15[i]           = -999.;
      RECOELE_e25max[i]        = -999.;
      RECOELE_e55[i]           = -999.;
      RECOELE_he[i]            = -999.;
      RECOELE_r9[i]            = -999.;
      // Particle flow
      RECOELE_mva[i] = -999.; 
      // Brem & Classifaction
      RECOELE_fbrem[i]  = -999.; 
      RECOELE_nbrems[i] = -999;
      RECOELE_Class[i]  = -999;
      RECOELE_fbrem_mode[i]=-999.;
      RECOELE_fbrem_mean[i]=-999.;
      // isolation
      RECOELE_EGMTRACKISO[i]=-999.; 
      RECOELE_EGMHCALISO[i]=-999.;
      RECOELE_EGMECALISO[i]=-999.;
      RECOELE_EGMX[i]=-999.;

      // PF isolation
      RECOELE_PFchAllPart[i]=-999.;
      RECOELE_PFchHad[i]=-999.;
      RECOELE_PFneuHad[i]=999.;
      RECOELE_PFphoton[i]=-999.;
      RECOELE_PFPUchAllPart[i]=-999.;
      RECOELE_PFX_dB[i]=-999.;
      RECOELE_PFX_rho[i]=-999.;

      // Electron regression
      RECOELE_regEnergy[i]=-999.;
      RECOELE_regEnergyError[i]=-999.;

      // TLE electrons
      RECOELE_TLE_ParentSC_X[i]=-999.;
      RECOELE_TLE_ParentSC_Y[i]=-999.;
      RECOELE_TLE_ParentSC_Z[i]=-999.;

      // IP
      RECOELE_IP[i]=-9999.;
      RECOELE_SIP[i]=-9999.;
      RECOELE_IPERROR[i]=-9999.;
      RECOELE_IP_KF[i]=-999.;
      RECOELE_SIP_KF[i]=-999.;
      RECOELE_IPERROR_KF[i]=-999.;
      RECOELE_SIP_GD[i]=-999.;
      RECOELE_SIP_GDEEEE[i]=-999.;
      RECOELE_SIP_Std[i]=-999.;
      RECOELE_SIP_StdEEEE[i]=-999.;
      RECOELE_SIP_Kin[i]=-999.;
      RECOELE_SIP_KinEEEE[i]=-999.;
      RECOELE_STIP[i]=-999.;
      RECOELE_TIP[i]=-999.;
      RECOELE_TIPERROR[i]=-999.;
      RECOELE_SLIP[i]=-999.;
      RECOELE_LIP[i]=-999.;
      RECOELE_LIPERROR[i]=-999.;
      ele_sclRawE[i]=-999. ;
      ele_sclX[i]=-999.; 
      ele_sclY[i]=-999.; 
      ele_sclZ[i]=-999.;
      ele_seedSubdet1[i]=-999.;
      ele_seedDphi1[i]=-999.; 
      ele_seedDrz1[i]=-999.;
      ele_seedSubdet2[i]=-999.;
      ele_seedDphi2[i]=-999.; 
      ele_seedDrz2[i]=-999.;
      ele_eidVeryLoose[i]=-999.; 
      ele_eidLoose[i]=-999.; 
      ele_eidMedium[i]=-999.; 
      ele_eidTight[i]=-999. ;
      ele_eidHZZVeryLoose[i]=-999.; 
      ele_eidHZZLoose[i]=-999.; 
      ele_eidHZZMedium[i]=-999.; 
      ele_eidHZZTight[i]=-999. ;
      RECOELE_mvaTrigV0[i]=-999.;
      RECOELE_mvaNonTrigV0[i]=-999.;
    
      // Conversion
      ConvMapDist[i]=-999.;
      ConvMapDcot[i]=-999.;
      
      // Muon block
      RECOMU_PT_MuHLTMatch[i]=-999.;
      RECOMU_ETA_MuHLTMatch[i]=-999.;
      RECOMU_PHI_MuHLTMatch[i]=-999.;
      RECOELE_PT_EleHLTMatch[i]=-999.;
      RECOELE_ETA_EleHLTMatch[i]=-999.;
      RECOELE_PHI_EleHLTMatch[i]=-999.;

      RECOMU_isPFMu[i]=false;
      RECOMU_isGlobalMu[i]=false;
      RECOMU_isStandAloneMu[i]=false;
      RECOMU_isTrackerMu[i]=false;
      RECOMU_isCaloMu[i]=false;
      RECOMU_isTrackerHighPtMu[i]=false;   
      RECOMU_isME0Muon[i]=false;

      RECOMU_E[i]=-999.;
      RECOMU_PT[i]=-999.;
      RECOMU_P[i]=-999.;
      RECOMU_PHI[i]=-999.;
      RECOMU_ETA[i]=-999.;
      RECOMU_THETA[i]=-999.;
      RECOMU_MASS[i]=-999.;      
      RECOMU_CHARGE[i]=-999.;
      
      for (int j=0; j<3; j++){
	for (int k=0;k<3;k++){
	  RECOMU_COV[i][j][k]=-999.;
	  RECOELE_COV[i][j][k]=-999.;
	}
      }

      RECOMU_TRACKISO[i]=-999.; 
      RECOMU_TRACKISO_SUMPT[i]=-999.; 
      RECOMU_ECALISO[i]=-999.;
      RECOMU_HCALISO[i]=-999.;
      RECOMU_X[i]=-999.;

      RECOMU_PFchHad[i]=-999.;
      RECOMU_PFneuHad[i]=999.;
      RECOMU_PFphoton[i]=-999.;
      RECOMU_PFX_dB[i]=-999.;
      RECOMU_PFX_rho[i]=-999.;

      RECOMU_IP[i]=-9999.;
      RECOMU_SIP[i]=-9999.;
      RECOMU_IPERROR[i]=-9999.;
      RECOMU_IP_KF[i]=-999.;
      RECOMU_SIP_KF[i]=-999.;
      RECOMU_IPERROR_KF[i]=-999.;

      RECOMU_SIP_GD[i]=-999.;
      RECOMU_SIP_GDMMMM[i]=-999.;
      RECOMU_SIP_Std[i]=-999.;
      RECOMU_SIP_StdMMMM[i]=-999.;
      RECOMU_SIP_Kin[i]=-999.;
      RECOMU_SIP_KinMMMM[i]=-999.;
      RECOMU_STIP[i]=-999.;
      RECOMU_TIP[i]=-999.;
      RECOMU_TIPERROR[i]=-999.;
      RECOMU_SLIP[i]=-999.;
      RECOMU_LIP[i]=-999.;
      RECOMU_LIPERROR[i]=-999.;
      
 

      RECOMU_numberOfMatches[i]=-999;
      RECOMU_numberOfMatchedStations[i]=-999.;
      RECOMU_caloCompatibility[i]=-999.;
      RECOMU_segmentCompatibility[i]=-999.;
      RECOMU_glbmuPromptTight[i]=false;
      
      RECOMU_mutrkPT[i]=-999.;
      RECOMU_mutrkPTError[i]=-999.;
      RECOMU_mutrkDxy[i]=-999.;
      RECOMU_mutrkDxyError[i]=-999.;
      RECOMU_mutrkDxyB[i]=-999.;
      RECOMU_mutrkDz[i]=-999.;
      RECOMU_mutrkDzError[i]=-999.;
      RECOMU_mutrkDzB[i]=-999.;
      RECOMU_mutrkChi2PerNdof[i]=-999.;
      RECOMU_mutrktrackerLayersWithMeasurement[i]=-999.;

      RECOMU_muInnertrkDxy[i]=-999.;
      RECOMU_muInnertrkDxyError[i]=-999.;
      RECOMU_muInnertrkDxyB[i]=-999.;
      RECOMU_muInnertrkDz[i]=-999.;
      RECOMU_muInnertrkDzError[i]=-999.;
      RECOMU_muInnertrkDzB[i]=-999.;
      RECOMU_muInnertrkChi2PerNdof[i]=-999.;
      RECOMU_muInnertrktrackerLayersWithMeasurement[i]=-999.;
      RECOMU_muInnertrkCharge[i]=-999.;
      RECOMU_muInnertrkNHits[i]=-999.;
      RECOMU_muInnertrkNPixHits[i]=-999.;
      RECOMU_muInnertrkNStripHits[i]=-999.;
      RECOMU_muInnertrkPT[i]=-999.;
      RECOMU_muInnertrkPTError[i]=-999.;
      
      RECOMU_mubesttrkType[i]=-999;
      RECOMU_mubesttrkDxy[i]=-999.;
      RECOMU_mubesttrkDxyB[i]=-999.;
      RECOMU_mubesttrkDxyError[i]=-999.;
      RECOMU_mubesttrkDz[i]=-999.;
      RECOMU_mubesttrkDzB[i]=-999.;
      RECOMU_mubesttrkDzError[i]=-999.;
      RECOMU_mubesttrkPTError[i]=-999.;

      RECOMU_mutrkCharge[i]=-999.;
      RECOMU_mutrkNHits[i]=-999.;
      RECOMU_mutrkNPixHits[i]=-999.;
      RECOMU_mutrkNMuonHits[i]=-999.;
      RECOMU_mutrkNStripHits[i]=-999.;
      RECOMU_trkmuArbitration[i]=false;
      RECOMU_trkmu2DCompatibilityLoose[i]=false;
      RECOMU_trkmu2DCompatibilityTight[i]=false;
      RECOMU_trkmuOneStationLoose[i]=false;
      RECOMU_trkmuOneStationTight[i]=false;
      RECOMU_trkmuLastStationLoose[i]=false;
      RECOMU_trkmuLastStationTight[i]=false;
      RECOMU_trkmuOneStationAngLoose[i]=false;
      RECOMU_trkmuOneStationAngTight[i]=false;
      RECOMU_trkmuLastStationAngLoose[i]=false;
      RECOMU_trkmuLastStationAngTight[i]=false;
      RECOMU_trkmuLastStationOptimizedLowPtLoose[i]=false;
      RECOMU_trkmuLastStationOptimizedLowPtTight[i]=false;

    }
    

    for (int i=0; i<50;i++){      
      //Matching 

      // Muons
      RECOMU_MatchingMCTruth[i]=false;
      RECOMU_MatchingMCpT[i]=-999.;
      RECOMU_MatchingMCEta[i]=-999.;
      RECOMU_MatchingMCPhi[i]=-999.;
      
      // Electrons
      RECOELE_MatchingMCTruth[i]=false;
      RECOELE_MatchingMCpT[i]=-999.;
      RECOELE_MatchingMCEta[i]=-999.;
      RECOELE_MatchingMCPhi[i]=-999.;


      //Gamma:
      RECOPHOT_MatchingMCTruth[i]=false;
      RECOPHOT_MatchingMCpT[i]=-999.;
      RECOPHOT_MatchingMCEta[i]=-999.;
      RECOPHOT_MatchingMCPhi[i]=-999.;

      
      //zToMuMu:
      RECOzMuMu_MatchingMCTruth[i]=false;
      RECOzMuMu_MatchingMCpT[i]=-999.;
      RECOzMuMu_MatchingMCmass[i]=-999.;
      RECOzMuMu_MatchingMCEta[i]=-999.;
      RECOzMuMu_MatchingMCPhi[i]=-999.;
      
      //zToEE:
      RECOzEE_MatchingMCTruth[i]=false;
      RECOzEE_MatchingMCpT[i]=-999.;
      RECOzEE_MatchingMCmass[i]=-999.;
    }

     for (int i=0; i<100;i++){
      //HTozzToMMMM:
      RECOHzzMMMM_MatchingMCTruth[i]=false;
      RECOHzzMMMM_MatchingMCpT[i]=-999.;
      RECOHzzMMMM_MatchingMCmass[i]=-999.;
      RECOHzzMMMM_MatchingMCEta[i]=-999.;
      //HTozzToEEMM:
      RECOHzzEEMM_MatchingMCTruth[i]=false;
      RECOHzzEEMM_MatchingMCpT[i]=-999.;
      RECOHzzEEMM_MatchingMCmass[i]=-999.;
      RECOHzzEEMM_MatchingMCEta[i]=-999.;
      //HTozzToEEEE:
      RECOHzzEEEE_MatchingMCTruth[i]=false;
      RECOHzzEEEE_MatchingMCpT[i]=-999.;
      RECOHzzEEEE_MatchingMCmass[i]=-999.;
      RECOHzzEEEE_MatchingMCEta[i]=-999.;
    }


    for (int i=0; i<200;i++){
      RECO_TRACK_PT[i]=-999.;
      RECO_TRACK_ETA[i]=-999.;
      RECO_TRACK_PHI[i]=-999.;
      RECO_TRACK_CHI2[i]=-999.;
      RECO_TRACK_CHI2RED[i]=-999.;
      RECO_TRACK_CHI2PROB[i]=-999.;
      RECO_TRACK_NHITS[i]=0;
      RECO_TRACK_DXY[i]=-999.;
      RECO_TRACK_DXYERR[i]=-999.;
      RECO_TRACK_DZ[i]=-999.;
      RECO_TRACK_DZERR[i]=-999.;   			       					
    }

    for (int i=0; i<50;i++){
      tCHighEff_BTagJet_PT[i]=-999.;
//      tCHighPur_BTagJet_PT[i]=-999.;
      cSV_BTagJet_PT[i]=-999.;
      tCHighEff_BTagJet_ETA[i]=-999.;
 //     tCHighPur_BTagJet_ETA[i]=-999.;
      cSV_BTagJet_ETA[i]=-999.;
      tCHighEff_BTagJet_PHI[i]=-999.;
 //     tCHighPur_BTagJet_PHI[i]=-999.;
      cSV_BTagJet_PHI[i]=-999.;
      tCHighEff_BTagJet_DISCR[i]=-999.;
 //     tCHighPur_BTagJet_DISCR[i]=-999.;
      cSV_BTagJet_DISCR[i]=-999.;
      cSV_BTagJet_ET[i]=-999.;
    }
    
  }
  
  void fillPU(const edm::Event& iEvent){
      edm::Handle<vector<PileupSummaryInfo> > PupInfo;
      iEvent.getByToken(PileupSrc_, PupInfo);

      if(!PupInfo.isValid()) return;
      for( vector<PileupSummaryInfo>::const_iterator cand = PupInfo->begin();cand != PupInfo->end(); ++ cand ) { 
	std::cout << " Pileup Information: bunchXing, nvtx: " << cand->getBunchCrossing() << " " << cand->getPU_NumInteractions() << std::endl;
	if (cand->getBunchCrossing() == 0) num_PU_vertices=cand->getTrueNumInteractions();;
	// num_PU_vertices=cand->getPU_NumInteractions(); in-time,out-of-time pileup
	PU_BunchCrossing=cand->getBunchCrossing();
      }	
  }

  void EventsMCReWeighting(const edm::Event& iEvent){

    // get the weight                                                                                                                                                     
    MC_weighting=0.;
    float EventWeight = 1.0;
    edm::Handle<GenEventInfoProduct> gen_ev_info;
    iEvent.getByToken(generator_, gen_ev_info);
    if(!gen_ev_info.isValid()) return;
    EventWeight = gen_ev_info->weight();
    //std::cout<<"mc_weight = "<< gen_ev_info->weight() <<std::endl;
                                                                                                                                                                        
    float mc_weight = ( EventWeight > 0 ) ? 1 : -1;
    //std::cout<<"mc_weight = "<< mc_weight <<std::endl;                                                                                                                  
    MC_weighting=mc_weight;

  }


  void fillHLTFired(const edm::Event& iEvent){
    edm::Handle<vector<std::string> > HLTfired_;
    iEvent.getByLabel(HLTInfoFired,HLTfired_);

    vector<string> HLTimported;
    string tmpstring="";

    for (vector<string>::const_iterator cand=HLTfired_->begin(); cand!=HLTfired_->end(); ++cand){
      unsigned int i=cand-HLTfired_->begin();
      HLTimported.push_back(cand->c_str());
      string newstr=HLTimported.at(i) + ":" + tmpstring;
      tmpstring=newstr;
    }

    cout << "HLTFiredString= " << tmpstring.c_str() << endl;
    if (!tmpstring.empty()) sprintf(HLTPathsFired,tmpstring.c_str());
        

  }


  void triggermatching(const edm::Event& iEvent){

    cout << "Start Trigger matching for muon" << endl;
    // check HLTrigger/Configuration/python/HLT_GRun_cff.py

//AOD    edm::Handle<trigger::TriggerEvent> handleTriggerEvent;
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    iEvent.getByToken(triggerObjects_, triggerObjects );
//    const trigger::TriggerObjectCollection & toc(handleTriggerEvent->getObjects());
    size_t nMuHLT =0, nEleHLT=0;


    //const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);    


    std::vector<pat::TriggerObjectStandAlone>  HLTMuMatched,HLTEleMatched;
    std::vector<string> HLTMuMatchedNames,HLTEleMatchedNames;
   
//MiniAOD
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
      // obj.unpackPathNames(names);
       for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
           std::string fullname = obj.filterLabels()[h];
           std::string name;
           size_t p = fullname.find_first_of(':');
      if ( p != std::string::npos) {
        name = fullname.substr(0, p);
      }
      else {
        name = fullname;
      } 
     

          if (name == triggerFilter.c_str()) {
            HLTMuMatched.push_back(obj);
            HLTMuMatchedNames.push_back(name);
            cout << "Matching " << triggerFilter.c_str()  << endl;
            nMuHLT++;
          }
          if (name == triggerEleFilter.c_str()) {
            HLTEleMatched.push_back(obj);
            HLTEleMatchedNames.push_back(name);
            cout << "Matching " << triggerEleFilter.c_str()  << endl;
            nEleHLT++;
          }
        }
      }
    


    edm::Handle<edm::View<pat::Muon> > MuCandidates;
    iEvent.getByToken(muonTag_, MuCandidates);
    float maxDeltaR_=0.2;
    float maxDPtRel_=1.0;
    int nMuHLTMatch=0;
    for (edm::View<pat::Muon>::const_iterator iCand = MuCandidates->begin(); iCand != MuCandidates->end(); ++iCand){
      unsigned int i=iCand-MuCandidates->begin();
      cout << "Muon with pt= " << iCand->pt() << ": check trigger matching" << endl;
      if (IsMuMatchedToHLTMu(*iCand,  HLTMuMatched , HLTMuMatchedNames, maxDeltaR_, maxDPtRel_)==true){
	nMuHLTMatch++;
	cout << "Muon HLT Matched with pT= " << iCand->pt() << endl;
	RECOMU_PT_MuHLTMatch[i] =iCand->pt();
	RECOMU_ETA_MuHLTMatch[i]=iCand->eta();
	RECOMU_PHI_MuHLTMatch[i]=iCand->phi();
      }
    }

    cout << "N. Muons HLT Matched= " << nMuHLTMatch << " FiredString:" << HLTPathsFired << endl;
    RECO_nMuHLTMatch    = nMuHLTMatch;


    cout << "Start Trigger matching for electron" << endl;

    int nEleHLTMatch=0;
    edm::Handle<edm::View<pat::Electron> > EleCandidates;
    iEvent.getByToken(electronEgmTag_, EleCandidates);

    for (edm::View<pat::Electron>::const_iterator iCand = EleCandidates->begin(); iCand != EleCandidates->end(); ++iCand){

      unsigned int i=iCand-EleCandidates->begin();
      cout << "Electron with pt= " << iCand->pt() << ": check trigger matching" << endl;
      if (IsEleMatchedToHLTEle(*iCand,  HLTEleMatched , HLTEleMatchedNames, maxDeltaR_, maxDPtRel_)==true){
	    cout << "Electron HLT Matched with pT= " << iCand->pt() << endl;
	    nEleHLTMatch++;
	    RECOELE_PT_EleHLTMatch[i]=iCand->pt();
	    RECOELE_ETA_EleHLTMatch[i]=iCand->eta();
	    RECOELE_PHI_EleHLTMatch[i]=iCand->phi();
      }
    }

    cout << "N. Electrons HLT Matched= " << nEleHLTMatch << endl;

    RECO_nEleHLTMatch = nEleHLTMatch;


  }

  bool IsMuMatchedToHLTMu ( const pat::Muon &mu, std::vector<pat::TriggerObjectStandAlone> HLTMu , std::vector<string> HLTMuNames, double DR, double DPtRel ) {
    size_t dim =  HLTMu.size();
    size_t nPass=0;
    if (dim==0) return false;
    for (size_t k =0; k< dim; k++ ) {
      //cout << "HLT mu filter is= " << HLTMuNames[k].c_str() << " Delta R= " << deltaR(HLTMu[k], mu) << " Delta pT= " << fabs(HLTMu[k].pt() - mu.pt())/ HLTMu[k].pt() << endl;
      if (  (deltaR(HLTMu[k], mu) < DR)   && (fabs(HLTMu[k].pt() - mu.pt())/ HLTMu[k].pt()<DPtRel)){ 
	cout << "HLT mu filter is= " << HLTMuNames[k].c_str() << " Delta R= " << deltaR(HLTMu[k], mu) << " Delta pT= " << fabs(HLTMu[k].pt() - mu.pt())/ HLTMu[k].pt() << endl;
	nPass++ ;
      }
    }
    return (nPass>0);
  }

  bool IsEleMatchedToHLTEle ( const pat::Electron &ele, std::vector<pat::TriggerObjectStandAlone> HLTEle , std::vector<string> HLTEleNames, double DR, double DPtRel ) {
    size_t dim =  HLTEle.size();
    size_t nPass=0;
    if (dim==0) return false;
    for (size_t k =0; k< dim; k++ ) {
      //cout << "HLT ele filter is= " << HLTEleNames[k].c_str() << " Delta R= " << deltaR(HLTEle[k], ele) << " Delta pT= " << fabs(HLTEle[k].pt() - ele.pt())/ HLTEle[k].pt() << endl;
      if (  (deltaR(HLTEle[k], ele) < DR)   && (fabs(HLTEle[k].pt() - ele.pt())/ HLTEle[k].pt()<DPtRel)){ 
	cout << "HLT ele filter is= " << HLTEleNames[k].c_str() << " Delta R= " << deltaR(HLTEle[k], ele) << " Delta pT= " << fabs(HLTEle[k].pt() - ele.pt())/ HLTEle[k].pt() << endl;
	nPass++ ;
      }
    }
    return (nPass>0);
  }

  bool isTrackerHighPtMu (const reco::Muon &mu, math::XYZPoint primaryVertex){
    return( mu.numberOfMatchedStations() > 1 &&
            (mu.muonBestTrack()->ptError()/mu.muonBestTrack()->pt()) < 0.3
            && std::abs(mu.muonBestTrack()->dxy(primaryVertex)) < 0.2
            && std::abs(mu.muonBestTrack()->dz(primaryVertex)) < 0.5
            && mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0
            && mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
            );
  }

  
  // PDT
  std::string getParticleName(int id) const{
    const ParticleData * pd = pdt_->particle( id );
    if (!pd) {
      std::ostringstream ss;
      ss << "P" << id;
      return ss.str();
    } else
      return pd->name();
  }

  
  void fillgenjets(const edm::Event& iEvent){
    edm::Handle<reco::GenJetCollection> genjetHandle;
    iEvent.getByToken(genjetTag_,genjetHandle);
    int i=0;
    for ( GenJetCollection::const_iterator igen=genjetHandle->begin(); igen!=genjetHandle->end(); igen++) {
      if (i>99) break;
      MC_GENJET_PT[i]=igen->pt();
      MC_GENJET_ETA[i]=igen->eta();
      MC_GENJET_PHI[i]=igen->phi();
      i++;
    }
  }

  // GenParticles  
  void fillgenparticles(const edm::Event& iEvent, const edm::EventSetup &es){    
    
    // get gen particle candidates 
    edm::Handle<reco::GenParticleCollection> genCandidates;           
    iEvent.getByToken(genParticles_, genCandidates);

    es.getData( pdt_ );
    vector<float> leptonpt,leptoneta,leptonphi,leptontheta,leptonpdgid;
    
    for ( GenParticleCollection::const_iterator mcIter=genCandidates->begin(); mcIter!=genCandidates->end(); ++mcIter ) {
      // lepton stable      
   
      if ( ( abs(mcIter->pdgId())==11 || abs(mcIter->pdgId())==13 || abs(mcIter->pdgId())==15  ) && mcIter->status()==1 ){	
	
	std::cout << " \n Found lepton with Id= " << mcIter->pdgId() << "  in the final state with mother= " ; 
	leptonpt.push_back(mcIter->pt());
	leptoneta.push_back(mcIter->eta());
	leptonphi.push_back(mcIter->phi());
	leptontheta.push_back(mcIter->theta());
        leptonpdgid.push_back(float(mcIter->pdgId()));

	if (mcIter->numberOfMothers()>0){ 
	  std::cout << getParticleName(int(mcIter->mother(0)->pdgId())) << " << ";
	  if (mcIter->mother(0)->status()==3) continue;
	  if (mcIter->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->status()!=3 ){
	    std::cout << getParticleName(int(mcIter->mother(0)->mother(0)->pdgId())) << " << ";
	    if (mcIter->mother(0)->mother(0)->status()==3) continue;
	    if (mcIter->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->status()!=3 ){
	      std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
	      if (mcIter->mother(0)->mother(0)->mother(0)->status()==3) continue;
	      if (mcIter->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
		if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		  std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
		  if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		  if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && 
mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		    std::cout << getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << " << std::endl;
		    if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		  }
		}
	      }
	    }
	  }
	}

    
      }
    }

    sort(leptonpt.rbegin(),leptonpt.rend());
    
    if (leptonpt.size()>0) { 
      cout << "\n The 4 highest pt leptons in MCtruth are= ";
      for (unsigned int i=0;i<leptonpt.size();i++){
	if (i>3) continue;
	cout << leptonpt.at(i) << " ";
	MC_LEPT_PT[i]=leptonpt.at(i);
	MC_LEPT_ETA[i]=leptoneta.at(i);
	MC_LEPT_PHI[i]=leptonphi.at(i);
	MC_LEPT_THETA[i]=leptontheta.at(i);
        MC_LEPT_PDGID[i]=leptonpdgid.at(i);
      }
      cout << endl;
    }

    // Check if Z bosons are generated (specially for Z + jet MC)
    int i=0;
      //j=0;
    for ( GenParticleCollection::const_iterator mcIter=genCandidates->begin(); mcIter!=genCandidates->end(); ++mcIter ) {
      // Select the Z decay: status 3 
      if ( abs(mcIter->pdgId())==23 && mcIter->status()==3 ) {
	int j=0;
	std::cout << "\n Found generated Z with PDGId = " << mcIter->pdgId() << " and status = "<< mcIter->status() << " and " << mcIter->numberOfDaughters() << " daughts |";
	std::cout << " pt = "<<std::setw(12)<<mcIter->pt()<<" GeV/c | eta = "<<std::setw(12)<<mcIter->eta()<<" | phi = "<<std::setw(12)<<mcIter->phi();
	std::cout << "   ===> Filled Z Bo at ["<<i<<"]["<<j<<"]";
	MC_Z_PT[i][0]=mcIter->pt();       MC_Z_ETA[i][0]=mcIter->eta();   MC_Z_PHI[i][0]=mcIter->phi();
        MC_Z_THETA[i][0]=mcIter->theta(); MC_Z_MASS[i][0]=mcIter->mass(); MC_Z_PDGID[i][0]=mcIter->pdgId();  
	
	// Check daughters of Z decay: status 3
	reco::GenParticle::daughters dst3 = mcIter->daughterRefVector();
	for (reco::GenParticle::daughters::const_iterator it_dst3 = dst3.begin(), est3 = dst3.end(); it_dst3 != est3; ++it_dst3) {
	  // Select status 3 electron, muon or tau
	  if (((abs((**it_dst3).pdgId()) == 11) || (abs((**it_dst3).pdgId()) == 13) || (abs((**it_dst3).pdgId()) == 15)) && (**it_dst3).status() == 3) {
	    std::cout<<"\n |--> Z daughter Particle  "<<std::setw(18)<<"| id = "<<std::setw(5)<<(**it_dst3).pdgId()<<" | st = "<<std::setw(5)<<(**it_dst3).status()<<" | pt = ";
	    std::cout<<std::setw(12)<<(**it_dst3).pt()<<" GeV/c | eta = "<<std::setw(12)<<(**it_dst3).eta()<<" | phi = "<<std::setw(12)<<(**it_dst3).phi();
	    // daughters of the status 3 electron, muon or tau
	    reco::GenParticle::daughters dst1 = (**it_dst3).daughterRefVector();
	    bool flag_emu_found = false; 
	    // Ask for status 1 daughters as long as no status 1 daughters are found 
	    while (flag_emu_found == false) {
	      reco::GenParticle::daughters::const_iterator it_dst1 = dst1.begin(), e_dst1 = dst1.end();
	      for (it_dst1 = dst1.begin(); it_dst1 != e_dst1; ++it_dst1) {
		std::cout<<"\n |--> Z daughter Particle  "<<std::setw(18)<<"| id = "<<std::setw(5)<<(**it_dst1).pdgId()<<" | st = "<<std::setw(5)<<(**it_dst1).status()<<" | pt = ";
		std::cout<<std::setw(12)<<(**it_dst1).pt()<<" GeV/c | eta = "<<std::setw(12)<<(**it_dst1).eta()<<" | phi = "<<std::setw(12)<<(**it_dst1).phi();
		if (((abs((**it_dst1).pdgId()) == 11) || (abs((**it_dst1).pdgId()) == 13)) && (**it_dst1).status() == 1) {
		  flag_emu_found = true;
		  ++j;
		  std::cout<<"   ===> Filled E/MU at ["<<i<<"]["<<j<<"]";
		  MC_Z_PT[i][j]=(**it_dst1).pt();          MC_Z_ETA[i][j]=(**it_dst1).eta();      MC_Z_PHI[i][j]=(**it_dst1).phi();
		  MC_Z_THETA[i][j]=(**it_dst1).theta();    MC_Z_MASS[i][j]=(**it_dst1).mass();    MC_Z_PDGID[i][j]=(**it_dst1).pdgId();
		}
		else if (((abs((**it_dst1).pdgId()) == 11) || (abs((**it_dst1).pdgId()) == 13) || (abs((**it_dst1).pdgId()) == 15) || 
			  (abs((**it_dst1).pdgId()) == 24)) && (**it_dst1).status() == 2) {
		  // std::cout<<"   ===> Not Status 1 || Provided daughterRefVector";
		  dst1 = (**it_dst1).daughterRefVector();
		}
		else if (abs((**it_dst1).pdgId()) > 22 ) { // exotic or hadronic tau detected (tau --> nu W)
		  // std::cout<<"   ===> Hadronic tau || Quit Loop";
		  flag_emu_found = true; // quit loop || no e mu candidate found
		}
		else if ((abs((**it_dst1).pdgId()) == 22) && ((**it_dst1).status() == 1)) { // FSR Photon
		  int ph_j = 0; if (flag_emu_found == true) ph_j = j+2; else ph_j = j+3;
		  if ( (**it_dst1).pt() > MC_Z_PT[i][ph_j] ) { // Save only highest pt photon
		    std::cout<<"   ===> Filled PHOT at ["<<i<<"]["<<ph_j<<"]";
		    MC_Z_PT[i][ph_j]=(**it_dst1).pt();          MC_Z_ETA[i][ph_j]=(**it_dst1).eta();      MC_Z_PHI[i][ph_j]=(**it_dst1).phi();
		    MC_Z_THETA[i][ph_j]=(**it_dst1).theta();    MC_Z_MASS[i][ph_j]=(**it_dst1).mass();    MC_Z_PDGID[i][ph_j]=(**it_dst1).pdgId();
		  }
		  else {
		    std::cout<<"   ===> NOT Filled, pt = "<<MC_Z_PT[i][ph_j]<<" GeV/c";
		  }
		}
		else {} // other particle than e or mu or photon
	      }
	    }
	  }
	}
	++i;
	std::cout<<"\n-----------------------------------------------------------------"<<std::endl;
      }
    }
    std::cout<<"\n"<<std::endl;




    // get 4l candidates
    i=0; 
    //j=0;
    edm::Handle<edm::View<Candidate> >  fourlCandidates;
    iEvent.getByToken(fourgenleptons_, fourlCandidates);
    for (edm::View<Candidate>::const_iterator mcIter=fourlCandidates->begin(); mcIter!=fourlCandidates->end(); ++mcIter ) {
      if (i>49 ) continue;
      cout << " MC 4l Mass= " << mcIter->mass()
	   << " Charge= " 
	   << mcIter->daughter(0)->daughter(0)->charge() << " " 
	   << mcIter->daughter(0)->daughter(1)->charge() << " "
	   << mcIter->daughter(0)->daughter(2)->charge() << " " 
	   << mcIter->daughter(1)->charge() << " " 
	   << endl;
      MC_fourl_MASS[i][0]=mcIter->p4().mass();
      MC_fourl_PT[i][0]=mcIter->p4().pt();
      MC_fourl_PDGID[i][0]=mcIter->pdgId();

      int ii=0; // l=0;
      for (unsigned j = 0; j < mcIter->numberOfDaughters(); ++j ) {
	//cout << "j= " << j << " " << abs(mcIter->daughter(j)->pdgId()) << endl;
	if (j==0){
	  for (unsigned k = 0; k < mcIter->daughter(j)->numberOfDaughters(); ++k ) {
	    //cout << abs(mcIter->daughter(j)->daughter(k)->pdgId()) << endl;
	    if ( abs(mcIter->daughter(j)->daughter(k)->pdgId())==13 || 
		 abs(mcIter->daughter(j)->daughter(k)->pdgId())==15 || 
		 abs(mcIter->daughter(j)->daughter(k)->pdgId())==11){
	      
	      //cout << "k+1= " << k+1 << endl;
	      MC_fourl_MASS[i][k+1]=mcIter->daughter(j)->daughter(k)->p4().mass();
	      MC_fourl_PT[i][k+1]=mcIter->daughter(j)->daughter(k)->p4().pt();
	      MC_fourl_PDGID[i][k+1]=mcIter->daughter(j)->daughter(k)->pdgId();
	      
	    }
	  }
	}
	
	if (j==1) {
	  //cout << "k+1= " << mcIter->daughter(0)->numberOfDaughters()+1 << endl;
	  MC_fourl_MASS[i][mcIter->daughter(0)->numberOfDaughters()+1]=mcIter->daughter(j)->p4().mass();
	  MC_fourl_PT[i][mcIter->daughter(0)->numberOfDaughters()+1]=mcIter->daughter(j)->p4().pt();
	  MC_fourl_PDGID[i][mcIter->daughter(0)->numberOfDaughters()+1]=mcIter->daughter(j)->pdgId();
	}

	ii++;
	
      }

      i++;
    }
    
    // get ZZ candidates
    i =0;
    edm::Handle<edm::View<Candidate> >  ZZCandidates;
    iEvent.getByToken(digenZ_, ZZCandidates);
    
    for (edm::View<Candidate>::const_iterator mcIterZZ=ZZCandidates->begin(); mcIterZZ!=ZZCandidates->end(); ++mcIterZZ ) {
      if (i>3 ) continue;
      cout << "MC ZZ Mass= " << mcIterZZ->p4().mass() 
	   << " and pT= " << mcIterZZ->p4().pt()  
	   << endl;
      
      
      MC_ZZ_MASS[i][0]   = mcIterZZ->p4().mass();
      MC_ZZ_PT[i][0]     = mcIterZZ->p4().pt();
      MC_ZZ_ETA[i][0]    = mcIterZZ->p4().eta();
      MC_ZZ_PHI[i][0]    = mcIterZZ->p4().phi();
      MC_ZZ_THETA[i][0]  = mcIterZZ->p4().theta();
      MC_ZZ_PDGID[i][0]  = mcIterZZ->pdgId();
      
      int ii=0,l=0;
      for (unsigned j = 0; j < mcIterZZ->numberOfDaughters(); ++j ) {
	// cout << " j= " << j << " " << abs(mcIterZZ->daughter(j)->pdgId()) << endl;
	
	MC_ZZ_MASS[i][j+1] = mcIterZZ->daughter(j)->p4().mass();
	MC_ZZ_PT[i][j+1]   = mcIterZZ->daughter(j)->p4().pt();
	MC_ZZ_ETA[i][j+1]  = mcIterZZ->daughter(j)->p4().eta();
	MC_ZZ_PHI[i][j+1]  = mcIterZZ->daughter(j)->p4().phi();
	MC_ZZ_THETA[i][j+1]= mcIterZZ->daughter(j)->p4().theta();
	MC_ZZ_PDGID[i][j+1]= mcIterZZ->daughter(j)->pdgId();
	
	//cout << mcIterZZ->daughter(j)->numberOfDaughters()<< endl;
	
	int kk=0;
	for (unsigned k = 0; k < mcIterZZ->daughter(j)->numberOfDaughters(); ++k ) {
	  // cout << " k= " << k << abs(mcIterZZ->daughter(j)->daughter(k)->pdgId()) << endl;
	  if ( abs(mcIterZZ->daughter(j)->daughter(k)->pdgId())==13 || 
	       abs(mcIterZZ->daughter(j)->daughter(k)->pdgId())==15 || 
	       abs(mcIterZZ->daughter(j)->daughter(k)->pdgId())==11){
	    l=ii+j+kk+3; 	      
	    
	    //cout << " l= " << l << " " << abs(mcIterZZ->daughter(j)->daughter(k)->pdgId()) << endl;
	    MC_ZZ_MASS[i][l] = mcIterZZ->daughter(j)->daughter(k)->p4().mass();
	    MC_ZZ_PT[i][l]   = mcIterZZ->daughter(j)->daughter(k)->p4().pt();
	    MC_ZZ_ETA[i][l]  = mcIterZZ->daughter(j)->daughter(k)->p4().eta();
	    MC_ZZ_PHI[i][l]  = mcIterZZ->daughter(j)->daughter(k)->p4().phi();
	    MC_ZZ_THETA[i][l]= mcIterZZ->daughter(j)->daughter(k)->p4().theta();
	    MC_ZZ_PDGID[i][l]= mcIterZZ->daughter(j)->daughter(k)->pdgId();
	    kk++;
	  }
	}
	ii++;	    	  
      }
    }
      
    
  }
  

  // MC Higgs 
  void fillmc(const edm::Event& iEvent){
    edm::Handle<edm::View<Candidate> > Candidates;
    iEvent.getByToken(MCcollName, Candidates);
    cout << "running fillmc" << endl;
    cout << "size= " << Candidates->end()-Candidates->begin() << endl;
    for( edm::View<Candidate>::const_iterator cand = Candidates->begin();cand != Candidates->end(); ++ cand ) { 
      cout << "Filling MC variables" << endl;
      const reco::Candidate& theParticle = (*cand);
      SetMCValues(theParticle,0);   
      int i=0,l=0;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	if (cand->daughter(j)->pdgId()==23){	  
	  const reco::Candidate& theParticle = (*cand->daughter(j));
	  SetMCValues(theParticle,j+1);
	  for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	    if ( abs(cand->daughter(j)->daughter(k)->pdgId())==13 || 
		 abs(cand->daughter(j)->daughter(k)->pdgId())==15 || 
		 abs(cand->daughter(j)->daughter(k)->pdgId())==11){
	      l=i+j+k+3; 	      
	      const reco::Candidate& theParticle = (*cand->daughter(j)->daughter(k));
	      SetMCValues(theParticle,l);
	    }
	  }
	  i++;
	}
      }
    }    
    //fillMCCP(iEvent);
  }

  
  
  struct SortCandByDecreasingPt {
    bool operator()( const reco::Candidate &c1, const reco::Candidate &c2) const {
      return c1.pt() > c2.pt();
    }
  };
  

  
  void fillAdditionalRECO(const edm::Event& iEvent){
   
    //Matching ZtoMuMu:
    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchZMuMu;
    iEvent.getByToken(goodZtoMuMuMCMatch_, GenParticlesMatchZMuMu);
    if (GenParticlesMatchZMuMu.isValid() ){
      cout << "The matched map collection has size= " <<  GenParticlesMatchZMuMu->size() << endl;
    }
    
    //Matching ZtoEE:
    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchZEE;
    iEvent.getByToken(goodZtoEEMCMatch_, GenParticlesMatchZEE);
    if (GenParticlesMatchZEE.isValid() ){
      cout << "The matched map collection has size= " <<  GenParticlesMatchZEE->size() << endl;
    }

    // di-leptons OS
    leptonscands_Z0->clear();
    leptonscands_Z1->clear();
   
    
    cout << "RECOcollNameZ size " << RECOcollNameZ.size() << endl;
    cout << "RECOcollNameZ string" << RECOcollNameZ.at(0) << endl;
    for (unsigned int i=0; i<RECOcollNameZ.size(); i++) {  
      edm::Handle<edm::View<Candidate> > CandidatesZ;
      if(i==0) iEvent.getByToken(zToMuMu, CandidatesZ); 
      else if(i==1) iEvent.getByToken(zToEE, CandidatesZ);
      int kk=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesZ->begin();cand != CandidatesZ->end(); ++ cand ) { 
	if (kk>49) continue;
	if (i==0) RECO_ZMM_MASS[kk]=cand->p4().mass();
	if (i==0) RECO_ZMM_PT[0][kk]=cand->p4().pt();
	if (i==0) RECO_ZMM_ETA[0][kk]=cand->p4().eta();
	if (i==0) RECO_ZMM_PHI[0][kk]=cand->p4().phi();
	if (i==1) RECO_ZEE_MASS[kk]=cand->p4().mass();
	if (i==1) RECO_ZEE_PT[0][kk]=cand->p4().pt();
	if (i==1) RECO_ZEE_ETA[0][kk]=cand->p4().eta();
	if (i==1) RECO_ZEE_PHI[0][kk]=cand->p4().phi();
	cout << "di-lepton candidate of type=" << RECOcollNameZ.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl; 
	  if (i==0) {
	    RECO_ZMM_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZMM_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZMM_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Z0->push_back( cand->daughter(j)->clone());
	  }
	  if (i==1) {
	    RECO_ZEE_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZEE_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZEE_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Z1->push_back( cand->daughter(j)->clone());
	  }
	}
	
	if (i==0 && fillMCTruth==true){
	  // Matching ZtoMuMu
	  edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesZ,kk);
	  edm::Ref<std::vector<reco::GenParticle> > genrefZMuMu = (*GenParticlesMatchZMuMu)[Ref]; 
	  if (!genrefZMuMu.isNull()){
	    cout << "Gen Z with pT= " << genrefZMuMu->p4().pt() << " and mass="<< genrefZMuMu->p4().mass() << endl;	  
	    RECOzMuMu_MatchingMCTruth[kk]= true;
	    RECOzMuMu_MatchingMCpT[kk]= genrefZMuMu->p4().pt();
	    RECOzMuMu_MatchingMCmass[kk]= genrefZMuMu->p4().mass();
	    RECOzMuMu_MatchingMCEta[kk]= genrefZMuMu->p4().eta();
	    RECOzMuMu_MatchingMCPhi[kk]= genrefZMuMu->p4().phi();
	  }
	}
	if (i==1 && fillMCTruth==true){
	  // Matching ZtoEE
	  edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesZ,kk);
	  edm::Ref<std::vector<reco::GenParticle> > genrefZEE = (*GenParticlesMatchZEE)[Ref]; 
	  if (!genrefZEE.isNull()){
	    cout << "Gen Z with pT= " << genrefZEE->p4().pt() << " and mass="<< genrefZEE->p4().mass() << endl;	  
	    RECOzEE_MatchingMCTruth[kk]= true;
	    RECOzEE_MatchingMCpT[kk]= genrefZEE->p4().pt();
	    RECOzEE_MatchingMCmass[kk]= genrefZEE->p4().mass();
	    RECOzEE_MatchingMCEta[kk]= genrefZEE->p4().eta();
	    RECOzEE_MatchingMCPhi[kk]= genrefZEE->p4().phi();
	  }
	}


	kk++;
      }
    }

    // di-leptons SS and cross-leptons
    leptonscands_Zss0->clear();
    leptonscands_Zss1->clear();
    leptonscands_Zcross->clear();


    cout << "RECOcollNameZss size " << RECOcollNameZss.size() << endl; 
    for (unsigned int i=0; i<RECOcollNameZss.size(); i++) {  
      edm::Handle<edm::View<Candidate> > CandidatesZss;
      if(i==0)   iEvent.getByToken(zToMuMussmerge, CandidatesZss);
      if(i==1)   iEvent.getByToken(zToEEssmerge, CandidatesZss);
      if(i==2)   iEvent.getByToken(zToCrossLeptons, CandidatesZss); 
      int kk=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesZss->begin();cand != CandidatesZss->end(); ++ cand ) { 
	if (kk>49) continue;
	if (i==0)  RECO_ZMMss_MASS[kk]=cand->p4().mass();
	if (i==0)  RECO_ZMMss_PT[0][kk]=cand->p4().pt();
	if (i==0)  RECO_ZMMss_ETA[0][kk]=cand->p4().eta();
	if (i==0)  RECO_ZMMss_PHI[0][kk]=cand->p4().phi();

	if (i==1)  RECO_ZEEss_MASS[kk]=cand->p4().mass();
	if (i==1)  RECO_ZEEss_PT[0][kk]=cand->p4().pt();
	if (i==1)  RECO_ZEEss_ETA[0][kk]=cand->p4().eta();
	if (i==1)  RECO_ZEEss_PHI[0][kk]=cand->p4().phi();

	if (i==2)  RECO_ZEM_MASS[kk]=cand->p4().mass();
	if (i==2)  RECO_ZEM_PT[0][kk]=cand->p4().pt();
	if (i==2)  RECO_ZEM_ETA[0][kk]=cand->p4().eta();
	if (i==2)  RECO_ZEM_PHI[0][kk]=cand->p4().phi();

	cout << "di-lepton candidate of type=" << RECOcollNameZss.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl; 
	  if (i==0) {
	    RECO_ZMMss_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZMMss_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZMMss_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Zss0->push_back( cand->daughter(j)->clone());
	  }
	  if (i==1) {
	    RECO_ZEEss_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZEEss_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZEEss_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Zss1->push_back( cand->daughter(j)->clone());
	  }
	  if (i==2) {
	    RECO_ZEM_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZEM_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZEM_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Zcross->push_back( cand->daughter(j)->clone());	  
	  }
	}
	kk++;
      }
    }


    // di-leptons ALL
    leptonscands_DiLep->clear();

    edm::Handle<edm::View<Candidate> > CandidatesDiLep;
    iEvent.getByToken(RECOcollNameDiLep_, CandidatesDiLep); 
    int kkk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesDiLep->begin();cand != CandidatesDiLep->end(); ++ cand ) { 
      if (kkk>49) continue;
      RECO_DiLep_MASS[kkk]=cand->p4().mass();
      RECO_DiLep_PT[0][kkk]=cand->p4().pt();
      RECO_DiLep_ETA[0][kkk]=cand->p4().eta();
      RECO_DiLep_PHI[0][kkk]=cand->p4().phi();
      
      cout << "di-lepton candidate of type=" << RECOcollNameDiLep.label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl;
	RECO_DiLep_PT[j+1][kkk]=cand->daughter(j)->p4().pt();
	RECO_DiLep_ETA[j+1][kkk]=cand->daughter(j)->p4().eta();
	RECO_DiLep_PHI[j+1][kkk]=cand->daughter(j)->p4().phi();
	leptonscands_DiLep->push_back( cand->daughter(j)->clone());
      }
      kkk++;
    }



    // MuMuMuMu
    int i=1; 
    leptonscands_MMMM->clear();
    edm::Handle<edm::View<Candidate> > CandidatesMMMM;
    iEvent.getByToken(RECOcollNameMMMM_, CandidatesMMMM);

    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchHMMMM;
    iEvent.getByToken(goodHiggsTozzToMMMMMCMatch_, GenParticlesMatchHMMMM);
    if (GenParticlesMatchHMMMM.isValid()){
      cout << "The matched map collection has size= " <<  GenParticlesMatchHMMMM->size() << endl;
    }

    int kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesMMMM->begin();cand != CandidatesMMMM->end(); ++ cand ) { 
      if (kk>99) break;
      RECO_MMMM_MASS[i-1][kk]=cand->p4().mass();
      RECO_MMMM_PT[i-1][kk]=cand->p4().pt();
      RECO_MMMM_ETA[i-1][kk]=cand->p4().eta();
      RECO_MMMM_PHI[i-1][kk]=cand->p4().phi();
      int l=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_MMMM_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_MMMM_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_MMMM_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_MMMM_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_MMMM_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_MMMM_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_MMMM_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_MMMM_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  leptonscands_MMMM->push_back( cand->daughter(j)->daughter(k)->clone());
	  //cout << "index" << i+j+k+l+2 <<endl;
	} 
	l++;
      }


      // Matching HiggsToMMMM
      if (fillMCTruth==true){
	edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesMMMM,kk);
	edm::Ref<std::vector<reco::GenParticle> > genrefHzzMMMM = (*GenParticlesMatchHMMMM)[Ref]; 
	if (!genrefHzzMMMM.isNull()){
	  cout << "Gen Z with pT= " << genrefHzzMMMM->p4().pt() << " and mass="<< genrefHzzMMMM->p4().mass() << endl;	  
	  RECOHzzMMMM_MatchingMCTruth[kk]= true;
	  RECOHzzMMMM_MatchingMCpT[kk]= genrefHzzMMMM->p4().pt();
	  RECOHzzMMMM_MatchingMCmass[kk]= genrefHzzMMMM->p4().mass();
	  RECOHzzMMMM_MatchingMCEta[kk]= genrefHzzMMMM->p4().eta();
	  RECOHzzMMMM_MatchingMCPhi[kk]= genrefHzzMMMM->p4().phi();
	}
      }

      kk++;
    }

    
    // EEEE
    i=1;
    leptonscands_EEEE->clear();
    edm::Handle<edm::View<Candidate> > CandidatesEEEE;
    iEvent.getByToken(RECOcollNameEEEE, CandidatesEEEE);

    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchHEEEE;
    iEvent.getByToken(goodHiggsTozzToEEEEMCMatch_, GenParticlesMatchHEEEE);
    if (GenParticlesMatchHEEEE.isValid()){
      cout << "The matched map collection has size= " <<  GenParticlesMatchHEEEE->size() << endl;
    }

    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEEE->begin();cand != CandidatesEEEE->end(); ++ cand ) {
      if (kk>99) break;
      RECO_EEEE_MASS[i-1][kk]=cand->p4().mass();
      RECO_EEEE_PT[i-1][kk]=cand->p4().pt();
      RECO_EEEE_ETA[i-1][kk]=cand->p4().eta();
      RECO_EEEE_PHI[i-1][kk]=cand->p4().phi();
      int l=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_EEEE_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_EEEE_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_EEEE_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_EEEE_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_EEEE_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_EEEE_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_EEEE_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_EEEE_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  leptonscands_EEEE->push_back( cand->daughter(j)->daughter(k)->clone());
	  //cout << "index" << i+j+k+l+2 <<endl;
	}
	l++;
      }

      // Matching HiggsToEEEE
      if (fillMCTruth==true){
	edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesEEEE,kk);
	edm::Ref<std::vector<reco::GenParticle> > genrefHzzEEEE = (*GenParticlesMatchHEEEE)[Ref]; 
	if (!genrefHzzEEEE.isNull()){
	  cout << "Gen Z with pT= " << genrefHzzEEEE->p4().pt() << " and mass="<< genrefHzzEEEE->p4().mass() << endl;	  
	  RECOHzzEEEE_MatchingMCTruth[kk]= true;
	  RECOHzzEEEE_MatchingMCpT[kk]= genrefHzzEEEE->p4().pt();
	  RECOHzzEEEE_MatchingMCmass[kk]= genrefHzzEEEE->p4().mass();
	  RECOHzzEEEE_MatchingMCEta[kk]= genrefHzzEEEE->p4().eta();
	  RECOHzzEEEE_MatchingMCPhi[kk]= genrefHzzEEEE->p4().phi();
	}
      }


      kk++;
    }


    // EEMM
    i=1;
    leptonscands_EEMM->clear();
    edm::Handle<edm::View<Candidate> > CandidatesEEMM;
    iEvent.getByToken(RECOcollNameEEMM, CandidatesEEMM);

    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchHEEMM;
    iEvent.getByToken(goodHiggsTozzToEEMMMCMatch_, GenParticlesMatchHEEMM);
    if (GenParticlesMatchHEEMM.isValid()){
      cout << "The matched map collection has size= " <<  GenParticlesMatchHEEMM->size() << endl;
    }

    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEMM->begin();cand != CandidatesEEMM->end(); ++ cand ) {
      if (kk>99) break;
      RECO_EEMM_MASS[i-1][kk]=cand->p4().mass();
      RECO_EEMM_PT[i-1][kk]=cand->p4().pt();
      RECO_EEMM_ETA[i-1][kk]=cand->p4().eta();
      RECO_EEMM_PHI[i-1][kk]=cand->p4().phi();
      int l=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_EEMM_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_EEMM_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_EEMM_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_EEMM_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_EEMM_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_EEMM_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_EEMM_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_EEMM_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();

	  leptonscands_EEMM->push_back( cand->daughter(j)->daughter(k)->clone());
	  //cout << "index" << i+j+k+l+2 <<endl;
	}
	l++;
      }

      // Matching HiggsToEEMM
      if (fillMCTruth==true){
	edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesEEMM,kk);
	edm::Ref<std::vector<reco::GenParticle> > genrefHzzEEMM = (*GenParticlesMatchHEEMM)[Ref]; 
	if (!genrefHzzEEMM.isNull()){
	  cout << "Gen Z with pT= " << genrefHzzEEMM->p4().pt() << " and mass="<< genrefHzzEEMM->p4().mass() << endl;	  
	  RECOHzzEEMM_MatchingMCTruth[kk]= true;
	  RECOHzzEEMM_MatchingMCpT[kk]= genrefHzzEEMM->p4().pt();
	  RECOHzzEEMM_MatchingMCmass[kk]= genrefHzzEEMM->p4().mass();
	  RECOHzzEEMM_MatchingMCEta[kk]= genrefHzzEEMM->p4().eta();
	  RECOHzzEEMM_MatchingMCPhi[kk]= genrefHzzEEMM->p4().phi();
	}
      }

      kk++;
    }
    
/*   
    // tri-leptons
    leptonscands_LLL0->clear();
    leptonscands_LLL1->clear();
    leptonscands_LLL2->clear();
    leptonscands_LLL3->clear();

    for (unsigned int i=0; i<RECOcollNameLLL.size(); i++) {  
      edm::Handle<edm::View<Candidate> > CandidatesLLL;
      if(i==0) iEvent.getByToken(triLeptonsMuMuMu, CandidatesLLL); 
      if(i==1) iEvent.getByToken(triLeptonsMuMuE, CandidatesLLL);
      if(i==2) iEvent.getByToken(triLeptonsMuEE, CandidatesLLL);
      if(i==3) iEvent.getByToken(triLeptonsEEE, CandidatesLLL);
      int k=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesLLL->begin();cand != CandidatesLLL->end(); ++ cand ) { 
	if (k>49) break;
	if (i==0) RECO_LLL0_MASS[k]=cand->p4().mass();
	if (i==0) RECO_LLL0_PT[0][k]=cand->p4().pt();
	if (i==1) RECO_LLL1_MASS[k]=cand->p4().mass();
	if (i==1) RECO_LLL1_PT[0][k]=cand->p4().pt();
	if (i==2) RECO_LLL2_MASS[k]=cand->p4().mass();
	if (i==2) RECO_LLL2_PT[0][k]=cand->p4().pt();
	if (i==3) RECO_LLL3_MASS[k]=cand->p4().mass();
	if (i==4) RECO_LLL3_PT[0][k]=cand->p4().pt();
	cout << "Tri-lepton candidate of type=" << RECOcollNameLLL.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  if (i==0) {
	    RECO_LLL0_PT[j+1][k]=cand->daughter(j)->p4().pt();
	    leptonscands_LLL0->push_back( cand->daughter(j)->clone());
	  }
	  if (i==1) {
	    RECO_LLL1_PT[j+1][k]=cand->daughter(j)->p4().pt();
	    leptonscands_LLL1->push_back( cand->daughter(j)->clone());
	  }
	  if (i==2) {
	    RECO_LLL2_PT[j+1][k]=cand->daughter(j)->p4().pt();
	    leptonscands_LLL2->push_back( cand->daughter(j)->clone()); 
	  }
	  if (i==3) {
	    RECO_LLL3_PT[j+1][k]=cand->daughter(j)->p4().pt();
	    leptonscands_LLL3->push_back( cand->daughter(j)->clone());
	  }
	}
	k++;
      }
    }
*/    
/*
    // 4-leptons SS
    leptonscands_LLLLss0->clear();
    leptonscands_LLLLss1->clear();
    leptonscands_LLLLss2->clear();
    
    for (unsigned int i=0; i<RECOcollNameLLLLss.size(); i++) {  
      edm::Handle<edm::View<Candidate> > CandidatesLLLLss;
      if(i==0) iEvent.getByToken(quadLeptons4Mu, CandidatesLLLLss);
      if(i==1) iEvent.getByToken(quadLeptons2Mu2E, CandidatesLLLLss);
      if(i==2) iEvent.getByToken(quadLeptons4E, CandidatesLLLLss);
      int kk=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesLLLLss->begin();cand != CandidatesLLLLss->end(); ++ cand ) { 
	if (kk>49) break;
	if (i==0) RECO_LLLL0ss_MASS[kk]=cand->p4().mass();
	if (i==0) RECO_LLLL0ss_PT[0][kk]=cand->p4().pt();
	if (i==1) RECO_LLLL1ss_MASS[kk]=cand->p4().mass();
	if (i==1) RECO_LLLL1ss_PT[0][kk]=cand->p4().pt();
	if (i==2) RECO_LLLL2ss_MASS[kk]=cand->p4().mass();
	if (i==2) RECO_LLLL2ss_PT[0][kk]=cand->p4().pt();
	cout << "4-lepton ss candidate of type=" << RECOcollNameLLLLss.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	    if (i==0) {
	      if (j==0) RECO_LLLL0ss_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      if (j==1) RECO_LLLL0ss_PT[j+k+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLLss0->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	    if (i==1) {
	      if (j==0) RECO_LLLL1ss_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      if (j==1) RECO_LLLL1ss_PT[j+k+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLLss1->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	    if (i==2) {
	      if (j==0) RECO_LLLL2ss_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      if (j==1) RECO_LLLL2ss_PT[j+k+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLLss2->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	  }
	}
	kk++;
      }      
    }
*/
/*
    // 4-leptons 3l+l
    leptonscands_LLLl0->clear();
    leptonscands_LLLl1->clear(); 
    
    for (unsigned int i=0; i<RECOcollNameLLLl.size(); i++){ 
      edm::Handle<edm::View<Candidate> > CandidatesLLLl; 
      if(i==0) iEvent.getByToken(quadLeptons3Mu1E, CandidatesLLLl);
      if(i==1) iEvent.getByToken(quadLeptons3E1Mu, CandidatesLLLl);  
      int kk=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesLLLl->begin();cand != CandidatesLLLl->end(); ++ cand ) { 
	if (kk>49) continue;
	if (i==0) RECO_LLLl0_MASS[kk]=cand->p4().mass();
	if (i==0) RECO_LLLl0_PT[0][kk]=cand->p4().pt();
	if (i==1) RECO_LLLl1_MASS[kk]=cand->p4().mass();
	if (i==1) RECO_LLLl1_PT[0][kk]=cand->p4().pt();
	cout << "4-lepton from 3l+l candidate of type=" << RECOcollNameLLLl.at(i).label()
	     << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  //cout << "j= " << j << endl;
	  if ( i==0 && j==1) {
	    RECO_LLLl0_PT[j+3][kk]=cand->daughter(j)->p4().pt();
	    leptonscands_LLLl0->push_back( cand->daughter(j)->clone());
	  }
	  if ( i==1 && j==1) {
	    RECO_LLLl1_PT[j+3][kk]=cand->daughter(j)->p4().pt(); 
	    leptonscands_LLLl1->push_back( cand->daughter(j)->clone());
	  }
	  for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) { 
	    //cout << "k= " << k << endl;
	    if (i==0) {
	      if (j==0) RECO_LLLl0_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLl0->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	    if (i==1) {
	      if (j==0) RECO_LLLl1_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLl1->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	    
	  } 
	} 
	kk++;
      }
    }


    // LLLL merged (no flavour, no charge)

    i=1;
    leptonscands_LLLL->clear();
    edm::Handle<edm::View<Candidate> > CandidatesLLLL;
    iEvent.getByToken(RECOcollNameLLLL, CandidatesLLLL);
    kk=0;    
    for( edm::View<Candidate>::const_iterator cand = CandidatesLLLL->begin();cand != CandidatesLLLL->end(); ++ cand ) {
      if (kk>49) break;
      cout << "4lepton (any flavour and charge) candidate of type= allLLLL"  << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
      
      RECO_LLLL_MASS[0][kk]=cand->p4().mass();
      RECO_LLLL_PT[0][kk]=cand->p4().pt();
      RECO_LLLL_ETA[0][kk]=cand->p4().eta();
      RECO_LLLL_PHI[0][kk]=cand->p4().phi();
      int l=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_LLLL_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_LLLL_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_LLLL_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_LLLL_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_LLLL_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_LLLL_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_LLLL_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_LLLL_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  
	  leptonscands_LLLL->push_back( cand->daughter(j)->daughter(k)->clone());
	  //cout << "index" << i+j+k+l+2 <<endl;
	}
	l++;
      }
      kk++;
    }
*/    
  }


  void SetMCValues(const reco::Candidate& cand, int nMC){
    MC_E[nMC]     = cand.p4().energy();
    MC_PT[nMC]    = cand.p4().pt();
    MC_ETA[nMC]   = cand.p4().eta();
    MC_THETA[nMC] = cand.p4().theta();
    MC_PHI[nMC]   = cand.p4().phi();
    MC_MASS[nMC]  = cand.p4().mass();
    MC_PDGID[nMC] = cand.pdgId();
  }

 
 
  
  bool match(double mass, double pt, int charge, const reco::CandidateCollection *c1Coll){
    
    bool found=false;
    for( reco::CandidateCollection::const_iterator pp = c1Coll->begin();pp != c1Coll->end(); ++ pp ) {
      if ((fabs(pp->p4().mass()-mass)  <0.001 ) &&
	  (fabs(pp->p4().pt()  -pt)    <0.001 ) &&
	  (fabs(pp->charge()   -charge)<0.001 )  ){
	found=true;
	//std::cout << "Found lepton in the higgs-like leptons collection" << std::endl;
      }
    }
    return found;
  }
  
  bool matchParticle(double mass, double pt, int charge, const reco::Candidate *c1){
    
    bool found=false;
    if ((fabs(c1->p4().mass()-mass)  <0.001 ) &&
	(fabs(c1->p4().pt()  -pt)    <0.001 ) &&
	(fabs(c1->charge()   -charge)<0.001 )  ){
      found=true;
      //std::cout << "Found lepton in the higgs-like leptons collection" << std::endl;
    }    
    return found;
  }
  
  
  void fillElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup){
   
 
    // Supercluster collection
    edm::Handle<reco::SuperClusterCollection> clusters;
    iEvent.getByLabel(clusterCollectionTag_,clusters);
    

    // EG isolation
    //edm::Handle<reco::GsfElectronRefVector> EleRefs;
    //iEvent.getByLabel(electronEgmTag_, EleRefs);
    edm::Handle<edm::View<pat::Electron> > EleRefs;
    iEvent.getByToken(electronEgmTag_, EleRefs);



    //Electron ID MVA Trig and Non Trig
//    edm::Handle<edm::ValueMap<float> >  mvaTrigV0map;
//    iEvent.getByToken(mvaTrigV0MapTag_, mvaTrigV0map);
    edm::Handle<edm::ValueMap<float> >  mvaNonTrigV0map;
    iEvent.getByToken(mvaNonTrigV0MapTag_, mvaNonTrigV0map);
    
    // Particle Flow Isolation
/*
    edm::Handle<edm::ValueMap<double> > isoPFChargedAllelemap;
    iEvent.getByToken(electronPFIsoValueChargedAllTag_, isoPFChargedAllelemap);
    
    edm::Handle<edm::ValueMap<double> > isoPFChargedelemap;
    iEvent.getByToken(electronPFIsoValueChargedTag_, isoPFChargedelemap);
    
    edm::Handle<edm::ValueMap<double> > isoPFNeutralelemap;
    iEvent.getByToken(electronPFIsoValueNeutralTag_, isoPFNeutralelemap);
    
    edm::Handle<edm::ValueMap<double> > isoPFGammaelemap;
    iEvent.getByToken(electronPFIsoValueGammaTag_, isoPFGammaelemap);
    
    edm::Handle<edm::ValueMap<double> > isoPFPUelemap;
    iEvent.getByToken(electronPFIsoValuePUTag_, isoPFPUelemap);
*/
    // electron regression
    //edm::Handle<edm::ValueMap<double> > eleRegressionEnergymap;
    //iEvent.getByLabel(eleRegressionEnergyTag_, eleRegressionEnergymap);

    //edm::Handle<edm::ValueMap<double> > eleRegressionEnergyErrormap;
    //iEvent.getByLabel(eleRegressionEnergyErrorTag_, eleRegressionEnergyErrormap);
    cout << "test 3" << endl;
    // Vertexing
    //3D DA
    edm::Handle<edm::View<pat::Electron> > VertEleCandidates;
    iEvent.getByLabel(electronTag_Vert, VertEleCandidates);
    
    edm::Handle<edm::ValueMap<float> > vertexelemap;
    iEvent.getByToken(electronMapTag_Vert, vertexelemap);
    
    edm::Handle<edm::ValueMap<float> > vertexelemapvalue;
    iEvent.getByToken(electronMapTag_VertValue, vertexelemapvalue);
    
    edm::Handle<edm::ValueMap<float> > vertexelemaperror;
    iEvent.getByToken(electronMapTag_VertError, vertexelemaperror);

    // 3D KF
    edm::Handle<edm::ValueMap<float> > vertexelemapKF;
    iEvent.getByToken(electronMapTag_VertKF, vertexelemapKF);
    
    edm::Handle<edm::ValueMap<float> > vertexelemapvalueKF;
    iEvent.getByToken(electronMapTag_VertValueKF, vertexelemapvalueKF);
    
    edm::Handle<edm::ValueMap<float> > vertexelemaperrorKF;
    iEvent.getByToken(electronMapTag_VertErrorKF, vertexelemaperrorKF);

    cout << "test 4" << endl;
    // STIP SLIP
    edm::Handle<edm::ValueMap<float> > stipelemap;
    iEvent.getByToken(electronSTIPMapTag_Vert, stipelemap);
    
    edm::Handle<edm::ValueMap<float> > slipelemap;
    iEvent.getByToken(electronSLIPMapTag_Vert, slipelemap);
    
    edm::Handle<edm::ValueMap<float> > stipelemapvalue;
    iEvent.getByToken(electronSTIPMapTag_VertValue, stipelemapvalue);
    
    edm::Handle<edm::ValueMap<float> > slipelemapvalue;
    iEvent.getByToken(electronSLIPMapTag_VertValue, slipelemapvalue);
    
    edm::Handle<edm::ValueMap<float> > stipelemaperror;
    iEvent.getByToken(electronSTIPMapTag_VertError, stipelemaperror);
    
    edm::Handle<edm::ValueMap<float> > slipelemaperror;
    iEvent.getByToken(electronSLIPMapTag_VertError, slipelemaperror);		
   
    cout << "test 5" << endl; 
    // Conversion Finder
//    edm::Handle<edm::ValueMap<float> > conversiondistmap;
//    iEvent.getByToken(ConvMapDistTag_, conversiondistmap);

//    edm::Handle<edm::ValueMap<float> > conversiondcotmap;
//    iEvent.getByToken(ConvMapDcotTag_, conversiondcotmap);


    // primary vertex
   

 
    Vertex primVertex;
    math::XYZPoint pVertex(0., 0., 0.);
    bool pvfound = (PV.size() != 0);
    cout << "pvfound=" << pvfound << endl;
    if(pvfound){       
      for(reco::VertexCollection::const_iterator it=PV.begin() ; it!=PV.end() ; ++it){
	if(!it->isFake() && it->ndof() > 4 && fabs(it->position().z()) <= 24 && fabs(it->position().rho()) <= 2){ 	  	  
	  pVertex = math::XYZPoint(it->position().x(), it->position().y(), it->position().z());
	  cout << "P vertex position used for computing dxy and dz for electron and muons is (x,y,z)= " << 
	    pVertex.x() << " " <<
	    pVertex.y() << " " <<
	    pVertex.z() << endl;
	  break;
	}
      }
    } 
    
    int index=0;
    RECO_NELE=EleRefs->size();
   
/*
    // Matching
    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchEle;
    iEvent.getByToken(goodElectronMCMatch_, GenParticlesMatchEle);
    edm::Handle<reco::CandidateCollection > CollEle;
    iEvent.getByToken(myElectrons_, CollEle);

    if (GenParticlesMatchEle.isValid()){
      cout << endl<< "Electrons:"<<endl<<"The reco collection to be matched has size= " <<  CollEle->size() << endl;
      cout << "The matched map collection has size= " <<  GenParticlesMatchEle->size() << endl;
    }*/
    //



    for (edm::View<pat::Electron>::const_iterator cand = EleRefs->begin(); 
	 cand != EleRefs->end(); ++cand) {
      
      if(index>99) break;

      edm::Ref<edm::View<pat::Electron> > eletrackref(EleRefs,index);
      edm::Ref<edm::View<pat::Electron> > eletrackrefv(VertEleCandidates,index);
  
     
      // Global variables
      RECOELE_isEcalDriven[index]    = cand->ecalDrivenSeed();
      RECOELE_isTrackerDriven[index] = cand->trackerDrivenSeed();

      std::cout << "\n Electron in the event: "
		<< "  isEcalDriven="    << RECOELE_isEcalDriven[index]  
		<< "  isTrackerDriven=" << RECOELE_isTrackerDriven[index] 
		<< std::endl;
      
      // kinematic
      RECOELE_E[index]       = cand->p4().energy();
      RECOELE_PT[index]      = cand->p4().pt();


      RECOELE_P[index]=sqrt(cand->p4().px()*cand->p4().px()+cand->p4().py()*cand->p4().py()+cand->p4().pz()*cand->p4().pz());
      RECOELE_ETA[index]     = cand->p4().eta();
      RECOELE_THETA[index]   = cand->p4().theta();
      RECOELE_PHI[index]     = cand->p4().phi();
      RECOELE_MASS[index]    = cand->p4().mass();
      RECOELE_CHARGE[index]  = cand->charge();


      std::cout << "--kinematic:" 		
		<< "  pT="     << RECOELE_PT[index]  
		<< "  E="      << RECOELE_E[index]  
		<< "  p="      << RECOELE_P[index] 
		<< "  eta="    << RECOELE_ETA[index] 
		<< "  theta="  << RECOELE_THETA[index] 
		<< "  phi="    << RECOELE_PHI[index] 
		<< "  mass="   << RECOELE_MASS[index] 
		<< "  charge=" << RECOELE_CHARGE[index] 
		<< std::endl;

      
      // SuperCluster
      reco::SuperClusterRef sclRef = cand->superCluster();
      math::XYZPoint sclPos = cand->superClusterPosition();

      //if(!cand->ecalDrivenSeed() && cand->trackerDrivenSeed())
      //clRef = cand->pflowSuperCluster();
      double R  = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
      double Rt = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
      ele_sclRawE[index]   = sclRef->rawEnergy() ;
      RECOELE_scl_E[index]   = sclRef->energy() ;
      RECOELE_scl_Et[index]  = sclRef->energy()*(Rt/R) ;
      RECOELE_scl_Eta[index] = sclRef->eta() ;
      RECOELE_scl_Phi[index] = sclRef->phi() ;
      ele_sclX[index]  = sclPos.X();
      ele_sclY[index] =  sclPos.Y();
      ele_sclZ[index] =  sclPos.Z();
			
      std::cout << "--supercluster properties:"
	<< "  E_sc="    << RECOELE_scl_E[index]
	<< "  Et_sc="   << RECOELE_scl_Et[index]
	<< "  Eeta_sc=" << RECOELE_scl_Eta[index]
	<< "  phi_sc="  << RECOELE_scl_Phi[index]
	<< std::endl;

      //Covariance matrix
      TMatrixDSym bigCov;
      double dp = 0.;
      if (cand->ecalDriven()) {
	dp = cand->p4Error(pat::Electron::P4_COMBINATION);	
	RECOELE_PTError[index]=dp;
      }
      
      else {
	// Parametrization from Claude Charlot, 
	// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/CJLST/ZZAnalysis/AnalysisStep/src/ZZMassErrors.cc?revision=1.2&view=markup
#if CMSSW_VERSION<500
	double ecalEnergy = cand->ecalEnergy() ;
#else
	double ecalEnergy = cand->correctedEcalEnergy() ;
#endif
	double err2 = 0.0;
	if (cand->isEB()) {
	  err2 += (5.24e-02*5.24e-02)/ecalEnergy;  
	  err2 += (2.01e-01*2.01e-01)/(ecalEnergy*ecalEnergy);
	  err2 += 1.00e-02*1.00e-02;
	} else if (cand->isEE()) {
	  err2 += (1.46e-01*1.46e-01)/ecalEnergy;  
	  err2 += (9.21e-01*9.21e-01)/(ecalEnergy*ecalEnergy);
	  err2 += 1.94e-03*1.94e-03;
	}
	dp = ecalEnergy * sqrt(err2);
	RECOELE_PTError[index]=dp;
      }
      // In order to produce a 3x3 matrix, we need a jacobian from (p) to (px,py,pz), i.e.
      //            [ Px/P  ]
      //  C_(3x3) = [ Py/P  ] * sigma^2(P) * [ Px/P Py/P Pz/P  ]
      //            [ Pz/P  ]
      AlgebraicMatrix31 ptop3;
      ptop3(0,0) = cand->px()/cand->p();
      ptop3(1,0) = cand->py()/cand->p();
      ptop3(2,0) = cand->pz()/cand->p();
      AlgebraicSymMatrix33 mat = ROOT::Math::Similarity(ptop3, AlgebraicSymMatrix11(dp*dp) );
      for (int i = 0; i < 3; ++i) { 
	for (int j = 0; j < 3; ++j) {
	  RECOELE_COV[index][i][j]=mat[i][j];
	}
      }
      



      // Egamma isolation
      //RECOELE_EGMTRACKISO[index]=(*egmisoTkelemap)[eletrackref]/eletrackref->pt();
      // RECOELE_EGMECALISO[index]=(*egmisoEcalelemap)[eletrackref]/eletrackref->pt();
      // RECOELE_EGMHCALISO[index]=(*egmisoHcalelemap)[eletrackref]/eletrackref->pt();
      // RECOELE_EGMX[index]=( (*egmisoTkelemap)[eletrackref] + (*egmisoEcalelemap)[eletrackref] + (*egmisoHcalelemap)[eletrackref])/eletrackref->pt();

      // temporary solution for the problem of rechits
      RECOELE_EGMECALISO[index]=(eletrackref->dr03EcalRecHitSumEt())/eletrackref->pt();
      RECOELE_EGMHCALISO[index]=(eletrackref->dr03HcalTowerSumEt())/eletrackref->pt();
      RECOELE_EGMX[index]=RECOELE_EGMTRACKISO[index]+RECOELE_EGMECALISO[index]+RECOELE_EGMHCALISO[index];
      
      std::cout << "--EG isolation: electron"
		<< "  X="     << RECOELE_EGMX[index]
		<< "  Track=" << RECOELE_EGMTRACKISO[index]
		<< "  Ecal= " << RECOELE_EGMECALISO[index]
		<< "  Hcal=" << RECOELE_EGMHCALISO[index]
		<< std::endl;


      // PF isolation  (R unknown, TODO)
      RECOELE_PFchHad[index]  = (cand->pfIsolationVariables().sumChargedHadronPt);
      RECOELE_PFneuHad[index] = (cand->pfIsolationVariables().sumNeutralHadronEt);
      RECOELE_PFphoton[index] = (cand->pfIsolationVariables().sumPhotonEt);
      RECOELE_PFPUchAllPart[index]= (cand->pfIsolationVariables().sumPUPt);
      RECOELE_PFchAllPart[index]= (cand->pfIsolationVariables().sumChargedParticlePt);

      //RECOELE_PFchAllPart[index]      =  (*isoPFChargedAllelemap)[eletrackref];
      //RECOELE_PFchHad[index]          =  (*isoPFChargedelemap)[eletrackref];
      //RECOELE_PFneuHad[index]         =  (*isoPFNeutralelemap)[eletrackref];
      //RECOELE_PFphoton[index]         =  (*isoPFGammaelemap)[eletrackref];
     // RECOELE_PFPUchAllPart[index]    =  (*isoPFPUelemap)[eletrackref];

      RECOELE_PFX_dB[index]           =  (RECOELE_PFchHad[index] + max(RECOELE_PFphoton[index]+RECOELE_PFneuHad[index]-0.5*RECOELE_PFPUchAllPart[index],0.))/eletrackref->pt();

      float EffectiveArea=0.;
      if (use2011EA){
	if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.18;
	if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.20;
	if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.15;
	if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.19;
	if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.21;
	if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.22;
	if (fabs(RECOELE_scl_Eta[index]) >= 2.4 )                                       EffectiveArea = 0.29;
      }
      else { 
	// 7_4_X use eta 
	//if (fabs(RECOELE_ETA[index]) >= 0.0   && fabs(RECOELE_ETA[index]) < 0.8 )   EffectiveArea = 0.1830;
	//if (fabs(RECOELE_ETA[index]) >= 0.8   && fabs(RECOELE_ETA[index]) < 1.3 )   EffectiveArea = 0.1734;
	//if (fabs(RECOELE_ETA[index]) >= 1.3   && fabs(RECOELE_ETA[index]) < 2.0 )   EffectiveArea = 0.1077;
	//if (fabs(RECOELE_ETA[index]) >= 2.0   && fabs(RECOELE_ETA[index]) < 2.2 )   EffectiveArea = 0.1565;
	//if (fabs(RECOELE_ETA[index]) >= 2.2 )                                       EffectiveArea = 0.2680;

	// 7_6_X use eta supercluster   
//	if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.1752;
//        if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.1862;
//        if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.1411;
//        if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.1534;
//        if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.1903;
//        if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.2243;
//        if (fabs(RECOELE_scl_Eta[index]) >= 2.4   && fabs(RECOELE_scl_Eta[index]) < 5.0  )  EffectiveArea = 0.2687;
        // 8_0_X  
        if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.1752;
        if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.1862;
        if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.1411;
        if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.1534;
        if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.1903;
        if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.2243;
        if (fabs(RECOELE_scl_Eta[index]) >= 2.4   && fabs(RECOELE_scl_Eta[index]) < 5.0  )  EffectiveArea = 0.2687;
      }


      //RECOELE_PFX_rho[index]=(RECOELE_PFchHad[index]+
      //		      max( (RECOELE_PFneuHad[index]+RECOELE_PFphoton[index]-max(RHO_ele,0.0)*EffectiveArea),0.0) )/cand->p4().pt(); 

      RECOELE_PFX_rho[index]= (cand->pfIsolationVariables().sumChargedHadronPt+
			       std::max(
					cand->pfIsolationVariables().sumPhotonEt+
					cand->pfIsolationVariables().sumNeutralHadronEt-
					max(RHO_ele,0.0)*EffectiveArea,
					0.0)
			       )/cand->p4().pt();

      RECOELE_PFchHad[index]  = (cand->pfIsolationVariables().sumChargedHadronPt);
      RECOELE_PFneuHad[index] = (cand->pfIsolationVariables().sumNeutralHadronEt); 
      RECOELE_PFphoton[index] = (cand->pfIsolationVariables().sumPhotonEt); 
      //RECOELE_PFsumPUPt[index]= (cand->pfIsolationVariables().sumPUPt); 
      RECOELE_PFPUchAllPart[index]=(cand->pfIsolationVariables().sumPUPt);   

      /*       RECOELE_PFX_dB[index]=(RECOELE_PFchHad[index]+max(0.,RECOELE_PFneuHad[index]+RECOELE_PFphoton[index]-0.5*RECOELE_PFsumPUPt[index]))/cand->p4().pt(); */      
      /*       float EffectiveArea=0.; */
      /*       if (fabs(RECOELE_ETA[index]) >= 0.0 && fabs(RECOELE_ETA[index]) < 1.0 ) EffectiveArea = 0.132; */
      /*       if (fabs(RECOELE_ETA[index]) >= 1.0 && fabs(RECOELE_ETA[index]) < 1.5 ) EffectiveArea = 0.120; */
      /*       if (fabs(RECOELE_ETA[index]) >= 1.5 && fabs(RECOELE_ETA[index]) < 2.0 ) EffectiveArea = 0.114; */
      /*       if (fabs(RECOELE_ETA[index]) >= 2.0 && fabs(RECOELE_ETA[index]) < 2.2 ) EffectiveArea = 0.139; */
      /*       if (fabs(RECOELE_ETA[index]) >= 2.2 && fabs(RECOELE_ETA[index]) < 2.3 ) EffectiveArea = 0.168; */
      /*       if (fabs(RECOELE_ETA[index]) >= 2.3 )                                     EffectiveArea = 0.189; */
      
      /*       RECOELE_PFX_rho[index]=(RECOELE_PFchHad[index]+  */
      /* 				max( (RECOELE_PFneuHad[index]+RECOELE_PFphoton[index]-max(RHO_ele,0.0)*EffectiveArea),0.0) )/cand->p4().pt(); */
      
      std::cout << "--PF isolation"
	        << "  chAllPart="<< RECOELE_PFchAllPart[index] 
                << "  chHad="    << RECOELE_PFchHad[index] 
                << "  neuHad="   << RECOELE_PFneuHad[index] 
                << "  photon="   << RECOELE_PFphoton[index]
                << "  PUchAllPart="  << RECOELE_PFPUchAllPart[index]
                << "  PFX_dB="  << RECOELE_PFX_dB[index] 
                << "  PFX_rho=" << RECOELE_PFX_rho[index]
                << std::endl;
      
      // Electron Regression
      // RECOELE_regEnergy[index]=(*eleRegressionEnergymap)[eletrackrefv];
      // RECOELE_regEnergyError[index]=(*eleRegressionEnergyErrormap)[eletrackrefv];

/*qier
      if (!cand->parentSuperCluster().isNull()){
        RECOELE_TLE_ParentSC_X[index]=cand->parentSuperCluster()->x();
        RECOELE_TLE_ParentSC_Y[index]=cand->parentSuperCluster()->y();
        RECOELE_TLE_ParentSC_Z[index]=cand->parentSuperCluster()->z();
      }

      cout << "Parent SC x,y,x= " <<RECOELE_TLE_ParentSC_X[index] << " " << RECOELE_TLE_ParentSC_Y[index] << " " << RECOELE_TLE_ParentSC_Z[index] << endl;
*/
      // Vertexing DA
      RECOELE_SIP[index]=(*vertexelemap)[eletrackrefv];
      RECOELE_IP[index]=(*vertexelemapvalue)[eletrackrefv];
      RECOELE_IPERROR[index]=(*vertexelemaperror)[eletrackrefv];

      // KF
//      RECOELE_SIP_KF[index]=(*vertexelemapKF)[eletrackrefv];
//      RECOELE_IP_KF[index]=(*vertexelemapvalueKF)[eletrackrefv];
//      RECOELE_IPERROR_KF[index]=(*vertexelemaperrorKF)[eletrackrefv];


      //RECOELE_SIP_GD[index]=(*vertexelemapGD)[eletrackrefv];
      //if (decaychannel=="4e" || decaychannel=="2e2mu" ) RECOELE_SIP_GDEEEE[index]=(*vertexelemapGDEEEE)[eletrackrefv];
      //RECOELE_SIP_Std[index]=(*vertexelemapStd)[eletrackrefv]; 
      //if (decaychannel=="4e" || decaychannel=="2e2mu" ) RECOELE_SIP_StdEEEE[index]=(*vertexelemapStdEEEE)[eletrackrefv]; 
      //RECOELE_SIP_Kin[index]=(*vertexelemapKin)[eletrackrefv]; 
      //if (decaychannel=="4e" || decaychannel=="2e2mu" ) RECOELE_SIP_KinEEEE[index]=(*vertexelemapKinEEEE)[eletrackrefv]; 
      
/*
      RECOELE_STIP[index]=(*stipelemap)[eletrackrefv];
      RECOELE_SLIP[index]=(*slipelemap)[eletrackrefv];
      RECOELE_TIP[index]=(*stipelemapvalue)[eletrackrefv] ;
      RECOELE_LIP[index]=(*slipelemapvalue)[eletrackrefv];
      RECOELE_TIPERROR[index]=(*stipelemaperror)[eletrackrefv] ;
      RECOELE_LIPERROR[index]=(*slipelemaperror)[eletrackrefv];
*/            
      std::cout << "--vertexing: electron" 
		<< "  Sign_3DIP=" << RECOELE_SIP[index]
		<< "  3DIP="      << RECOELE_IP[index]
		<< "  3DIPerr="   << RECOELE_IPERROR[index] 	     
		<< "  Sign_TIP="  << RECOELE_STIP[index]  
		<< "  TIP="       << RECOELE_TIP[index] 
		<< "  TIPerror="  << RECOELE_TIPERROR[index] 
		<< "  Sign_LIP="  << RECOELE_SLIP[index]
		<< "  LIP="       << RECOELE_LIP[index] 
		<< "  LIPerror="  << RECOELE_LIPERROR[index]
		<< std::endl;


      // GsfTrack
      RECOELE_gsftrack_NPixHits[index]  = cand->gsfTrack()->hitPattern().numberOfValidPixelHits();
      RECOELE_gsftrack_NStripHits[index]= cand->gsfTrack()->hitPattern().numberOfValidStripHits();
      RECOELE_gsftrack_chi2[index]      = cand->gsfTrack()->normalizedChi2();
      RECOELE_gsftrack_dxyB[index]      = cand->gsfTrack()->dxy(bs.position()) ;
      RECOELE_gsftrack_dxy[index]       = cand->gsfTrack()->dxy(pVertex);
      RECOELE_gsftrack_dxyError[index]  = cand->gsfTrack()->dxyError();
      RECOELE_gsftrack_dzB[index]       = cand->gsfTrack()->dz(bs.position());
      RECOELE_gsftrack_dz[index]        = cand->gsfTrack()->dz(pVertex);
      RECOELE_gsftrack_dzError[index]   = cand->gsfTrack()->dzError();

      //Conversion variables
      RECOELE_gsftrack_losthits[index]  = cand->gsfTrack()->numberOfLostHits();
      RECOELE_gsftrack_validhits[index] = cand->gsfTrack()->numberOfValidHits();
      RECOELE_gsftrack_expected_inner_hits[index] = cand->gsfTrack()->hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS);

      std::cout << "--gfstrack properties: " 
	<< "  nPixhits=" << RECOELE_gsftrack_NPixHits[index]
	<< "  nStriphits=" << RECOELE_gsftrack_NStripHits[index]
	<< "  chi2="     << RECOELE_gsftrack_chi2[index] 
	<< "  dxy="      << RECOELE_gsftrack_dxy[index] 
	<< "  dxyB="     << RECOELE_gsftrack_dxyB[index] 
	<< "  dxyError=" << RECOELE_gsftrack_dxyError[index] 
	<< "  dz="       << RECOELE_gsftrack_dz[index]
	<< "  dzB="      << RECOELE_gsftrack_dzB[index]
        << "  dzError="  << RECOELE_gsftrack_dzError[index] 
	<< "  losthits=" << RECOELE_gsftrack_losthits[index] 
	<< "  validhits="<< RECOELE_gsftrack_validhits[index] 
        << "  innerhits="<< RECOELE_gsftrack_expected_inner_hits[index] 
	<< std::endl;

      // Track-Cluster matching attributes
      RECOELE_ep[index]             = cand->eSuperClusterOverP(); 
      RECOELE_eSeedp[index]         = cand->eSeedClusterOverP();     
      RECOELE_eSeedpout[index]      = cand->eSeedClusterOverPout(); 
      RECOELE_eElepout[index]       = cand->eEleClusterOverPout();        
      
      RECOELE_deltaEtaIn [index]  = cand->deltaEtaSuperClusterTrackAtVtx(); 
      RECOELE_deltaEtaSeed[index] = cand->deltaEtaSeedClusterTrackAtCalo(); 
      RECOELE_deltaEtaEle[index]  = cand->deltaEtaEleClusterTrackAtCalo();  
      RECOELE_deltaPhiIn[index]   = cand->deltaPhiSuperClusterTrackAtVtx();
      RECOELE_deltaPhiSeed[index] = cand->deltaPhiSeedClusterTrackAtCalo(); 
      RECOELE_deltaPhiEle[index]  = cand->deltaPhiEleClusterTrackAtCalo() ;  

      std::cout << "--track-cluster matching: " 
	<< "  eSC_p="        << RECOELE_ep[index] 
	<< "  eSeed_p="      << RECOELE_eSeedp[index] 
	<< "  eSeed_pout="   << RECOELE_eSeedpout[index] 
	<< "  eEle_pout="    << RECOELE_eElepout[index] 
	<< "  deltaEtaIn="   << RECOELE_deltaEtaIn[index] 
        << "  deltaEtaSeed=" << RECOELE_deltaEtaSeed[index] 
        << "  deltaEtaEle="  << RECOELE_deltaEtaEle[index] 
        << "  deltaPhiIn="   << RECOELE_deltaPhiIn[index] 
        << "  deltaPhiSeed=" << RECOELE_deltaPhiSeed[index] 
	<< "  deltaPhiEle="  << RECOELE_deltaPhiEle[index] 
	<< std::endl;
      
       // Fiducial flags
      if (cand->isEB()) RECOELE_isbarrel[index] = 1 ; 
      else  RECOELE_isbarrel[index] = 0 ;
      if (cand->isEE()) RECOELE_isendcap[index] = 1 ; 
      else  RECOELE_isendcap[index] = 0 ;
      if (cand->isEBEtaGap())  RECOELE_isEBetaGap[index]  = 1 ;  
      if (cand->isEBPhiGap())  RECOELE_isEBphiGap[index]  = 1 ;  
      if (cand->isEEDeeGap())  RECOELE_isEEdeeGap[index]  = 1 ;  
      if (cand->isEERingGap()) RECOELE_isEEringGap[index] = 1 ;
      if (cand->isGap())       RECOELE_isGap[index] = 1 ;

      std::cout << "--fiducial flags: " 
	<< "  isEB="        << RECOELE_isbarrel[index] 
	<< "  isEBetagap="  << RECOELE_isEBetaGap[index] 
	<< "  isEBphigap="  << RECOELE_isEBphiGap[index] 
	<< "  isEE="        << RECOELE_isendcap[index] 
        << "  isEEdeegap="  << RECOELE_isEEdeeGap[index] 
	<< "  isEEringgap=" << RECOELE_isEEringGap[index] 
        << "  isGap="       << RECOELE_isGap[index]
	<< std::endl;
      

      // Shower shape
      RECOELE_sigmaIetaIeta[index] = cand->sigmaIetaIeta() ; 
      RECOELE_sigmaEtaEta[index]   = cand->sigmaEtaEta() ;
      RECOELE_e15[index]           = cand->e1x5() ;
      RECOELE_e25max[index]        = cand->e2x5Max() ;
      RECOELE_e55[index]           = cand->e5x5() ;
      RECOELE_he[index]            = cand->hadronicOverEm() ;
      // RECOELE_r9[index]            = cand->r9() ;

      std::cout << "--shower shape:" 
	<< "  sigmaIetaIeta=" << RECOELE_sigmaIetaIeta[index] 
	<< "  sigmaEtaEta=" << RECOELE_sigmaEtaEta[index]  
	<< "  e15=" << RECOELE_e15[index] 
	<< "  e25max=" << RECOELE_e25max[index] 
	<< "  e55=" << RECOELE_e55[index]  
	<< "  he=" << RECOELE_he[index] 
	<< "  r9=" << RECOELE_r9[index] 
	<< std::endl;
      
      // also variables on hcal
      // Isolation
      // Particle flow
      //RECOELE_mva[index] = cand->mva();
      std::cout << "--PF: mva = " << RECOELE_mva[index] << std::endl;

      // Brem & Classifaction
      RECOELE_fbrem[index]         = cand->fbrem() ;
      RECOELE_nbrems[index]        = cand->numberOfBrems();
      RECOELE_Class[index]         = cand->classification() ;
      if (useRECOformat) RECOELE_fbrem_mean[index]=1. - cand->gsfTrack()->outerMomentum().R()/cand->gsfTrack()->innerMomentum().R();
      RECOELE_fbrem_mode[index]=cand->fbrem();

      std::cout << "--brem & classification: fbrem/nbrems/Class/fbrem_mean/fbrem_mode = "
		<< "  fbrem="      << RECOELE_fbrem[index]
		<< "  nbrems="     << RECOELE_nbrems[index]
        	<< "  class="      << RECOELE_Class[index]
		<< "  fbrem_mean=" << RECOELE_fbrem_mean[index]
		<< "  fbrem_mode=" << RECOELE_fbrem_mode[index]
		<< std::endl;

      // Corrections
      /*        RECOELE_isEcalEnergyCorrected[index] = cand->isEcalEnergyCorrected(); */
              RECOELE_ecalEnergy[index]            = cand->ecalEnergy(); 
      /*        RECOELE_ecalEnergyError[index]       = cand->ecalEnergyError(); */
      /*        RECOELE_isMomentumCorrected[index]   = cand->isMomentumCorrected(); */
      /*        RECOELE_trackMomentumError[index]    = cand->trackMomentumError(); */
      /*        RECOELE_electronMomentumError[index] = cand->electronMomentumError(); */
      

    
      // Seed Collection
      if (useRECOformat) {
	edm::RefToBase<TrajectorySeed> seed = cand->gsfTrack()->extra()->seedRef();
	reco::ElectronSeedRef MyS = seed.castTo<reco::ElectronSeedRef>();
	ele_seedSubdet2[index] = int(MyS->subDet2());
	if(fabs(MyS->dPhi2()) < 100.) ele_seedDphi2[index] = double(MyS->dPhi2());
	if(fabs(MyS->dRz2()) < 100.)  ele_seedDrz2[index]  = double(MyS->dRz2());
	
	ele_seedSubdet1[index] = int(MyS->subDet1());
	if(fabs(MyS->dPhi1()) < 100.) ele_seedDphi1[index] = double(MyS->dPhi1());
	if(fabs(MyS->dRz1()) < 100.)  ele_seedDrz1[index]  = double(MyS->dRz1());
      }
      //
      edm::Handle<edm::View<pat::Electron> > elecoll;
      iEvent.getByToken(mvaElectronTag_, elecoll);
      int indextris=0;
      for (edm::View<pat::Electron>::const_iterator candid = elecoll->begin(); candid!=elecoll->end(); ++candid) {
	edm::Ref<edm::View<pat::Electron> > eletrackrefa(elecoll,indextris);
	
	//cout << "Elebis" << RECOELE_PT[index] << " " << candid->p4().pt() << " " << RECOELE_sigmaIetaIeta[index] << " " << candid->sigmaIetaIeta() << " " << fabs(RECOELE_sigmaIetaIeta[index]-candid->sigmaIetaIeta()) << endl;       

	if ( RECOELE_sigmaEtaEta[index]==candid->sigmaEtaEta() && RECOELE_sigmaIetaIeta[index]==candid->sigmaIetaIeta() ){
	  //cout << "Trovato electtrone= " << RECOELE_PT[index] << " " << candid->p4().pt() << endl;
          //qier test 
           cout << "test11" << endl;
//	   RECOELE_mvaTrigV0[index]=(*mvaTrigV0map)[eletrackrefa];
	   RECOELE_mvaNonTrigV0[index]=(*mvaNonTrigV0map)[eletrackrefa];
	   
	   std::cout << "BDT MVA eleID flag = "
		   //  << RECOELE_mvaTrigV0[index] << " "
		     << RECOELE_mvaNonTrigV0[index] << " "
		     << std::endl;
	   
	}
      	indextris++;
      }



      // Conversion Finder
/*
      ConvMapDist[index]=(*conversiondistmap)[eletrackref];
      ConvMapDcot[index]=(*conversiondcotmap)[eletrackref];
      std::cout << "Conversion finder = "
		<< ConvMapDist[index] << " " 
		<< ConvMapDcot[index] << " " 
		<< std::endl;
*/
/*     
      // Matching
      if (fillMCTruth==true){
	int i=0;
	for ( reco::CandidateCollection::const_iterator hIter=CollEle->begin(); hIter!= CollEle->end(); ++hIter ){
	  //cout << "Reco Electron with pT= " << hIter->pt() << " " << RECOELE_PT[index] << " " << fabs(hIter->pt()-RECOELE_PT[index]) << " and mass="<< hIter->mass()<< endl;
	  if (fabs(hIter->pt()-RECOELE_PT[index])<0.01){
	    i=hIter-(CollEle->begin());
	    CandidateRef Ref( CollEle, i );
	    edm::Ref<std::vector<reco::GenParticle> > genrefEle = (*GenParticlesMatchEle)[Ref];
	    if (!genrefEle.isNull()){
	      cout << "GenElectron with pT= " << genrefEle->p4().pt() << " and mass="<< genrefEle->p4().mass()<< endl;
	      RECOELE_MatchingMCTruth[i]= true;
	      RECOELE_MatchingMCpT[i]= genrefEle->p4().pt();
	      RECOELE_MatchingMCEta[i]= genrefEle->p4().eta();
	      RECOELE_MatchingMCPhi[i]= genrefEle->p4().phi();	    
	    } 
	    else {
	      cout << "There is no a reference to a genElectron" << endl;
	    }
	  }   
	}    
      }
*/
      index ++;
    }    
  }
  

  void fillMuons(const edm::Event& iEvent,const edm::EventSetup& iSetup){

    // Get the B-field
    //edm::ESHandle<MagneticField> magfield_;
    iSetup.get<IdealMagneticFieldRecord>().get( magfield_ );        


    // Muons
    edm::Handle<edm::View<pat::Muon> > MuCandidates;
    iEvent.getByToken(muonTag_, MuCandidates);

    edm::Handle<edm::ValueMap<float> > corrpterrormumap;
    iEvent.getByToken(muonCorrPtErrorMapTag_,corrpterrormumap);
    
    // Particle Flow Isolation
/*
    edm::Handle<edm::ValueMap<double> > isoPFChargedAllmumap;
    iEvent.getByToken(muonPFIsoValueChargedAllTag_, isoPFChargedAllmumap);
    
    edm::Handle<edm::ValueMap<double> > isoPFChargedmumap;
    iEvent.getByToken(muonPFIsoValueChargedTag_, isoPFChargedmumap);
    
    edm::Handle<edm::ValueMap<double> > isoPFNeutralmumap;
    iEvent.getByToken(muonPFIsoValueNeutralTag_, isoPFNeutralmumap);
    
    edm::Handle<edm::ValueMap<double> > isoPFGammamumap;
    iEvent.getByToken(muonPFIsoValueGammaTag_, isoPFGammamumap);
    
    edm::Handle<edm::ValueMap<double> > isoPFPUmumap;
    iEvent.getByToken(muonPFIsoValuePUTag_, isoPFPUmumap);
*/

    // Vertexing
    // 3D
    edm::Handle<edm::View<pat::Muon> > VertMuCandidates;
    iEvent.getByLabel(muonTag_Vert, VertMuCandidates);
    
    // 3D w.r.t primary vertex DA
    edm::Handle<edm::ValueMap<float> > vertexmumap;
    iEvent.getByToken(muonMapTag_Vert, vertexmumap);
    
    edm::Handle<edm::ValueMap<float> > vertexmumapvalue;
    iEvent.getByToken(muonMapTag_VertValue, vertexmumapvalue);
    
    edm::Handle<edm::ValueMap<float> > vertexmumaperror;
    iEvent.getByToken(muonMapTag_VertError, vertexmumaperror);

    // 3D w.r.t primary vertex KF
    edm::Handle<edm::ValueMap<float> > vertexmumapKF;
//    iEvent.getByToken(muonMapTag_VertKF, vertexmumapKF);
    
    edm::Handle<edm::ValueMap<float> > vertexmumapvalueKF;
//    iEvent.getByToken(muonMapTag_VertValueKF, vertexmumapvalueKF);
    
    edm::Handle<edm::ValueMap<float> > vertexmumaperrorKF;
//    iEvent.getByToken(muonMapTag_VertErrorKF, vertexmumaperrorKF);


    // STIP SLIP
    edm::Handle<edm::ValueMap<float> > stipmumap;
//    iEvent.getByToken(muonSTIPMapTag_Vert, stipmumap);
    
    edm::Handle<edm::ValueMap<float> > slipmumap;
//    iEvent.getByToken(muonSLIPMapTag_Vert, slipmumap);
    
    edm::Handle<edm::ValueMap<float> > stipmumapvalue;
//    iEvent.getByToken(muonSTIPMapTag_VertValue, stipmumapvalue);
    
    edm::Handle<edm::ValueMap<float> > slipmumapvalue;
//    iEvent.getByToken(muonSLIPMapTag_VertValue, slipmumapvalue);
    
    edm::Handle<edm::ValueMap<float> > stipmumaperror;
//    iEvent.getByToken(muonSTIPMapTag_VertError, stipmumaperror);
    
    edm::Handle<edm::ValueMap<float> > slipmumaperror;
//    iEvent.getByToken(muonSLIPMapTag_VertError, slipmumaperror);
        
    // primary vertex
 
    Vertex primVertex;
    math::XYZPoint pVertex(0., 0., 0.);
    bool pvfound = (PV.size() != 0);
    cout << "pvfound=" << pvfound << endl;
    if(pvfound){       
      for(reco::VertexCollection::const_iterator it=PV.begin() ; it!=PV.end() ; ++it){
	if(!it->isFake() && it->ndof() > 4 && fabs(it->position().z()) <= 24 && fabs(it->position().rho()) <= 2){ 	  	  
	  pVertex = math::XYZPoint(it->position().x(), it->position().y(), it->position().z());
	  cout << "P vertex position used for computing dxy and dz for electron and muons is (x,y,z)= " << 
	    pVertex.x() << " " <<
	    pVertex.y() << " " <<
	    pVertex.z() << endl;
	  break;
	}
      }
    }
    
    //pVertex = math::XYZPoint(                                                                                                                                                    
    //			     PV.at(2).position().x(),                                                                                                                       
    //			     PV.at(2).position().y(),                                                                                                                               
    //			     PV.at(2).position().z());                            

    int indexbis=0;
    RECO_NMU=MuCandidates->size();
    cout << "Event has muon numbers of " << RECO_NMU << endl;

/*
    // Matching
    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchMu;
    iEvent.getByToken(goodMuonMCMatch_, GenParticlesMatchMu);
    edm::Handle<reco::CandidateCollection > CollMu;
    iEvent.getByToken(myMuons_, CollMu);

    if (GenParticlesMatchMu.isValid()){
      cout << endl<< "Muons:"<<endl<<"The reco collection to be matched has size= " <<  CollMu->size() << endl;
      cout << "The matched map collection has size= " <<  GenParticlesMatchMu->size() << endl;
    }
*/
    //


    for (edm::View<pat::Muon>::const_iterator cand = MuCandidates->begin(); cand != MuCandidates->end(); ++cand) {
            
      if(indexbis>99) break;
      
      edm::Ref<edm::View<pat::Muon> > mutrackref(MuCandidates,indexbis); 
      edm::Ref<edm::View<pat::Muon> > mutrackrefv(VertMuCandidates,indexbis); 
     
      RECOMU_isPFMu[indexbis]=cand->isPFMuon();

      RECOMU_isGlobalMu[indexbis]=cand->isGlobalMuon();
      RECOMU_isStandAloneMu[indexbis]=cand->isStandAloneMuon();
      RECOMU_isTrackerMu[indexbis]=cand->isTrackerMuon();
      RECOMU_isCaloMu[indexbis]=cand->isCaloMuon();
      RECOMU_isTrackerHighPtMu[indexbis]=isTrackerHighPtMu(*cand,pVertex);   

      std::cout << "\n Muon in the event: "
	        <<   "  isPF=" << RECOMU_isPFMu[indexbis]
		<<   "  isGB=" << RECOMU_isGlobalMu[indexbis]
		<< "   isSTA=" << RECOMU_isStandAloneMu[indexbis]
		<<   "  isTM=" << RECOMU_isTrackerMu[indexbis]
		<< "  isCalo=" << RECOMU_isCaloMu[indexbis]
		<< std::endl;

      // Kinematic of muon
      RECOMU_E[indexbis]=cand->p4().energy();
      RECOMU_PT[indexbis]=cand->p4().pt();
      RECOMU_P[indexbis]=sqrt(cand->p4().px()*cand->p4().px()+cand->p4().py()*cand->p4().py()+cand->p4().pz()*cand->p4().pz());
      RECOMU_ETA[indexbis]=cand->p4().eta();
      RECOMU_THETA[indexbis]=cand->p4().theta();
      RECOMU_PHI[indexbis]=cand->p4().phi();
      RECOMU_MASS[indexbis]=cand->p4().mass();
      RECOMU_CHARGE[indexbis]=cand->charge();

      std::cout << "--kinematic:"
		<< "  pT="     << RECOMU_PT[indexbis]
	        << "  E="      << RECOMU_E[indexbis]
		<< "  p="      << RECOMU_P[indexbis]
		<< "  eta="    << RECOMU_ETA[indexbis]
		<< "  theta="  << RECOMU_THETA[indexbis]
		<< "  phi="    << RECOMU_PHI[indexbis]
		<< "  mass="   << RECOMU_MASS[indexbis]
		<< "  charge=" << RECOMU_CHARGE[indexbis]
		<< std::endl;
	
      // Covariance matrix
      
      if(cand->innerTrack().isAvailable()){

	GlobalTrajectoryParameters gp(GlobalPoint(cand->innerTrack()->vx(), cand->innerTrack()->vy(),  cand->innerTrack()->vz()),
				      GlobalVector(cand->innerTrack()->px(),cand->innerTrack()->py(),cand->innerTrack()->pz()),
				      cand->innerTrack()->charge(),
				      magfield_.product());
	JacobianCurvilinearToCartesian curv2cart(gp);
	CartesianTrajectoryError cartErr= ROOT::Math::Similarity(curv2cart.jacobian(), cand->innerTrack()->covariance());
	const AlgebraicSymMatrix66 mat = cartErr.matrix();
	for (int i = 0; i < 3; ++i) { 
	  for (int j = 0; j < 3; ++j) { 
	    //bigCov(i,j) = mat[i+3][j+3]; 
	    //	  cout << "index= " << indexbis << "i= " << i << "j= " << j <<endl;
	    RECOMU_COV[indexbis][i][j]=mat[i+3][j+3]; 
	  }  
	}  
	
      }

      // Isolation
      //      RECOMU_TRACKISO[indexbis]=(*isoTkmumap)[mutrackref]/cand->p4().pt();
      RECOMU_TRACKISO_SUMPT[indexbis]=(cand->isolationR03().sumPt)/cand->p4().pt();
      //RECOMU_ECALISO[indexbis]=(*isoEcalmumap)[mutrackref]/cand->p4().pt();
      //RECOMU_HCALISO[indexbis]=(*isoHcalmumap)[mutrackref]/cand->p4().pt();
      //RECOMU_X[indexbis]=[(*isoTkmumap)[mutrackref]+(*isoEcalmumap)[mutrackref]+(*isoHcalmumap)[mutrackref]]/cand->p4().pt();
      // temporary solution for reducedrechit problem
      RECOMU_ECALISO[indexbis]=(cand->isolationR03().emEt)/cand->p4().pt();
      RECOMU_HCALISO[indexbis]=(cand->isolationR03().hadEt)/cand->p4().pt();
      RECOMU_X[indexbis]=RECOMU_TRACKISO[indexbis]+RECOMU_ECALISO[indexbis]+RECOMU_HCALISO[indexbis];

      RECOMU_PFchHad[indexbis]  = (cand->pfIsolationR03().sumChargedHadronPt);
      RECOMU_PFneuHad[indexbis] = (cand->pfIsolationR03().sumNeutralHadronEt);
      RECOMU_PFphoton[indexbis] = (cand->pfIsolationR03().sumPhotonEt);
      RECOMU_PFPUchAllPart[indexbis]= (cand->pfIsolationR03().sumPUPt);

//      RECOMU_PFchHad[indexbis]          =  (*isoPFChargedmumap)[mutrackref];
//      RECOMU_PFneuHad[indexbis]         =  (*isoPFNeutralmumap)[mutrackref];
//      RECOMU_PFphoton[indexbis]         =  (*isoPFGammamumap)[mutrackref];
//      RECOMU_PFPUchAllPart[indexbis]    =  (*isoPFPUmumap)[mutrackref];

      RECOMU_PFX_dB[indexbis]=(RECOMU_PFchHad[indexbis]+max(0.,RECOMU_PFneuHad[indexbis]+RECOMU_PFphoton[indexbis]-0.5*RECOMU_PFPUchAllPart[indexbis]))/cand->p4().pt();

      float EffectiveArea=0.;
      if (use2011EA){
	if (fabs(RECOMU_ETA[indexbis]) >= 0.0 && fabs(RECOMU_ETA[indexbis]) < 1.0 ) EffectiveArea = 0.132;
	if (fabs(RECOMU_ETA[indexbis]) >= 1.0 && fabs(RECOMU_ETA[indexbis]) < 1.5 ) EffectiveArea = 0.120;
	if (fabs(RECOMU_ETA[indexbis]) >= 1.5 && fabs(RECOMU_ETA[indexbis]) < 2.0 ) EffectiveArea = 0.114;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.0 && fabs(RECOMU_ETA[indexbis]) < 2.2 ) EffectiveArea = 0.139;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.2 && fabs(RECOMU_ETA[indexbis]) < 2.3 ) EffectiveArea = 0.168;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.3 )                                     EffectiveArea = 0.189;
      }
      else {       
	if (fabs(RECOMU_ETA[indexbis]) >= 0.0 && fabs(RECOMU_ETA[indexbis]) < 1.0 ) EffectiveArea = 0.674;
	if (fabs(RECOMU_ETA[indexbis]) >= 1.0 && fabs(RECOMU_ETA[indexbis]) < 1.5 ) EffectiveArea = 0.565;
	if (fabs(RECOMU_ETA[indexbis]) >= 1.5 && fabs(RECOMU_ETA[indexbis]) < 2.0 ) EffectiveArea = 0.442;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.0 && fabs(RECOMU_ETA[indexbis]) < 2.2 ) EffectiveArea = 0.515;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.2 && fabs(RECOMU_ETA[indexbis]) < 2.3 ) EffectiveArea = 0.821;
	if (fabs(RECOMU_ETA[indexbis]) >= 2.3 )                                     EffectiveArea = 0.660;
      }
       
      RECOMU_PFX_rho[indexbis]=(RECOMU_PFchHad[indexbis]+max( (RECOMU_PFneuHad[indexbis]+RECOMU_PFphoton[indexbis]-max(RHO_mu,0.0)*(EffectiveArea)),0.0) 
)/double(cand->p4().pt());


      std::cout << "--isolation: muon"
		<< "  X="       << RECOMU_X[indexbis]
		<< "  Tracker=" << RECOMU_TRACKISO[indexbis] <<  "  no veto=" << RECOMU_TRACKISO_SUMPT[indexbis]
		<< "  Ecal="    << RECOMU_ECALISO[indexbis]
		<< "  Hcal="    << RECOMU_HCALISO[indexbis]
	        << endl;

      std::cout << "--PF isolation"
	        << "  chHad="   << RECOMU_PFchHad[indexbis] 
	        << "  neuHad="  << RECOMU_PFneuHad[indexbis]
	        << "  photon="  << RECOMU_PFphoton[indexbis] 
	        << "  sumPUPt=" << RECOMU_PFPUchAllPart[indexbis]
                << "  PFX_dB="  << RECOMU_PFX_dB[indexbis] 
                << "  PFX_rho=" << RECOMU_PFX_rho[indexbis]
		<< std::endl;

      // Vertexing

      //Float_t mySIP = cand->userFloat("SIP");
      //Float_t mySIP = std::abs(cand->dB(pat::Muon::PV3D))/cand->edB(pat::Muon::PV3D);;
      //cout << "mySIP= " <<  mySIP << endl;
      RECOMU_SIP[indexbis]=(*vertexmumap)[mutrackrefv];
      RECOMU_IP[indexbis]=(*vertexmumapvalue)[mutrackrefv];
      RECOMU_IPERROR[indexbis]=(*vertexmumaperror)[mutrackrefv];
/*
      RECOMU_SIP_KF[indexbis]=(*vertexmumapKF)[mutrackrefv];
      RECOMU_IP_KF[indexbis]=(*vertexmumapvalueKF)[mutrackrefv];
      RECOMU_IPERROR_KF[indexbis]=(*vertexmumaperrorKF)[mutrackrefv];
*/      

      //RECOMU_SIP_GD[indexbis]=(*vertexmumapGD)[mutrackrefv];
      //if (decaychannel=="4mu" || decaychannel=="2e2mu" ) RECOMU_SIP_GDMMMM[indexbis]=(*vertexmumapGDMMMM)[mutrackrefv];
      //      RECOMU_SIP_Std[indexbis]=(*vertexmumapStd)[mutrackrefv];
      //if (decaychannel=="4mu" || decaychannel=="2e2mu" ) RECOMU_SIP_StdMMMM[indexbis]=(*vertexmumapStdMMMM)[mutrackrefv];
      //RECOMU_SIP_Kin[indexbis]=(*vertexmumapKin)[mutrackrefv];
      //if (decaychannel=="4mu" || decaychannel=="2e2mu" ) RECOMU_SIP_KinMMMM[indexbis]=(*vertexmumapKinMMMM)[mutrackrefv];

/*
      RECOMU_STIP[indexbis]=(*stipmumap)[mutrackrefv];
      RECOMU_SLIP[indexbis]=(*slipmumap)[mutrackrefv];
      RECOMU_TIP[indexbis]=(*stipmumapvalue)[mutrackrefv];
      RECOMU_LIP[indexbis]=(*slipmumapvalue)[mutrackrefv];
      RECOMU_TIPERROR[indexbis]=(*stipmumaperror)[mutrackrefv];
      RECOMU_LIPERROR[indexbis]=(*slipmumaperror)[mutrackrefv];
*/      

      std::cout << "--vertexing: muon"
		<< "  Sign3DIP=" << RECOMU_SIP[indexbis]
		<< "  3DIP="     << RECOMU_IP[indexbis]
		<< "  3DIPerr="  << RECOMU_IPERROR[indexbis]
		<< "  SignTIP="  << RECOMU_STIP[indexbis]
		<< "  TIP="      << RECOMU_TIP[indexbis]
		<< "  TIPerr="   << RECOMU_TIPERROR[indexbis]
		<< "  SignLIP="  << RECOMU_SLIP[indexbis]
		<< "  LIP="      << RECOMU_LIP[indexbis]
		<< "  LIPerror=" << RECOMU_LIPERROR[indexbis]
		<< std::endl;
      
      // Other properties
      RECOMU_numberOfMatches[indexbis]=cand->numberOfMatches();
      RECOMU_numberOfMatchedStations[indexbis]=cand->numberOfMatchedStations();
      RECOMU_caloCompatibility[indexbis]=cand->caloCompatibility();
      RECOMU_segmentCompatibility[indexbis]=(muon::segmentCompatibility( (*cand)));
      RECOMU_glbmuPromptTight[indexbis]=(muon::isGoodMuon( (*cand),muon::GlobalMuonPromptTight));


      std::cout	<< "--other properties:"
		<< "  n.matches="   << RECOMU_numberOfMatches[indexbis]
		<< "  caloComp="    << RECOMU_caloCompatibility[indexbis]
		<< "  segmentComp=" << RECOMU_segmentCompatibility[indexbis]
	        << "  glbmuPromptTight=" << RECOMU_glbmuPromptTight[indexbis]
		<< std::endl;

  
      // Track properties
      if(cand->muonBestTrack().isAvailable()){
	RECOMU_mubesttrkType[indexbis]=cand->muonBestTrackType();
	RECOMU_mubesttrkDxy[indexbis]=cand->muonBestTrack()->dxy(pVertex);
	RECOMU_mubesttrkDxyB[indexbis]=cand->muonBestTrack()->dxy(bs.position());
	RECOMU_mubesttrkDxyError[indexbis]=cand->muonBestTrack()->dxyError();
	RECOMU_mubesttrkDz[indexbis]=cand->muonBestTrack()->dz(pVertex);
	RECOMU_mubesttrkDzB[indexbis]=cand->muonBestTrack()->dz(bs.position());
	RECOMU_mubesttrkDzError[indexbis]=cand->muonBestTrack()->dzError();
	//RECOMU_mubesttrkPTError[indexbis]=(*corrpterrormumap)[mutrackref];;
      }

      if(cand->globalTrack().isAvailable()){
	RECOMU_mutrkPT[indexbis]=cand->globalTrack()->pt();
	RECOMU_mutrkPTError[indexbis]=cand->globalTrack()->ptError();
	//RECOMU_mutrkPTError[indexbis]=cand->bestTrack()->ptError(); // bestTrack
	RECOMU_mutrkDxy[indexbis]=cand->globalTrack()->dxy(pVertex);
	RECOMU_mutrkDxyError[indexbis]=cand->globalTrack()->dxyError();
	RECOMU_mutrkDxyB[indexbis]=cand->globalTrack()->dxy(bs.position()) ;
	RECOMU_mutrkDz[indexbis]=cand->globalTrack()->dz(pVertex);
	RECOMU_mutrkDzError[indexbis]=cand->globalTrack()->dzError();
	RECOMU_mutrkDzB[indexbis]=cand->globalTrack()->dz(bs.position());
	RECOMU_mutrkChi2PerNdof[indexbis]=cand->globalTrack()->normalizedChi2();
	RECOMU_mutrkCharge[indexbis]=cand->globalTrack()->charge();
 	RECOMU_mutrkNHits[indexbis]=cand->globalTrack()->numberOfValidHits(); 
	RECOMU_mutrkNPixHits[indexbis]=cand->globalTrack()->hitPattern().numberOfValidPixelHits();
        RECOMU_mutrkNStripHits[indexbis]=cand->globalTrack()->hitPattern().numberOfValidStripHits();
	RECOMU_mutrkNMuonHits[indexbis]=cand->globalTrack()->hitPattern().numberOfValidMuonHits(); 
	RECOMU_mutrktrackerLayersWithMeasurement[indexbis]=cand->globalTrack()->hitPattern().trackerLayersWithMeasurement(); 

	RECOMU_muInnertrkDxy[indexbis]=cand->innerTrack()->dxy(pVertex);
	RECOMU_muInnertrkDxyError[indexbis]=cand->innerTrack()->dxyError();
	RECOMU_muInnertrkDxyB[indexbis]=cand->innerTrack()->dxy(bs.position()) ;
	RECOMU_muInnertrkDz[indexbis]=cand->innerTrack()->dz(pVertex);
	RECOMU_muInnertrkDzError[indexbis]=cand->innerTrack()->dzError();
	RECOMU_muInnertrkDzB[indexbis]=cand->innerTrack()->dz(bs.position());
	RECOMU_muInnertrkChi2PerNdof[indexbis]=cand->innerTrack()->normalizedChi2();
	RECOMU_muInnertrktrackerLayersWithMeasurement[indexbis]=cand->innerTrack()->hitPattern().trackerLayersWithMeasurement(); 
	RECOMU_muInnertrkPT[indexbis]=cand->innerTrack()->pt();	
	//RECOMU_muInnertrkPTError[indexbis]=cand->innerTrack()->ptError();
	RECOMU_muInnertrkPTError[indexbis]=cand->bestTrack()->ptError(); // Besttrack

	RECOMU_muInnertrkCharge[indexbis]=cand->innerTrack()->charge();
 	RECOMU_muInnertrkNHits[indexbis]=cand->innerTrack()->numberOfValidHits(); 
	RECOMU_muInnertrkNPixHits[indexbis]=cand->innerTrack()->hitPattern().numberOfValidPixelHits();
        RECOMU_muInnertrkNStripHits[indexbis]=cand->innerTrack()->hitPattern().numberOfValidStripHits();

      }
      else if(cand->innerTrack().isAvailable()){
	RECOMU_muInnertrkDxy[indexbis]=cand->innerTrack()->dxy(pVertex);
	RECOMU_muInnertrkDxyError[indexbis]=cand->innerTrack()->dxyError();
	RECOMU_muInnertrkDxyB[indexbis]=cand->innerTrack()->dxy(bs.position()) ;
	RECOMU_muInnertrkDz[indexbis]=cand->innerTrack()->dz(pVertex);
	RECOMU_muInnertrkDzError[indexbis]=cand->innerTrack()->dzError();
	RECOMU_muInnertrkDzB[indexbis]=cand->innerTrack()->dz(bs.position());
	RECOMU_muInnertrkChi2PerNdof[indexbis]=cand->innerTrack()->normalizedChi2();
	RECOMU_muInnertrktrackerLayersWithMeasurement[indexbis]=cand->innerTrack()->hitPattern().trackerLayersWithMeasurement(); 
	RECOMU_muInnertrkPT[indexbis]=cand->innerTrack()->pt();	
	//RECOMU_muInnertrkPTError[indexbis]=cand->innerTrack()->ptError();
	RECOMU_muInnertrkPTError[indexbis]=cand->bestTrack()->ptError(); // Besttrack
	RECOMU_muInnertrkCharge[indexbis]=cand->innerTrack()->charge();
 	RECOMU_muInnertrkNHits[indexbis]=cand->innerTrack()->numberOfValidHits(); 
	RECOMU_muInnertrkNPixHits[indexbis]=cand->innerTrack()->hitPattern().numberOfValidPixelHits();
        RECOMU_muInnertrkNStripHits[indexbis]=cand->innerTrack()->hitPattern().numberOfValidStripHits();
      }

      if(cand->globalTrack().isAvailable() || cand->innerTrack().isAvailable() ){
	std::cout << "--muon track properties: "
	          << "  pt="        << RECOMU_mutrkPT[indexbis] 
	          << "  bestTrackType= " << RECOMU_mubesttrkType[indexbis]
		  << "  dxy="       << RECOMU_mubesttrkDxy[indexbis]
		  << "  dxyError="  << RECOMU_mubesttrkDxyError[indexbis]
		  << "  dxyB="      << RECOMU_mubesttrkDxyB[indexbis] 
		  << "  dz="        << RECOMU_mubesttrkDz[indexbis] 
		  << "  dzError="   << RECOMU_mubesttrkDzError[indexbis]
		  << "  dzB="       << RECOMU_mubesttrkDzB[indexbis]
		  << "  PtError="   << RECOMU_mubesttrkPTError[indexbis]
		  << "  chi2_nodf=" << RECOMU_mutrkChi2PerNdof[indexbis]
		  << "  charge="    << RECOMU_mutrkCharge[indexbis]
		  << "  nhits="     << RECOMU_mutrkNHits[indexbis] 
		  << "  nPixhits="   << RECOMU_mutrkNPixHits[indexbis]
	          << "  nStriphits=" << RECOMU_mutrkNStripHits[indexbis]
	          << "  nMuonhits="  << RECOMU_mutrkNMuonHits[indexbis]
		  << std::endl;
	
	
	// Tracker muon properties
	// RECOMU_trkmuArbitration[indexbis]=(muon::isGoodMuon( (*cand),muon::TrackerMuonArbitrated));
	RECOMU_trkmuArbitration[indexbis]=(muon::segmentCompatibility((*cand),pat::Muon::SegmentAndTrackArbitration));
	RECOMU_trkmu2DCompatibilityLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TM2DCompatibilityLoose));
	RECOMU_trkmu2DCompatibilityTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TM2DCompatibilityTight));
	RECOMU_trkmuOneStationLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationLoose));
	RECOMU_trkmuOneStationTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationTight));
	RECOMU_trkmuLastStationLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationLoose));
	RECOMU_trkmuLastStationTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationTight));
	RECOMU_trkmuOneStationAngLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationAngLoose));
	RECOMU_trkmuOneStationAngTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationAngTight));
	RECOMU_trkmuLastStationAngLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationAngLoose));
	RECOMU_trkmuLastStationAngTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationAngTight));
	RECOMU_trkmuLastStationOptimizedLowPtLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationOptimizedLowPtLoose));
	RECOMU_trkmuLastStationOptimizedLowPtTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationOptimizedLowPtTight));
	
	std::cout << "--tracker muon properties:"
		  << "  arbitration="              << RECOMU_trkmuArbitration[indexbis]
		  << "  2DCompLoose="              << RECOMU_trkmu2DCompatibilityLoose[indexbis]
		  << "  2DCompTight="              << RECOMU_trkmu2DCompatibilityTight[indexbis]
		  << "  1StationLoose="            << RECOMU_trkmuOneStationLoose[indexbis]
		  << "  1StationTight="            << RECOMU_trkmuOneStationTight[indexbis]
		  << "  LastStationLoose="         << RECOMU_trkmuLastStationLoose[indexbis]
		  << "  LastStationTight="         << RECOMU_trkmuLastStationTight[indexbis]
		  << "  1StationAngLoose="         << RECOMU_trkmuOneStationAngLoose[indexbis]
		  << "  1StationAngTight="         << RECOMU_trkmuOneStationAngTight[indexbis]
		  << "  LastStationAngLoose="      << RECOMU_trkmuLastStationAngLoose[indexbis]
		  << "  LastStationAngTight="      << RECOMU_trkmuLastStationAngTight[indexbis]
		  << "  LastStationOptLowptLoose=" << RECOMU_trkmuLastStationOptimizedLowPtLoose[indexbis]
		  << "  LastStationOptLowptTight=" << RECOMU_trkmuLastStationOptimizedLowPtTight[indexbis]
		  << std::endl;
      }
  
/*    
      // Matching
      if (fillMCTruth==true){
	int i=0;
	for ( reco::CandidateCollection::const_iterator hIter=CollMu->begin(); hIter!= CollMu->end(); ++hIter ){
	  //cout << "Reco Muon with pT= " << hIter->pt() << " and mass="<< hIter->mass()<< endl;
	  if (fabs(hIter->pt()-RECOMU_PT[indexbis])<0.01){
	    i=hIter-(CollMu->begin());
	    CandidateRef Ref( CollMu, i );
	    edm::Ref<std::vector<reco::GenParticle> > genrefMu = (*GenParticlesMatchMu)[Ref];
	    if (!genrefMu.isNull()){
	      cout << "GenMuon with pT= " << genrefMu->p4().pt() << " and mass="<< genrefMu->p4().mass()<< endl;
	      RECOMU_MatchingMCTruth[i]= true;
	      RECOMU_MatchingMCpT[i]= genrefMu->p4().pt();
	      RECOMU_MatchingMCEta[i]= genrefMu->p4().eta();
	      RECOMU_MatchingMCPhi[i]= genrefMu->p4().phi();	    
	    } 
	    else {
	      cout << "There is no reference to a genMuon" << endl;
	    }
	  }   
	}   
      }  
*/ 
      indexbis++;
    }
  }

  void fillPhotons(const edm::Event& iEvent){
    // Photons
    //edm::Handle<edm::View<reco::Candidate> > photons;
    edm::Handle<edm::View<pat::Photon> > photons;
    iEvent.getByToken(photonsTag_, photons);
    RECO_NPHOT=photons->size();
/*
    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchPhot;
    iEvent.getByToken(goodGammaMCMatch_, GenParticlesMatchPhot);
    edm::Handle<reco::CandidateCollection > CollPhot;
    iEvent.getByToken(myGammas_, CollPhot);
    bool ismyGammas=false;
    if (CollPhot.isValid()) ismyGammas=true;
*/    
    int iphot=0;
    for (edm::View<pat::Photon>::const_iterator cand = photons->begin(); cand != photons->end(); ++cand) {
      if (iphot>19) break;
      RECOPHOT_PT[iphot]=cand->pt();
      RECOPHOT_ETA[iphot]=cand->eta();
      RECOPHOT_PHI[iphot]=cand->phi();
      RECOPHOT_THETA[iphot]=cand->theta();
/*
      if (!cand->parentSuperCluster().isNull()){
        RECOPHOT_TLE_ParentSC_X[iphot]=cand->parentSuperCluster()->x();
        RECOPHOT_TLE_ParentSC_Y[iphot]=cand->parentSuperCluster()->y();
        RECOPHOT_TLE_ParentSC_Z[iphot]=cand->parentSuperCluster()->z();
      }

      cout << "Parent SC x,y,x= " <<RECOPHOT_TLE_ParentSC_X[iphot] << " " << RECOPHOT_TLE_ParentSC_Y[iphot] << " " << RECOPHOT_TLE_ParentSC_Z[iphot] << endl;
*/

      cout << "Reco Photon with pT= " << RECOPHOT_PT[iphot] << " eta= " << RECOPHOT_ETA[iphot] << " phi= " << RECOPHOT_PHI[iphot] << endl;
/*
      // Matching
      int i=0;
      if (ismyGammas){
	for ( reco::CandidateCollection::const_iterator hIter=CollPhot->begin(); hIter!= CollPhot->end(); ++hIter ){
	  //cout << "Reco Photon with pT= " << hIter->pt() << " and mass="<< hIter->mass()<< endl;
	  if (fabs(hIter->pt()-RECOPHOT_PT[iphot])<0.01){
	    i=hIter-(CollPhot->begin());
	    CandidateRef Ref( CollPhot, i );
	    edm::Ref<std::vector<reco::GenParticle> > genrefPhot = (*GenParticlesMatchPhot)[Ref];
	    if (!genrefPhot.isNull()){
	      cout << "GenMuon with pT= " << genrefPhot->p4().pt() << " and mass="<< genrefPhot->p4().mass()<< endl;
	      RECOPHOT_MatchingMCTruth[i]= true;
	      RECOPHOT_MatchingMCpT[i]= genrefPhot->p4().pt();
	      RECOPHOT_MatchingMCEta[i]= genrefPhot->p4().eta();
	      RECOPHOT_MatchingMCPhi[i]= genrefPhot->p4().phi();	    
	    } 
	  }   
	}
      }  */ 
      iphot++;
    }

    // PF ISR photon

    // Photons
//    edm::Handle<edm::View<pat::PackedCandidate> > pfphotons;
//    iEvent.getByToken(pfphotonsTag_, pfphotons);

/*
    edm::Handle<edm::ValueMap<double> > isoPFChargedAllphmap; 
    iEvent.getByToken(photonPFIsoValueChargedAllTag_, isoPFChargedAllphmap); 
    
    edm::Handle<edm::ValueMap<double> > isoPFChargedphmap; 
    iEvent.getByToken(photonPFIsoValueChargedTag_, isoPFChargedphmap); 
    
    edm::Handle<edm::ValueMap<double> > isoPFNeutralphmap; 
    iEvent.getByToken(photonPFIsoValueNeutralTag_, isoPFNeutralphmap); 
    
    edm::Handle<edm::ValueMap<double> > isoPFGammaphmap; 
    iEvent.getByToken(photonPFIsoValueGammaTag_, isoPFGammaphmap); 
    
    edm::Handle<edm::ValueMap<double> > isoPFPUphmap; 
    iEvent.getByToken(photonPFIsoValuePUTag_, isoPFPUphmap); 
*/

    iphot=0;
    RECO_NPFPHOT=photons->size();
    cout << "There are " << photons->size() << " photons" << endl;
    for (edm::View<pat::Photon>::const_iterator cand = photons->begin(); cand != photons->end(); ++cand) {

      if (iphot>19) break;

      edm::Ref<edm::View<pat::Photon> > phtrackref(photons,iphot); 

      RECOPFPHOT_PT[iphot]=cand->pt();

      // error on pt of photon taken from #include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"
      double C;
      double S;
      double N;
      if(TMath::Abs(cand->eta())<1.48){C=0.35/100; S=5.51/100; N=98./1000.;}
      else{C=0; S=12.8/100; N=440./1000.;} 
      double  perr = TMath::Sqrt(C*C*cand->p4().e()*cand->p4().e() + S*S*cand->p4().e() + N*N);
      double pterr = perr*cand->pt()/cand->p(); 
      RECOPFPHOT_PTError[iphot]=float(pterr);
      
      cout << "PF Photon : pT= " << RECOPFPHOT_PT[iphot] << " pTerr= " << RECOPFPHOT_PTError[iphot]<< endl;

      RECOPFPHOT_ETA[iphot]=cand->eta();
      RECOPFPHOT_PHI[iphot]=cand->phi();
      RECOPFPHOT_THETA[iphot]=cand->theta();
/*      
      RECOPFPHOT_PFchAllPart[iphot]      =  (*isoPFChargedAllphmap)[phtrackref]; 
      RECOPFPHOT_PFchHad[iphot]          =  (*isoPFChargedphmap)[phtrackref]; 
      RECOPFPHOT_PFneuHad[iphot]         =  (*isoPFNeutralphmap)[phtrackref]; 
      RECOPFPHOT_PFphoton[iphot]         =  (*isoPFGammaphmap)[phtrackref]; 
      RECOPFPHOT_PFPUchAllPart[iphot]    =  (*isoPFPUphmap)[phtrackref]; 
*/
//TODO      RECOPFPHOT_PFchAllPart[iphot]      =  (cand->PflowIsolationVariables().sumChargedParticlePt); 
      RECOPFPHOT_PFchHad[iphot]          =  (cand->chargedHadronIso());
      RECOPFPHOT_PFneuHad[iphot]         =  (cand->neutralHadronIso()); 
      RECOPFPHOT_PFphoton[iphot]         =  (cand->photonIso()); 
      RECOPFPHOT_PFPUchAllPart[iphot]    =  (cand->puChargedHadronIso());


      
      
      RECOPFPHOT_PFX_rho[iphot]=( 
 				 RECOPFPHOT_PFchHad[iphot]+ 
 				 RECOPFPHOT_PFneuHad[iphot]+ 
 				 RECOPFPHOT_PFphoton[iphot]+ 
 				 RECOPFPHOT_PFPUchAllPart[iphot] 
				  )/double(cand->p4().pt()); 
      
      iphot++;
      
    }
    
  }

void fillTracks(const edm::Event& iEvent){
    // Tracks
  using namespace edm; using namespace std; using namespace reco;
  edm::Handle<edm::View<pat::PackedCandidate>> cands;
  iEvent.getByToken(pfTag_,cands);
 
  int countk=0;
  RECO_NTRACK=0;
  
//  for(edm::View<pat::PackedCandidate>::const_iterator i=cands->begin(); i!=cands->end(); i++){
  for(unsigned int i=0;i<cands->size();i++){
     if (countk>199) break;
     const pat::PackedCandidate & c = (*cands)[i];
     if(!(c.charge() != 0 && c.numberOfHits()> 0)) continue;
     RECO_NTRACK++;
    RECO_TRACK_PT[countk]=c.pseudoTrack().pt();
    RECO_TRACK_ETA[countk]=c.pseudoTrack().eta();
    RECO_TRACK_PHI[countk]=c.pseudoTrack().phi();
    RECO_TRACK_CHI2[countk]=c.pseudoTrack().chi2();
    RECO_TRACK_CHI2RED[countk]=c.pseudoTrack().normalizedChi2();
    //RECO_TRACK_CHI2PROB=TMath::Prob(i->chi2(),i->ndof());
    RECO_TRACK_CHI2PROB[countk]=ChiSquaredProbability(c.pseudoTrack().chi2(),c.pseudoTrack().ndof());
    RECO_TRACK_NHITS[countk]=c.pseudoTrack().numberOfValidHits();
    RECO_TRACK_DXY[countk]=c.pseudoTrack().dxy();
    RECO_TRACK_DXYERR[countk]=c.pseudoTrack().dxyError();
    RECO_TRACK_DZ[countk]=c.pseudoTrack().dz();
    RECO_TRACK_DZERR[countk]=c.pseudoTrack().dzError();
    countk++;
  }
 cout << "Number of Tracks in the event= " << RECO_NTRACK << endl;
}
        

  
 /* 
  double SetMET(const edm::Event& iEvent, edm::InputTag myTag_ ){
    double met=-999.;
    edm::Handle<reco::METCollection> metHandle;
    iEvent.getByLabel(myTag_,metHandle);
    for ( pat::METCollection::const_iterator iMet=metHandle->begin(); iMet!=metHandle->end(); iMet++) {
      met = iMet->pt();
    }
    return met;
  }
  
  double SetCaloMET(const edm::Event& iEvent, edm::InputTag myTag_ ){
    double met=-999.;
    edm::Handle<reco::CaloMETCollection> metHandle;
    iEvent.getByLabel(myTag_,metHandle);
    for ( CaloMETCollection::const_iterator iMet=metHandle->begin(); iMet!=metHandle->end(); iMet++) {
      met = iMet->pt();
    }
    return met;
  }
  
*/
  void fillMET(const edm::Event& iEvent){
    //TCMET
    //tcmet = SetMET(iEvent,trackermetTag_);
    //CALOMET
    //calomet = SetCaloMET(iEvent,calometTag_);
/*     calometnohf = SetCaloMET(iEvent,calometnohfTag_); */

/*     if (useAdditionalMET_==true){ */
/*       calometho = SetCaloMET(iEvent,calomethoTag_); */
/*       calometopt = SetCaloMET(iEvent,calometoptTag_); */
/*       calometoptnohf =  SetCaloMET(iEvent,calometoptnohfTag_); */
/*       calometoptnohfho =  SetCaloMET(iEvent,calometoptnohfhoTag_); */
/*       calometoptho = SetCaloMET(iEvent,calomethoTag_); */
/*       calometnohfho = SetCaloMET(iEvent,calometnohfhoTag_); */
/*     } */

    //Type I Muon Correction MET
    //cormetmuons = SetCaloMET(iEvent,cormetMuTag_);
    
    //PFMET
    std::cout<<"met uncer"<< pat::MET::METUncertainty(12) <<endl;
    edm::Handle<pat::METCollection> pfmetHandle;
    iEvent.getByToken(pfmetTag_,pfmetHandle);
    for ( pat::METCollection::const_iterator i=pfmetHandle->begin(); i!=pfmetHandle->end(); i++) {
      cormetmuons = i->pt();
      pfmet       = i->uncorPt();     
      pfmet_x     = i->uncorPx();
      pfmet_y     = i->uncorPy();
      pfmet_phi   = i->uncorPhi();
      pfmet_theta = i->uncorP3().theta();
      //std::cout<<"met uncer"<< i->MET::METUncertainty(0) <<endl;
    }




    std::cout << "MET:"
	      << "  tcMET="            << tcmet 
	      << "  caloMET="          << calomet 
/* 	      << "  caloMETopt="       << calometopt */
/* 	      << "  caloMEToptnohf="   << calometoptnohf */
/* 	      << "  caloMEToptnohfho=" << calometoptnohfho */
/* 	      << "  caloMEToptho="     << calometoptho */
/* 	      << "  caloMETnohf="      << calometnohf */
/* 	      << "  caloMETnohfho="    << calometnohfho */
/* 	      << "  caloMETho="        << calometho */
	      << "  TypeIMucorrMET="   << cormetmuons
	      << "  pfMET="            << pfmet
	      << std::endl;
    
  }

  void fillGD2e2mu(const edm::Event& iEvent){
    edm::Handle<vector<double> > GeomD;
    iEvent.getByLabel(ftsigma_Vert, GeomD);
    int jjj=0;
    for (vector<double>::const_iterator cand=GeomD->begin(); cand!=GeomD->end(); ++cand){
      ftsigma[jjj]=(*cand);
      jjj++;
    }

   
    edm::Handle<vector<double> > GeomDlag;
    iEvent.getByLabel(ftsigmalag_Vert, GeomDlag);
    jjj=0;
    for (vector<double>::const_iterator cand=GeomDlag->begin(); cand!=GeomDlag->end(); ++cand){
      ftsigmalag[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdX_;
    iEvent.getByLabel(gdX_Vert, gdX_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdX_->begin(); cand!=gdX_->end(); ++cand){
      gdX[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagX_;
    iEvent.getByLabel(gdlagX_Vert, gdlagX_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagX_->begin(); cand!=gdlagX_->end(); ++cand){
      gdlagX[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdY_;
    iEvent.getByLabel(gdY_Vert, gdY_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdY_->begin(); cand!=gdY_->end(); ++cand){
      gdY[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagY_;
    iEvent.getByLabel(gdlagY_Vert, gdlagY_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagY_->begin(); cand!=gdlagY_->end(); ++cand){
      gdlagY[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdZ_;
    iEvent.getByLabel(gdZ_Vert, gdZ_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdZ_->begin(); cand!=gdZ_->end(); ++cand){
      gdZ[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagZ_;
    iEvent.getByLabel(gdlagZ_Vert, gdlagZ_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagZ_->begin(); cand!=gdlagZ_->end(); ++cand){
      gdlagZ[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagProb_;
    iEvent.getByLabel(gdlagProb_Vert, gdlagProb_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagProb_->begin(); cand!=gdlagProb_->end(); ++cand){
      gdlagProb[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagNdof_;
    iEvent.getByLabel(gdlagNdof_Vert, gdlagNdof_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagNdof_->begin(); cand!=gdlagNdof_->end(); ++cand){
      gdlagNdof[jjj]=(*cand);
      jjj++;
    }

  }

  void fillGD4mu(const edm::Event& iEvent){
    edm::Handle<vector<double> > GeomD_MMMM;
    iEvent.getByLabel(ftsigma_VertMMMM, GeomD_MMMM);
    int jjj=0;
    for (vector<double>::const_iterator cand=GeomD_MMMM->begin(); cand!=GeomD_MMMM->end(); ++cand){
      ftsigmaMMMM[jjj]=(*cand);
      jjj++;
    }
  
    edm::Handle<vector<double> > GeomDlag_MMMM;
    iEvent.getByLabel(ftsigmalag_VertMMMM, GeomDlag_MMMM);
    jjj=0;
    for (vector<double>::const_iterator cand=GeomDlag_MMMM->begin(); cand!=GeomDlag_MMMM->end(); ++cand){
      ftsigmalagMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdXMMMM_;
    iEvent.getByLabel(gdX_VertMMMM, gdXMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdXMMMM_->begin(); cand!=gdXMMMM_->end(); ++cand){
      gdXMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagXMMMM_;
    iEvent.getByLabel(gdlagX_VertMMMM, gdlagXMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagXMMMM_->begin(); cand!=gdlagXMMMM_->end(); ++cand){
      gdlagXMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdYMMMM_;
    iEvent.getByLabel(gdY_VertMMMM, gdYMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdYMMMM_->begin(); cand!=gdYMMMM_->end(); ++cand){
      gdYMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagYMMMM_;
    iEvent.getByLabel(gdlagY_VertMMMM, gdlagYMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagYMMMM_->begin(); cand!=gdlagYMMMM_->end(); ++cand){
      gdlagYMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdZMMMM_;
    iEvent.getByLabel(gdZ_VertMMMM, gdZMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdZMMMM_->begin(); cand!=gdZMMMM_->end(); ++cand){
      gdZMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagZMMMM_;
    iEvent.getByLabel(gdlagZ_VertMMMM, gdlagZMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagZMMMM_->begin(); cand!=gdlagZMMMM_->end(); ++cand){
      gdlagZMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagProbMMMM_;
    iEvent.getByLabel(gdlagProb_VertMMMM, gdlagProbMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagProbMMMM_->begin(); cand!=gdlagProbMMMM_->end(); ++cand){
      gdlagProbMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagNdofMMMM_;
    iEvent.getByLabel(gdlagNdof_VertMMMM, gdlagNdofMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagNdofMMMM_->begin(); cand!=gdlagNdofMMMM_->end(); ++cand){
      gdlagNdofMMMM[jjj]=(*cand);
      jjj++;
    }

  }

  void fillGD4e(const edm::Event& iEvent){
    edm::Handle<vector<double> > GeomD_EEEE;
    iEvent.getByLabel(ftsigma_VertEEEE, GeomD_EEEE);
    int jjj=0;
    for (vector<double>::const_iterator cand=GeomD_EEEE->begin(); cand!=GeomD_EEEE->end(); ++cand){
      ftsigmaEEEE[jjj]=(*cand);
      jjj++;
    }
  
    edm::Handle<vector<double> > GeomDlag_EEEE;
    iEvent.getByLabel(ftsigmalag_VertEEEE, GeomDlag_EEEE);
    jjj=0;
    for (vector<double>::const_iterator cand=GeomDlag_EEEE->begin(); cand!=GeomDlag_EEEE->end(); ++cand){
      ftsigmalagEEEE[jjj]=(*cand);
      jjj++;
    }


    edm::Handle<vector<double> > gdXEEEE_;
    iEvent.getByLabel(gdX_VertEEEE, gdXEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdXEEEE_->begin(); cand!=gdXEEEE_->end(); ++cand){
      gdXEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagXEEEE_;
    iEvent.getByLabel(gdlagX_VertEEEE, gdlagXEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagXEEEE_->begin(); cand!=gdlagXEEEE_->end(); ++cand){
      gdlagXEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdYEEEE_;
    iEvent.getByLabel(gdY_VertEEEE, gdYEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdYEEEE_->begin(); cand!=gdYEEEE_->end(); ++cand){
      gdYEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagYEEEE_;
    iEvent.getByLabel(gdlagY_VertEEEE, gdlagYEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagYEEEE_->begin(); cand!=gdlagYEEEE_->end(); ++cand){
      gdlagYEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdZEEEE_;
    iEvent.getByLabel(gdZ_VertEEEE, gdZEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdZEEEE_->begin(); cand!=gdZEEEE_->end(); ++cand){
      gdZEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagZEEEE_;
    iEvent.getByLabel(gdlagZ_VertEEEE, gdlagZEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagZEEEE_->begin(); cand!=gdlagZEEEE_->end(); ++cand){
      gdlagZEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagProbEEEE_;
    iEvent.getByLabel(gdlagProb_VertEEEE, gdlagProbEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagProbEEEE_->begin(); cand!=gdlagProbEEEE_->end(); ++cand){
      gdlagProbEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagNdofEEEE_;
    iEvent.getByLabel(gdlagNdof_VertEEEE, gdlagNdofEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagNdofEEEE_->begin(); cand!=gdlagNdofEEEE_->end(); ++cand){
      gdlagNdofEEEE[jjj]=(*cand);
      jjj++;
    }

  }

  void fillConstraintVtx2e2mu(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtx_;
    iEvent.getByLabel(StandardFitVertex, StandardFitVtx_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtx_->begin(); cand!=StandardFitVtx_->end(); ++cand){
      if (jjj > 99) break;
      StdFitVertexX[jjj]=cand->position().x();
      StdFitVertexY[jjj]=cand->position().y();
      StdFitVertexZ[jjj]=cand->position().z();
      StdFitVertexChi2r[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProb[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit: " 
	   <<  StdFitVertexX[jjj] << " " 
	   <<  StdFitVertexY[jjj] << " " 
	   <<  StdFitVertexZ[jjj] << " " 
	   <<  StdFitVertexChi2r[jjj] << " " 
	   <<  StdFitVertexProb[jjj] << endl;

      
      // Refitted tracks     
      bool hasRefittedTracks = cand->hasRefittedTracks();
      if(hasRefittedTracks){	
	vector<reco::Track> refit_tks= cand->refittedTracks(); 	
	for(unsigned int i=0; i< refit_tks.size(); i++){	 
	  if (i > 3) break;
	  cout << "Track momentum Refit=" << refit_tks[i].pt() << endl;
	  StdFitVertexTrack_PT[i][jjj] =refit_tks[i].pt() ;
	  StdFitVertexTrack_ETA[i][jjj]=refit_tks[i].eta() ;
	  StdFitVertexTrack_PHI[i][jjj]=refit_tks[i].phi();
	}
      }
                   
      jjj++;
    }

    

    edm::Handle<edm::View<Candidate> > CandidatesEEMM;
    iEvent.getByToken(RECOcollNameEEMM, CandidatesEEMM);
    edm::Handle<edm::ValueMap<float> > refmassmap;
    iEvent.getByLabel(RefittedMass, refmassmap);

    int kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEMM->begin();cand != CandidatesEEMM->end(); ++ cand ) {
      if (kk > 99) break;
      edm::Ref<edm::View<Candidate> > Ref(CandidatesEEMM,kk);
      cout << "Original 2e2mu mass is= " << cand->p4().mass() << " Refitted mass is= " << (*refmassmap)[Ref] << endl;
      RECO_EEMM_MASS_REFIT[kk]=(*refmassmap)[Ref];
      kk++;
    }


    edm::Handle<reco::VertexCollection> KinematicFitVtx_;
    iEvent.getByLabel(KinematicFitVertex, KinematicFitVtx_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=KinematicFitVtx_->begin(); cand!=KinematicFitVtx_->end(); ++cand){
      if (jjj > 99) break;
      KinFitVertexX[jjj]=cand->position().x();
      KinFitVertexY[jjj]=cand->position().y();
      KinFitVertexZ[jjj]=cand->position().z();
      KinFitVertexChi2r[jjj]=cand->chi2()/cand->ndof();
      KinFitVertexProb[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Kin Fit: " 
	   <<  KinFitVertexX[jjj] << " " 
	   <<  KinFitVertexY[jjj] << " " 
	   <<  KinFitVertexZ[jjj] << " " 
	   <<  KinFitVertexChi2r[jjj] << " " 
	   <<  KinFitVertexProb[jjj] << endl;
      jjj++;
    }
  }

  void fillConstraintVtx4mu(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtx_;
    iEvent.getByLabel(StandardFitVertexMMMM, StandardFitVtx_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtx_->begin(); cand!=StandardFitVtx_->end(); ++cand){
      if (jjj > 99) break;
      StdFitVertexXMMMM[jjj]=cand->position().x();
      StdFitVertexYMMMM[jjj]=cand->position().y();
      StdFitVertexZMMMM[jjj]=cand->position().z();
      StdFitVertexChi2rMMMM[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbMMMM[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit MMMM: " 
	   <<  StdFitVertexXMMMM[jjj] << " " 
	   <<  StdFitVertexYMMMM[jjj] << " " 
	   <<  StdFitVertexZMMMM[jjj] << " " 
	   <<  StdFitVertexChi2rMMMM[jjj] << " " 
	   <<  StdFitVertexProbMMMM[jjj] << endl;

      // Refitted tracks     
      bool hasRefittedTracks = cand->hasRefittedTracks();
      if(hasRefittedTracks){	
	vector<reco::Track> refit_tks= cand->refittedTracks(); 	
	for(unsigned int i=0; i< refit_tks.size(); i++){
	  if (i > 3) break;
	  cout << "Track momentum Refit=" << refit_tks[i].pt() << endl;
	  StdFitVertexTrackMMMM_PT[i][jjj] =refit_tks[i].pt() ;
	  StdFitVertexTrackMMMM_ETA[i][jjj]=refit_tks[i].eta() ;
	  StdFitVertexTrackMMMM_PHI[i][jjj]=refit_tks[i].phi();
	}
      }

      jjj++;

    }


    edm::Handle<edm::View<Candidate> > CandidatesMMMM;
    iEvent.getByToken(RECOcollNameMMMM_, CandidatesMMMM);
    edm::Handle<edm::ValueMap<float> > refmassmap;
    iEvent.getByLabel(RefittedMassMMMM, refmassmap);

    int kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesMMMM->begin();cand != CandidatesMMMM->end(); ++ cand ) {
      if (kk > 99) break;
      edm::Ref<edm::View<Candidate> > Ref(CandidatesMMMM,kk);
      cout << "Original 4mu mass is= " << cand->p4().mass() << " Refitted mass is= " << (*refmassmap)[Ref] << endl;
      RECO_MMMM_MASS_REFIT[kk]=(*refmassmap)[Ref];
      kk++;
    }



    edm::Handle<reco::VertexCollection> KinematicFitVtx_;
    iEvent.getByLabel(KinematicFitVertexMMMM, KinematicFitVtx_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=KinematicFitVtx_->begin(); cand!=KinematicFitVtx_->end(); ++cand){
      if (jjj > 99) break;
      KinFitVertexXMMMM[jjj]=cand->position().x();
      KinFitVertexYMMMM[jjj]=cand->position().y();
      KinFitVertexZMMMM[jjj]=cand->position().z();
      KinFitVertexChi2rMMMM[jjj]=cand->chi2()/cand->ndof();
      KinFitVertexProbMMMM[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Kin Fit MMMM: " 
	   <<  KinFitVertexXMMMM[jjj] << " " 
	   <<  KinFitVertexYMMMM[jjj] << " " 
	   <<  KinFitVertexZMMMM[jjj] << " " 
	   <<  KinFitVertexChi2rMMMM[jjj] << " " 
	   <<  KinFitVertexProbMMMM[jjj] << endl;
      jjj++;
    }
  }

  void fillConstraintVtx4e(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtx_;
    iEvent.getByLabel(StandardFitVertexEEEE, StandardFitVtx_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtx_->begin(); cand!=StandardFitVtx_->end(); ++cand){
      if (jjj > 99) break;
      StdFitVertexXEEEE[jjj]=cand->position().x();
      StdFitVertexYEEEE[jjj]=cand->position().y();
      StdFitVertexZEEEE[jjj]=cand->position().z();
      StdFitVertexChi2rEEEE[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbEEEE[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit EEEE: " 
	   <<  StdFitVertexXEEEE[jjj] << " " 
	   <<  StdFitVertexYEEEE[jjj] << " " 
	   <<  StdFitVertexZEEEE[jjj] << " " 
	   <<  StdFitVertexChi2rEEEE[jjj] << " " 
	   <<  StdFitVertexProbEEEE[jjj] << endl;

      // Refitted tracks     
      bool hasRefittedTracks = cand->hasRefittedTracks();
      if(hasRefittedTracks){	
	vector<reco::Track> refit_tks= cand->refittedTracks(); 	
	for(unsigned int i=0; i< refit_tks.size(); i++){
	  if (i > 3) break;
	  cout << "Track momentum Refit=" << refit_tks[i].pt() << endl;
	  StdFitVertexTrackEEEE_PT[i][jjj] =refit_tks[i].pt() ;
	  StdFitVertexTrackEEEE_ETA[i][jjj]=refit_tks[i].eta() ;
	  StdFitVertexTrackEEEE_PHI[i][jjj]=refit_tks[i].phi();
	}
      }

      jjj++;      

    }



    edm::Handle<edm::View<Candidate> > CandidatesEEEE;
    iEvent.getByToken(RECOcollNameEEEE, CandidatesEEEE);
    edm::Handle<edm::ValueMap<float> > refmassmap;
    iEvent.getByLabel(RefittedMassEEEE, refmassmap);

    int kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEEE->begin();cand != CandidatesEEEE->end(); ++ cand ) {
      if (kk > 99) break;
      edm::Ref<edm::View<Candidate> > Ref(CandidatesEEEE,kk);
      cout << "Original 4e mass is= " << cand->p4().mass() << " Refitted mass is= " << (*refmassmap)[Ref] << endl;
      RECO_EEEE_MASS_REFIT[kk]=(*refmassmap)[Ref];
      kk++;
    }


    edm::Handle<reco::VertexCollection> KinematicFitVtx_;
    iEvent.getByLabel(KinematicFitVertexEEEE, KinematicFitVtx_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=KinematicFitVtx_->begin(); cand!=KinematicFitVtx_->end(); ++cand){
      if (jjj > 99) break;
      KinFitVertexXEEEE[jjj]=cand->position().x();
      KinFitVertexYEEEE[jjj]=cand->position().y();
      KinFitVertexZEEEE[jjj]=cand->position().z();
      KinFitVertexChi2rEEEE[jjj]=cand->chi2()/cand->ndof();
      KinFitVertexProbEEEE[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Kin Fit EEEE: " 
	   <<  KinFitVertexXEEEE[jjj] << " " 
	   <<  KinFitVertexYEEEE[jjj] << " " 
	   <<  KinFitVertexZEEEE[jjj] << " " 
	   <<  KinFitVertexChi2rEEEE[jjj] << " " 
	   <<  KinFitVertexProbEEEE[jjj] << endl;
      jjj++;
    }
  }


  void fillConstraintVtxDiLeptons(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtxDiLep_;
    iEvent.getByLabel(StandardFitVertexDiLep, StandardFitVtxDiLep_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxDiLep_->begin(); cand!=StandardFitVtxDiLep_->end(); ++cand){
      StdFitVertexChi2rDiLep[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbDiLep[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit DiLeptons: " 
	   <<  StdFitVertexChi2rDiLep[jjj] << " " 
	   <<  StdFitVertexProbDiLep[jjj] << endl;
      jjj++;
    }
  }


  void fillConstraintVtxTriLeptons(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtxMMM_;
    iEvent.getByLabel(StandardFitVertexMMM, StandardFitVtxMMM_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxMMM_->begin(); cand!=StandardFitVtxMMM_->end(); ++cand){
      if (jjj > 39) break;
      StdFitVertexChi2rMMM[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbMMM[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit MMM: " 
	   <<  StdFitVertexChi2rMMM[jjj] << " " 
	   <<  StdFitVertexProbMMM[jjj] << endl;
      jjj++;
    }

    edm::Handle<reco::VertexCollection> StandardFitVtxMME_;
    iEvent.getByLabel(StandardFitVertexMME, StandardFitVtxMME_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxMME_->begin(); cand!=StandardFitVtxMME_->end(); ++cand){
      if (jjj > 19) break;
      StdFitVertexChi2rMME[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbMME[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit MME: " 
	   <<  StdFitVertexChi2rMME[jjj] << " " 
	   <<  StdFitVertexProbMME[jjj] << endl;
      jjj++;
    }

    edm::Handle<reco::VertexCollection> StandardFitVtxEEE_;
    iEvent.getByLabel(StandardFitVertexEEE, StandardFitVtxEEE_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxEEE_->begin(); cand!=StandardFitVtxEEE_->end(); ++cand){
      if (jjj > 19) break;
      StdFitVertexChi2rEEE[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbEEE[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit EEE: " 
	   <<  StdFitVertexChi2rEEE[jjj] << " " 
	   <<  StdFitVertexProbEEE[jjj] << endl;
      jjj++;
    }

    edm::Handle<reco::VertexCollection> StandardFitVtxMEE_;
    iEvent.getByLabel(StandardFitVertexMEE, StandardFitVtxMEE_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxMEE_->begin(); cand!=StandardFitVtxMEE_->end(); ++cand){
      if (jjj > 19) break;
      StdFitVertexChi2rMEE[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbMEE[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit MEE: " 
	   <<  StdFitVertexChi2rMEE[jjj] << " " 
	   <<  StdFitVertexProbMEE[jjj] << endl;
      jjj++;
    }

    
  }

  
  		      
  void filljets(const edm::Event& iEvent){
    edm::Handle<pat::JetCollection> pfjets,pfjetsmva;

   iEvent.getByToken(jetsTag_, pfjets);
   iEvent.getByToken(jetsMVATag_, pfjetsmva);
/*    
    if (fillMCTruth == 1){
      cout << "Using Jet energy collection for MC " << endl;

      iEvent.getByToken(jetsTag_, pfjets);
      iEvent.getByToken(jetsMVATag_, pfjetsmva);
      iEvent.getByToken(PuJetMvaMCfullDiscr_,puJetIdMVAMC);
    }
    else{
      iEvent.getByToken(jetsDataTag_, pfjets);
      iEvent.getByToken(jetsMVATag_, pfjetsmva);
      iEvent.getByToken(PuJetMvaDatafullDiscr_,puJetIdMVAData);      
      
    }
*/

    RECO_PFJET_N = pfjets->size();
    cout << "Number of PFJets in the event= " << RECO_PFJET_N << endl;
    
    int index_jets = 0;
    
    for ( pat::JetCollection::const_iterator i=pfjets->begin(); i!=pfjets->end(); i++) {  
      if (index_jets>99) continue;
      
      edm::Ref<pat::JetCollection> pfjetref(pfjets,index_jets);      
      edm::Ref<pat::JetCollection> pfjetrefmva(pfjetsmva,index_jets);
      
      float mva = 0.;
      int pupass = 1;
  
//      mva = i->userFloat("vtxPx");    
      mva = i->userFloat("pileupJetId:fullDiscriminant");
/*
      if (fillMCTruth == 1){
	mva = (*puJetIdMVAMC)[pfjetrefmva];
      }
      else{
	mva = (*puJetIdMVAData)[pfjetrefmva];
      }
*/      
      //      New Selection
      if(i->pt()>20.){
	if(fabs(i->eta())>3.){
	  if(mva<=-0.45) pupass=0;
	}else if(fabs(i->eta())>2.75){
	  if(mva<=-0.55) pupass=0;
	}else if(fabs(i->eta())>2.5){
	  if(mva<=-0.6) pupass=0;
	}else if(mva<=-0.63) pupass=0;
      }else{
	if(fabs(i->eta())>3.){
	  if(mva<=-0.95) pupass=0;
	}else if(fabs(i->eta())>2.75){
	  if(mva<=-0.94) pupass=0;
	}else if(fabs(i->eta())>2.5){
	  if(mva<=-0.96) pupass=0;
	}else if(mva<=-0.95) pupass=0;
      }
        
      
      RECO_PFJET_CHARGE[index_jets] = i->charge();
      RECO_PFJET_ET[index_jets]     = i->et();
      RECO_PFJET_PT[index_jets]     = i->pt();
      RECO_PFJET_ETA[index_jets]    = i->eta();
      RECO_PFJET_PHI[index_jets]    = i->phi();
      RECO_PFJET_PUID[index_jets]     = pupass;
      RECO_PFJET_PUID_MVA[index_jets] = mva;
      
      cout 
	<< "PF Jet with ET= " << RECO_PFJET_ET[index_jets]   
	<< " PT="   << RECO_PFJET_PT[index_jets]   
	<< " ETA="  << RECO_PFJET_ETA[index_jets]  
	<< " PHI="  << RECO_PFJET_PHI[index_jets]  
	<< " PUID=" << RECO_PFJET_PUID[index_jets] 
	<< " PUID_MVA=" << RECO_PFJET_PUID_MVA[index_jets]
       //qier test
        << " Uncorrected Pt=" << i->correctedP4("Uncorrected").Pt()
        << " L1FastJet Pt=" << i->correctedP4("L1FastJet").Pt()
        << " L2Relative Pt=" << i->correctedP4("L2Relative").Pt()
        << " L3Absolute Pt=" << i->correctedP4("L3Absolute").Pt()
	<< endl;
      
      index_jets++;
    } // for loop on PFJets jets
    
    edm::Handle<double> rhoHandle;

    iEvent.getByToken(rhojetsTag_,rhoHandle); 
    if (rhoHandle.isValid() ) {
      RHO_mu=*rhoHandle;
      cout << "RHO mu fastjet= " << RHO_mu << endl; 
    }
    else {
      cout << "Not valid RHO mu collection" << endl;
    }
    
    iEvent.getByToken(rhojetsTag_,rhoHandle); 
    if (rhoHandle.isValid() ) {
      RHO_ele=*rhoHandle;
      cout << "RHO ele fastjet= " << RHO_ele << endl; 
    }
    else {
      cout << "Not valid RHO ele collection" << endl;
    } 
    
  }
  


  void fillBTagging(const edm::Event& iEvent){

    // trackCountingHighEffBJetTags
    edm::Handle<pat::JetCollection> bTagHandle;
    iEvent.getByToken(jetsTag_ ,bTagHandle);
    int l=0;
    for (pat::JetCollection::const_iterator btagIter=bTagHandle->begin(); btagIter!=bTagHandle->end();++btagIter) {
      if(l>=49) continue;
      double discrCSV1 = btagIter->bDiscriminator(tCHighEff_bTag_);
	cout<<" Jet "<< l
	    <<" has b tag discriminator trackCountingHighEffBJetTags = "<< discrCSV1 
	    << " and jet Pt = "<<btagIter->pt()<<endl;      
	tCHighEff_BTagJet_PT[l]=btagIter->pt();
	tCHighEff_BTagJet_ETA[l]=btagIter->eta();
	tCHighEff_BTagJet_PHI[l]=btagIter->phi();
	tCHighEff_BTagJet_DISCR[l]=discrCSV1;
      

    // trackCountingHighPurBJetTags
     double discrCSV2 = btagIter->bDiscriminator(tCHighPur_bTag_);
	cout<<" Jet "<< l
	    <<" has b tag discriminator trackCountingHighPurBJetTags = "<< discrCSV2
	    << " and jet Pt = "<<btagIter->pt()<<endl;      
	tCHighPur_BTagJet_PT[l]=btagIter->pt();
	tCHighPur_BTagJet_ETA[l]=btagIter->eta();
	tCHighPur_BTagJet_PHI[l]=btagIter->phi();
	tCHighPur_BTagJet_DISCR[l]=discrCSV2;

    // combinedSecondaryVertexBJetTags
    double discrCSV3 = btagIter->bDiscriminator(cSV_bTag_);
	cout<<" Jet "<< l
	    <<" has b tag discriminator combinedSecondaryVertexBJetTags = "<< discrCSV3
	    << " and jet Pt = "<<btagIter->pt()<<endl;      
	cSV_BTagJet_PT[l]=btagIter->pt();
	cSV_BTagJet_ETA[l]=btagIter->eta();
	cSV_BTagJet_PHI[l]=btagIter->phi();
	cSV_BTagJet_ET[l]=btagIter->et();
	cSV_BTagJet_DISCR[l]=discrCSV3;
      l++;
    }

  }

  
  
  void fillP3Covariance(const reco::PFCandidate &c, TMatrixDSym &bigCov, int offset) const {
    double dp = PFEnergyResolution().getEnergyResolutionEm(c.energy(), c.eta());
    // In order to produce a 3x3 matrix, we need a jacobian from (p) to (px,py,pz), i.e.
    //            [ Px/P  ]                
    //  C_(3x3) = [ Py/P  ] * sigma^2(P) * [ Px/P Py/P Pz/P  ]
    //            [ Pz/P  ]                
    AlgebraicMatrix31 ptop3;
    ptop3(0,0) = c.px()/c.p();
    ptop3(1,0) = c.py()/c.p();
    ptop3(2,0) = c.pz()/c.p();
    AlgebraicSymMatrix33 mat = ROOT::Math::Similarity(ptop3, AlgebraicSymMatrix11(dp*dp) );
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) {
	bigCov(offset+i,offset+j) = mat(i,j);
      } } 
  }
  
  
  
 private:
 
  std::vector<std::string> module_to_search;
  std::string par_to_search;
  
  // ROOT definition
  TFile *theFile_ ;
  TTree *theTree_ ;
  
  bool useRECOformat;
  std::string inputfileName;

  // Input tags
  std::string decaychannel;
  std::string flaginst;
  std::vector<std::string> flagtags;
  std::string rootFileName;

  // PU
  bool fillPUinfo;
  int num_PU_vertices;
  int PU_BunchCrossing; // bunch crossing for the PU event added
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PileupSrc_;     

  // MC weight
  float MC_weighting;
  edm::EDGetTokenT<GenEventInfoProduct> generator_;

  // HLT
  bool fillHLTinfo;
  edm::InputTag HLTInfoFired;
  std::string HLTAnalysisinst;
  std::vector<edm::InputTag> flagHLTnames; 

  edm::EDGetTokenT<trigger::TriggerEvent> triggerEvent;     

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;


  edm::InputTag triggerMatchObjectEle;
  edm::EDGetTokenT<edm::Association<std::vector<pat::TriggerObjectStandAlone> > > triggerMatchObject;
  std::string triggerFilter,triggerEleFilter,triggerHLTcollection;
  
  
  // SkimEarlyData
  std::string SkimEarlyDataAnalysisinst;
  std::vector<edm::InputTag> flagSkimEarlyDatanames; 
 
  // MC truth
  bool fillMCTruth;
  edm::EDGetTokenT<edm::View<reco::Candidate> > MCcollName;
  // GenParticles;
  edm::EDGetTokenT<vector<reco::GenParticle> > genParticles_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > fourgenleptons_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > digenZ_;

	
  // RECO
  bool useAdditionalRECO;
  bool use2011EA;
  std::vector<edm::InputTag> RECOcollNameBest2e2mu,RECOcollNameBest4e,RECOcollNameBest4mu,
    RECOcollNameBestRestFrame2e2mu,RECOcollNameBestRestFrame4mu,RECOcollNameBestRestFrame4e;
   edm::EDGetTokenT<edm::View<reco::Candidate> > zToMuMu, zToEE,
    zToMuMussmerge, zToEEssmerge, zToCrossLeptons,
    RECOcollNameDiLep_,RECOcollNameMMMM_,RECOcollNameEEEE,RECOcollNameEEMM,
    quadLeptons4Mu, quadLeptons2Mu2E, quadLeptons4E,
    triLeptonsMuMuMu, triLeptonsMuMuE, triLeptonsMuEE, triLeptonsEEE,
    quadLeptons3Mu1E, quadLeptons3E1Mu, quadLeptonsSSOSele,
    quadLeptonsSSOSmu, quadLeptonsSSOSelemu, quadLeptonsSSOSmuele,
    RECOcollNameLLLL;
    std::vector<edm::InputTag> RECOcollNameZ,RECOcollNameZss,
 //   RECOcollNameMMMM,RECOcollNameEEEE,RECOcollNameEEMM,
    RECOcollNameLLLLss,RECOcollNameLLL,RECOcollNameLLLl,RECOcollNameLLLLssos;
  edm::InputTag RECOcollNameDiLep;//RECOcollNameLLLL,RECOcollNameDiLep;
  
  // electron and muon tags
  bool useBestCandidate;
  edm::InputTag BestCandidatesLeptonsTag_;
  edm::InputTag clusterCollectionTag_,gsftrackCollection_;
  edm::EDGetTokenT<edm::View<pat::Muon> > muonPFTag_;
  edm::EDGetTokenT<edm::View<pat::Electron> > electronEgmTag_;
  edm::EDGetTokenT<edm::View<pat::Muon> > muonTag_;

  edm::EDGetTokenT<edm::ValueMap<float> > muonCorrPtErrorMapTag_;
    
  edm::InputTag electronEgmTkMapTag_;
  edm::InputTag electronEgmEcalMapTag_;
  edm::InputTag electronEgmHcalMapTag_;
  edm::EDGetTokenT<edm::View<pat::Electron> > mvaElectronTag_;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaTrigV0MapTag_,mvaNonTrigV0MapTag_;
/*
  edm::EDGetTokenT<edm::ValueMap<double> >
    electronPFIsoValueChargedAllTag_,
    electronPFIsoValueChargedTag_,
    electronPFIsoValueNeutralTag_,
    electronPFIsoValueGammaTag_,
    electronPFIsoValuePUTag_;   


  edm::InputTag eleRegressionEnergyErrorTag_,eleRegressionEnergyTag_;

  edm::EDGetTokenT<edm::ValueMap<double> >
    muonPFIsoValueChargedAllTag_,
    muonPFIsoValueChargedTag_,
    muonPFIsoValueNeutralTag_,
    muonPFIsoValueGammaTag_,
    muonPFIsoValuePUTag_;   


  edm::EDGetTokenT<edm::ValueMap<double> >
    photonPFIsoValueChargedAllTag_,
    photonPFIsoValueChargedTag_,
    photonPFIsoValueNeutralTag_,
    photonPFIsoValueGammaTag_,
    photonPFIsoValuePUTag_;   

*/


  // vertexing 3D 
  edm::InputTag electronTag_Vert,muonTag_Vert;
  edm::EDGetTokenT<edm::ValueMap<float> > electronMapTag_Vert,electronMapTag_VertKF,muonMapTag_Vert,muonMapTag_VertKF,
    electronMapTag_VertValue,electronMapTag_VertValueKF,muonMapTag_VertValue,muonMapTag_VertValueKF,
    electronMapTag_VertError,electronMapTag_VertErrorKF,muonMapTag_VertError,muonMapTag_VertErrorKF;

  edm::InputTag electronMapTag_VertGD,muonMapTag_VertGD;
  edm::InputTag electronMapTag_VertGDEEEE,muonMapTag_VertGDMMMM;
  edm::InputTag electronMapTag_VertStd,muonMapTag_VertStd;
  edm::InputTag electronMapTag_VertStdEEEE,muonMapTag_VertStdMMMM;
  edm::InputTag electronMapTag_VertKin,muonMapTag_VertKin;
  edm::InputTag electronMapTag_VertKinEEEE,muonMapTag_VertKinMMMM;

  // vertexing STIP / SLIP
  edm::EDGetTokenT<edm::ValueMap<float> >  electronSTIPMapTag_Vert,muonSTIPMapTag_Vert,
    electronSLIPMapTag_Vert,muonSLIPMapTag_Vert,
    electronSTIPMapTag_VertValue,muonSTIPMapTag_VertValue,
    electronSLIPMapTag_VertValue,muonSLIPMapTag_VertValue,
    electronSTIPMapTag_VertError,muonSTIPMapTag_VertError,
    electronSLIPMapTag_VertError,muonSLIPMapTag_VertError;
  
  // vertexing GD
  edm::InputTag ftsigma_Vert,ftsigmalag_Vert,
                ftsigma_VertMMMM,ftsigmalag_VertMMMM,
                ftsigma_VertEEEE,ftsigmalag_VertEEEE;
  edm::InputTag 
    gdX_Vert,gdY_Vert,gdZ_Vert,
    gdX_VertMMMM,gdY_VertMMMM,gdZ_VertMMMM,
    gdX_VertEEEE,gdY_VertEEEE,gdZ_VertEEEE;
   edm::InputTag 
    gdlagX_Vert,gdlagY_Vert,gdlagZ_Vert,gdlagProb_Vert,gdlagNdof_Vert,
     gdlagX_VertMMMM,gdlagY_VertMMMM,gdlagZ_VertMMMM,gdlagProb_VertMMMM,gdlagNdof_VertMMMM,
     gdlagX_VertEEEE,gdlagY_VertEEEE,gdlagZ_VertEEEE,gdlagProb_VertEEEE,gdlagNdof_VertEEEE;
  //

   // ConstraintVtx
   edm::InputTag 
     StandardFitVertex,StandardFitVertexMMMM,StandardFitVertexEEEE,
     KinematicFitVertex,KinematicFitVertexMMMM,KinematicFitVertexEEEE,
     RefittedMass,RefittedMassMMMM,RefittedMassEEEE;

   edm::InputTag 
     StandardFitVertexEEE,StandardFitVertexMMM,StandardFitVertexMEE,StandardFitVertexMME,StandardFitVertexDiLep;

  //electronID
  std::vector<edm::InputTag> eleIDTag_;
  edm::InputTag EleID_VeryLooseTag_ ;
  edm::InputTag EleID_LooseTag_  ;
  edm::InputTag EleID_MediumTag_  ;
  edm::InputTag EleID_TightTag_ ;

  edm::InputTag EleID_HZZVeryLooseTag_ ;
  edm::InputTag EleID_HZZLooseTag_  ;
  edm::InputTag EleID_HZZMediumTag_  ;
  edm::InputTag EleID_HZZTightTag_ ;
  
  edm::InputTag allelectronsColl, allmuonsColl;
  edm::InputTag isoVarTagElectronsCal,isoVarTag,isoVarTagElectronsTracker,isoVarTagElectronsX;
  edm::InputTag isoVarTagMuonsHCalIso,isoVarTagMuonsECalIso,
    isoVarTagMuonsCalIso,isoVarTagMuonsTracker,isoVarTagMuonsX;
  
  edm::InputTag theECALIsoDepositLabel;    //EM calorimeter Isolation deposit label
  edm::InputTag theHCALIsoDepositLabel;    //Hadron calorimeter Isolation deposit label
  edm::InputTag theHOCALIsoDepositLabel;   //Outer calorimeter Isolation deposit label
  edm::InputTag theTrackerIsoDepositLabel; //Tracker Isolation deposit label 
  
  // Photon, Tracks, Jets, Vertices
//  edm::EDGetTokenT<vector<reco::Track> > tracksTag_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > pfTag_;
  edm::EDGetTokenT<edm::View<pat::Photon> > photonsTag_;
  vector<reco::Vertex> PV;
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesTag_;
  edm::EDGetTokenT<double> rhojetsTag_;
  edm::EDGetTokenT<vector<pat::Jet> > jetsTag_,jetsMVATag_;
//  edm::EDGetTokenT<edm::ValueMap<float> >  PuJetMvaMCfullDiscr_,PuJetMvaMCfullId_;
//  edm::EDGetTokenT<edm::ValueMap<float> >  PuJetMvaDatafullDiscr_,PuJetMvaDatafullId_;

  // genJet
  edm::EDGetTokenT<vector<reco::GenJet> >  genjetTag_;
  
  // MET
  edm::InputTag trackermetTag_;
  edm::EDGetTokenT<vector<pat::MET> >  pfmetTag_;
  edm::EDGetTokenT<vector<reco::GenMET> > genmetTag_;
  edm::InputTag calometTag_,calometoptTag_,calometoptnohfTag_,calometoptnohfhoTag_;
  edm::InputTag calometopthoTag_,calometnohfTag_,calometnohfhoTag_,calomethoTag_;
  bool useAdditionalMET_;
  edm::InputTag htmetic5Tag_,htmetkt4Tag_,htmetkt6Tag_,htmetsc5Tag_,htmetsc7Tag_;
  edm::InputTag jescormetic5Tag_,jescormetkt4Tag_,jescormetkt6Tag_,jescormetsc5Tag_,jescormetsc7Tag_;
  edm::InputTag cormetMuTag_;

  // Conversion
//  edm::EDGetTokenT<edm::ValueMap<float> > ConvMapDistTag_,ConvMapDcotTag_;

  // Matching
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > > goodElectronMCMatch_;
  edm::EDGetTokenT<reco::CandidateCollection> myElectrons_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > > goodMuonMCMatch_;
  edm::EDGetTokenT<reco::CandidateCollection> myMuons_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > > goodGammaMCMatch_;
  edm::EDGetTokenT<reco::CandidateCollection> myGammas_;

  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodZtoMuMuMCMatch_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodZtoEEMCMatch_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodHiggsTozzToEEMMMCMatch_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodHiggsTozzToMMMMMCMatch_;
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > >goodHiggsTozzToEEEEMCMatch_;
  // Beam Spot
  edm::EDGetTokenT<reco::BeamSpot> offlineBeamSpot_;

  // bTagging
//  edm::EDGetTokenT<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,std::vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference> >  
  std::string tCHighEff_bTag_, tCHighPur_bTag_,cSV_bTag_;

//  const std::vector<std::string> bDiscriminators_;

  // counters
  int irun, ievt,ils,nevt;
  float Avginstlumi;

  // HLT
  int RECO_nMuHLTMatch,RECO_nEleHLTMatch;
  float RECOMU_PT_MuHLTMatch[100],RECOMU_ETA_MuHLTMatch[100],RECOMU_PHI_MuHLTMatch[100];
  float RECOELE_PT_EleHLTMatch[100],RECOELE_ETA_EleHLTMatch[100],RECOELE_PHI_EleHLTMatch[100];
  
  char HLTPathsFired[20000];

 
  // tmp candidate collections
  reco::CandidateCollection *leptonscands2e2mu_;
  reco::CandidateCollection *leptonscands2e2murf_;
  reco::CandidateCollection *leptonscands4mu_;
  reco::CandidateCollection *leptonscands4murf_;
  reco::CandidateCollection *leptonscands4e_;
  reco::CandidateCollection *leptonscands4erf_;
  
  reco::CandidateCollection *leptonscands_Z0;
  reco::CandidateCollection *leptonscands_Z1;
  reco::CandidateCollection *leptonscands_Zss0;
  reco::CandidateCollection *leptonscands_Zss1;
  reco::CandidateCollection *leptonscands_Zcross;
  reco::CandidateCollection *leptonscands_DiLep;
  reco::CandidateCollection *leptonscands_MMMM;
  reco::CandidateCollection *leptonscands_EEEE;
  reco::CandidateCollection *leptonscands_EEMM;
  reco::CandidateCollection *leptonscands_LLL0;
  reco::CandidateCollection *leptonscands_LLL1;
  reco::CandidateCollection *leptonscands_LLL2;
  reco::CandidateCollection *leptonscands_LLL3;
  reco::CandidateCollection *leptonscands_LLLLss0;
  reco::CandidateCollection *leptonscands_LLLLss1;
  reco::CandidateCollection *leptonscands_LLLLss2;
  reco::CandidateCollection *leptonscands_LLLl0;
  reco::CandidateCollection *leptonscands_LLLl1;
  reco::CandidateCollection *leptonscands_LLLL;


  // MC info
  edm::ESHandle<ParticleDataTable>  pdt_;
  
  // MC truth
  float MC_E[7],MC_PT[7],MC_ETA[7],MC_THETA[7],MC_PHI[7],MC_MASS[7],MC_PDGID[7];

  double MCRF_cosTheta1_spin[10], MCRF_cosTheta2_spin[10], MCRF_cosThetaStar_spin[10], MCRF_Phi_spin[10], 
    MCRF_Phi1_spin[10], MCRF_Phi2_spin[10], MCRF_phi1RF_spin[10], MCRF_phi2RF_spin[10], MCRF_MELA[10];
  
  float MC_LEPT_PT[4],MC_LEPT_ETA[4],MC_LEPT_PHI[4],MC_LEPT_THETA[4],MC_LEPT_PDGID[4];
  float MC_Z_MASS[2][5],MC_Z_PT[2][5],MC_Z_ETA[2][5],MC_Z_PHI[2][5],MC_Z_THETA[2][5],MC_Z_PDGID[2][5];

  float MC_fourl_MASS[100][5],MC_fourl_PT[100][5],MC_fourl_PDGID[100][5];
  float MC_ZZ_MASS[4][7],MC_ZZ_PT[4][7],MC_ZZ_ETA[4][7],MC_ZZ_PHI[4][7],MC_ZZ_THETA[4][7],MC_ZZ_PDGID[4][7];

  // RECO collection
 
  
  // RECORF
    
  double RECORF_2e2mu_cosTheta1_spin[100], RECORF_2e2mu_cosTheta2_spin[100], RECORF_2e2mu_cosThetaStar_spin[100], RECORF_2e2mu_Phi_spin[100], 
    RECORF_2e2mu_Phi1_spin[100], RECORF_2e2mu_Phi2_spin[100], RECORF_2e2mu_phi1RF_spin[100], RECORF_2e2mu_phi2RF_spin[100],RECORF_2e2mu_MELA[100];
  
    
  double RECORF_4e_cosTheta1_spin[100], RECORF_4e_cosTheta2_spin[100], RECORF_4e_cosThetaStar_spin[100], RECORF_4e_Phi_spin[100],RECORF_4e_MELA[100], 
    RECORF_4e_Phi1_spin[100], RECORF_4e_Phi2_spin[100], RECORF_4e_phi1RF_spin[100], RECORF_4e_phi2RF_spin[100],RECORF_4mu_MELA[100];
  
    
  double RECORF_4mu_cosTheta1_spin[100], RECORF_4mu_cosTheta2_spin[100], RECORF_4mu_cosThetaStar_spin[100], RECORF_4mu_Phi_spin[100], 
    RECORF_4mu_Phi1_spin[100], RECORF_4mu_Phi2_spin[100], RECORF_4mu_phi1RF_spin[100], RECORF_4mu_phi2RF_spin[100];
  

  int leptonflavor;

  // RECO additional
  float 
    RECO_ZMM_MASS[50],RECO_ZMM_PT[3][50],RECO_ZMM_ETA[3][50],RECO_ZMM_PHI[3][50],
    RECO_ZEE_MASS[50],RECO_ZEE_PT[3][50],RECO_ZEE_ETA[3][50],RECO_ZEE_PHI[3][50],
    RECO_ZMMss_MASS[50],RECO_ZMMss_PT[3][50],RECO_ZMMss_ETA[3][50],RECO_ZMMss_PHI[3][50],
    RECO_ZEEss_MASS[50],RECO_ZEEss_PT[3][50],RECO_ZEEss_ETA[3][50],RECO_ZEEss_PHI[3][50],
    RECO_ZEM_MASS[50],RECO_ZEM_PT[3][50],RECO_ZEM_ETA[3][50],RECO_ZEM_PHI[3][50],
    RECO_DiLep_MASS[50],RECO_DiLep_PT[3][50],RECO_DiLep_ETA[3][50],RECO_DiLep_PHI[3][50];
  

  float 
    RECO_EEMM_MASS[7][100],RECO_MMMM_MASS[7][100],RECO_EEEE_MASS[7][100],
    RECO_EEMM_PT[7][100],  RECO_MMMM_PT[7][100],  RECO_EEEE_PT[7][100],
    RECO_EEMM_ETA[7][100],RECO_MMMM_ETA[7][100],RECO_EEEE_ETA[7][100],
    RECO_EEMM_PHI[7][100],  RECO_MMMM_PHI[7][100],  RECO_EEEE_PHI[7][100],
    RECO_MMMM_MASS_REFIT[100],RECO_EEMM_MASS_REFIT[100],RECO_EEEE_MASS_REFIT[100],
    RECO_LLLL_MASS[7][50],RECO_LLLL_PT[7][50],RECO_LLLL_ETA[7][50],RECO_LLLL_PHI[7][50];

  float	RECO_LLL0_MASS[50],RECO_LLL1_MASS[50],RECO_LLL2_MASS[50],RECO_LLL3_MASS[50];
  float	RECO_LLL0_PT[4][50],RECO_LLL1_PT[4][50],RECO_LLL2_PT[4][50],RECO_LLL3_PT[4][50];

  float	RECO_LLLl0_MASS[20],RECO_LLLl1_MASS[20];
  float	RECO_LLLl0_PT[5][20],RECO_LLLl1_PT[5][20];

  float 
    RECO_LLLL0ss_MASS[20],RECO_LLLL0ss_PT[5][20],
    RECO_LLLL1ss_MASS[20],RECO_LLLL1ss_PT[5][20],
    RECO_LLLL2ss_MASS[20],RECO_LLLL2ss_PT[5][20];
  
  // RECO electrons
  edm::ESHandle<CaloGeometry> theCaloGeom_;  
  float RECOELE_E[100],RECOELE_PT[100],RECOELE_PTError[100],RECOELE_P[100],RECOELE_ETA[100],RECOELE_THETA[100],RECOELE_PHI[100],RECOELE_MASS[100];
  float RECOELE_CHARGE[100];

  bool RECOELE_isEcalDriven[100], RECOELE_isTrackerDriven[100];
  float 
    RECOELE_gsftrack_NPixHits[100],
    RECOELE_gsftrack_NStripHits[100],
    RECOELE_gsftrack_chi2[100], 
    RECOELE_gsftrack_dxyB[100], RECOELE_gsftrack_dxy[100], RECOELE_gsftrack_dxyError[100],
    RECOELE_gsftrack_dzB[100], RECOELE_gsftrack_dz[100],RECOELE_gsftrack_dzError[100];
  int RECOELE_gsftrack_losthits[100],RECOELE_gsftrack_validhits[100],RECOELE_gsftrack_expected_inner_hits[100]; 
  float
    RECOELE_scl_E[100],RECOELE_scl_Et[100],RECOELE_scl_Eta[100],RECOELE_scl_Phi[100]; 
  //float RECOELE_eSuperClusterOverP[100], RECOELE_eSeedClusterOverPout[100], RECOELE_deltaEtaSuperClusterTrackAtVtx[100], RECOELE_deltaPhiSuperClusterTrackAtVtx[100];
  float 
    RECOELE_ep[100], RECOELE_eSeedp[100], RECOELE_eSeedpout[100], RECOELE_eElepout[100],
    RECOELE_deltaEtaIn[100],RECOELE_deltaEtaSeed[100],RECOELE_deltaEtaEle[100],RECOELE_deltaPhiIn[100],
    RECOELE_deltaPhiSeed[100],RECOELE_deltaPhiEle[100],RECOELE_ecalEnergy[100];
  int RECOELE_isbarrel[100], RECOELE_isendcap[100], RECOELE_isEBetaGap[100], RECOELE_isEBphiGap[100], RECOELE_isEEdeeGap[100], RECOELE_isEEringGap[100], RECOELE_isGap[100]; 
  float RECOELE_sigmaIetaIeta[100], RECOELE_sigmaEtaEta[100], RECOELE_e15[100], RECOELE_e25max[100], RECOELE_e55[100], RECOELE_he[100], RECOELE_r9[100];
  float RECOELE_mva[100], RECOELE_fbrem[100],RECOELE_fbrem_mean[100],RECOELE_fbrem_mode[100];
  int RECOELE_nbrems[100], RECOELE_Class[100];
  
  float RECOELE_EGMTRACKISO[100],RECOELE_EGMECALISO[100],RECOELE_EGMHCALISO[100],RECOELE_EGMX[100],
    RECOELE_IP[100],RECOELE_SIP[100],RECOELE_IPERROR[100],
    RECOELE_IP_KF[100],RECOELE_SIP_KF[100],RECOELE_IPERROR_KF[100],
    RECOELE_STIP[100],RECOELE_TIP[100],RECOELE_TIPERROR[100],
    RECOELE_SLIP[100],RECOELE_LIP[100],RECOELE_LIPERROR[100],
    RECOELE_SIP_GD[100], RECOELE_SIP_GDEEEE[100],
    RECOELE_SIP_Std[100], RECOELE_SIP_StdEEEE[100], 
    RECOELE_SIP_Kin[100], RECOELE_SIP_KinEEEE[100];

  double RECOELE_PFchAllPart[100],RECOELE_PFchHad[100],RECOELE_PFneuHad[100],RECOELE_PFphoton[100],
    RECOELE_PFPUchAllPart[100],RECOELE_PFX_dB[100],RECOELE_PFX_rho[100],RECOELE_PF_RingsIsoMVA[100];

  double RECOELE_regEnergy[100],RECOELE_regEnergyError[100];

  float RECOELE_TLE_ParentSC_X[100],RECOELE_TLE_ParentSC_Y[100],RECOELE_TLE_ParentSC_Z[100];
  
  int RECOELE_EEEE_MATCHED[100],RECOELE_EEMM_MATCHED[100],RECOELE_ZEE_MATCHED[100],RECOELE_ZssEE_MATCHED[10],RECOELE_ZEM_MATCHED[100],
    RECOELE_LLL0_MATCHED[100],RECOELE_LLL1_MATCHED[100],RECOELE_LLL2_MATCHED[100],RECOELE_LLL3_MATCHED[100],
    RECOELE_LLLLss0_MATCHED[100],RECOELE_LLLLss1_MATCHED[100],RECOELE_LLLLss2_MATCHED[100],
    RECOELE_LLLl0_MATCHED[100],RECOELE_LLLl1_MATCHED[100],RECOELE_LLLL_MATCHED[100];
  double RECOELE_mvaTrigV0[100],RECOELE_mvaNonTrigV0[100];
  
  double ele_sclRawE[100] ;
  double ele_sclX[100], ele_sclY[100], ele_sclZ[100];
  int ele_seedSubdet1[100];
  double ele_seedDphi1[100], ele_seedDrz1[100];
  int ele_seedSubdet2[100];
  double ele_seedDphi2[100], ele_seedDrz2[100];
  double ele_eidVeryLoose[100], ele_eidLoose[100], ele_eidMedium[100], ele_eidTight[100] ;
  double ele_eidHZZVeryLoose[100], ele_eidHZZLoose[100], ele_eidHZZMedium[100], ele_eidHZZTight[100] ;
  double RECOELE_COV[100][3][3];

  // RECO muons
  bool RECOMU_isPFMu[100],RECOMU_isGlobalMu[100],RECOMU_isStandAloneMu[100],RECOMU_isTrackerMu[100],RECOMU_isCaloMu[100],RECOMU_isTrackerHighPtMu[100],RECOMU_isME0Muon[100];
  float RECOMU_E[100],RECOMU_PT[100],RECOMU_P[100],RECOMU_ETA[100],RECOMU_THETA[100],RECOMU_PHI[100],RECOMU_MASS[100],RECOMU_CHARGE[100];
  double RECOMU_COV[100][3][3];

  float 
    RECOMU_TRACKISO[100],RECOMU_TRACKISO_SUMPT[100],RECOMU_ECALISO[100],RECOMU_HCALISO[100], RECOMU_X[100],   
    RECOMU_IP[100],RECOMU_SIP[100],RECOMU_IPERROR[100],
    RECOMU_IP_KF[100],RECOMU_SIP_KF[100],RECOMU_IPERROR_KF[100],
    RECOMU_STIP[100],RECOMU_TIP[100],RECOMU_TIPERROR[100],
    RECOMU_SLIP[100],RECOMU_LIP[100],RECOMU_LIPERROR[100],
    RECOMU_SIP_GD[100], RECOMU_SIP_GDMMMM[100],
    RECOMU_SIP_Std[100], RECOMU_SIP_StdMMMM[100], 
    RECOMU_SIP_Kin[100], RECOMU_SIP_KinMMMM[100];

  double 
    RECOMU_PFchHad[100],RECOMU_PFneuHad[100],RECOMU_PFphoton[100],RECOMU_PFsumPUPt[100],RECOMU_PFX_dB[100],RECOMU_PFPUchAllPart[100],RECOMU_PFchAllPart[100],RECOMU_PFX_rho[100],
    RECOMU_PFchHad42[100],RECOMU_PFneuHad42[100],RECOMU_PFphoton42[100],RECOMU_PFPUchAllPart42[100],RECOMU_PFchAllPart42[100],
    RECOMU_PF_RingsIsoMVA[100],RECOMU_PF_RingsIDMVA[100];

  float
    RECOMU_caloCompatibility[100],RECOMU_segmentCompatibility[100];
  bool RECOMU_glbmuPromptTight[100];
 
  int RECOMU_MMMM_MATCHED[100],RECOMU_EEMM_MATCHED[100],
      RECOMU_ZMM_MATCHED[100],RECOMU_ZssMM_MATCHED[100],RECOMU_ZEM_MATCHED[100],
      RECOMU_LLL0_MATCHED[100],RECOMU_LLL1_MATCHED[100],RECOMU_LLL2_MATCHED[100],RECOMU_LLL3_MATCHED[100],
    RECOMU_LLLLss0_MATCHED[100],RECOMU_LLLLss1_MATCHED[100],RECOMU_LLLLss2_MATCHED[100],
    RECOMU_LLLl0_MATCHED[100],RECOMU_LLLl1_MATCHED[100],RECOMU_LLLL_MATCHED[100];
  int RECOMU_numberOfMatches[100],RECOMU_numberOfMatchedStations[100];
  int RECOMU_mubesttrkType[100];
  
  float 
    RECOMU_mutrkPT[100],RECOMU_mutrkPTError[100],
    RECOMU_mutrkDxy[100],RECOMU_mutrkDxyError[100],RECOMU_mutrkDxyB[100],
    RECOMU_mutrkDz[100],RECOMU_mutrkDzError[100],RECOMU_mutrkDzB[100],
    RECOMU_mutrkChi2PerNdof[100],
    RECOMU_mutrktrackerLayersWithMeasurement[100],
    RECOMU_muInnertrkDxy[100],RECOMU_muInnertrkDxyError[100],RECOMU_muInnertrkDxyB[100],
    RECOMU_muInnertrkDz[100],RECOMU_muInnertrkDzError[100],RECOMU_muInnertrkDzB[100],
    RECOMU_mubesttrkDxy[100],RECOMU_mubesttrkDxyError[100],RECOMU_mubesttrkDxyB[100],
    RECOMU_mubesttrkDz[100],RECOMU_mubesttrkDzError[100],RECOMU_mubesttrkDzB[100],
    RECOMU_mubesttrkPTError[100],
    RECOMU_muInnertrkChi2PerNdof[100],
    RECOMU_muInnertrktrackerLayersWithMeasurement[100],RECOMU_muInnertrkPT[100],RECOMU_muInnertrkPTError[100],
    RECOMU_muInnertrkCharge[100],RECOMU_muInnertrkNHits[100],RECOMU_muInnertrkNPixHits[100],RECOMU_muInnertrkNStripHits[100],
    RECOMU_mutrkCharge[100],RECOMU_mutrkNHits[100],RECOMU_mutrkNPixHits[100],RECOMU_mutrkNStripHits[100],RECOMU_mutrkNMuonHits[100];
  bool RECOMU_trkmuArbitration[100],RECOMU_trkmu2DCompatibilityLoose[100],RECOMU_trkmu2DCompatibilityTight[100];
  bool RECOMU_trkmuOneStationLoose[100],RECOMU_trkmuOneStationTight[100];
  bool RECOMU_trkmuLastStationLoose[100],RECOMU_trkmuLastStationTight[100];
  bool RECOMU_trkmuLastStationAngLoose[100],RECOMU_trkmuLastStationAngTight[100];
  bool RECOMU_trkmuOneStationAngLoose[100],RECOMU_trkmuOneStationAngTight[100];
  bool RECOMU_trkmuLastStationOptimizedLowPtLoose[100],RECOMU_trkmuLastStationOptimizedLowPtTight[100];
  
   // Photons
  float RECOPHOT_PT[20],RECOPHOT_ETA[20],RECOPHOT_PHI[20],RECOPHOT_THETA[20],RECOPHOT_TLE_ParentSC_X[20],RECOPHOT_TLE_ParentSC_Y[20],RECOPHOT_TLE_ParentSC_Z[20];
  float RECOPFPHOT_PT[20],RECOPFPHOT_PTError[20],RECOPFPHOT_ETA[20],RECOPFPHOT_PHI[20],RECOPFPHOT_THETA[20];
  double RECOPFPHOT_PFchAllPart[20],RECOPFPHOT_PFchHad[20],RECOPFPHOT_PFneuHad[20],RECOPFPHOT_PFphoton[20],
    RECOPFPHOT_PFPUchAllPart[20],RECOPFPHOT_PFX_rho[20];
  

  // Vertexing
  double ftsigma[100],ftsigmalag[100],ftsigmaMMMM[100],ftsigmalagMMMM[100],ftsigmaEEEE[100],ftsigmalagEEEE[100];
  double gdX[100],gdY[100],gdZ[100],gdXMMMM[100],gdYMMMM[100],gdZMMMM[100],gdXEEEE[100],gdYEEEE[100],gdZEEEE[100];
  double gdlagX[100],gdlagY[100],gdlagZ[100],gdlagProb[100],gdlagNdof[100],
    gdlagXMMMM[100],gdlagYMMMM[100],gdlagZMMMM[100],gdlagProbMMMM[100],gdlagNdofMMMM[100],
    gdlagXEEEE[100],gdlagYEEEE[100],gdlagZEEEE[100],gdlagProbEEEE[100],gdlagNdofEEEE[100];

  //ConstraintFit 4l
  double  StdFitVertexX[100], StdFitVertexY[100], StdFitVertexZ[100], StdFitVertexChi2r[100], StdFitVertexProb[100];
  double  KinFitVertexX[100], KinFitVertexY[100], KinFitVertexZ[100], KinFitVertexChi2r[100], KinFitVertexProb[100];
  double  StdFitVertexXMMMM[100], StdFitVertexYMMMM[100], StdFitVertexZMMMM[100], StdFitVertexChi2rMMMM[100], StdFitVertexProbMMMM[100];
  double  KinFitVertexXMMMM[100], KinFitVertexYMMMM[100], KinFitVertexZMMMM[100], KinFitVertexChi2rMMMM[100], KinFitVertexProbMMMM[100];
  double  StdFitVertexXEEEE[100], StdFitVertexYEEEE[100], StdFitVertexZEEEE[100], StdFitVertexChi2rEEEE[100], StdFitVertexProbEEEE[100];
  double  KinFitVertexXEEEE[100], KinFitVertexYEEEE[100], KinFitVertexZEEEE[100], KinFitVertexChi2rEEEE[100], KinFitVertexProbEEEE[100];

  float StdFitVertexTrack_PT[4][100],StdFitVertexTrack_ETA[4][100],StdFitVertexTrack_PHI[4][100],
    StdFitVertexTrackMMMM_PT[4][100],StdFitVertexTrackMMMM_ETA[4][100],StdFitVertexTrackMMMM_PHI[4][100],
    StdFitVertexTrackEEEE_PT[4][100],StdFitVertexTrackEEEE_ETA[4][100],StdFitVertexTrackEEEE_PHI[4][100];

   //ConstraintFit 2l
  double  StdFitVertexChi2rDiLep[40], StdFitVertexProbDiLep[40];

  //ConstraintFit 3l
  double  StdFitVertexChi2rMMM[100], StdFitVertexProbMMM[100];
  double  StdFitVertexChi2rMME[100], StdFitVertexProbMME[100];
  double  StdFitVertexChi2rEEE[100], StdFitVertexProbEEE[100];
  double  StdFitVertexChi2rMEE[100], StdFitVertexProbMEE[100];
   
  
  //Muons Matching
  bool RECOMU_MatchingMCTruth[100];
  float RECOMU_MatchingMCpT[100];
  float RECOMU_MatchingMCEta[100];
  float RECOMU_MatchingMCPhi[100];
  
  //Electrons:
  bool RECOELE_MatchingMCTruth[100];
  float RECOELE_MatchingMCpT[100];
  float RECOELE_MatchingMCEta[100];
  float RECOELE_MatchingMCPhi[100];
  //Gamma:
  bool RECOPHOT_MatchingMCTruth[50];
  float RECOPHOT_MatchingMCpT[50];
  float RECOPHOT_MatchingMCEta[50];
  float RECOPHOT_MatchingMCPhi[50];

  //zToMuMu:
  bool RECOzMuMu_MatchingMCTruth[50];
  float RECOzMuMu_MatchingMCpT[50];
  float RECOzMuMu_MatchingMCmass[50];
  float RECOzMuMu_MatchingMCEta[50];
  float RECOzMuMu_MatchingMCPhi[50];

  //zToEE:
  bool RECOzEE_MatchingMCTruth[50];
  float RECOzEE_MatchingMCpT[50];
  float RECOzEE_MatchingMCmass[50];
  float RECOzEE_MatchingMCEta[50];
  float RECOzEE_MatchingMCPhi[50];

  //HtoZtoEEEE:
  bool RECOHzzEEEE_MatchingMCTruth[100];
  float RECOHzzEEEE_MatchingMCpT[100];
  float RECOHzzEEEE_MatchingMCmass[100];
  float RECOHzzEEEE_MatchingMCEta[100];
  float RECOHzzEEEE_MatchingMCPhi[100];

  //HtoZtoMMMM:
  bool RECOHzzMMMM_MatchingMCTruth[100];
  float RECOHzzMMMM_MatchingMCpT[100];
  float RECOHzzMMMM_MatchingMCmass[100];
  float RECOHzzMMMM_MatchingMCEta[100];
  float RECOHzzMMMM_MatchingMCPhi[100];

  //HtoZtoEEMM:
  bool RECOHzzEEMM_MatchingMCTruth[100];
  float RECOHzzEEMM_MatchingMCpT[100];
  float RECOHzzEEMM_MatchingMCmass[100];
  float RECOHzzEEMM_MatchingMCEta[100];
  float RECOHzzEEMM_MatchingMCPhi[100];
  
 

  // RECO counters
  int RECO_NMU, RECO_NELE, RECO_NTRACK, RECO_NPHOT, RECO_NPFPHOT,RECO_NJET, RECO_NVTX;
  float RECO_TRACK_PT[200], RECO_TRACK_ETA[200], RECO_TRACK_PHI[200],
    RECO_TRACK_CHI2[200],RECO_TRACK_CHI2RED[200],RECO_TRACK_CHI2PROB[200], 
    RECO_TRACK_DXY[200],RECO_TRACK_DXYERR[200], 
    RECO_TRACK_DZ[200],RECO_TRACK_DZERR[200];
  int RECO_TRACK_NHITS[200];
  
  // Primary Vertices
  float RECO_VERTEX_x[15], RECO_VERTEX_y[15], RECO_VERTEX_z[15],RECO_VERTEX_ndof[15],RECO_VERTEX_chi2[15],RECO_VERTEXPROB[15],RECO_VERTEX_TRACK_PT[15][100];
  bool RECO_VERTEX_isValid[15];
  int RECO_VERTEX_ntracks[15];
  
  // RECO JETS
  int RECO_PFJET_N, RECO_PFJET_CHARGE[200],RECO_PFJET_PUID[200];
  float RECO_PFJET_ET[200], RECO_PFJET_PT[200], RECO_PFJET_ETA[200], RECO_PFJET_PHI[200],RECO_PFJET_PUID_MVA[200];
  double RHO,RHO_ele,RHO_mu;

  // GenJET
  float MC_GENJET_PT[100], MC_GENJET_ETA[100], MC_GENJET_PHI[100];

  // RECO MET
  float genmet;
  float calomet;
    //calometopt,calometoptnohf,calometoptnohfho,calometoptho,calometnohf,calometnohfho,calometho;       
  float pfmet,pfmet_x,pfmet_y,pfmet_phi,pfmet_theta;

    //htmetic5,htmetkt4,htmetkt6,htmetsc5,htmetsc7;        
    float tcmet;
    //jescormetic5,jescormetkt4,jescormetkt6,jescormetsc5,jescormetsc7;    
  float cormetmuons;  
  
  // Beam Spot
  double BeamSpot_X,BeamSpot_Y,BeamSpot_Z;	
  BeamSpot bs;
  
 
  
  float tCHighEff_BTagJet_PT[50],
    tCHighPur_BTagJet_PT[50],
    cSV_BTagJet_PT[50];
  float tCHighEff_BTagJet_ETA[50],
    tCHighPur_BTagJet_ETA[50],
    cSV_BTagJet_ETA[50];
  float tCHighEff_BTagJet_PHI[50],
    tCHighPur_BTagJet_PHI[50],
    cSV_BTagJet_PHI[50];
  float tCHighEff_BTagJet_DISCR[50],
    tCHighPur_BTagJet_DISCR[50],
    cSV_BTagJet_DISCR[50];
  float cSV_BTagJet_ET[50];

  float ConvMapDist[100],ConvMapDcot[100];

  // MVA Ring Isolation
  //MuonMVAEstimator *fMuonIsoMVA;

  // Magnetic Field
  edm::ESHandle<MagneticField> magfield_;

};

#endif



