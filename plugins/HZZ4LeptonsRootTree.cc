/** \class  HZZ4LeptonsRootTree
 *
 *  Root Tree for H->ZZ->4l analysis.
 *
 *  Author: N. De Filippis - Politecnico and INFN Bari
 *
 *  Modified (from array to vector format) by Sherif Elgammal / Nicol De Filippis / Giorgia Miniello
 *  11/3/2016    
 *  CMSSW_7_6_X    
 */

#include "HZZ4LeptonsRootTree.h"

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
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

// Trigger
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
#include "DataFormats/Math/interface/deltaR.h"

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
//#include "Muon/MuonAnalysisTools/interface/MuonMVAEstimator.h"
//#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"

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


class MultiTrajectoryStateMode ;
class EgammaTowerIsolation ;


// Class to create TTree variables
#include <TFile.h> 
#include <TTree.h> 

//==================== Maps ================================
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Nicola
// code for accessing registry info
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/Registry.h"

#include "FWCore/Framework/interface/ConstProductRegistry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Provenance/interface/ProcessHistoryRegistry.h"
//nicola


// Chi2 Prob
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

// Run
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "FWCore/Framework/interface/Run.h"

// Maps
#include "DataFormats/Common/interface/ValueMap.h"


// Namespaces
using namespace edm;
using namespace reco;
using namespace std;
using namespace pat;
//====================== end a part for photons ==========================
#include <vector>
#include <set>
#include <stdio.h>
#include "TFile.h"
#include <math.h>
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include <string>

//========================================================================
HZZ4LeptonsRootTree::HZZ4LeptonsRootTree( const edm::ParameterSet& pset )
//========================================================================
{
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
  generator_                = consumes<GenEventInfoProduct >(pset.getParameter<edm::InputTag>("Generator"));

  // Get HLT flags
  fillHLTinfo               = pset.getUntrackedParameter<bool>("fillHLTinfo");
  HLTInfoFired              = pset.getParameter<edm::InputTag>("HLTInfoFired");
  HLTAnalysisinst           = pset.getParameter<string>("HLTAnalysisinst");
  flagHLTnames              = pset.getParameter<vtag>("flagHLTnames");
  // Get HLT matching
  triggerEvent              = consumes<trigger::TriggerEvent >(pset.getParameter<edm::InputTag>("triggerEvent"));
  triggerFilter             = pset.getParameter<std::string>("triggerFilter");
  triggerEleFilter          = pset.getParameter<std::string>("triggerEleFilter");
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
  RECOcollNameDiLep         = pset.getParameter<edm::InputTag>("RECOcollNameDiLep");
  RECOcollNameMMMM          = pset.getParameter<vtag>("RECOcollNameMMMM");
  RECOcollNameEEEE          = pset.getParameter<vtag>("RECOcollNameEEEE");
  RECOcollNameEEMM          = pset.getParameter<vtag>("RECOcollNameEEMM");
  RECOcollNameLLLLss        = pset.getParameter<vtag>("RECOcollNameLLLLss");
  RECOcollNameLLLLssos      = pset.getParameter<vtag>("RECOcollNameLLLLssos");
  RECOcollNameLLL           = pset.getParameter<vtag>("RECOcollNameLLL");
  RECOcollNameLLLl          = pset.getParameter<vtag>("RECOcollNameLLLl");
  RECOcollNameLLLL          = pset.getParameter<edm::InputTag>("RECOcollNameLLLL");

  // electrons and muons tags
  use2011EA                 = pset.getUntrackedParameter<bool>("use2011EA");
  muonTag_                  = consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("MuonsLabel"));
  muonCorrPtErrorMapTag_    = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("MuonsCorrPtErrorMapLabel"));

  muonPFTag_                  = consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("PFMuonsLabel"));
  muonPFIsoValueChargedAllTag_= consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValueChargedAll"));
  muonPFIsoValueChargedTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValueCharged"));
  muonPFIsoValueNeutralTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValueNeutral"));
  muonPFIsoValueGammaTag_     = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValueGamma"));
  muonPFIsoValuePUTag_        = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("MuonPFIsoValuePU"));
      
  clusterCollectionTag_    = pset.getParameter<edm::InputTag>("SuperClustersLabel");
  gsftrackCollection_      = pset.getParameter<edm::InputTag>("GsfTracksElectronsLabel");
    
  electronEgmTag_          = consumes<edm::View<reco::GsfElectron> >(pset.getParameter<edm::InputTag>("ElectronsEgmLabel"));
  electronEgmTkMapTag_     = pset.getParameter<edm::InputTag>("ElectronsEgmTkMapLabel");
  electronEgmEcalMapTag_   = pset.getParameter<edm::InputTag>("ElectronsEgmEcalMapLabel");
  electronEgmHcalMapTag_   = pset.getParameter<edm::InputTag>("ElectronsEgmHcalMapLabel");
    
  mvaElectronTag_          = consumes<edm::View<reco::GsfElectron> >(pset.getParameter<edm::InputTag>("mvaElectronTag"));
  mvaTrigV0MapTag_         = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("mvaTrigV0MapTag"));
  mvaNonTrigV0MapTag_      = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("mvaNonTrigV0MapTag"));
    

  electronPFIsoValueChargedAllTag_= consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueChargedAll"));
  electronPFIsoValueChargedTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueCharged"));
  electronPFIsoValueNeutralTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueNeutral"));
  electronPFIsoValueGammaTag_     = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueGamma"));
  electronPFIsoValuePUTag_        = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValuePU"));
    
  eleRegressionEnergyErrorTag_= pset.getParameter<edm::InputTag>("eleRegressionEnergyErrorLabel");
  eleRegressionEnergyTag_     = pset.getParameter<edm::InputTag>("eleRegressionEnergyLabel");
        

  // PF photons
  pfphotonsTag_                 = consumes<edm::View<reco::PFCandidate>>(pset.getParameter<edm::InputTag>("PFPhotonsLabel"));
  photonPFIsoValueChargedAllTag_= consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValueChargedAll"));
  photonPFIsoValueChargedTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValueCharged"));
  photonPFIsoValueNeutralTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValueNeutral"));
  photonPFIsoValueGammaTag_     = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValueGamma"));
  photonPFIsoValuePUTag_        = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("PhotonPFIsoValuePU"));

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
  eleIDTag_                = pset.getParameter<vtag>("eleIDLabel");

  // Other objets
  photonsTag_              = consumes<edm::View<reco::Photon> >(pset.getParameter<edm::InputTag>("PhotonsLabel"));
  tracksTag_               = consumes<vector<reco::Track> >(pset.getParameter<edm::InputTag>("TracksLabel"));
  jetsTag_                 = consumes<vector<reco::PFJet> >(pset.getParameter<edm::InputTag>("JetsLabel"));
  jetsDataTag_             = consumes<vector<reco::PFJet> >(pset.getParameter<edm::InputTag>("JetsDataLabel"));
  jetsMVATag_              = consumes<vector<reco::PFJet> >(pset.getParameter<edm::InputTag>("JetsMVALabel"));
  PuJetMvaMCfullDiscr_     = consumes<edm::ValueMap<float> > (pset.getParameter<edm::InputTag>("PuJetMvaMCfullDiscrLabel"));
  PuJetMvaMCfullId_        = consumes<edm::ValueMap<float> > (pset.getParameter<edm::InputTag>("PuJetMvaMCfullIdLabel"));
  PuJetMvaDatafullDiscr_   = consumes<edm::ValueMap<float> > (pset.getParameter<edm::InputTag>("PuJetMvaDatafullDiscrLabel"));
  PuJetMvaDatafullId_      = consumes<edm::ValueMap<float> > (pset.getParameter<edm::InputTag>("PuJetMvaDatafullIdLabel"));
  rhojetsTag_              = consumes<double>(pset.getParameter<edm::InputTag>("RhoJetsLabel"));
  verticesTag_             = consumes<std::vector<reco::Vertex> >(pset.getParameter<edm::InputTag>("VerticesLabel"));
  
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
  
  // PF MET
  pfmetTag_               = consumes<vector<reco::PFMET> >(pset.getParameter<edm::InputTag>("PfMETLabel"));
  
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

   // btagging
  tCHighEff_bTag_         = consumes<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,std::vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference> >(pset.getParameter<edm::InputTag>("tCHighEff_bTagLabel"));
  tCHighPur_bTag_         = consumes<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,std::vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference> >(pset.getParameter<edm::InputTag>("tCHighPur_bTagLabel"));
  cSV_bTag_               = consumes<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,std::vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference> >(pset.getParameter<edm::InputTag>("cSV_bTagLabel"));
  
  
  // Conversion finder
  ConvMapDistTag_       = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ConvMapDist"));
  ConvMapDcotTag_       = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ConvMapDcot"));
    
  // Matching
  goodElectronMCMatch_  = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodElectronMCMatch"));
  myElectrons_          = consumes<reco::CandidateCollection>(pset.getParameter<edm::InputTag>("myElectrons"));
  goodMuonMCMatch_      = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodMuonMCMatch"));
  myMuons_              = consumes<reco::CandidateCollection>(pset.getParameter<edm::InputTag>("myMuons"));
  goodGammaMCMatch_     = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodGammaMCMatch"));
  myGammas_             = consumes<reco::CandidateCollection>(pset.getParameter<edm::InputTag>("myGammas"));

  // Beam Spot
  offlineBeamSpot_      = consumes<reco::BeamSpot>(pset.getParameter<edm::InputTag>("offlineBeamSpot"));

// get names from module parameters, then derive slot numbers
  n_                  = 0;
  firstevent_         = true;  
  //produces<vector<std::string> >();

  // PG and FRC 06-07-11
  debug	=	pset.getUntrackedParameter<bool> ("debug", false);
  
  //===============================================================================================
  //Create the root file
  theFile_ = new TFile(rootFileName.c_str(), "RECREATE");
  //theFile_->cd();

  //outputFile_ = pset.getParameter<std::string>("outputFile");
  //rootFile_   = TFile::Open(outputFile_.c_str(),"RECREATE"); // open output file to store histograms

  // Counter of number of analyzed events
  nevt=0;
 
}

 
//========================================================================
HZZ4LeptonsRootTree::~HZZ4LeptonsRootTree()
//========================================================================    
{
  delete rootFile_;
}

//========================================================================
void HZZ4LeptonsRootTree::beginJob() {
  //======================================================================== 
  // go to *OUR* rootfile and book histograms
  //rootFile_->cd();
  //========================================================
  mytree = new TTree("HZZ4LeptonsAnalysis", "HZZ4Leptons Analysis Tree");

  //=============================================================
  //
  //           Create Branchs for Nb of event,run,lumi
  //
  //=============================================================
  /*TotalNbEvents = 0;
  NbEventsPass  = 0;
  Run       = 0;
  Event     = 0;
  lumi      = 0;
  bunch     = 0;
  mytree->Branch("NbEventsPassTrigger",&NbEventsPassTrigger);
  mytree->Branch("NbEventsPassTriggerandPVcond",&NbEventsPassTriggerandPVcond);
  mytree->Branch("event_runNo",  &Run,   "event_runNo/I");
  mytree->Branch("event_evtNo",  &Event, "event_evtNo/I");
  mytree->Branch("event_lumi",   &lumi,  "event_lumi/I");
  mytree->Branch("event_bunch",  &bunch, "event_bunch/I");*/

  // Run event
  mytree->Branch("Run",&irun,"irun/i");
  mytree->Branch("Event",&ievt,"ievt/i");
  mytree->Branch("LumiSection",&ils,"ils/i");
  mytree->Branch("Avginstlumi",&Avginstlumi,"Avginstlumi/F");
  //=============================================================
  //                   
  //           Create Branches for Vertices variables
  //
  //=============================================================
  mytree->Branch( "RECO_NVTX", &RECO_NVTX, "RECO_NVTX/I");
  mytree->Branch( "RECO_VERTEX_x", &RECO_VERTEX_x);
  mytree->Branch( "RECO_VERTEX_y", &RECO_VERTEX_y);
  mytree->Branch( "RECO_VERTEX_z", &RECO_VERTEX_z);
  mytree->Branch( "RECO_VERTEX_ndof", &RECO_VERTEX_ndof);
  mytree->Branch( "RECO_VERTEX_chi2", &RECO_VERTEX_chi2);
  mytree->Branch( "RECO_VERTEX_ntracks", &RECO_VERTEX_ntracks);
  mytree->Branch( "RECO_VERTEXPROB", &RECO_VERTEXPROB);
  mytree->Branch( "RECO_VERTEX_isValid", &RECO_VERTEX_isValid);
  mytree->Branch( "RECO_VERTEX_TRACK_PT", &RECO_VERTEX_TRACK_PT);

  //=============================================================
  //                   
  //           Create Branches for Tracks variables
  //
  //=============================================================
  mytree->Branch( "RECO_NTRACK", &RECO_NTRACK, "RECO_NTRACK/I");
  mytree->Branch( "RECO_TRACK_PT", &RECO_TRACK_PT);
  mytree->Branch( "RECO_TRACK_ETA", &RECO_TRACK_ETA);
  mytree->Branch( "RECO_TRACK_PHI", &RECO_TRACK_PHI);
  mytree->Branch( "RECO_TRACK_CHI2", &RECO_TRACK_CHI2);
  mytree->Branch( "RECO_TRACK_CHI2RED", &RECO_TRACK_CHI2RED);
  mytree->Branch( "RECO_TRACK_CHI2PROB", &RECO_TRACK_CHI2PROB);
  mytree->Branch( "RECO_TRACK_NHITS", &RECO_TRACK_NHITS);
  mytree->Branch( "RECO_TRACK_DXY", &RECO_TRACK_DXY);
  mytree->Branch( "RECO_TRACK_DXYERR", &RECO_TRACK_DXYERR);
  mytree->Branch( "RECO_TRACK_DZ", &RECO_TRACK_DZ);
  mytree->Branch( "RECO_TRACK_DZERR", &RECO_TRACK_DZERR);
  //=============================================================
  //                   
  //           Create Branches for MET reco
  //
  //=============================================================
  mytree->Branch("MC_GENMET", &genmet, "MC_GENMET/F");
  //=============================================================
  //                   
  //           Create Branches for Muons match HLT variables
  //
  //=============================================================
  // HLT 
  mytree->Branch("HLTPathsFired",HLTPathsFired,"HLTPathsFired/C");
  mytree->Branch("RECO_nMuHLTMatch",&RECO_nMuHLTMatch,"RECO_nMuHLTMatch/I");
  mytree->Branch("RECO_nEleHLTMatch",&RECO_nEleHLTMatch,"RECO_nEleHLTMatch/I");
  mytree->Branch("RECOMU_PT_MuHLTMatch",&RECOMU_PT_MuHLTMatch);
  mytree->Branch("RECOMU_ETA_MuHLTMatch",&RECOMU_ETA_MuHLTMatch);
  mytree->Branch("RECOMU_PHI_MuHLTMatch",&RECOMU_PHI_MuHLTMatch);
  mytree->Branch("RECOELE_PT_EleHLTMatch",&RECOELE_PT_EleHLTMatch);
  mytree->Branch("RECOELE_ETA_EleHLTMatch",&RECOELE_ETA_EleHLTMatch);
  mytree->Branch("RECOELE_PHI_EleHLTMatch",&RECOELE_PHI_EleHLTMatch);
  //=============================================================  
  //
  //           Create Branchs for gen particles variables
  //
  //=============================================================
  mytree->Branch("MC_LEPT_PT",&MC_LEPT_PT);
  mytree->Branch("MC_LEPT_ETA",&MC_LEPT_ETA);
  mytree->Branch("MC_LEPT_PHI",&MC_LEPT_PHI);
  mytree->Branch("MC_LEPT_THETA",&MC_LEPT_THETA);
  mytree->Branch("MC_LEPT_PDGID",&MC_LEPT_PDGID);
  mytree->Branch("MC_Z_PT",&MC_Z_PT);
  mytree->Branch("MC_Z_ETA",&MC_Z_ETA);
  mytree->Branch("MC_Z_PHI",&MC_Z_PHI);
  mytree->Branch("MC_Z_THETA",&MC_Z_THETA);
  mytree->Branch("MC_Z_MASS",&MC_Z_MASS);
  mytree->Branch("MC_Z_PDGID",&MC_Z_PDGID);
  mytree->Branch("MC_fourl_MASS",&MC_fourl_MASS);
  mytree->Branch("MC_fourl_PT",&MC_fourl_PT);
  mytree->Branch("MC_fourl_PDGID",&MC_fourl_PDGID);
  mytree->Branch("MC_ZZ_MASS",&MC_ZZ_MASS);
  mytree->Branch("MC_ZZ_PT",&MC_ZZ_PT);
  mytree->Branch("MC_ZZ_ETA",&MC_ZZ_ETA);
  mytree->Branch("MC_ZZ_PHI",&MC_ZZ_PHI);
  mytree->Branch("MC_ZZ_THETA",&MC_ZZ_THETA);
  mytree->Branch("MC_ZZ_PDGID",&MC_ZZ_PDGID);
  mytree->Branch("MC_E",&MC_E);
  mytree->Branch("MC_PT",&MC_PT);
  mytree->Branch("MC_ETA",&MC_ETA);
  mytree->Branch("MC_THETA",&MC_THETA);
  mytree->Branch("MC_PHI",&MC_PHI);
  mytree->Branch("MC_MASS",&MC_MASS);
  mytree->Branch("MC_PDGID",&MC_PDGID);
  //=============================================================
  //
  //           Create Branchs for Muons variables
  //
  //=============================================================
  mytree->Branch("RECO_NMU",&RECO_NMU, "RECO_NMU/I");
  mytree->Branch("RECOMU_isPFMu",&RECOMU_isPFMu);
  mytree->Branch("RECOMU_isGlobalMu",&RECOMU_isGlobalMu);
  mytree->Branch("RECOMU_isStandAloneMu",&RECOMU_isStandAloneMu);
  mytree->Branch("RECOMU_isTrackerMu",&RECOMU_isTrackerMu);
  mytree->Branch("RECOMU_isCaloMu",&RECOMU_isCaloMu);
  // Kinematic of muon
  mytree->Branch("RECOMU_E",&RECOMU_E);
  mytree->Branch("RECOMU_PT",&RECOMU_PT);
  mytree->Branch("RECOMU_P",&RECOMU_P);
  mytree->Branch("RECOMU_ETA",&RECOMU_ETA);
  mytree->Branch("RECOMU_THETA",&RECOMU_THETA);
  mytree->Branch("RECOMU_PHI",&RECOMU_PHI);
  mytree->Branch("RECOMU_MASS",&RECOMU_MASS);
  mytree->Branch("RECOMU_CHARGE",&RECOMU_CHARGE);
  // Covariance matrix
  mytree->Branch("RECOMU_COV",&RECOMU_COV);
  // Isolation
  mytree->Branch("RECOMU_TRACKISO",&RECOMU_TRACKISO);
  mytree->Branch("RECOMU_TRACKISO_SUMPT",&RECOMU_TRACKISO_SUMPT);
  //temporary solution for reducedrechit problem
  mytree->Branch("RECOMU_ECALISO",&RECOMU_ECALISO);
  mytree->Branch("RECOMU_HCALISO",&RECOMU_HCALISO);
  mytree->Branch("RECOMU_X",&RECOMU_X);
  mytree->Branch("RECOMU_PFchHad",&RECOMU_PFchHad);
  mytree->Branch("RECOMU_PFneuHad",&RECOMU_PFneuHad);
  mytree->Branch("RECOMU_PFphoton",&RECOMU_PFphoton);
  mytree->Branch("RECOMU_PFPUchAllPart",&RECOMU_PFPUchAllPart);
  // Vertexing
  mytree->Branch("RECOMU_SIP",&RECOMU_SIP);
  mytree->Branch("RECOMU_IP",&RECOMU_IP);
  mytree->Branch("RECOMU_IPERROR",&RECOMU_IPERROR);
  mytree->Branch("RECOMU_SIP_KF",&RECOMU_SIP_KF);
  mytree->Branch("RECOMU_IP_KF",&RECOMU_IP_KF);
  mytree->Branch("RECOMU_IPERROR_KF",&RECOMU_IPERROR_KF);
  mytree->Branch("RECOMU_PFX_dB",&RECOMU_PFX_dB);
  mytree->Branch("RECOMU_PFX_rho",&RECOMU_PFX_rho);
  mytree->Branch("RECOMU_STIP",&RECOMU_STIP);
  mytree->Branch("RECOMU_SLIP",&RECOMU_SLIP);
  mytree->Branch("RECOMU_TIP",&RECOMU_TIP);
  mytree->Branch("RECOMU_LIP",&RECOMU_LIP);
  mytree->Branch("RECOMU_TIPERROR",&RECOMU_TIPERROR);
  mytree->Branch("RECOMU_LIPERROR",&RECOMU_LIPERROR);
  // Other properties
  mytree->Branch("RECOMU_numberOfMatches",&RECOMU_numberOfMatches);
  mytree->Branch("RECOMU_numberOfMatchedStations",&RECOMU_numberOfMatchedStations);
  mytree->Branch("RECOMU_caloCompatibility",&RECOMU_caloCompatibility);
  mytree->Branch("RECOMU_segmentCompatibility",&RECOMU_segmentCompatibility);
  mytree->Branch("RECOMU_glbmuPromptTight",&RECOMU_glbmuPromptTight);
  // Track properties
  mytree->Branch("RECOMU_mubesttrkType",&RECOMU_mubesttrkType);
  mytree->Branch("RECOMU_mubesttrkDxy",&RECOMU_mubesttrkDxy);
  mytree->Branch("RECOMU_mubesttrkDxyB",&RECOMU_mubesttrkDxyB);
  mytree->Branch("RECOMU_mubesttrkDxyError",&RECOMU_mubesttrkDxyError);
  mytree->Branch("RECOMU_mubesttrkDz",&RECOMU_mubesttrkDz);
  mytree->Branch("RECOMU_mubesttrkDzB",&RECOMU_mubesttrkDzB);
  mytree->Branch("RECOMU_mubesttrkDzError",&RECOMU_mubesttrkDzError);
  mytree->Branch("RECOMU_mubesttrkPTError",&RECOMU_mubesttrkPTError);
  mytree->Branch("RECOMU_mutrkPT",&RECOMU_mutrkPT);
  mytree->Branch("RECOMU_mutrkPTError",&RECOMU_mutrkPTError);
  //mytree->Branch("RECOMU_mutrkPTError",&RECOMU_mutrkPTError);
  mytree->Branch("RECOMU_mutrkDxy",&RECOMU_mutrkDxy);
  mytree->Branch("RECOMU_mutrkDxyError",&RECOMU_mutrkDxyError);
  mytree->Branch("RECOMU_mutrkDxyB",&RECOMU_mutrkDxyB);
  mytree->Branch("RECOMU_mutrkDz",&RECOMU_mutrkDz);
  mytree->Branch("RECOMU_mutrkDzError",&RECOMU_mutrkDzError);
  mytree->Branch("RECOMU_mutrkDzB",&RECOMU_mutrkDzB);
  mytree->Branch("RECOMU_mutrkChi2PerNdof",&RECOMU_mutrkChi2PerNdof);
  mytree->Branch("RECOMU_mutrkCharge",&RECOMU_mutrkCharge);
  mytree->Branch("RECOMU_mutrkNHits",&RECOMU_mutrkNHits);
  mytree->Branch("RECOMU_mutrkNPixHits",&RECOMU_mutrkNPixHits);
  mytree->Branch("RECOMU_mutrkNStripHits",&RECOMU_mutrkNStripHits);
  mytree->Branch("RECOMU_mutrkNMuonHits",&RECOMU_mutrkNMuonHits);
  mytree->Branch("RECOMU_mutrktrackerLayersWithMeasurement",&RECOMU_mutrktrackerLayersWithMeasurement);
  mytree->Branch("RECOMU_muInnertrkDxy",&RECOMU_muInnertrkDxy);
  mytree->Branch("RECOMU_muInnertrkDxyError",&RECOMU_muInnertrkDxyError);
  mytree->Branch("RECOMU_muInnertrkDxyB",&RECOMU_muInnertrkDxyB);
  mytree->Branch("RECOMU_muInnertrkDz",&RECOMU_muInnertrkDz);
  mytree->Branch("RECOMU_muInnertrkDzError",&RECOMU_muInnertrkDzError);
  mytree->Branch("RECOMU_muInnertrkDzB",&RECOMU_muInnertrkDzB);
  mytree->Branch("RECOMU_muInnertrkChi2PerNdof",&RECOMU_muInnertrkChi2PerNdof);
  mytree->Branch("RECOMU_muInnertrktrackerLayersWithMeasurement",&RECOMU_muInnertrktrackerLayersWithMeasurement);
  mytree->Branch("RECOMU_muInnertrkPT",&RECOMU_muInnertrkPT);
  //mytree->Branch("RECOMU_muInnertrkPTError",&RECOMU_muInnertrkPTError);
  mytree->Branch("RECOMU_muInnertrkPTError",&RECOMU_muInnertrkPTError);
  mytree->Branch("RECOMU_muInnertrkCharge",&RECOMU_muInnertrkCharge);
  mytree->Branch("RECOMU_muInnertrkNHits",&RECOMU_muInnertrkNHits);
  mytree->Branch("RECOMU_muInnertrkNPixHits",&RECOMU_muInnertrkNPixHits);
  mytree->Branch("RECOMU_muInnertrkNStripHits",&RECOMU_muInnertrkNStripHits);
  // Tracker muon properties
  // mytree->Branch("RECOMU_trkmuArbitration",&RECOMU_trkmuArbitration);
  mytree->Branch("RECOMU_trkmuArbitration",&RECOMU_trkmuArbitration);
  mytree->Branch("RECOMU_trkmu2DCompatibilityLoose",&RECOMU_trkmu2DCompatibilityLoose);
  mytree->Branch("RECOMU_trkmu2DCompatibilityTight",&RECOMU_trkmu2DCompatibilityTight);
  mytree->Branch("RECOMU_trkmuOneStationLoose",&RECOMU_trkmuOneStationLoose);
  mytree->Branch("RECOMU_trkmuOneStationTight",&RECOMU_trkmuOneStationTight);
  mytree->Branch("RECOMU_trkmuLastStationLoose",&RECOMU_trkmuLastStationLoose);
  mytree->Branch("RECOMU_trkmuLastStationTight",&RECOMU_trkmuLastStationTight);
  mytree->Branch("RECOMU_trkmuOneStationAngLoose",&RECOMU_trkmuOneStationAngLoose);
  mytree->Branch("RECOMU_trkmuOneStationAngTight",&RECOMU_trkmuOneStationAngTight);
  mytree->Branch("RECOMU_trkmuLastStationAngLoose",&RECOMU_trkmuLastStationAngLoose);
  mytree->Branch("RECOMU_trkmuLastStationAngTight",&RECOMU_trkmuLastStationAngTight);
  mytree->Branch("RECOMU_trkmuLastStationOptimizedLowPtLoose",&RECOMU_trkmuLastStationOptimizedLowPtLoose);
  mytree->Branch("RECOMU_trkmuLastStationOptimizedLowPtTight",&RECOMU_trkmuLastStationOptimizedLowPtTight);        
  mytree->Branch("RECOMU_MatchingMCTruth",&RECOMU_MatchingMCTruth); 
  mytree->Branch("RECOMU_MatchingMCpT",&RECOMU_MatchingMCpT); 
  mytree->Branch("RECOMU_MatchingMCEta",&RECOMU_MatchingMCEta); 
  mytree->Branch("RECOMU_MatchingMCPhi",&RECOMU_MatchingMCPhi); 

  //=============================================================
  //
  //  Create Branches for Electron Tree
  //
  //=============================================================
  mytree->Branch("RECO_NELE",&RECO_NELE, "RECO_NELE/I");
  mytree->Branch("RECOELE_isEcalDriven",&RECOELE_isEcalDriven);
  mytree->Branch("RECOELE_isTrackerDriven",&RECOELE_isTrackerDriven);
  mytree->Branch("RECOELE_CHARGE",&RECOELE_CHARGE);
  mytree->Branch("RECOELE_E",&RECOELE_E);
  mytree->Branch("RECOELE_PT",&RECOELE_PT);
  mytree->Branch("RECOELE_P",&RECOELE_P);
  mytree->Branch("RECOELE_ETA",&RECOELE_ETA);
  mytree->Branch("RECOELE_THETA",&RECOELE_THETA);
  mytree->Branch("RECOELE_PHI",&RECOELE_PHI);
  mytree->Branch("RECOELE_MASS",&RECOELE_MASS);
  mytree->Branch("ele_sclRawE",&ele_sclRawE);
  mytree->Branch("RECOELE_scl_E",&RECOELE_scl_E);
  mytree->Branch("RECOELE_scl_Et",&RECOELE_scl_Et);
  mytree->Branch("RECOELE_scl_Eta",&RECOELE_scl_Eta);
  mytree->Branch("RECOELE_scl_Phi",&RECOELE_scl_Phi);
  mytree->Branch("ele_sclX",&ele_sclX);
  mytree->Branch("ele_sclY",&ele_sclY);
  mytree->Branch("ele_sclZ",&ele_sclZ);
  mytree->Branch("RECOELE_PTError",&RECOELE_PTError);
  mytree->Branch("RECOELE_COV",&RECOELE_COV);
  mytree->Branch("RECOELE_EGMECALISO",&RECOELE_EGMECALISO);
  mytree->Branch("RECOELE_EGMHCALISO",&RECOELE_EGMHCALISO);
  mytree->Branch("RECOELE_EGMX",&RECOELE_EGMX);
  // PF isolation
  mytree->Branch("RECOELE_PFchAllPart",&RECOELE_PFchAllPart);
  mytree->Branch("RECOELE_PFchHad",&RECOELE_PFchHad);
  mytree->Branch("RECOELE_PFneuHad",&RECOELE_PFneuHad);
  mytree->Branch("RECOELE_PFphoton",&RECOELE_PFphoton);
  mytree->Branch("RECOELE_PFPUchAllPart",&RECOELE_PFPUchAllPart);
  mytree->Branch("RECOELE_PFX_dB",&RECOELE_PFX_dB);
  mytree->Branch("RECOELE_PFX_rho",&RECOELE_PFX_rho);
  // Vertexing DA
  mytree->Branch("RECOELE_SIP",&RECOELE_SIP);
  mytree->Branch("RECOELE_IP",&RECOELE_IP);
  mytree->Branch("RECOELE_IPERROR",&RECOELE_IPERROR);
  // KF
  mytree->Branch("RECOELE_SIP_KF",&RECOELE_SIP_KF);
  mytree->Branch("RECOELE_IP_KF",&RECOELE_IP_KF);
  mytree->Branch("RECOELE_IPERROR_KF",&RECOELE_IPERROR_KF);
  mytree->Branch("RECOELE_STIP",&RECOELE_STIP);
  mytree->Branch("RECOELE_SLIP",&RECOELE_SLIP);
  mytree->Branch("RECOELE_TIP",&RECOELE_TIP);
  mytree->Branch("RECOELE_LIP",&RECOELE_LIP);
  mytree->Branch("RECOELE_TIPERROR",&RECOELE_TIPERROR);
  mytree->Branch("RECOELE_LIPERROR",&RECOELE_LIPERROR);
  // GsfTrack
  mytree->Branch("RECOELE_gsftrack_NPixHits",&RECOELE_gsftrack_NPixHits);
  mytree->Branch("RECOELE_gsftrack_NStripHits",&RECOELE_gsftrack_NStripHits);
  mytree->Branch("RECOELE_gsftrack_chi2",&RECOELE_gsftrack_chi2);
  mytree->Branch("RECOELE_gsftrack_dxyB",&RECOELE_gsftrack_dxyB);
  mytree->Branch("RECOELE_gsftrack_dxy",&RECOELE_gsftrack_dxy);
  mytree->Branch("RECOELE_gsftrack_dxyError",&RECOELE_gsftrack_dxyError);
  mytree->Branch("RECOELE_gsftrack_dzB",&RECOELE_gsftrack_dzB);
  mytree->Branch("RECOELE_gsftrack_dz",&RECOELE_gsftrack_dz);
  mytree->Branch("RECOELE_gsftrack_dzError",&RECOELE_gsftrack_dzError);
  //Conversion variables
  mytree->Branch("RECOELE_gsftrack_losthits",&RECOELE_gsftrack_losthits);
  mytree->Branch("RECOELE_gsftrack_validhits",&RECOELE_gsftrack_validhits);
  mytree->Branch("RECOELE_gsftrack_expected_inner_hits",&RECOELE_gsftrack_expected_inner_hits);
  // Track-Cluster matching attributes
  mytree->Branch("RECOELE_ep",&RECOELE_ep);
  mytree->Branch("RECOELE_eSeedp",&RECOELE_eSeedp);
  mytree->Branch("RECOELE_eSeedpout",&RECOELE_eSeedpout);
  mytree->Branch("RECOELE_eElepout",&RECOELE_eElepout);
  mytree->Branch("RECOELE_deltaEtaIn",&RECOELE_deltaEtaIn);
  mytree->Branch("RECOELE_deltaEtaSeed",&RECOELE_deltaEtaSeed);
  mytree->Branch("RECOELE_deltaEtaEle",&RECOELE_deltaEtaEle);
  mytree->Branch("RECOELE_deltaPhiIn",&RECOELE_deltaPhiIn);
  mytree->Branch("RECOELE_deltaPhiSeed",&RECOELE_deltaPhiSeed);
  mytree->Branch("RECOELE_deltaPhiEle",&RECOELE_deltaPhiEle);
  // Fiducial flags
  mytree->Branch("RECOELE_isbarrel",&RECOELE_isbarrel);
  mytree->Branch("RECOELE_isendcap",&RECOELE_isendcap);
  mytree->Branch("RECOELE_isGap",&RECOELE_isGap);
  mytree->Branch("RECOELE_isEBetaGap",&RECOELE_isEBetaGap);
  mytree->Branch("RECOELE_isEBphiGap",&RECOELE_isEBphiGap);
  mytree->Branch("RECOELE_isEEdeeGap",&RECOELE_isEEdeeGap);
  mytree->Branch("RECOELE_isEEringGap",&RECOELE_isEEringGap);
  // Shower shape
  mytree->Branch("RECOELE_sigmaIetaIeta",&RECOELE_sigmaIetaIeta);
  mytree->Branch("RECOELE_sigmaEtaEta",&RECOELE_sigmaEtaEta);
  mytree->Branch("RECOELE_e15",&RECOELE_e15);
  mytree->Branch("RECOELE_e25max",&RECOELE_e25max);
  mytree->Branch("RECOELE_e55",&RECOELE_e55);
  mytree->Branch("RECOELE_he",&RECOELE_he);
  //mytree->Branch("RECOELE_r9",&RECOELE_r9);
  // Brem & Classifaction
  mytree->Branch("RECOELE_fbrem",&RECOELE_fbrem);
  mytree->Branch("RECOELE_nbrems",&RECOELE_nbrems);
  mytree->Branch("RECOELE_Class",&RECOELE_Class);
  mytree->Branch("RECOELE_fbrem_mean",&RECOELE_fbrem_mean);
  mytree->Branch("RECOELE_fbrem_mode",&RECOELE_fbrem_mode);
  // Corrections
  mytree->Branch("RECOELE_ecalEnergy",&RECOELE_ecalEnergy);
  // Seed Collection
  mytree->Branch("ele_seedSubdet2",&ele_seedSubdet2);
  mytree->Branch("ele_seedDphi2",&ele_seedDphi2);
  mytree->Branch("ele_seedDrz2",&ele_seedDrz2);
  mytree->Branch("ele_seedSubdet1",&ele_seedSubdet1);
  mytree->Branch("ele_seedDphi1",&ele_seedDphi1);
  mytree->Branch("ele_seedDrz1",&ele_seedDrz1);
  mytree->Branch("RECOELE_mvaTrigV0",&RECOELE_mvaTrigV0);
  mytree->Branch("RECOELE_mvaNonTrigV0",&RECOELE_mvaNonTrigV0);
  // Conversion Finder
  mytree->Branch("ConvMapDist",&ConvMapDist);
  mytree->Branch("ConvMapDcot",&ConvMapDcot);
  // Matching
  mytree->Branch("RECOELE_MatchingMCTruth",&RECOELE_MatchingMCTruth);
  mytree->Branch("RECOELE_MatchingMCpT",&RECOELE_MatchingMCpT);
  mytree->Branch("RECOELE_MatchingMCEta",&RECOELE_MatchingMCEta);
  mytree->Branch("RECOELE_MatchingMCPhi",&RECOELE_MatchingMCPhi);
  
  //=============================================================
  //
  //  Create vectors for Photons Tree
  //
  //=============================================================
  mytree->Branch("RECO_NPHOT",&RECO_NPHOT,"RECO_NPHOT/I");
  mytree->Branch("RECOPHOT_PT",&RECOPHOT_PT);
  mytree->Branch("RECOPHOT_ETA",&RECOPHOT_ETA);
  mytree->Branch("RECOPHOT_PHI",&RECOPHOT_PHI);
  mytree->Branch("RECOPHOT_THETA",&RECOPHOT_THETA);
  mytree->Branch("RECO_NPFPHOT",&RECO_NPFPHOT,"RECO_NPFPHOT/I");
  mytree->Branch("RECOPFPHOT_PT",&RECOPFPHOT_PT);
  mytree->Branch("RECOPFPHOT_PTError",&RECOPFPHOT_PTError);
  mytree->Branch("RECOPFPHOT_ETA",&RECOPFPHOT_ETA);
  mytree->Branch("RECOPFPHOT_PHI",&RECOPFPHOT_PHI);
  mytree->Branch("RECOPFPHOT_THETA",&RECOPFPHOT_THETA);
  mytree->Branch("RECOPFPHOT_PFchAllPart",&RECOPFPHOT_PFchAllPart);
  mytree->Branch("RECOPFPHOT_PFchHad",&RECOPFPHOT_PFchHad);
  mytree->Branch("RECOPFPHOT_PFneuHad",&RECOPFPHOT_PFneuHad);
  mytree->Branch("RECOPFPHOT_PFphoton",&RECOPFPHOT_PFphoton);
  mytree->Branch("RECOPFPHOT_PFPUchAllPart",&RECOPFPHOT_PFPUchAllPart);
  mytree->Branch("RECOPFPHOT_PFX_rho",&RECOPFPHOT_PFX_rho);
  mytree->Branch("RECOPHOT_MatchingMCTruth",&RECOPHOT_MatchingMCTruth);
  mytree->Branch("RECOPHOT_MatchingMCpT",&RECOPHOT_MatchingMCpT);
  mytree->Branch("RECOPHOT_MatchingMCEta",&RECOPHOT_MatchingMCEta);
  mytree->Branch("RECOPHOT_MatchingMCPhi",&RECOPHOT_MatchingMCPhi);

  //=============================================================
  //                   
  //           Create Branches for PF MET
  //
  //=============================================================
  mytree->Branch("PFMet_pt",&PFMet_pt,"PFMet_pt/F");
  mytree->Branch("PFMet_eta",&PFMet_eta,"PFMet_eta/F");
  mytree->Branch("PFMet_phi",&PFMet_phi,"PFMet_phi/F");
  mytree->Branch("PFMet_en",&PFMet_en,"PFMet_en/F");
  mytree->Branch("PFMet_px",&PFMet_px,"PFMet_px/F");
  mytree->Branch("PFMet_py",&PFMet_py,"PFMet_py/F");
  mytree->Branch("PFMet_pz",&PFMet_pz,"PFMet_pz/F");
  mytree->Branch("PFMet_sumEt",&PFMet_sumEt,"PFMet_sumEt/F");
  mytree->Branch("METSign",&METSign,"METSign/D");

  //=============================================================
  //
  //  Create vectors for BTagging Tree
  //
  //=============================================================
  mytree->Branch("tCHighEff_nb",&tCHighEff_nb);
  mytree->Branch("tCHighEff_BTagJet_PT",&tCHighEff_BTagJet_PT);
  mytree->Branch("tCHighEff_BTagJet_ETA",&tCHighEff_BTagJet_ETA);
  mytree->Branch("tCHighEff_BTagJet_PHI",&tCHighEff_BTagJet_PHI);
  mytree->Branch("tCHighEff_BTagJet_DISCR",&tCHighEff_BTagJet_DISCR);
  mytree->Branch("tCHighPur_nb",&tCHighPur_nb);
  mytree->Branch("tCHighPur_BTagJet_PT",&tCHighPur_BTagJet_PT);
  mytree->Branch("tCHighPur_BTagJet_ETA",&tCHighPur_BTagJet_ETA);
  mytree->Branch("tCHighPur_BTagJet_PHI",&tCHighPur_BTagJet_PHI);
  mytree->Branch("tCHighPur_BTagJet_DISCR",&tCHighPur_BTagJet_DISCR);
  mytree->Branch("cSV_nb",&cSV_nb);
  mytree->Branch("cSV_BTagJet_PT",&cSV_BTagJet_PT);
  mytree->Branch("cSV_BTagJet_ETA",&cSV_BTagJet_ETA);
  mytree->Branch("cSV_BTagJet_PHI",&cSV_BTagJet_PHI);
  mytree->Branch("cSV_BTagJet_DISCR",&cSV_BTagJet_DISCR);
  mytree->Branch("cSV_BTagJet_ET",&cSV_BTagJet_ET);
  //=============================================================
  //
  //  Create Branches for reco leptons Tree
  //
  //=============================================================
  // RECO additional block for reconstructed higgs, Z and their daughters
  mytree->Branch("RECO_ZMM_MASS",&RECO_ZMM_MASS);
  mytree->Branch("RECO_ZEE_MASS",&RECO_ZEE_MASS);
  mytree->Branch("RECO_DiLep_MASS",&RECO_DiLep_MASS);
  mytree->Branch("RECO_ZMM_PT",&RECO_ZMM_PT);
  mytree->Branch("RECO_ZEE_PT",&RECO_ZEE_PT);
  mytree->Branch("RECO_DiLep_PT",&RECO_DiLep_PT);
  mytree->Branch("RECO_ZMM_ETA",&RECO_ZMM_ETA);
  mytree->Branch("RECO_ZEE_ETA",&RECO_ZEE_ETA);
  mytree->Branch("RECO_DiLep_ETA",&RECO_DiLep_ETA);
  mytree->Branch("RECO_ZMM_PHI",&RECO_ZMM_PHI);
  mytree->Branch("RECO_ZEE_PHI",&RECO_ZEE_PHI);
  mytree->Branch("RECO_DiLep_PHI",&RECO_DiLep_PHI);                              
  mytree->Branch("RECO_ZMMss_MASS",&RECO_ZMMss_MASS);
  mytree->Branch("RECO_ZEEss_MASS",&RECO_ZEEss_MASS);
  mytree->Branch("RECO_ZEM_MASS",&RECO_ZEM_MASS);
  mytree->Branch("RECO_ZMMss_PT",&RECO_ZMMss_PT);
  mytree->Branch("RECO_ZEEss_PT",&RECO_ZEEss_PT);
  mytree->Branch("RECO_ZEM_PT",&RECO_ZEM_PT);
  mytree->Branch("RECO_ZMMss_ETA",&RECO_ZMMss_ETA);
  mytree->Branch("RECO_ZEEss_ETA",&RECO_ZEEss_ETA);
  mytree->Branch("RECO_ZEM_ETA",&RECO_ZEM_ETA);
  mytree->Branch("RECO_ZMMss_PHI",&RECO_ZMMss_PHI);
  mytree->Branch("RECO_ZEEss_PHI",&RECO_ZEEss_PHI);
  mytree->Branch("RECO_ZEM_PHI",&RECO_ZEM_PHI);    
  mytree->Branch("RECO_MMMM_MASS",&RECO_MMMM_MASS);
  mytree->Branch("RECO_MMMM_PT",&RECO_MMMM_PT);
  mytree->Branch("RECO_MMMM_ETA",&RECO_MMMM_ETA);
  mytree->Branch("RECO_MMMM_PHI",&RECO_MMMM_PHI);
  mytree->Branch("RECO_MMMM_MASS_REFIT",&RECO_MMMM_MASS_REFIT);
  mytree->Branch("RECO_EEEE_MASS",&RECO_EEEE_MASS);
  mytree->Branch("RECO_EEEE_PT",&RECO_EEEE_PT);
  mytree->Branch("RECO_EEEE_ETA",&RECO_EEEE_ETA);
  mytree->Branch("RECO_EEEE_PHI",&RECO_EEEE_PHI);
  mytree->Branch("RECO_EEEE_MASS_REFIT",&RECO_EEEE_MASS_REFIT);
  mytree->Branch("RECO_EEMM_MASS",&RECO_EEMM_MASS);
  mytree->Branch("RECO_EEMM_PT",&RECO_EEMM_PT);
  mytree->Branch("RECO_EEMM_ETA",&RECO_EEMM_ETA);
  mytree->Branch("RECO_EEMM_PHI",&RECO_EEMM_PHI);
  mytree->Branch("RECO_EEMM_MASS_REFIT",&RECO_EEMM_MASS_REFIT);
  mytree->Branch("RECO_LLL0_MASS",&RECO_LLL0_MASS);
  mytree->Branch("RECO_LLL1_MASS",&RECO_LLL1_MASS);
  mytree->Branch("RECO_LLL2_MASS",&RECO_LLL2_MASS);
  mytree->Branch("RECO_LLL3_MASS",&RECO_LLL3_MASS);
  mytree->Branch("RECO_LLL0_PT",&RECO_LLL0_PT);
  mytree->Branch("RECO_LLL1_PT",&RECO_LLL1_PT);
  mytree->Branch("RECO_LLL2_PT",&RECO_LLL2_PT);
  mytree->Branch("RECO_LLL3_PT",&RECO_LLL3_PT);
  mytree->Branch("RECO_LLLl0_MASS",&RECO_LLLl0_MASS);
  mytree->Branch("RECO_LLLl1_MASS",&RECO_LLLl1_MASS);
  mytree->Branch("RECO_LLLl0_PT",&RECO_LLLl0_PT);
  mytree->Branch("RECO_LLLl1_PT",&RECO_LLLl1_PT);
  mytree->Branch("RECO_LLLL0ss_MASS",&RECO_LLLL0ss_MASS);
  mytree->Branch("RECO_LLLL0ss_PT",&RECO_LLLL0ss_PT);
  mytree->Branch("RECO_LLLL1ss_MASS",&RECO_LLLL1ss_MASS);
  mytree->Branch("RECO_LLLL1ss_PT",&RECO_LLLL1ss_PT);
  mytree->Branch("RECO_LLLL2ss_MASS",&RECO_LLLL2ss_MASS);
  mytree->Branch("RECO_LLLL2ss_PT",&RECO_LLLL2ss_PT);
  //mytree->Branch("RECOcollNameLLLLssos_MASS",&RECOcollNameLLLLssos_MASS);
  //mytree->Branch("RECOcollNameLLLLssos_PT",&RECOcollNameLLLLssos_PT);
  mytree->Branch("RECO_LLLL_MASS",&RECO_LLLL_MASS);
  mytree->Branch("RECO_LLLL_PT",&RECO_LLLL_PT);
  mytree->Branch("RECO_LLLL_ETA",&RECO_LLLL_ETA);
  mytree->Branch("RECO_LLLL_PHI",&RECO_LLLL_PHI);
  // Matching ZtoMuMu
  mytree->Branch("RECOzMuMu_MatchingMCTruth",&RECOzMuMu_MatchingMCTruth);
  mytree->Branch("RECOzMuMu_MatchingMCpT",&RECOzMuMu_MatchingMCpT);
  mytree->Branch("RECOzMuMu_MatchingMCmass",&RECOzMuMu_MatchingMCmass);
  mytree->Branch("RECOzMuMu_MatchingMCEta",&RECOzMuMu_MatchingMCEta);
  mytree->Branch("RECOzMuMu_MatchingMCPhi",&RECOzMuMu_MatchingMCPhi);
  // Matching ZtoEE
  mytree->Branch("RECOzEE_MatchingMCTruth",&RECOzEE_MatchingMCTruth);
  mytree->Branch("RECOzEE_MatchingMCpT",&RECOzEE_MatchingMCpT);
  mytree->Branch("RECOzEE_MatchingMCmass",&RECOzEE_MatchingMCmass);
  mytree->Branch("RECOzEE_MatchingMCEta",&RECOzEE_MatchingMCEta);
  mytree->Branch("RECOzEE_MatchingMCPhi",&RECOzEE_MatchingMCPhi);
  // Matching HiggsToMMMM
  mytree->Branch("RECOHzzMMMM_MatchingMCTruth",&RECOHzzMMMM_MatchingMCTruth);
  mytree->Branch("RECOHzzMMMM_MatchingMCpT",&RECOHzzMMMM_MatchingMCpT);
  mytree->Branch("RECOHzzMMMM_MatchingMCmass",&RECOHzzMMMM_MatchingMCmass);
  mytree->Branch("RECOHzzMMMM_MatchingMCEta",&RECOHzzMMMM_MatchingMCEta);
  mytree->Branch("RECOHzzMMMM_MatchingMCPhi",&RECOHzzMMMM_MatchingMCPhi);
  // Matching HiggsToEEEE
  mytree->Branch("RECOHzzEEEE_MatchingMCTruth",&RECOHzzEEEE_MatchingMCTruth);
  mytree->Branch("RECOHzzEEEE_MatchingMCpT",&RECOHzzEEEE_MatchingMCpT);
  mytree->Branch("RECOHzzEEEE_MatchingMCmass",&RECOHzzEEEE_MatchingMCmass);
  mytree->Branch("RECOHzzEEEE_MatchingMCEta",&RECOHzzEEEE_MatchingMCEta);
  mytree->Branch("RECOHzzEEEE_MatchingMCPhi",&RECOHzzEEEE_MatchingMCPhi);
  // Matching HiggsToEEMM
  mytree->Branch("RECOHzzEEMM_MatchingMCTruth",&RECOHzzEEMM_MatchingMCTruth);
  mytree->Branch("RECOHzzEEMM_MatchingMCpT",&RECOHzzEEMM_MatchingMCpT);
  mytree->Branch("RECOHzzEEMM_MatchingMCmass",&RECOHzzEEMM_MatchingMCmass);
  mytree->Branch("RECOHzzEEMM_MatchingMCEta",&RECOHzzEEMM_MatchingMCEta);
  mytree->Branch("RECOHzzEEMM_MatchingMCPhi",&RECOHzzEEMM_MatchingMCPhi);
  //=============================================================
  //                   
  //           Create Branchs for PileUp tree
  //
  //=============================================================
  mytree->Branch("num_PU_vertices",&num_PU_vertices,"num_PU_vertices/I");
  mytree->Branch("PU_BunchCrossing",&PU_BunchCrossing,"PU_BunchCrossing/I");

  //=============================================================
  //                   
  //           Create Branch for Rho
  //
  //=============================================================
  mytree->Branch("Rho2",&Rho2);
  //=============================================================
  //                   
  //           Create Branch for events reweighting
  //
  //=============================================================
  mytree->Branch("MC_weighting",&MC_weighting);
  //=============================================================
  //                   
  //           Create Branchs for ConstraintVtx2e2mu tree
  //
  //=============================================================
  // ConstraintFit 4l
  mytree->Branch("StdFitVertexX",        &StdFitVertexX);
  mytree->Branch("StdFitVertexY",        &StdFitVertexY);
  mytree->Branch("StdFitVertexZ",        &StdFitVertexZ);
  mytree->Branch("StdFitVertexChi2r",    &StdFitVertexChi2r);
  mytree->Branch("StdFitVertexProb",     &StdFitVertexProb);
  mytree->Branch("StdFitVertexTrack_PT", &StdFitVertexTrack_PT);
  mytree->Branch("StdFitVertexTrack_ETA",&StdFitVertexTrack_ETA);
  mytree->Branch("StdFitVertexTrack_PHI",&StdFitVertexTrack_PHI);
  mytree->Branch("KinFitVertexX",        &KinFitVertexX);
  mytree->Branch("KinFitVertexY",        &KinFitVertexY);
  mytree->Branch("KinFitVertexZ",        &KinFitVertexZ);
  mytree->Branch("KinFitVertexChi2r",    &KinFitVertexChi2r);
  mytree->Branch("KinFitVertexProb",     &KinFitVertexProb);
  //Create Branchs for fillConstraintVtx4mu tree
  mytree->Branch("StdFitVertexXMMMM",        &StdFitVertexXMMMM);
  mytree->Branch("StdFitVertexYMMMM",        &StdFitVertexYMMMM);
  mytree->Branch("StdFitVertexZMMMM",        &StdFitVertexZMMMM);
  mytree->Branch("StdFitVertexChi2rMMMM",    &StdFitVertexChi2rMMMM);
  mytree->Branch("StdFitVertexProbMMMM",     &StdFitVertexProbMMMM);
  mytree->Branch("StdFitVertexTrackMMMM_PT", &StdFitVertexTrackMMMM_PT);
  mytree->Branch("StdFitVertexTrackMMMM_ETA",&StdFitVertexTrackMMMM_ETA);
  mytree->Branch("StdFitVertexTrackMMMM_PHI",&StdFitVertexTrackMMMM_PHI);
  mytree->Branch("KinFitVertexXMMMM",        &KinFitVertexXMMMM);
  mytree->Branch("KinFitVertexYMMMM",        &KinFitVertexYMMMM);
  mytree->Branch("KinFitVertexZMMMM",        &KinFitVertexZMMMM);
  mytree->Branch("KinFitVertexChi2rMMMM",    &KinFitVertexChi2rMMMM);
  mytree->Branch("KinFitVertexProbMMMM",     &KinFitVertexProbMMMM);
  //Create Branchs for fillConstraintVtx4e tree
  mytree->Branch("StdFitVertexXEEEE",        &StdFitVertexXEEEE);
  mytree->Branch("StdFitVertexYEEEE",        &StdFitVertexYEEEE);
  mytree->Branch("StdFitVertexZEEEE",        &StdFitVertexZEEEE);
  mytree->Branch("StdFitVertexChi2rEEEE",    &StdFitVertexChi2rEEEE);
  mytree->Branch("StdFitVertexProbEEEE",     &StdFitVertexProbEEEE);
  mytree->Branch("StdFitVertexTrackEEEE_PT", &StdFitVertexTrackEEEE_PT);
  mytree->Branch("StdFitVertexTrackEEEE_ETA",&StdFitVertexTrackEEEE_ETA);
  mytree->Branch("StdFitVertexTrackEEEE_PHI",&StdFitVertexTrackEEEE_PHI);
  mytree->Branch("KinFitVertexXEEEE",        &KinFitVertexXEEEE);
  mytree->Branch("KinFitVertexYEEEE",        &KinFitVertexYEEEE);
  mytree->Branch("KinFitVertexZEEEE",        &KinFitVertexZEEEE);
  mytree->Branch("KinFitVertexChi2rEEEE",    &KinFitVertexChi2rEEEE);
  mytree->Branch("KinFitVertexProbEEEE",     &KinFitVertexProbEEEE);
  // constrintFit Dileptons
  mytree->Branch("StdFitVertexChi2rDiLep",   &StdFitVertexChi2rDiLep);
  mytree->Branch("StdFitVertexProbDiLep",    &StdFitVertexProbDiLep);
  // constrintFit 3l
  mytree->Branch("StdFitVertexChi2rMMM",     &StdFitVertexChi2rMMM);
  mytree->Branch("StdFitVertexProbMMM",      &StdFitVertexProbMMM);
  mytree->Branch("StdFitVertexChi2rMME",     &StdFitVertexChi2rMME);
  mytree->Branch("StdFitVertexProbMME",      &StdFitVertexProbMME);
  mytree->Branch("StdFitVertexChi2rEEE",     &StdFitVertexChi2rEEE);
  mytree->Branch("StdFitVertexProbEEE",      &StdFitVertexProbEEE);
  mytree->Branch("StdFitVertexChi2rMEE",     &StdFitVertexChi2rMEE);
  mytree->Branch("StdFitVertexProbMEE",      &StdFitVertexProbMEE);
  //=============================================================
  //                   
  //           Create Branchs for Geom. Discri. tree
  //
  //=============================================================
  mytree->Branch("ftsigma",        &ftsigma);
  mytree->Branch("gdX",            &gdX);
  mytree->Branch("gdY",            &gdY);
  mytree->Branch("gdZ",            &gdZ);
  mytree->Branch("ftsigmalag",     &ftsigmalag);
  mytree->Branch("gdlagX",         &gdlagX);
  mytree->Branch("gdlagY",         &gdlagY);
  mytree->Branch("gdlagZ",         &gdlagZ);
  mytree->Branch("gdlagProb",      &gdlagProb);
  mytree->Branch("gdlagNdof",      &gdlagNdof);
  mytree->Branch("ftsigmaMMMM",    &ftsigmaMMMM);
  mytree->Branch("gdXMMMM",        &gdXMMMM);
  mytree->Branch("gdYMMMM",        &gdYMMMM);
  mytree->Branch("gdZMMMM",        &gdZMMMM);
  mytree->Branch("ftsigmalagMMMM", &ftsigmalagMMMM);
  mytree->Branch("gdlagXMMMM",     &gdlagXMMMM);
  mytree->Branch("gdlagYMMMM",     &gdlagYMMMM);
  mytree->Branch("gdlagZMMMM",     &gdlagZMMMM);
  mytree->Branch("gdlagProbMMMM",  &gdlagProbMMMM);
  mytree->Branch("gdlagNdofMMMM",  &gdlagNdofMMMM);
  mytree->Branch("ftsigmaEEEE",    &ftsigmaEEEE);
  mytree->Branch("gdXEEEE",        &gdXEEEE);
  mytree->Branch("gdYEEEE",        &gdYEEEE);
  mytree->Branch("gdZEEEE",        &gdZEEEE);
  mytree->Branch("ftsigmalagEEEE", &ftsigmalagEEEE);
  mytree->Branch("gdlagXEEEE",     &gdlagXEEEE);
  mytree->Branch("gdlagYEEEE",     &gdlagYEEEE);
  mytree->Branch("gdlagZEEEE",     &gdlagZEEEE);
  mytree->Branch("gdlagProbEEEE",  &gdlagProbEEEE);
  mytree->Branch("gdlagNdofEEEE",  &gdlagNdofEEEE);
  //=============================================================
  //                   
  //           Create Branchs for  RECORF block 2e2mu tree
  //
  //=============================================================
  mytree->Branch("RECORF_2e2mu_cosTheta1_spin",&RECORF_2e2mu_cosTheta1_spin);
  mytree->Branch("RECORF_2e2mu_cosTheta2_spin",&RECORF_2e2mu_cosTheta2_spin);
  mytree->Branch("RECORF_2e2mu_cosThetaStar_spin",&RECORF_2e2mu_cosThetaStar_spin);
  mytree->Branch("RECORF_2e2mu_Phi_spin",&RECORF_2e2mu_Phi_spin);
  mytree->Branch("RECORF_2e2mu_Phi1_spin",&RECORF_2e2mu_Phi1_spin);
  mytree->Branch("RECORF_2e2mu_Phi2_spin",&RECORF_2e2mu_Phi2_spin);
  mytree->Branch("RECORF_2e2mu_phi1RF_spin",&RECORF_2e2mu_phi1RF_spin);
  mytree->Branch("RECORF_2e2mu_phi2RF_spin",&RECORF_2e2mu_phi2RF_spin);
  mytree->Branch("RECORF_2e2mu_MELA",&RECORF_2e2mu_MELA);
  
  mytree->Branch("RECORF_4e_cosTheta1_spin",&RECORF_4e_cosTheta1_spin);
  mytree->Branch("RECORF_4e_cosTheta2_spin",&RECORF_4e_cosTheta2_spin);
  mytree->Branch("RECORF_4e_cosThetaStar_spin",&RECORF_4e_cosThetaStar_spin);
  mytree->Branch("RECORF_4e_Phi_spin",&RECORF_4e_Phi_spin);
  mytree->Branch("RECORF_4e_Phi1_spin",&RECORF_4e_Phi1_spin);
  mytree->Branch("RECORF_4e_Phi2_spin",&RECORF_4e_Phi2_spin);
  mytree->Branch("RECORF_4e_phi1RF_spin",&RECORF_4e_phi1RF_spin);
  mytree->Branch("RECORF_4e_phi2RF_spin",&RECORF_4e_phi2RF_spin);
  mytree->Branch("RECORF_4e_MELA",&RECORF_4e_MELA);
      
  mytree->Branch("RECORF_4mu_cosTheta1_spin",&RECORF_4mu_cosTheta1_spin);
  mytree->Branch("RECORF_4mu_cosTheta2_spin",&RECORF_4mu_cosTheta2_spin);
  mytree->Branch("RECORF_4mu_cosThetaStar_spin",&RECORF_4mu_cosThetaStar_spin);
  mytree->Branch("RECORF_4mu_Phi_spin",&RECORF_4mu_Phi_spin);
  mytree->Branch("RECORF_4mu_Phi1_spin",&RECORF_4mu_Phi1_spin);
  mytree->Branch("RECORF_4mu_Phi2_spin",&RECORF_4mu_Phi2_spin);
  mytree->Branch("RECORF_4mu_phi1RF_spin",&RECORF_4mu_phi1RF_spin);
  mytree->Branch("RECORF_4mu_phi2RF_spin",&RECORF_4mu_phi2RF_spin);
  mytree->Branch("RECORF_4mu_MELA",&RECORF_4mu_MELA);

  // MCCP variables
  mytree->Branch("MCRF_cosTheta1_spin",&MCRF_cosTheta1_spin);
  mytree->Branch("MCRF_cosTheta2_spin",&MCRF_cosTheta2_spin);
  mytree->Branch("MCRF_cosThetaStar_spin",&MCRF_cosThetaStar_spin);
  mytree->Branch("MCRF_Phi_spin",&MCRF_Phi_spin);
  mytree->Branch("MCRF_Phi1_spin",&MCRF_Phi1_spin);
  mytree->Branch("MCRF_Phi2_spin",&MCRF_Phi2_spin);
  mytree->Branch("MCRF_phi1RF_spin",&MCRF_phi1RF_spin);
  mytree->Branch("MCRF_phi2RF_spin",&MCRF_phi2RF_spin);
  mytree->Branch("MCRF_MELA",&MCRF_MELA);
  //=============================================================
  //                   
  //           Create Branchs for Beam Spot position tree
  //
  //=============================================================
  mytree->Branch("BeamSpot_X",  &BeamSpot_X, "BeamSpot_X/D");
  mytree->Branch("BeamSpot_Y",  &BeamSpot_Y, "BeamSpot_Y/D");
  mytree->Branch("BeamSpot_Z",  &BeamSpot_Z, "BeamSpot_Z/D");
  //=============================================================
  //
  //           Create Branchs for Pimary Vertice variables
  //
  //=============================================================
  /*mytree->Branch("nbPv",&nbPv);
  mytree->Branch("Nbdof",&Nbdof);
  mytree->Branch("PositionRho",&PositionRho);
  mytree->Branch("PositionX",&PositionX);
  mytree->Branch("PositionY",&PositionY);
  mytree->Branch("PositionZ",&PositionZ);*/
  //=============================================================
  //                   
  //           Create Branches for jets variables
  //=============================================================
  // PFJets
  mytree->Branch( "RECO_PFJET_N",   &RECO_PFJET_N,"RECO_PFJET_N/I");
  mytree->Branch( "RECO_PFJET_CHARGE",  &RECO_PFJET_CHARGE);
  mytree->Branch( "RECO_PFJET_ET",  &RECO_PFJET_ET);
  mytree->Branch( "RECO_PFJET_PT",  &RECO_PFJET_PT);
  mytree->Branch( "RECO_PFJET_ETA", &RECO_PFJET_ETA);
  mytree->Branch( "RECO_PFJET_PHI", &RECO_PFJET_PHI);
  mytree->Branch( "RECO_PFJET_PUID", &RECO_PFJET_PUID);
  mytree->Branch( "RECO_PFJET_PUID_MVA", &RECO_PFJET_PUID_MVA);
  mytree->Branch( "RHO_ele", &RHO_ele, "RHO_ele/D");
  mytree->Branch( "RHO_mu", &RHO_mu, "RHO_mu/D");
}

//========================================================================================
void HZZ4LeptonsRootTree::analyze( const edm::Event& evt, const edm::EventSetup& es ) {
//======================================================================================
  using namespace edm; // needed for all fwk related classe
  //==============================================
  //=        Begin of the main program           =
  //============================================== 
  // Incrementing counter of events
  nevt++;
  cout << "Dumping the information of the event in a ROOT tree: " << nevt << std::endl;
  IntialValues();
  // Dump Run, Event, LumiSection
  irun=evt.id().run();
  ievt=evt.id().event();
  ils=evt.luminosityBlock();
  cout << "Dumping the information of run=" << irun << "  event=" << ievt << "  lumisection=" << ils << std::endl;
  edm::Handle<LumiSummary> l;
  evt.getLuminosityBlock().getByLabel("lumiProducer", l); 
  // Check that there is something
  if (l.isValid()){
    Avginstlumi=l->avgInsDelLumi();
    cout << "Instataneous luminosity= " << Avginstlumi << endl;
  }

  //GIORGIA fillTracks(evt);
  
  if (fillPUinfo) fillPU(evt);
  EventsMCReWeighting(evt);

  // fill HLT block
  //fillHLTFired(evt);

  //trigger matching
  triggermatching(evt);

  // Get the MC Truth particles, H, ZZ and 4 leptons
  if ( fillMCTruth) {
    cout << "Filling MCtruth variables" << endl;
    fillgenparticles(evt,es);
    fillmc(evt);
  }

  // Vertices
  edm::Handle<reco::VertexCollection> recoPrimaryVertexCollection;
  evt.getByToken(verticesTag_,recoPrimaryVertexCollection);
  PV = *recoPrimaryVertexCollection;
  fillVertices(evt);
  
  
  // Fill RECO block in the rootple
  // PF Jets
  filljets(evt);
  
  if (useAdditionalRECO==true) {
    fillAdditionalRECO(evt);
  }
  
  
  //GENMET
  if (fillMCTruth) {
    edm::Handle<reco::GenMETCollection> genmetHandle;
    evt.getByToken(genmetTag_,genmetHandle);
    for ( GenMETCollection::const_iterator i=genmetHandle->begin(); i!=genmetHandle->end(); i++) {
      genmet = i->pt();
    }
  }

  
  // Filling electron and muons vectors
  fillMuons(evt,es);
  fillElectrons(evt,es);
  fillPhotons(evt);  

  //btagging
  fillBTagging(evt);

  //RECO MET
  fillMET(evt);

  //Beam Spot
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  evt.getByToken(offlineBeamSpot_,recoBeamSpotHandle) ;
  bs = *recoBeamSpotHandle;
  BeamSpot_X=bs.position().x();
  BeamSpot_Y=bs.position().y();
  BeamSpot_Z=bs.position().z();

  cout << "BeamSpot:"
    << "  bs_X=" << BeamSpot_X
    << "  bs_Y=" << BeamSpot_Y
    << "  bs_Z=" << BeamSpot_Z
    << endl;

  // Geometrical discriminant 
  // if (useAdditionalRECO==false){
  //  if (decaychannel=="2e2mu" ) fillGD2e2mu(evt);
  //  if (decaychannel=="4mu" )   fillGD4mu(evt);
  //  if (decaychannel=="4e" )    fillGD4e(evt);
  // }
  // else if (useAdditionalRECO==true) {
  //  fillGD2e2mu(evt);
  //  fillGD4mu(evt);
  //  fillGD4e(evt);
  // }

  // ConstraintVertexFit
  //if (useAdditionalRECO==false){
  //  if (decaychannel=="2e2mu" ) fillConstraintVtx2e2mu(evt);
  //  if (decaychannel=="4mu" )   fillConstraintVtx4mu(evt);
  //  if (decaychannel=="4e" )    fillConstraintVtx4e(evt);
  //}
  //else if (useAdditionalRECO==true) {
  //  fillConstraintVtx2e2mu(evt);
  //  fillConstraintVtx4mu(evt);
  //  fillConstraintVtx4e(evt);
  //}
  // fillConstraintVtxTriLeptons(evt);
  // fillConstraintVtxDiLeptons(evt);
  //fillRho(evt); 
  //==============================================
  //=        End of the main program             =
  //============================================== 
  mytree->Fill();
}

//========================================================================
void HZZ4LeptonsRootTree::endJob() {
  //========================================================================
  // go to *OUR* root file and store histograms
  //rootFile_->cd();
  //mytree->Write();
  //rootFile_->Close();

    theFile_->cd();
    mytree->Write();
    theFile_->Close();

}

//=============================================================
//
//         Method for finding good Primary Vertex
//
//=============================================================

//=============================================================
//
//            Method for Photons Tree
//
//=============================================================
void HZZ4LeptonsRootTree::fillPhotons(const edm::Event& iEvent){
  RECO_NPHOT=0;
  RECOPHOT_PT.clear();
  RECOPHOT_ETA.clear();
  RECOPHOT_PHI.clear();
  RECOPHOT_THETA.clear();
  RECO_NPFPHOT=0;
  RECOPFPHOT_PT.clear();
  RECOPFPHOT_PTError.clear();
  RECOPFPHOT_ETA.clear();
  RECOPFPHOT_PHI.clear();
  RECOPFPHOT_THETA.clear();
  RECOPFPHOT_PFchAllPart.clear();
  RECOPFPHOT_PFchHad.clear();
  RECOPFPHOT_PFneuHad.clear();
  RECOPFPHOT_PFphoton.clear();
  RECOPFPHOT_PFPUchAllPart.clear();
  RECOPFPHOT_PFX_rho.clear();
  RECOPHOT_MatchingMCTruth.clear();
  RECOPHOT_MatchingMCpT.clear();
  RECOPHOT_MatchingMCEta.clear();
  RECOPHOT_MatchingMCPhi.clear();
 
  // Photons
  edm::Handle<edm::View<reco::Photon> > photons;
  iEvent.getByToken(photonsTag_, photons);
  RECO_NPHOT=photons->size();
  edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchPhot;
  iEvent.getByToken(goodGammaMCMatch_, GenParticlesMatchPhot);
  edm::Handle<reco::CandidateCollection > CollPhot;
  iEvent.getByToken(myGammas_, CollPhot);
  bool ismyGammas=false;
  if (CollPhot.isValid()) ismyGammas=true;
  int iphot=0;
  
  for (edm::View<reco::Photon>::const_iterator cand = photons->begin(); cand != photons->end(); ++cand) {
    if (iphot>19) break;
    //RECO_NPHOT.push_back(iphot);
    RECOPHOT_PT.push_back(cand->pt());
    RECOPHOT_ETA.push_back(cand->eta());
    RECOPHOT_PHI.push_back(cand->phi());
    RECOPHOT_THETA.push_back(cand->theta());
    // Matching
    int i=0;
    if (ismyGammas){
      for ( reco::CandidateCollection::const_iterator hIter=CollPhot->begin(); hIter!= CollPhot->end(); ++hIter ){
	//cout << "Reco Photon with pT= " << hIter->pt() << " and mass="<< hIter->mass()<< endl;
	if (fabs(hIter->pt()-cand->pt())<0.01){
	  i=hIter-(CollPhot->begin());
	  CandidateRef Ref( CollPhot, i );
	  edm::Ref<std::vector<reco::GenParticle> > genrefPhot = (*GenParticlesMatchPhot)[Ref];
	  if (!genrefPhot.isNull()){
	    cout << "GenMuon with pT= " << genrefPhot->p4().pt() << " and mass="<< genrefPhot->p4().mass()<< endl;
	    RECOPHOT_MatchingMCTruth.push_back(true);
	    RECOPHOT_MatchingMCpT.push_back(genrefPhot->p4().pt());
	    RECOPHOT_MatchingMCEta.push_back(genrefPhot->p4().eta());
	    RECOPHOT_MatchingMCPhi.push_back(genrefPhot->p4().phi());	    
	  } 
	}   
      }
    }   
    iphot++;
  }
  
  // PF ISR photon
  
  // Photons
  edm::Handle<edm::View<reco::PFCandidate> > pfphotons;
  iEvent.getByToken(pfphotonsTag_, pfphotons);
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
  iphot=0;
  RECO_NPFPHOT=pfphotons->size();
  cout << "There are " << pfphotons->size() << " photons" << endl;
  for (edm::View<reco::PFCandidate>::const_iterator cand = pfphotons->begin(); cand != pfphotons->end(); ++cand) {
    edm::Ref<edm::View<reco::PFCandidate> > phtrackref(pfphotons,iphot); 
    RECOPFPHOT_PT.push_back(cand->pt());
    // error on pt of photon taken from #include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"
    double C;
    double S;
    double N;
    if(TMath::Abs(cand->eta())<1.48){C=0.35/100; S=5.51/100; N=98./1000.;}
    else{C=0; S=12.8/100; N=440./1000.;} 
    double  perr = TMath::Sqrt(C*C*cand->p4().e()*cand->p4().e() + S*S*cand->p4().e() + N*N);
    double pterr = perr*cand->pt()/cand->p(); 
    RECOPFPHOT_PTError.push_back(float(pterr));
    cout << "Photon : pT= " << RECOPFPHOT_PT[iphot] << " pTerr= " << RECOPFPHOT_PTError[iphot]<< endl;
    
    RECOPFPHOT_ETA.push_back(cand->eta());
    RECOPFPHOT_PHI.push_back(cand->phi());
    RECOPFPHOT_THETA.push_back(cand->theta());

    RECOPFPHOT_PFchAllPart.push_back((*isoPFChargedAllphmap)[phtrackref]); 
    RECOPFPHOT_PFchHad.push_back((*isoPFChargedphmap)[phtrackref]); 
    RECOPFPHOT_PFneuHad.push_back((*isoPFNeutralphmap)[phtrackref]); 
    RECOPFPHOT_PFphoton.push_back((*isoPFGammaphmap)[phtrackref]); 
    RECOPFPHOT_PFPUchAllPart.push_back((*isoPFPUphmap)[phtrackref]); 
    
    
    RECOPFPHOT_PFX_rho.push_back(( 
				  (*isoPFChargedphmap)[phtrackref]+ 
				  (*isoPFNeutralphmap)[phtrackref]+ 
                                  (*isoPFGammaphmap)[phtrackref]+ 
				  (*isoPFPUphmap)[phtrackref]
				   )/double(cand->p4().pt())); 
    iphot++;
  }
}


//=============================================================
//
//            Method for BTagging Tree
//
//=============================================================
void HZZ4LeptonsRootTree::fillBTagging(const edm::Event& iEvent)
{
  tCHighEff_nb.clear();
  tCHighEff_BTagJet_PT.clear();
  tCHighEff_BTagJet_ETA.clear();
  tCHighEff_BTagJet_PHI.clear();
  tCHighEff_BTagJet_DISCR.clear();
  tCHighPur_nb.clear();
  tCHighPur_BTagJet_PT.clear();
  tCHighPur_BTagJet_ETA.clear();
  tCHighPur_BTagJet_PHI.clear();
  tCHighPur_BTagJet_DISCR.clear();
  cSV_nb.clear();
  cSV_BTagJet_PT.clear();
  cSV_BTagJet_ETA.clear();
  cSV_BTagJet_PHI.clear();
  cSV_BTagJet_DISCR.clear();
  cSV_BTagJet_ET.clear();
  // trackCountingHighEffBJetTags
  edm::Handle<reco::JetTagCollection> bTagHandle;
  iEvent.getByToken(tCHighEff_bTag_,bTagHandle);
  for (reco::JetTagCollection::const_iterator btagIter=bTagHandle->begin(); btagIter!=bTagHandle->end();++btagIter) {
    iBtag++;
    if (btagIter->second>-99999. && btagIter->second<99999.){
      /*cout<<" Jet "<< iBtag
	  <<" has b tag discriminator trackCountingHighEffBJetTags = "<<btagIter->second
	  << " and jet Pt = "<<btagIter->first->pt()<<endl;*/     
      tCHighEff_nb.push_back(iBtag);
      tCHighEff_BTagJet_PT.push_back(btagIter->first->pt());
      tCHighEff_BTagJet_ETA.push_back(btagIter->first->eta());
      tCHighEff_BTagJet_PHI.push_back(btagIter->first->phi());
      tCHighEff_BTagJet_DISCR.push_back(btagIter->second);
    }
  }
  // trackCountingHighPurBJetTags
  edm::Handle<reco::JetTagCollection> bTagHandle_b;
  iEvent.getByToken(tCHighPur_bTag_,bTagHandle_b);
  for (reco::JetTagCollection::const_iterator btagIter=bTagHandle_b->begin(); btagIter!=bTagHandle_b->end();++btagIter) {
    iBtag++;
    if (btagIter->second>-99999. && btagIter->second<99999.){
      /*cout<<" Jet "<< iBtag
	  <<" has b tag discriminator trackCountingHighPurBJetTags = "<<btagIter->second
	  << " and jet Pt = "<<btagIter->first->pt()<<endl;*/      
      tCHighPur_nb.push_back(iBtag);
      tCHighPur_BTagJet_PT.push_back(btagIter->first->pt());
      tCHighPur_BTagJet_ETA.push_back(btagIter->first->eta());
      tCHighPur_BTagJet_PHI.push_back(btagIter->first->phi());
      tCHighPur_BTagJet_DISCR.push_back(btagIter->second);
    }
  }
  // combinedSecondaryVertexBJetTags
  edm::Handle<reco::JetTagCollection> bTagHandle_g;
  iEvent.getByToken(cSV_bTag_,bTagHandle_g);
  for (reco::JetTagCollection::const_iterator btagIter=bTagHandle_g->begin(); btagIter!=bTagHandle_g->end();++btagIter) {
    iBtag++;
    if (btagIter->second>-99999. && btagIter->second<99999.){
      /*cout<<" Jet "<< iBtag
	  <<" has b tag discriminator combinedSecondaryVertexBJetTags = "<<btagIter->second
	  << " and jet Pt = "<<btagIter->first->pt()<<endl;*/      
      cSV_nb.push_back(iBtag);
      cSV_BTagJet_PT.push_back(btagIter->first->pt());
      cSV_BTagJet_ETA.push_back(btagIter->first->eta());
      cSV_BTagJet_PHI.push_back(btagIter->first->phi());
      cSV_BTagJet_DISCR.push_back(btagIter->second);
      cSV_BTagJet_ET.push_back(btagIter->first->et());
    }
  }
}
  
//=============================================================
//
//            Method for Jets Tree
//
//=============================================================
void HZZ4LeptonsRootTree::filljets(const edm::Event& iEvent)
{
  RECO_PFJET_N=0;
  RECO_PFJET_CHARGE.clear();
  RECO_PFJET_ET.clear();
  RECO_PFJET_PT.clear();
  RECO_PFJET_ETA.clear();
  RECO_PFJET_PHI.clear();
  RECO_PFJET_PUID.clear();
  RECO_PFJET_PUID_MVA.clear();
  edm::Handle<reco::PFJetCollection> pfjets,pfjetsmva;
  edm::Handle<edm::ValueMap<float> > puJetIdMVAMC;
  edm::Handle<edm::ValueMap<float> > puJetIdMVAData;
 
  if (fillMCTruth == 1){
    iEvent.getByToken(jetsTag_, pfjets);
    iEvent.getByToken(jetsMVATag_, pfjetsmva);
    iEvent.getByToken(PuJetMvaMCfullDiscr_,puJetIdMVAMC);
  }
  else{
    iEvent.getByToken(jetsDataTag_, pfjets);
    iEvent.getByToken(jetsMVATag_, pfjetsmva);
    iEvent.getByToken(PuJetMvaDatafullDiscr_,puJetIdMVAData);      
    
  }
  RECO_PFJET_N=pfjets->size();
  cout << "Number of PFJets in the event= " << RECO_PFJET_N << endl;

  int index_jets = 0;

  for ( PFJetCollection::const_iterator i=pfjets->begin(); i!=pfjets->end(); i++) {
    
    edm::Ref<reco::PFJetCollection> pfjetref(pfjets,index_jets);
    edm::Ref<reco::PFJetCollection> pfjetrefmva(pfjetsmva,index_jets);
	    
    float mva = 0.;
    int pupass = 1;
    
    if (fillMCTruth == 1){
      mva = (*puJetIdMVAMC)[pfjetrefmva];
    }
    else{
      mva = (*puJetIdMVAData)[pfjetrefmva];
    }
      
    //      New Selection
    if(i->pt()>20.){
      if(i->eta()>3.){
	if(mva<=-0.45)pupass=0;
      }else if(i->eta()>2.75){
	if(mva<=-0.55)pupass=0;
      }else if(i->eta()>2.5){
	if(mva<=-0.6)pupass=0;
      }else if(mva<=-0.63)pupass=0;
    }else{
      if(i->eta()>3.){
	if(mva<=-0.95)pupass=0;
      }else if(i->eta()>2.75){
	if(mva<=-0.94)pupass=0;
      }else if(i->eta()>2.5){
	if(mva<=-0.96)pupass=0;
      }else if(mva<=-0.95)pupass=0;
    }
   
    
    RECO_PFJET_CHARGE.push_back(i->charge());
    RECO_PFJET_ET.push_back(i->et());
    RECO_PFJET_PT.push_back(i->pt());
    RECO_PFJET_ETA.push_back(i->eta());
    RECO_PFJET_PHI.push_back(i->phi());
    RECO_PFJET_PUID.push_back(pupass);
    RECO_PFJET_PUID_MVA.push_back(mva);
    
    cout 
      << "PF Jet with ET= " << RECO_PFJET_ET[index_jets]   
      << " PT="   << RECO_PFJET_PT[index_jets]   
      << " ETA="  << RECO_PFJET_ETA[index_jets]  
      << " PHI="  << RECO_PFJET_PHI[index_jets]  
      << " PUID=" << RECO_PFJET_PUID[index_jets] 
      << " PUID_MVA=" << RECO_PFJET_PUID_MVA[index_jets] 
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
    cout << "Not valid RHO collection" << endl;
  }

  iEvent.getByToken(rhojetsTag_,rhoHandle); 
  if (rhoHandle.isValid() ) {
    RHO_ele=*rhoHandle;
    cout << "RHO ele fastjet= " << RHO_ele << endl; 
  }
  else {
    cout << "Not valid RHO collection" << endl;
  }
}
  


float HZZ4LeptonsRootTree::delR(float eta1,float phi1,float eta2,float phi2){
  float mpi=3.14;
  float dp=std::abs(phi1-phi2);
  if (dp>mpi) dp-=float(2*mpi);
  return sqrt((eta1-eta2)*(eta1-eta2) + dp*dp);
}



//=============================================================
//
//                Method for PF MET Tree
//
//=============================================================
void HZZ4LeptonsRootTree::fillMET(const edm::Event& evt)
{
  edm::Handle<reco::PFMETCollection> pfmetcoll;
  evt.getByToken(pfmetTag_, pfmetcoll);
  if(!pfmetcoll.isValid()) return;
  const PFMETCollection *pfmetcol = pfmetcoll.product();
  const PFMET *pfmet;
  pfmet     = &(pfmetcol->front());
  PFMet_pt  = pfmet->pt();
  PFMet_eta = pfmet->eta();
  PFMet_phi = pfmet->phi();
  PFMet_en  = pfmet->energy();
  PFMet_px  = pfmet->px();
  PFMet_py  = pfmet->py();
  PFMet_pz  = pfmet->pz();
  //scalar sum of transverse energy over all objects
  PFMet_sumEt = pfmet->sumEt();
  //edm::Handle<double> metsighandle;
  //evt.getByLabel(theMETSignificance_, metsighandle);
  //METSign=*metsighandle;

}
//=============================================================
//
//                Method for PU Tree
//
//=============================================================
void HZZ4LeptonsRootTree::fillPU(const edm::Event& iEvent)
{
  edm::Handle<vector<PileupSummaryInfo> > PupInfo;
  iEvent.getByToken(PileupSrc_, PupInfo);
  if(!PupInfo.isValid()) return;
  for( vector<PileupSummaryInfo>::const_iterator cand = PupInfo->begin();cand != PupInfo->end(); ++ cand ) {
    if (cand->getBunchCrossing() == 0) num_PU_vertices=cand->getTrueNumInteractions();
    // num_PU_vertices = cand->getPU_NumInteractions();
    PU_BunchCrossing = cand->getBunchCrossing();
  }
}
//=============================================================                                                                                                           //                                                                                                                                                                        //          madgraph MC samples reweighing                                                                                                                                //                                                                                                                                                                        //=============================================================          
void HZZ4LeptonsRootTree::EventsMCReWeighting(const edm::Event& iEvent){
  MC_weighting.clear();
  float EventWeight = 1.0;
  edm::Handle<GenEventInfoProduct> gen_ev_info;
  iEvent.getByToken(generator_, gen_ev_info);
  if(!gen_ev_info.isValid()) return;
  EventWeight = gen_ev_info->weight();
  //std::cout<<"mc_weight = "<< gen_ev_info->weight() <<std::endl;                                                                                                                                       
  float mc_weight = ( EventWeight > 0 ) ? 1 : -1;
  //std::cout<<"mc_weight = "<< mc_weight <<std::endl;                                                                                                                                                   
  MC_weighting.push_back(mc_weight);
}


void HZZ4LeptonsRootTree::fillHLTFired(const edm::Event& iEvent)
{
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

//=============================================================
//
//                Method for trigger matching Tree
//
//=============================================================
void HZZ4LeptonsRootTree::triggermatching(const edm::Event& iEvent)
{
  RECOMU_PT_MuHLTMatch.clear();
  RECOELE_PT_EleHLTMatch.clear();
  cout << "Start Trigger matching for muon" << endl;
  // check HLTrigger/Configuration/python/HLT_GRun_cff.py
  
  edm::Handle<trigger::TriggerEvent> handleTriggerEvent;
  iEvent.getByToken(triggerEvent, handleTriggerEvent );
  const trigger::TriggerObjectCollection & toc(handleTriggerEvent->getObjects());
 
  size_t nMuHLT=0, nEleHLT=0;
  std::vector<reco::Particle>  HLTMuMatched, HLTEleMatched;
  std::vector<string> HLTMuMatchedNames,HLTEleMatchedNames;
  
  for ( size_t ia = 0; ia < handleTriggerEvent->sizeFilters(); ++ ia) {
    std::string fullname = handleTriggerEvent->filterTag(ia).encode();
    //std::cout<< "Trigger fullname::== " << fullname<< std::endl;
    std::string name;
    size_t p = fullname.find_first_of(':');
    if ( p != std::string::npos) {
      name = fullname.substr(0, p);
    }
    else {
      name = fullname;
    }
    //std::cout<< "name::== " << name<< std::endl;
    if ( &toc !=0 ) {
      const trigger::Keys & k = handleTriggerEvent->filterKeys(ia);
      for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
	// looking at all the single muon l3 trigger present, for example hltSingleMu15L3Filtered15.....
	// cout << "name=" << name << endl;
	// cout << "Trigger Filter " << triggerFilter.c_str()  << endl;
	//if (name.find("hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q")<100 && name.find("Mu") <100 && toc[*ki].pt()>40. && (fabs(toc[*ki].eta())<2.1)) {
	//  HLTMuMatched.push_back(toc[*ki].particle());
	//  HLTMuMatchedNames.push_back(name);
	//  nMuHLT++;
	//}
	if (name == triggerFilter.c_str()) { 
	  HLTMuMatched.push_back(toc[*ki].particle());
	  HLTMuMatchedNames.push_back(name); 
	  cout << "Matching " << triggerFilter.c_str()  << endl; 
	  nMuHLT++; 
	}
	if (name == triggerEleFilter.c_str()) { 
	  HLTEleMatched.push_back(toc[*ki].particle());
	  HLTEleMatchedNames.push_back(name); 
	  cout << "Matching " << triggerEleFilter.c_str()  << endl; 
	  nEleHLT++; 
	}
      }
    }
  }

  
  // Based on processing
  edm::Handle<edm::View<reco::Muon> > MuCandidates;
  iEvent.getByToken(muonTag_, MuCandidates);

  float maxDeltaR_=0.2;
  float maxDPtRel_=1.0;
  int nMuHLTMatch=0;
  
  for (edm::View<reco::Muon>::const_iterator iCand = MuCandidates->begin(); iCand != MuCandidates->end(); ++iCand){
    cout << "Muon with pt= " << iCand->pt() << ": check trigger matching" << endl;
    if (IsMuMatchedToHLTMu(*iCand,  HLTMuMatched , HLTMuMatchedNames, maxDeltaR_, maxDPtRel_)==true){
      nMuHLTMatch++;
      cout << "Muon HLT Matched with pT= " << iCand->pt() << endl;
      RECOMU_PT_MuHLTMatch.push_back(iCand->pt());
      RECOMU_ETA_MuHLTMatch.push_back(iCand->eta());
      RECOMU_PHI_MuHLTMatch.push_back(iCand->phi());
    }
  }
  cout << "N. Muons HLT Matched= " << nMuHLTMatch << " FiredString:" << HLTPathsFired << endl;
  RECO_nMuHLTMatch    = nMuHLTMatch;

  cout << "Start Trigger matching for electron" << endl;
  int nEleHLTMatch=0;
  edm::Handle<edm::View<reco::GsfElectron> > EleCandidates;
  iEvent.getByToken(electronEgmTag_, EleCandidates);

  for (edm::View<reco::GsfElectron>::const_iterator iCand = EleCandidates->begin(); iCand != EleCandidates->end(); ++iCand){
    cout << "Electron with pt= " << iCand->pt() << ": check trigger matching" << endl;
    if (IsEleMatchedToHLTEle(*iCand,  HLTEleMatched , HLTEleMatchedNames, maxDeltaR_, maxDPtRel_)==true){
      cout << "Electron HLT Matched with pT= " << iCand->pt() << endl;
      nEleHLTMatch++;
      RECOELE_PT_EleHLTMatch.push_back(iCand->pt());
      RECOELE_ETA_EleHLTMatch.push_back(iCand->eta());
      RECOELE_PHI_EleHLTMatch.push_back(iCand->phi());
    }
  }
  
  cout << "N. Electrons HLT Matched= " << nEleHLTMatch << endl;
  RECO_nEleHLTMatch = nEleHLTMatch;

  
}

bool HZZ4LeptonsRootTree::IsMuMatchedToHLTMu ( const reco::Muon &mu, std::vector<reco::Particle> HLTMu , std::vector<string> HLTMuNames, double DR, double DPtRel ) {
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

bool HZZ4LeptonsRootTree::IsEleMatchedToHLTEle ( const reco::GsfElectron &ele, std::vector<reco::Particle> HLTEle , std::vector<string> HLTEleNames, double DR, double DPtRel ) {
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
  

//=============================================================
//
//            Method for Genrated Particles Tree
//
//=============================================================  
// PDT
const std::string HZZ4LeptonsRootTree::getParticleName(int id){
  const ParticleData * pd = pdt_->particle( id );
  if (!pd) {
    std::ostringstream ss;
    ss << "P" << id;
    return ss.str();
  } else
    return pd->name();
}
  
// GenParticles  
void HZZ4LeptonsRootTree::fillgenparticles(const edm::Event& iEvent, const edm::EventSetup &es)
{
  MC_LEPT_PT.clear();
  MC_LEPT_ETA.clear();
  MC_LEPT_PHI.clear();
  MC_LEPT_THETA.clear();
  MC_LEPT_PDGID.clear();
  MC_Z_PT.clear();
  MC_Z_ETA.clear();
  MC_Z_PHI.clear();
  MC_Z_THETA.clear();
  MC_Z_MASS.clear();
  MC_Z_PDGID.clear();
  MC_fourl_MASS.clear();
  MC_fourl_PT.clear();
  MC_fourl_PDGID.clear();
  MC_ZZ_MASS.clear();
  MC_ZZ_PT.clear();
  MC_ZZ_ETA.clear();
  MC_ZZ_PHI.clear();
  MC_ZZ_THETA.clear();
  MC_ZZ_PDGID.clear();
 
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
      MC_LEPT_PT.push_back(leptonpt.at(i));
      MC_LEPT_ETA.push_back(leptoneta.at(i));
      MC_LEPT_PHI.push_back(leptonphi.at(i));
      MC_LEPT_THETA.push_back(leptontheta.at(i));
      MC_LEPT_PDGID.push_back(leptonpdgid.at(i));
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
      MC_Z_PT.push_back(mcIter->pt());       
      MC_Z_ETA.push_back(mcIter->eta());   
      MC_Z_PHI.push_back(mcIter->phi());
      MC_Z_THETA.push_back(mcIter->theta()); 
      MC_Z_MASS.push_back(mcIter->mass()); 
      MC_Z_PDGID.push_back(mcIter->pdgId());  
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
		MC_Z_PT.push_back((**it_dst1).pt());          
		MC_Z_ETA.push_back((**it_dst1).eta());      
		MC_Z_PHI.push_back((**it_dst1).phi());
		MC_Z_THETA.push_back((**it_dst1).theta());    
		MC_Z_MASS.push_back((**it_dst1).mass());    
		MC_Z_PDGID.push_back((**it_dst1).pdgId());
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
		if ( (**it_dst1).pt() > mcIter->pt() ) { // Save only highest pt photon
		  std::cout<<"   ===> Filled PHOT at ["<<i<<"]["<<ph_j<<"]";
		  MC_Z_PT.push_back((**it_dst1).pt());          
		  MC_Z_ETA.push_back((**it_dst1).eta());      
		  MC_Z_PHI.push_back((**it_dst1).phi());
		  MC_Z_THETA.push_back((**it_dst1).theta());    
		  MC_Z_MASS.push_back((**it_dst1).mass());    
		  MC_Z_PDGID.push_back((**it_dst1).pdgId());
		}
		else {
		  //std::cout<<"   ===> NOT Filled, pt = "<<MC_Z_PT[i][ph_j]<<" GeV/c";
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
    cout << " MC 4l Mass= " << mcIter->mass()
	 << " Charge= " 
	 << mcIter->daughter(0)->daughter(0)->charge() << " " 
	 << mcIter->daughter(0)->daughter(1)->charge() << " "
	 << mcIter->daughter(0)->daughter(2)->charge() << " " 
	 << mcIter->daughter(1)->charge() << " " 
	 << endl;
    MC_fourl_MASS.push_back(mcIter->p4().mass());
    MC_fourl_PT.push_back(mcIter->p4().pt());
    MC_fourl_PDGID.push_back(mcIter->pdgId());
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
	    MC_fourl_MASS.push_back(mcIter->daughter(j)->daughter(k)->p4().mass());
	    MC_fourl_PT.push_back(mcIter->daughter(j)->daughter(k)->p4().pt());
	    MC_fourl_PDGID.push_back(mcIter->daughter(j)->daughter(k)->pdgId());
	  }
	}
      }
      if (j==1) {
	//cout << "k+1= " << mcIter->daughter(0)->numberOfDaughters()+1 << endl;
	MC_fourl_MASS.push_back(mcIter->daughter(j)->p4().mass());
	MC_fourl_PT.push_back(mcIter->daughter(j)->p4().pt());
	MC_fourl_PDGID.push_back(mcIter->daughter(j)->pdgId());
      }
      ii++;
    }
    i++;
  }
  // get ZZ candidates
  i =0;
  kk=0;
  ii=0;
  l=0;
  edm::Handle<edm::View<Candidate> >  ZZCandidates;
  iEvent.getByToken(digenZ_, ZZCandidates);
  for (edm::View<Candidate>::const_iterator mcIterZZ=ZZCandidates->begin(); mcIterZZ!=ZZCandidates->end(); ++mcIterZZ ) {
    if (i>3 ) continue;
    cout << "MC ZZ Mass= " << mcIterZZ->p4().mass() 
	 << " and pT= " << mcIterZZ->p4().pt()  
	 << endl;
    MC_ZZ_MASS.push_back(mcIterZZ->p4().mass());
    MC_ZZ_PT.push_back(mcIterZZ->p4().pt());
    MC_ZZ_ETA.push_back(mcIterZZ->p4().eta());
    MC_ZZ_PHI.push_back(mcIterZZ->p4().phi());
    MC_ZZ_THETA.push_back(mcIterZZ->p4().theta());
    MC_ZZ_PDGID.push_back(mcIterZZ->pdgId());
    //int ii=0,l=0;
    for (unsigned j = 0; j < mcIterZZ->numberOfDaughters(); ++j ) {
      MC_ZZ_MASS.push_back(mcIterZZ->daughter(j)->p4().mass());
      MC_ZZ_PT.push_back(mcIterZZ->daughter(j)->p4().pt());
      MC_ZZ_ETA.push_back(mcIterZZ->daughter(j)->p4().eta());
      MC_ZZ_PHI.push_back(mcIterZZ->daughter(j)->p4().phi());
      MC_ZZ_THETA.push_back(mcIterZZ->daughter(j)->p4().theta());
      MC_ZZ_PDGID.push_back(mcIterZZ->daughter(j)->pdgId());
      for (unsigned k = 0; k < mcIterZZ->daughter(j)->numberOfDaughters(); ++k ) {
	if ( abs(mcIterZZ->daughter(j)->daughter(k)->pdgId())==13 || 
	     abs(mcIterZZ->daughter(j)->daughter(k)->pdgId())==15 || 
	     abs(mcIterZZ->daughter(j)->daughter(k)->pdgId())==11){
	  l=ii+j+kk+3; 	      
	  MC_ZZ_MASS.push_back(mcIterZZ->daughter(j)->daughter(k)->p4().mass());
	  MC_ZZ_PT.push_back(mcIterZZ->daughter(j)->daughter(k)->p4().pt());
	  MC_ZZ_ETA.push_back(mcIterZZ->daughter(j)->daughter(k)->p4().eta());
	  MC_ZZ_PHI.push_back(mcIterZZ->daughter(j)->daughter(k)->p4().phi());
	  MC_ZZ_THETA.push_back(mcIterZZ->daughter(j)->daughter(k)->p4().theta());
	  MC_ZZ_PDGID.push_back(mcIterZZ->daughter(j)->daughter(k)->pdgId());
	  kk++;
	}
      }
      ii++;	    	  
    }
  }
}

// MC Higgs 
void HZZ4LeptonsRootTree::fillmc(const edm::Event& iEvent)
{
  edm::Handle<edm::View<Candidate> > Candidates;
  iEvent.getByToken(MCcollName, Candidates);
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

void HZZ4LeptonsRootTree::SetMCValues(const reco::Candidate& cand, int nMC)
{
  MC_E.clear();
  MC_PT.clear();
  MC_ETA.clear();
  MC_THETA.clear();
  MC_PHI.clear();
  MC_MASS.clear();
  MC_PDGID.clear();
  //MC_E[nMC]     = cand.p4().energy();
  MC_E.push_back(cand.p4().energy());
  MC_PT.push_back(cand.p4().pt());
  MC_ETA.push_back(cand.p4().eta());
  MC_THETA.push_back(cand.p4().theta());
  MC_PHI.push_back(cand.p4().phi());
  MC_MASS.push_back(cand.p4().mass());
  MC_PDGID.push_back(cand.pdgId());
}



bool HZZ4LeptonsRootTree::match(double mass, double pt, int charge, 
			      const reco::CandidateCollection *c1Coll)
{
  bool found=false;
  for( CandidateCollection::const_iterator pp = c1Coll->begin();pp != c1Coll->end(); ++ pp ) {
    if ((abs(pp->p4().mass()-mass) <0.001 ) &&
	(abs(pp->p4().pt() -pt) <0.001 ) &&
	(abs(pp->charge() -charge)<0.001 ) ){
      found=true;
      cout << "Found lepton in the best leptons collection" << endl;
    }
  }
  return found;
}

bool HZZ4LeptonsRootTree::matchParticle(double mass, double pt, int charge, const reco::Candidate *c1)
{
  bool found=false;
  if ((fabs(c1->p4().mass()-mass)  <0.001 ) &&
      (fabs(c1->p4().pt()  -pt)    <0.001 ) &&
      (fabs(c1->charge()   -charge)<0.001 )  ){
    found=true;
    //std::cout << "Found lepton in the higgs-like leptons collection" << std::endl;
  }    
  return found;
}
  







//=============================================================
//
//                rho for isolation
//
//=============================================================
/*
void HZZ4LeptonsRootTree::fillRho(const edm::Event& evt){
  Rho2.clear();
  edm::Handle<double> rhoHandle;
  rhojetsTag_=edm::InputTag("kt6PFJetsCentralNeutral:rho");      
  iEvent.getByLabel(rhojetsTag_,rhoHandle); 
  if (rhoHandle.isValid() ) {
    RHO_mu=*rhoHandle;
    cout << "RHO mu fastjet= " << RHO_mu << endl; 
  }
      


  edm::Handle<double> rhoHandle;
  evt.getByLabel(rhoIsoInputTag_, rhoHandle);
  if(!rhoHandle.isValid()) return;
  //float RhoIsoValue = *rhoHandle;
  double RhoIsoValue = *(rhoHandle.product());
  Rho2.push_back(RhoIsoValue);
}
*/
//=============================================================
//
//     Method for initializing values for the variables
//
//=============================================================
void HZZ4LeptonsRootTree::IntialValues()
{
  value3_   = 0;
  iphot     = 0;
  iBtag     = 0;
}
//=============================================================
//
//                Method for Reco Muon Tree
//
//=============================================================
void HZZ4LeptonsRootTree::fillMuons(const edm::Event& iEvent,const edm::EventSetup& iSetup)
{ 
  RECO_NMU=0;;
  RECOMU_isPFMu.clear();
  RECOMU_isGlobalMu.clear();
  RECOMU_isStandAloneMu.clear();
  RECOMU_isTrackerMu.clear();
  RECOMU_isCaloMu.clear();
  // Kinematic of muon
  RECOMU_E.clear();
  RECOMU_PT.clear();
  RECOMU_P.clear();
  RECOMU_ETA.clear();
  RECOMU_THETA.clear();
  RECOMU_PHI.clear();
  RECOMU_MASS.clear();
  RECOMU_CHARGE.clear();
  // Covariance matrix
  RECOMU_COV.clear();
  // Isolation
  RECOMU_TRACKISO.clear();
  RECOMU_TRACKISO_SUMPT.clear();
  //temporary solution for reducedrechit problem
  RECOMU_ECALISO.clear();
  RECOMU_HCALISO.clear();
  RECOMU_X.clear();
  RECOMU_PFchHad.clear();
  RECOMU_PFneuHad.clear();
  RECOMU_PFphoton.clear();
  RECOMU_PFPUchAllPart.clear();
  // Vertexing
  RECOMU_SIP.clear();
  RECOMU_IP.clear();
  RECOMU_IPERROR.clear();
  RECOMU_SIP_KF.clear();
  RECOMU_IP_KF.clear();
  RECOMU_IPERROR_KF.clear();
  RECOMU_PFX_dB.clear();
  RECOMU_PFX_rho.clear();
  RECOMU_STIP.clear();
  RECOMU_SLIP.clear();
  RECOMU_TIP.clear();
  RECOMU_LIP.clear();
  RECOMU_TIPERROR.clear();
  RECOMU_LIPERROR.clear();
  // Other properties
  RECOMU_numberOfMatches.clear();
  RECOMU_numberOfMatchedStations.clear();
  RECOMU_caloCompatibility.clear();
  RECOMU_segmentCompatibility.clear();
  RECOMU_glbmuPromptTight.clear();
  // Track properties
  RECOMU_mubesttrkType.clear();
  RECOMU_mubesttrkDxy.clear();
  RECOMU_mubesttrkDxyB.clear();
  RECOMU_mubesttrkDxyError.clear();
  RECOMU_mubesttrkDz.clear();
  RECOMU_mubesttrkDzB.clear();
  RECOMU_mubesttrkDzError.clear();
  RECOMU_mubesttrkPTError.clear();
  
  RECOMU_mutrkPT.clear();
  RECOMU_mutrkPTError.clear();
  //RECOMU_mutrkPTError.clear();
  RECOMU_mutrkDxy.clear();
  RECOMU_mutrkDxyError.clear();
  RECOMU_mutrkDxyB.clear();
  RECOMU_mutrkDz.clear();
  RECOMU_mutrkDzError.clear();
  RECOMU_mutrkDzB.clear();
  RECOMU_mutrkChi2PerNdof.clear();
  RECOMU_mutrkCharge.clear();
  RECOMU_mutrkNHits.clear();
  RECOMU_mutrkNPixHits.clear();
  RECOMU_mutrkNStripHits.clear();
  RECOMU_mutrkNMuonHits.clear();
  RECOMU_mutrktrackerLayersWithMeasurement.clear();
  
  RECOMU_muInnertrkDxy.clear();
  RECOMU_muInnertrkDxyError.clear();
  RECOMU_muInnertrkDxyB.clear();
  RECOMU_muInnertrkDz.clear();
  RECOMU_muInnertrkDzError.clear();
  RECOMU_muInnertrkDzB.clear();
  RECOMU_muInnertrkChi2PerNdof.clear();
  RECOMU_muInnertrktrackerLayersWithMeasurement.clear();
  RECOMU_muInnertrkPT.clear();
  //RECOMU_muInnertrkPTError.clear();
  RECOMU_muInnertrkPTError.clear();
  RECOMU_muInnertrkCharge.clear();
  RECOMU_muInnertrkNHits.clear();
  RECOMU_muInnertrkNPixHits.clear();
  RECOMU_muInnertrkNStripHits.clear();
  
  // Tracker muon properties
  // RECOMU_trkmuArbitration.clear();
  RECOMU_trkmuArbitration.clear();
  RECOMU_trkmu2DCompatibilityLoose.clear();
  RECOMU_trkmu2DCompatibilityTight.clear();
  RECOMU_trkmuOneStationLoose.clear();
  RECOMU_trkmuOneStationTight.clear();
  RECOMU_trkmuLastStationLoose.clear();
  RECOMU_trkmuLastStationTight.clear();
  RECOMU_trkmuOneStationAngLoose.clear();
  RECOMU_trkmuOneStationAngTight.clear();
  RECOMU_trkmuLastStationAngLoose.clear();
  RECOMU_trkmuLastStationAngTight.clear();
  RECOMU_trkmuLastStationOptimizedLowPtLoose.clear();
  RECOMU_trkmuLastStationOptimizedLowPtTight.clear();

  RECOMU_MatchingMCTruth.clear();
  RECOMU_MatchingMCpT.clear();
  RECOMU_MatchingMCEta.clear();
  RECOMU_MatchingMCPhi.clear();

  // Magnetic Field
  edm::ESHandle<MagneticField> magfield_;
  iSetup.get<IdealMagneticFieldRecord>().get( magfield_ );        
  // Muons
  edm::Handle<edm::View<reco::Muon> > MuCandidates;
  iEvent.getByToken(muonTag_, MuCandidates);

  edm::Handle<edm::ValueMap<float> > corrpterrormumap;
  iEvent.getByToken(muonCorrPtErrorMapTag_,corrpterrormumap);
  
  // Vertexing
  // 3D
  edm::Handle<edm::View<reco::Muon> > VertMuCandidates;
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
  iEvent.getByToken(muonMapTag_VertKF, vertexmumapKF);
  
  edm::Handle<edm::ValueMap<float> > vertexmumapvalueKF;
  iEvent.getByToken(muonMapTag_VertValueKF, vertexmumapvalueKF);
  
  edm::Handle<edm::ValueMap<float> > vertexmumaperrorKF;
  iEvent.getByToken(muonMapTag_VertErrorKF, vertexmumaperrorKF);
  
  // Particle Flow Isolation
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
  
  // STIP SLIP
  edm::Handle<edm::ValueMap<float> > stipmumap;
  iEvent.getByToken(muonSTIPMapTag_Vert, stipmumap);
  
  edm::Handle<edm::ValueMap<float> > slipmumap;
  iEvent.getByToken(muonSLIPMapTag_Vert, slipmumap);
  
  edm::Handle<edm::ValueMap<float> > stipmumapvalue;
  iEvent.getByToken(muonSTIPMapTag_VertValue, stipmumapvalue);
  
  edm::Handle<edm::ValueMap<float> > slipmumapvalue;
  iEvent.getByToken(muonSLIPMapTag_VertValue, slipmumapvalue);
  
  edm::Handle<edm::ValueMap<float> > stipmumaperror;
  iEvent.getByToken(muonSTIPMapTag_VertError, stipmumaperror);
  
  edm::Handle<edm::ValueMap<float> > slipmumaperror;
  iEvent.getByToken(muonSLIPMapTag_VertError, slipmumaperror);
  
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

  // Matching
  edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchMu;
  iEvent.getByToken(goodMuonMCMatch_, GenParticlesMatchMu);
  edm::Handle<reco::CandidateCollection > CollMu;
  iEvent.getByToken(myMuons_, CollMu);
  if (GenParticlesMatchMu.isValid()){
    cout << endl<< "Muons:"<<endl<<"The reco collection to be matched has size= " <<  CollMu->size() << endl;
    cout << "The matched map collection has size= " <<  GenParticlesMatchMu->size() << endl;
  }
  //
  int indexbis=0;
  cout<<"SIZE MuCandidates :"<<MuCandidates->size()<<endl; 
  RECO_NMU=MuCandidates->size();
  
  for (edm::View<reco::Muon>::const_iterator cand = MuCandidates->begin();
       cand != MuCandidates->end(); ++cand) {

    edm::Ref<edm::View<reco::Muon> > mutrackref(MuCandidates,indexbis);
    edm::Ref<edm::View<reco::Muon> > mutrackrefv(VertMuCandidates,indexbis); 

    RECOMU_isPFMu.push_back(cand->isPFMuon());
    cout<<"RECOMU_isPFmu SIZE: "<<RECOMU_isPFMu.size()<<endl;
    RECOMU_isGlobalMu.push_back(cand->isGlobalMuon());
    cout<<"RECOMU_isGlobal SIZE: "<<RECOMU_isGlobalMu.size()<<endl;
    RECOMU_isStandAloneMu.push_back(cand->isStandAloneMuon());
    cout<<"RECOMU_isSTAMU SIZE: "<<RECOMU_isStandAloneMu.size()<<endl;
    RECOMU_isTrackerMu.push_back(cand->isTrackerMuon());
    cout<<"RECOMU_IS_trkMU SIZE: "<<RECOMU_isTrackerMu.size()<<endl;
    RECOMU_isCaloMu.push_back(cand->isCaloMuon());
    cout<<"RECOMU_IS_CaloMU SIZE: "<<RECOMU_isCaloMu.size()<<endl;
    
    /* std::cout << "\n Muon in the event: "
	<<   "  isPF=" << RECOMU_isPFMu[indexbis]
	<<   "  isGB=" << RECOMU_isGlobalMu[indexbis]
	<< "   isSTA=" << RECOMU_isStandAloneMu[indexbis]
	<<   "  isTM=" << RECOMU_isTrackerMu[indexbis]
	<< "  isCalo=" << RECOMU_isCaloMu[indexbis]
	<< std::endl;*/
    // Kinematic of muon
    RECOMU_E.push_back(cand->p4().energy());
    RECOMU_PT.push_back(cand->p4().pt());
    RECOMU_P.push_back(sqrt(cand->p4().px()*cand->p4().px()+cand->p4().py()*cand->p4().py()+cand->p4().pz()*cand->p4().pz()));
    RECOMU_ETA.push_back(cand->p4().eta());
    RECOMU_THETA.push_back(cand->p4().theta());
    RECOMU_PHI.push_back(cand->p4().phi());
    RECOMU_MASS.push_back(cand->p4().mass());
    RECOMU_CHARGE.push_back(cand->charge());
    /*std::cout << "--kinematic:"
      << "  pT="     << RECOMU_PT[indexbis]
      << "  E="      << RECOMU_E[indexbis]
      << "  p="      << RECOMU_P[indexbis]
      << "  eta="    << RECOMU_ETA[indexbis]
      << "  theta="  << RECOMU_THETA[indexbis]
      << "  phi="    << RECOMU_PHI[indexbis]
      << "  mass="   << RECOMU_MASS[indexbis]
      << "  charge=" << RECOMU_CHARGE[indexbis]
      << std::endl;*/
    
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
	  //RECOMU_COV[indexbis][i][j]=mat[i+3][j+3]; 
	  RECOMU_COV.push_back(mat[i+3][j+3]); 
	}  
      }  
    }

    // Isolation
    //RECOMU_TRACKISO.push_back((*isoTkmumap)[mutrackref]/cand->p4().pt());
    RECOMU_TRACKISO_SUMPT.push_back((cand->isolationR03().sumPt)/cand->p4().pt());
    //RECOMU_ECALISO.push_back((*isoEcalmumap)[mutrackref]/cand->p4().pt());
    //RECOMU_HCALISO.push_back((*isoHcalmumap)[mutrackref]/cand->p4().pt());
    //RECOMU_X.push_back([(*isoTkmumap)[mutrackref]+(*isoEcalmumap)[mutrackref]+(*isoHcalmumap)[mutrackref]]/cand->p4().pt());
    //temporary solution for reducedrechit problem
    RECOMU_ECALISO.push_back((cand->isolationR03().emEt)/cand->p4().pt());
    RECOMU_HCALISO.push_back((cand->isolationR03().hadEt)/cand->p4().pt());
    RECOMU_X.push_back((cand->isolationR03().sumPt)/cand->p4().pt()+(cand->isolationR03().emEt)/cand->p4().pt()+(cand->isolationR03().hadEt)/cand->p4().pt());
    //RECOMU_PFchHad.push_back((cand->pfIsolationR04().sumChargedHadronPt));
    //RECOMU_PFneuHad.push_back((cand->pfIsolationR04().sumNeutralHadronEt));
    //RECOMU_PFphoton.push_back((cand->pfIsolationR04().sumPhotonEt));
    //RECOMU_PFPUchAllPart.push_back((cand->pfIsolationR04().sumPUPt));
    RECOMU_PFchHad.push_back((*isoPFChargedmumap)[mutrackref]);
    RECOMU_PFneuHad.push_back((*isoPFNeutralmumap)[mutrackref]);
    RECOMU_PFphoton.push_back((*isoPFGammamumap)[mutrackref]);
    RECOMU_PFPUchAllPart.push_back((*isoPFPUmumap)[mutrackref]);
    // Vertexing
    RECOMU_SIP.push_back((*vertexmumap)[mutrackrefv]);
    RECOMU_IP.push_back((*vertexmumapvalue)[mutrackrefv]);
    RECOMU_IPERROR.push_back((*vertexmumaperror)[mutrackrefv]);
    
    RECOMU_SIP_KF.push_back((*vertexmumapKF)[mutrackrefv]);
    RECOMU_IP_KF.push_back((*vertexmumapvalueKF)[mutrackrefv]);
    RECOMU_IPERROR_KF.push_back((*vertexmumaperrorKF)[mutrackrefv]);
    
    /* begin Ahha*/
    //RECOMU_SIP_GD.push_back((*vertexmumapGD)[mutrackrefv];
    //if (decaychannel=="4mu" || decaychannel=="2e2mu" ) RECOMU_SIP_GDMMMM.push_back((*vertexmumapGDMMMM)[mutrackrefv];
    //      RECOMU_SIP_Std.push_back((*vertexmumapStd)[mutrackrefv];
    //if (decaychannel=="4mu" || decaychannel=="2e2mu" ) RECOMU_SIP_StdMMMM.push_back((*vertexmumapStdMMMM)[mutrackrefv];
    //RECOMU_SIP_Kin.push_back((*vertexmumapKin)[mutrackrefv];
    //if (decaychannel=="4mu" || decaychannel=="2e2mu" ) RECOMU_SIP_KinMMMM.push_back((*vertexmumapKinMMMM)[mutrackrefv];
    
    
    RECOMU_STIP.push_back((*stipmumap)[mutrackrefv]);
    RECOMU_SLIP.push_back((*slipmumap)[mutrackrefv]);
    RECOMU_TIP.push_back((*stipmumapvalue)[mutrackrefv]);
    RECOMU_LIP.push_back((*slipmumapvalue)[mutrackrefv]);
    RECOMU_TIPERROR.push_back((*stipmumaperror)[mutrackrefv]);
    RECOMU_LIPERROR.push_back((*slipmumaperror)[mutrackrefv]);
    
    
    /*std::cout << "--vertexing: muon"
      << "  Sign3DIP=" << RECOMU_SIP[indexbis]
      << "  3DIP="     << RECOMU_IP[indexbis]
      << "  3DIPerr="  << RECOMU_IPERROR[indexbis]
      << "  SignTIP="  << RECOMU_STIP[indexbis]
      << "  TIP="      << RECOMU_TIP[indexbis]
      << "  TIPerr="   << RECOMU_TIPERROR[indexbis]
      << "  SignLIP="  << RECOMU_SLIP[indexbis]
      << "  LIP="      << RECOMU_LIP[indexbis]
      << "  LIPerror=" << RECOMU_LIPERROR[indexbis]
      << std::endl;*/

    RECOMU_PFX_dB.push_back((cand->pfIsolationR04().sumChargedHadronPt + max(0.0,cand->pfIsolationR04().sumNeutralHadronEt+cand->pfIsolationR04().sumPhotonEt-0.5*cand->pfIsolationR04().sumPUPt))/cand->p4().pt());

    EffectiveArea=0.;
    if (use2011EA){
      if (fabs(cand->p4().eta()) >= 0.0 && fabs(cand->p4().eta()) < 1.0 ) EffectiveArea = 0.132;
      if (fabs(cand->p4().eta()) >= 1.0 && fabs(cand->p4().eta()) < 1.5 ) EffectiveArea = 0.120;
      if (fabs(cand->p4().eta()) >= 1.5 && fabs(cand->p4().eta()) < 2.0 ) EffectiveArea = 0.114;
      if (fabs(cand->p4().eta()) >= 2.0 && fabs(cand->p4().eta()) < 2.2 ) EffectiveArea = 0.139;
      if (fabs(cand->p4().eta()) >= 2.2 && fabs(cand->p4().eta()) < 2.3 ) EffectiveArea = 0.168;
      if (fabs(cand->p4().eta()) >= 2.3 )                                 EffectiveArea = 0.189;
    }

    else {  // 7_4_X use eta      
      if (fabs(cand->p4().eta()) >= 0.0 && fabs(cand->p4().eta()) < 1.0 ) EffectiveArea = 0.674;
      if (fabs(cand->p4().eta()) >= 1.0 && fabs(cand->p4().eta()) < 1.5 ) EffectiveArea = 0.565;
      if (fabs(cand->p4().eta()) >= 1.5 && fabs(cand->p4().eta()) < 2.0 ) EffectiveArea = 0.442;
      if (fabs(cand->p4().eta()) >= 2.0 && fabs(cand->p4().eta()) < 2.2 ) EffectiveArea = 0.515;
      if (fabs(cand->p4().eta()) >= 2.2 && fabs(cand->p4().eta()) < 2.3 ) EffectiveArea = 0.821;
      if (fabs(cand->p4().eta()) >= 2.3 )                                 EffectiveArea = 0.660;
    }

    
    //RECOMU_PFX_rho.push_back(rhoIso);
    RECOMU_PFX_rho.push_back((cand->pfIsolationR04().sumChargedHadronPt + max( (cand->pfIsolationR04().sumNeutralHadronEt + cand->pfIsolationR04().sumPhotonEt - max(RHO_mu,0.0)*(EffectiveArea)),0.0))/double(cand->p4().pt()));

    
    // Other properties
    RECOMU_numberOfMatches.push_back(cand->numberOfMatches());
    RECOMU_numberOfMatchedStations.push_back(cand->numberOfMatchedStations());
    RECOMU_caloCompatibility.push_back(cand->caloCompatibility());
    RECOMU_segmentCompatibility.push_back((muon::segmentCompatibility( (*cand))));
    RECOMU_glbmuPromptTight.push_back((muon::isGoodMuon( (*cand),muon::GlobalMuonPromptTight)));
    
    
    /*std::cout	<< "--other properties:"
      << "  n.matches="   << RECOMU_numberOfMatches[indexbis]
      << "  caloComp="    << RECOMU_caloCompatibility[indexbis]
      << "  segmentComp=" << RECOMU_segmentCompatibility[indexbis]
      << "  glbmuPromptTight=" << RECOMU_glbmuPromptTight[indexbis]
      << std::endl;*/
    
    
    // Track properties
    if(cand->muonBestTrack().isAvailable()){
      RECOMU_mubesttrkType.push_back(cand->muonBestTrackType());
      RECOMU_mubesttrkDxy.push_back(cand->muonBestTrack()->dxy(pVertex));
      RECOMU_mubesttrkDxyB.push_back(cand->muonBestTrack()->dxy(bs.position()));
      RECOMU_mubesttrkDxyError.push_back(cand->muonBestTrack()->dxyError());
      RECOMU_mubesttrkDz.push_back(cand->muonBestTrack()->dz(pVertex));
      RECOMU_mubesttrkDzB.push_back(cand->muonBestTrack()->dz(bs.position()));
      RECOMU_mubesttrkDzError.push_back(cand->muonBestTrack()->dzError());
      //RECOMU_mubesttrkPTError.push_back((*corrpterrormumap)[mutrackref]);
    }
    
    if(cand->globalTrack().isAvailable()){
      RECOMU_mutrkPT.push_back(cand->globalTrack()->pt());
      RECOMU_mutrkPTError.push_back(cand->globalTrack()->ptError());
      //RECOMU_mutrkPTError.push_back(cand->bestTrack()->ptError()); // bestTrack
      RECOMU_mutrkDxy.push_back(cand->globalTrack()->dxy(pVertex));
      RECOMU_mutrkDxyError.push_back(cand->globalTrack()->dxyError());
      RECOMU_mutrkDxyB.push_back(cand->globalTrack()->dxy(bs.position()));
      RECOMU_mutrkDz.push_back(cand->globalTrack()->dz(pVertex));
      RECOMU_mutrkDzError.push_back(cand->globalTrack()->dzError());
      RECOMU_mutrkDzB.push_back(cand->globalTrack()->dz(bs.position()));
      RECOMU_mutrkChi2PerNdof.push_back(cand->globalTrack()->normalizedChi2());
      RECOMU_mutrkCharge.push_back(cand->globalTrack()->charge());
      RECOMU_mutrkNHits.push_back(cand->globalTrack()->numberOfValidHits()); 
      RECOMU_mutrkNPixHits.push_back(cand->globalTrack()->hitPattern().numberOfValidPixelHits());
      RECOMU_mutrkNStripHits.push_back(cand->globalTrack()->hitPattern().numberOfValidStripHits());
      RECOMU_mutrkNMuonHits.push_back(cand->globalTrack()->hitPattern().numberOfValidMuonHits()); 
      RECOMU_mutrktrackerLayersWithMeasurement.push_back(cand->globalTrack()->hitPattern().trackerLayersWithMeasurement()); 
      
      RECOMU_muInnertrkDxy.push_back(cand->innerTrack()->dxy(pVertex));
      RECOMU_muInnertrkDxyError.push_back(cand->innerTrack()->dxyError());
      RECOMU_muInnertrkDxyB.push_back(cand->innerTrack()->dxy(bs.position()));
      RECOMU_muInnertrkDz.push_back(cand->innerTrack()->dz(pVertex));
      RECOMU_muInnertrkDzError.push_back(cand->innerTrack()->dzError());
      RECOMU_muInnertrkDzB.push_back(cand->innerTrack()->dz(bs.position()));
      RECOMU_muInnertrkChi2PerNdof.push_back(cand->innerTrack()->normalizedChi2());
      RECOMU_muInnertrktrackerLayersWithMeasurement.push_back(cand->innerTrack()->hitPattern().trackerLayersWithMeasurement()); 
      RECOMU_muInnertrkPT.push_back(cand->innerTrack()->pt());	
      //RECOMU_muInnertrkPTError.push_back(cand->innerTrack()->ptError());
      RECOMU_muInnertrkPTError.push_back(cand->bestTrack()->ptError()); // Besttrack

      RECOMU_muInnertrkCharge.push_back(cand->innerTrack()->charge());
      RECOMU_muInnertrkNHits.push_back(cand->innerTrack()->numberOfValidHits()); 
      RECOMU_muInnertrkNPixHits.push_back(cand->innerTrack()->hitPattern().numberOfValidPixelHits());
      RECOMU_muInnertrkNStripHits.push_back(cand->innerTrack()->hitPattern().numberOfValidStripHits());

    }
    else if(cand->innerTrack().isAvailable()){
      RECOMU_muInnertrkDxy.push_back(cand->innerTrack()->dxy(pVertex));
      RECOMU_muInnertrkDxyError.push_back(cand->innerTrack()->dxyError());
      RECOMU_muInnertrkDxyB.push_back(cand->innerTrack()->dxy(bs.position()));
      RECOMU_muInnertrkDz.push_back(cand->innerTrack()->dz(pVertex));
      RECOMU_muInnertrkDzError.push_back(cand->innerTrack()->dzError());
      RECOMU_muInnertrkDzB.push_back(cand->innerTrack()->dz(bs.position()));
      RECOMU_muInnertrkChi2PerNdof.push_back(cand->innerTrack()->normalizedChi2());
      RECOMU_muInnertrktrackerLayersWithMeasurement.push_back(cand->innerTrack()->hitPattern().trackerLayersWithMeasurement()); 
      RECOMU_muInnertrkPT.push_back(cand->innerTrack()->pt());	
      //RECOMU_muInnertrkPTError.push_back(cand->innerTrack()->ptError());
      RECOMU_muInnertrkPTError.push_back(cand->bestTrack()->ptError()); // Besttrack
      RECOMU_muInnertrkCharge.push_back(cand->innerTrack()->charge());
      RECOMU_muInnertrkNHits.push_back(cand->innerTrack()->numberOfValidHits()); 
      RECOMU_muInnertrkNPixHits.push_back(cand->innerTrack()->hitPattern().numberOfValidPixelHits());
      RECOMU_muInnertrkNStripHits.push_back(cand->innerTrack()->hitPattern().numberOfValidStripHits());
    }

    if(cand->globalTrack().isAvailable() || cand->innerTrack().isAvailable() ){
      /*std::cout << "--muon track properties: "
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
		  << std::endl;*/
	
	
      // Tracker muon properties
      // RECOMU_trkmuArbitration.push_back((muon::isGoodMuon( (*cand),muon::TrackerMuonArbitrated)));
      RECOMU_trkmuArbitration.push_back((muon::segmentCompatibility((*cand),reco::Muon::SegmentAndTrackArbitration)));
      RECOMU_trkmu2DCompatibilityLoose.push_back((muon::isGoodMuon( (*cand),muon::TM2DCompatibilityLoose)));
      RECOMU_trkmu2DCompatibilityTight.push_back((muon::isGoodMuon( (*cand),muon::TM2DCompatibilityTight)));
      RECOMU_trkmuOneStationLoose.push_back((muon::isGoodMuon( (*cand),muon::TMOneStationLoose)));
      RECOMU_trkmuOneStationTight.push_back((muon::isGoodMuon( (*cand),muon::TMOneStationTight)));
      RECOMU_trkmuLastStationLoose.push_back((muon::isGoodMuon( (*cand),muon::TMLastStationLoose)));
      RECOMU_trkmuLastStationTight.push_back((muon::isGoodMuon( (*cand),muon::TMLastStationTight)));
      RECOMU_trkmuOneStationAngLoose.push_back((muon::isGoodMuon( (*cand),muon::TMOneStationAngLoose)));
      RECOMU_trkmuOneStationAngTight.push_back((muon::isGoodMuon( (*cand),muon::TMOneStationAngTight)));
      RECOMU_trkmuLastStationAngLoose.push_back((muon::isGoodMuon( (*cand),muon::TMLastStationAngLoose)));
      RECOMU_trkmuLastStationAngTight.push_back((muon::isGoodMuon( (*cand),muon::TMLastStationAngTight)));
      RECOMU_trkmuLastStationOptimizedLowPtLoose.push_back((muon::isGoodMuon( (*cand),muon::TMLastStationOptimizedLowPtLoose)));
      RECOMU_trkmuLastStationOptimizedLowPtTight.push_back((muon::isGoodMuon( (*cand),muon::TMLastStationOptimizedLowPtTight)));
      /*std::cout << "--tracker muon properties:"
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
	<< std::endl;*/
    }
      
    // Matching
    if (fillMCTruth==true){
      int i=0;
      for ( reco::CandidateCollection::const_iterator hIter=CollMu->begin(); hIter!= CollMu->end(); ++hIter ){
	//cout << "Reco Muon with pT= " << hIter->pt() << " and mass="<< hIter->mass()<< endl;
	if (fabs(hIter->pt()-cand->p4().pt())<0.01){
	  i=hIter-(CollMu->begin());
	  CandidateRef Ref( CollMu, i );
	  edm::Ref<std::vector<reco::GenParticle> > genrefMu = (*GenParticlesMatchMu)[Ref];
	  if (!genrefMu.isNull()){
	    cout << "GenMuon with pT= " << genrefMu->p4().pt() << " and mass="<< genrefMu->p4().mass()<< endl;
	    RECOMU_MatchingMCTruth.push_back(true);
	    RECOMU_MatchingMCpT.push_back(genrefMu->p4().pt());
	    RECOMU_MatchingMCEta.push_back(genrefMu->p4().eta());
	    RECOMU_MatchingMCPhi.push_back(genrefMu->p4().phi());	    
	  } 
	  else {
	    cout << "There is no reference to a genMuon" << endl;
	  }
	}   
      }   
    }   


      /*end Ahha */
    indexbis++;
  }
}


//void fillGD2e2mu(const edm::Event& iEvent){
//=============================================================
//
//                Method for Reco GSF Ele Tree
//
//=============================================================
void HZZ4LeptonsRootTree::fillElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  RECO_NELE=0;
  RECOELE_isEcalDriven.clear();
  RECOELE_isTrackerDriven.clear();
  RECOELE_CHARGE.clear();
  RECOELE_E.clear();
  RECOELE_PT.clear();
  RECOELE_P.clear();
  RECOELE_ETA.clear();
  RECOELE_THETA.clear();
  RECOELE_PHI.clear();
  RECOELE_MASS.clear();
  ele_sclRawE.clear();
  RECOELE_scl_E.clear();
  RECOELE_scl_Et.clear();
  RECOELE_scl_Eta.clear();
  RECOELE_scl_Phi.clear();
  ele_sclX.clear();
  ele_sclY.clear();
  ele_sclZ.clear();
  RECOELE_PTError.clear();
  RECOELE_COV.clear();
  RECOELE_EGMECALISO.clear();
  RECOELE_EGMHCALISO.clear();
  RECOELE_EGMX.clear();
  // PF isolation
  RECOELE_PFchAllPart.clear();
  RECOELE_PFchHad.clear();
  RECOELE_PFneuHad.clear();
  RECOELE_PFphoton.clear();
  RECOELE_PFPUchAllPart.clear();
  RECOELE_PFX_dB.clear();
  RECOELE_PFX_rho.clear();
  // Vertexing DA
  RECOELE_SIP.clear();
  RECOELE_IP.clear();
  RECOELE_IPERROR.clear();
  // KF
  RECOELE_SIP_KF.clear();
  RECOELE_IP_KF.clear();
  RECOELE_IPERROR_KF.clear();
  RECOELE_STIP.clear();
  RECOELE_SLIP.clear();
  RECOELE_TIP.clear();
  RECOELE_LIP.clear();
  RECOELE_TIPERROR.clear();
  RECOELE_LIPERROR.clear();
  // GsfTrack
  RECOELE_gsftrack_NPixHits.clear();
  RECOELE_gsftrack_NStripHits.clear();
  RECOELE_gsftrack_chi2.clear();
  RECOELE_gsftrack_dxyB.clear();
  RECOELE_gsftrack_dxy.clear();
  RECOELE_gsftrack_dxyError.clear();
  RECOELE_gsftrack_dzB.clear();
  RECOELE_gsftrack_dz.clear();
  RECOELE_gsftrack_dzError.clear();
  //Conversion variables
  RECOELE_gsftrack_losthits.clear();
  RECOELE_gsftrack_validhits.clear();
  RECOELE_gsftrack_expected_inner_hits.clear();
  // Track-Cluster matching attributes
  RECOELE_ep.clear();
  RECOELE_eSeedp.clear();
  RECOELE_eSeedpout.clear();
  RECOELE_eElepout.clear();
  RECOELE_deltaEtaIn.clear();
  RECOELE_deltaEtaSeed.clear();
  RECOELE_deltaEtaEle.clear();
  RECOELE_deltaPhiIn.clear();
  RECOELE_deltaPhiSeed.clear();
  RECOELE_deltaPhiEle.clear();
  // Fiducial flags
  RECOELE_isbarrel.clear();
  RECOELE_isendcap.clear();
  RECOELE_isGap.clear();
  RECOELE_isEBetaGap.clear();
  RECOELE_isEBphiGap.clear();
  RECOELE_isEEdeeGap.clear();
  RECOELE_isEEringGap.clear();
  // Shower shape
  RECOELE_sigmaIetaIeta.clear();
  RECOELE_sigmaEtaEta.clear();
  RECOELE_e15.clear();
  RECOELE_e25max.clear();
  RECOELE_e55.clear();
  RECOELE_he.clear();
  //RECOELE_r9.clear();
  // Brem & Classifaction
  RECOELE_fbrem.clear();
  RECOELE_nbrems.clear();
  RECOELE_Class.clear();
  RECOELE_fbrem_mean.clear();
  RECOELE_fbrem_mode.clear();
  // Corrections
  RECOELE_ecalEnergy.clear();
  // Seed Collection
  ele_seedSubdet2.clear();
  ele_seedDphi2.clear();
  ele_seedDrz2.clear();
  ele_seedSubdet1.clear();
  ele_seedDphi1.clear();
  ele_seedDrz1.clear();
  RECOELE_mvaTrigV0.clear();
  RECOELE_mvaNonTrigV0.clear();
  // Conversion Finder
  ConvMapDist.clear();
  ConvMapDcot.clear();
  // Matching
  RECOELE_MatchingMCTruth.clear();
  RECOELE_MatchingMCpT.clear();
  RECOELE_MatchingMCEta.clear();
  RECOELE_MatchingMCPhi.clear();
  
  // Supercluster collection
  edm::Handle<reco::SuperClusterCollection> clusters;
  iEvent.getByLabel(clusterCollectionTag_,clusters);
  
  // EG isolation
  //edm::Handle<reco::GsfElectronRefVector> EleRefs;
  //iEvent.getByLabel(electronEgmTag_, EleRefs);
  edm::Handle<edm::View<reco::GsfElectron> > EleRefs;
  iEvent.getByToken(electronEgmTag_, EleRefs);
  
  //    edm::Handle<edm::ValueMap<double> > egmisoTkelemap;
  //iEvent.getByLabel(electronEgmTkMapTag_, egmisoTkelemap);
  
  //    edm::Handle<edm::ValueMap<double> > egmisoEcalelemap;
  //iEvent.getByLabel(electronEgmEcalMapTag_, egmisoEcalelemap);
  
  //edm::Handle<edm::ValueMap<double> > egmisoHcalelemap;
  //iEvent.getByLabel(electronEgmHcalMapTag_, egmisoHcalelemap);
  
  //Electron ID MVA Trig and Non Trig
  edm::Handle<edm::ValueMap<float> >  mvaTrigV0map;
  iEvent.getByToken(mvaTrigV0MapTag_, mvaTrigV0map);
  edm::Handle<edm::ValueMap<float> >  mvaNonTrigV0map;
  iEvent.getByToken(mvaNonTrigV0MapTag_, mvaNonTrigV0map);
  
  // Particle Flow Isolation
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
  
  // electron regression
  //edm::Handle<edm::ValueMap<double> > eleRegressionEnergymap;
  //iEvent.getByLabel(eleRegressionEnergyTag_, eleRegressionEnergymap);
  
  //edm::Handle<edm::ValueMap<double> > eleRegressionEnergyErrormap;
  //iEvent.getByLabel(eleRegressionEnergyErrorTag_, eleRegressionEnergyErrormap);
  
  // Vertexing
  //3D DA
  edm::Handle<edm::View<reco::GsfElectron> > VertEleCandidates;
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
  
  // Conversion Finder
  edm::Handle<edm::ValueMap<float> > conversiondistmap;
  iEvent.getByToken(ConvMapDistTag_, conversiondistmap);
  
  edm::Handle<edm::ValueMap<float> > conversiondcotmap;
  iEvent.getByToken(ConvMapDcotTag_, conversiondcotmap);
  
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
  // Matching
  edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchEle;
  iEvent.getByToken(goodElectronMCMatch_, GenParticlesMatchEle);
  edm::Handle<reco::CandidateCollection > CollEle;
  iEvent.getByToken(myElectrons_, CollEle);
  if (GenParticlesMatchEle.isValid()){
    cout << endl<< "Electrons:"<<endl<<"The reco collection to be matched has size= " <<  CollEle->size() << endl;
    cout << "The matched map collection has size= " <<  GenParticlesMatchEle->size() << endl;
  }
  //
  for (edm::View<reco::GsfElectron>::const_iterator cand = EleRefs->begin(); 
       cand != EleRefs->end(); ++cand) {
    edm::Ref<edm::View<reco::GsfElectron> > eletrackref(EleRefs,index);
    edm::Ref<edm::View<reco::GsfElectron> > eletrackrefv(VertEleCandidates,index);
    // Global variables
    RECOELE_isEcalDriven.push_back(cand->ecalDrivenSeed());
    RECOELE_isTrackerDriven.push_back(cand->trackerDrivenSeed());
    // kinematic
    RECOELE_E.push_back(cand->p4().energy());
    RECOELE_PT.push_back(cand->p4().pt());
    RECOELE_P.push_back(sqrt(cand->p4().px()*cand->p4().px()+cand->p4().py()*cand->p4().py()+cand->p4().pz()*cand->p4().pz()));
    RECOELE_ETA.push_back(cand->p4().eta());
    RECOELE_THETA.push_back(cand->p4().theta());
    RECOELE_PHI.push_back(cand->p4().phi());
    RECOELE_MASS.push_back(cand->p4().mass());
    RECOELE_CHARGE.push_back(cand->charge());
    // SuperCluster
    reco::SuperClusterRef sclRef = cand->superCluster();
    math::XYZPoint sclPos = cand->superClusterPosition();
    //if(!cand->ecalDrivenSeed() && cand->trackerDrivenSeed())
    //clRef = cand->pflowSuperCluster();
    double R  = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
    double Rt = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
    ele_sclRawE.push_back(sclRef->rawEnergy()) ;
    RECOELE_scl_E.push_back(sclRef->energy()) ;
    RECOELE_scl_Et.push_back(sclRef->energy()*(Rt/R)) ;
    RECOELE_scl_Eta.push_back(sclRef->eta()) ;
    RECOELE_scl_Phi.push_back(sclRef->phi()) ;
    ele_sclX.push_back(sclPos.X());
    ele_sclY.push_back(sclPos.Y());
    ele_sclZ.push_back(sclPos.Z());
    //Covariance matrix
    TMatrixDSym bigCov;
    double dp = 0.;
    if (cand->ecalDriven()) {
      dp = cand->p4Error(reco::GsfElectron::P4_COMBINATION);	
      RECOELE_PTError.push_back(dp);
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
      RECOELE_PTError.push_back(dp);
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
	RECOELE_COV.push_back(mat[i][j]);
      }
    }
    
    // Egamma isolation 
    //RECOELE_EGMTRACKISO.push_back((*egmisoTkelemap)[eletrackref]/eletrackref->pt());
    // RECOELE_EGMECALISO.push_back((*egmisoEcalelemap)[eletrackref]/eletrackref->pt());
    // RECOELE_EGMHCALISO.push_back((*egmisoHcalelemap)[eletrackref]/eletrackref->pt());
    // RECOELE_EGMX.push_back(( (*egmisoTkelemap)[eletrackref] + (*egmisoEcalelemap)[eletrackref] + (*egmisoHcalelemap)[eletrackref])/eletrackref->pt());
    
    // temporary solution for the problem of rechits
    RECOELE_EGMECALISO.push_back((eletrackref->dr03EcalRecHitSumEt())/eletrackref->pt());
    RECOELE_EGMHCALISO.push_back((eletrackref->dr03HcalTowerSumEt())/eletrackref->pt());
    RECOELE_EGMX.push_back(eletrackref->dr03TkSumPt()/eletrackref->pt()+
			   (eletrackref->dr03EcalRecHitSumEt())/eletrackref->pt()+
			   (eletrackref->dr03HcalTowerSumEt())/eletrackref->pt());
    
    // PF isolation
    //RECOELE_PFchAllPart.push_back((*isoPFChargedAllelemap)[eletrackref]);
    //RECOELE_PFchHad.push_back((*isoPFChargedelemap)[eletrackref]);
    //RECOELE_PFneuHad.push_back((*isoPFNeutralelemap)[eletrackref]);
    //RECOELE_PFphoton.push_back((*isoPFGammaelemap)[eletrackref]);
    //RECOELE_PFPUchAllPart.push_back((*isoPFPUelemap)[eletrackref]);
    RECOELE_PFX_dB.push_back(((*isoPFChargedelemap)[eletrackref] + max((*isoPFGammaelemap)[eletrackref]+(*isoPFNeutralelemap)[eletrackref]-0.5*(*isoPFPUelemap)[eletrackref],0.0))/eletrackref->pt());  

    
    EffectiveArea=0.0;
    if (use2011EA){
      if (fabs(sclRef->eta()) >= 0.0   && fabs(sclRef->eta()) < 1.0 )   EffectiveArea = 0.18;
      if (fabs(sclRef->eta()) >= 1.0   && fabs(sclRef->eta()) < 1.479 ) EffectiveArea = 0.20;
      if (fabs(sclRef->eta()) >= 1.479 && fabs(sclRef->eta()) < 2.0 )   EffectiveArea = 0.15;
      if (fabs(sclRef->eta()) >= 2.0   && fabs(sclRef->eta()) < 2.2 )   EffectiveArea = 0.19;
      if (fabs(sclRef->eta()) >= 2.2   && fabs(sclRef->eta()) < 2.3 )   EffectiveArea = 0.21;
      if (fabs(sclRef->eta()) >= 2.3   && fabs(sclRef->eta()) < 2.4 )   EffectiveArea = 0.22;
      if (fabs(sclRef->eta()) >= 2.4 )                                  EffectiveArea = 0.29;
    }
    else {
      // 7_4_X use eta 
      // if (fabs(cand->p4().eta()) >= 0.0   && fabs(cand->p4().eta()) < 0.8 )   EffectiveArea = 0.1830;
      // if (fabs(cand->p4().eta()) >= 0.8   && fabs(cand->p4().eta()) < 1.3 )   EffectiveArea = 0.1734;
      // if (fabs(cand->p4().eta()) >= 1.3   && fabs(cand->p4().eta()) < 2.0 )   EffectiveArea = 0.1077;
      // if (fabs(cand->p4().eta()) >= 2.0   && fabs(cand->p4().eta()) < 2.2 )   EffectiveArea = 0.1565;
      //if (fabs(cand->p4().eta()) >= 2.2 )                                     EffectiveArea = 0.2680;

      // 7_6_X use eta supercluster   
        if (fabs(sclRef->eta()) >= 0.0   && fabs(sclRef->eta()) < 1.0 )   EffectiveArea = 0.1752;
        if (fabs(sclRef->eta()) >= 1.0   && fabs(sclRef->eta()) < 1.479 ) EffectiveArea = 0.1862;
        if (fabs(sclRef->eta()) >= 1.479 && fabs(sclRef->eta()) < 2.0 )   EffectiveArea = 0.1411;
        if (fabs(sclRef->eta()) >= 2.0   && fabs(sclRef->eta()) < 2.2 )   EffectiveArea = 0.1534;
        if (fabs(sclRef->eta()) >= 2.2   && fabs(sclRef->eta()) < 2.3 )   EffectiveArea = 0.1903;
        if (fabs(sclRef->eta()) >= 2.3   && fabs(sclRef->eta()) < 2.4 )   EffectiveArea = 0.2243;
        if (fabs(sclRef->eta()) >= 2.4   && fabs(sclRef->eta()) < 5.0  )  EffectiveArea = 0.2687;
    }

    //RECOELE_PFX_rho.push_back(((*isoPFChargedelemap)[eletrackref]+ max( ((*isoPFNeutralelemap)[eletrackref]+(*isoPFGammaelemap)[eletrackref]- max(RHO_ele,0.0)*EffectiveArea),0.0))/cand->p4().pt()); 

    RECOELE_PFX_rho.push_back((cand->pfIsolationVariables().sumChargedHadronPt+
			       std::max(
					cand->pfIsolationVariables().sumPhotonEt+
					cand->pfIsolationVariables().sumNeutralHadronEt-
					max(RHO_ele,0.0)*EffectiveArea,
					0.0)
			       )/cand->p4().pt());
    
    RECOELE_PFchHad.push_back(cand->pfIsolationVariables().sumChargedHadronPt);
    RECOELE_PFneuHad.push_back(cand->pfIsolationVariables().sumNeutralHadronEt);
    RECOELE_PFphoton.push_back(cand->pfIsolationVariables().sumPhotonEt);
    //RECOELE_PFsumPUPt.push_back(cand->pfIsolationVariables().sumPUPt);                                                                                                             
    RECOELE_PFPUchAllPart.push_back(cand->pfIsolationVariables().sumPUPt);
       
    // Vertexing DA
    RECOELE_SIP.push_back((*vertexelemap)[eletrackrefv]);
    RECOELE_IP.push_back((*vertexelemapvalue)[eletrackrefv]);
    RECOELE_IPERROR.push_back((*vertexelemaperror)[eletrackrefv]);
    
    // KF
    RECOELE_SIP_KF.push_back((*vertexelemapKF)[eletrackrefv]);
    RECOELE_IP_KF.push_back((*vertexelemapvalueKF)[eletrackrefv]);
    RECOELE_IPERROR_KF.push_back((*vertexelemaperrorKF)[eletrackrefv]);
    
    
    //RECOELE_SIP_GD.push_back((*vertexelemapGD)[eletrackrefv]);
    //if (decaychannel=="4e" || decaychannel=="2e2mu" ) RECOELE_SIP_GDEEEE.push_back((*vertexelemapGDEEEE)[eletrackrefv]);
    //RECOELE_SIP_Std.push_back((*vertexelemapStd)[eletrackrefv]; 
    //if (decaychannel=="4e" || decaychannel=="2e2mu" ) RECOELE_SIP_StdEEEE.push_back((*vertexelemapStdEEEE)[eletrackrefv]); 
    //RECOELE_SIP_Kin.push_back((*vertexelemapKin)[eletrackrefv]; 
    //if (decaychannel=="4e" || decaychannel=="2e2mu" ) RECOELE_SIP_KinEEEE.push_back((*vertexelemapKinEEEE)[eletrackrefv]); 
    
    
    RECOELE_STIP.push_back((*stipelemap)[eletrackrefv]);
    RECOELE_SLIP.push_back((*slipelemap)[eletrackrefv]);
    RECOELE_TIP.push_back((*stipelemapvalue)[eletrackrefv]);
    RECOELE_LIP.push_back((*slipelemapvalue)[eletrackrefv]);
    RECOELE_TIPERROR.push_back((*stipelemaperror)[eletrackrefv]);
    RECOELE_LIPERROR.push_back((*slipelemaperror)[eletrackrefv]);
    // GsfTrack
    RECOELE_gsftrack_NPixHits.push_back(cand->gsfTrack()->hitPattern().numberOfValidPixelHits());
    RECOELE_gsftrack_NStripHits.push_back(cand->gsfTrack()->hitPattern().numberOfValidStripHits());
    RECOELE_gsftrack_chi2.push_back(cand->gsfTrack()->normalizedChi2());
    RECOELE_gsftrack_dxyB.push_back(cand->gsfTrack()->dxy(bs.position()));
    RECOELE_gsftrack_dxy.push_back(cand->gsfTrack()->dxy(pVertex));
    RECOELE_gsftrack_dxyError.push_back(cand->gsfTrack()->dxyError());
    RECOELE_gsftrack_dzB.push_back(cand->gsfTrack()->dz(bs.position()));
    RECOELE_gsftrack_dz.push_back(cand->gsfTrack()->dz(pVertex));
    RECOELE_gsftrack_dzError.push_back(cand->gsfTrack()->dzError());
    
    //Conversion variables
    RECOELE_gsftrack_losthits.push_back(cand->gsfTrack()->numberOfLostHits());
    RECOELE_gsftrack_validhits.push_back(cand->gsfTrack()->numberOfValidHits());
    RECOELE_gsftrack_expected_inner_hits.push_back(cand->gsfTrack()->hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS));
    /*std::cout << "--gfstrack properties: " 
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
      << std::endl;*/

    // Track-Cluster matching attributes
    RECOELE_ep.push_back(cand->eSuperClusterOverP()); 
    RECOELE_eSeedp.push_back(cand->eSeedClusterOverP());     
    RECOELE_eSeedpout.push_back(cand->eSeedClusterOverPout()); 
    RECOELE_eElepout.push_back(cand->eEleClusterOverPout());        
    RECOELE_deltaEtaIn.push_back(cand->deltaEtaSuperClusterTrackAtVtx()); 
    RECOELE_deltaEtaSeed.push_back(cand->deltaEtaSeedClusterTrackAtCalo()); 
    RECOELE_deltaEtaEle.push_back(cand->deltaEtaEleClusterTrackAtCalo());  
    RECOELE_deltaPhiIn.push_back(cand->deltaPhiSuperClusterTrackAtVtx());
    RECOELE_deltaPhiSeed.push_back(cand->deltaPhiSeedClusterTrackAtCalo()); 
    RECOELE_deltaPhiEle.push_back(cand->deltaPhiEleClusterTrackAtCalo());  
    /*std::cout << "--track-cluster matching: " 
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
      << std::endl;*/
    
    // Fiducial flags
    if (cand->isEB()) RECOELE_isbarrel.push_back(1); 
    else  RECOELE_isbarrel.push_back(0);
    if (cand->isEE()) RECOELE_isendcap.push_back(1) ; 
    else  RECOELE_isendcap.push_back(0);
    if (cand->isGap())       RECOELE_isGap.push_back(1);
    if (cand->isEBEtaGap())  RECOELE_isEBetaGap.push_back(1);  
    if (cand->isEBPhiGap())  RECOELE_isEBphiGap.push_back(1);  
    if (cand->isEEDeeGap())  RECOELE_isEEdeeGap.push_back(1);  
    if (cand->isEERingGap()) RECOELE_isEEringGap.push_back(1);
    /*std::cout << "--fiducial flags: " 
      << "  isEB="        << RECOELE_isbarrel[index] 
      << "  isEBetagap="  << RECOELE_isEBetaGap[index] 
      << "  isEBphigap="  << RECOELE_isEBphiGap[index] 
      << "  isEE="        << RECOELE_isendcap[index] 
      << "  isEEdeegap="  << RECOELE_isEEdeeGap[index] 
      << "  isEEringgap=" << RECOELE_isEEringGap[index] 
      << std::endl;*/
    // Shower shape
    RECOELE_sigmaIetaIeta.push_back(cand->sigmaIetaIeta()); 
    RECOELE_sigmaEtaEta.push_back(cand->sigmaEtaEta());
    RECOELE_e15.push_back(cand->e1x5());
    RECOELE_e25max.push_back(cand->e2x5Max());
    RECOELE_e55.push_back(cand->e5x5());
    RECOELE_he.push_back(cand->hadronicOverEm());
    // RECOELE_r9.push_back(cand->r9());
    /*std::cout << "--shower shape:" 
      << "  sigmaIetaIeta=" << RECOELE_sigmaIetaIeta[index] 
      << "  sigmaEtaEta=" << RECOELE_sigmaEtaEta[index]  
      << "  e15=" << RECOELE_e15[index] 
      << "  e25max=" << RECOELE_e25max[index] 
      << "  e55=" << RECOELE_e55[index]  
      << "  he=" << RECOELE_he[index] 
      << "  r9=" << RECOELE_r9[index] 
      << std::endl;*/
      
    // also variables on hcal
    // Isolation
    // Particle flow
    //RECOELE_mva.push_back(cand->mva();
    //std::cout << "--PF: mva = " << RECOELE_mva[index] << std::endl;
    // Brem & Classifaction
    RECOELE_fbrem.push_back(cand->fbrem());
    RECOELE_nbrems.push_back(cand->numberOfBrems());
    RECOELE_Class.push_back(cand->classification());
    if (useRECOformat) RECOELE_fbrem_mean.push_back(1. - cand->gsfTrack()->outerMomentum().R()/cand->gsfTrack()->innerMomentum().R());
    RECOELE_fbrem_mode.push_back(cand->fbrem());
    /*std::cout << "--brem & classification: fbrem/nbrems/Class/fbrem_mean/fbrem_mode = "
      << "  fbrem="      << RECOELE_fbrem[index]
      << "  nbrems="     << RECOELE_nbrems[index]
      << "  class="      << RECOELE_Class[index]
      << "  fbrem_mean=" << RECOELE_fbrem_mean[index]
      << "  fbrem_mode=" << RECOELE_fbrem_mode[index]
      << std::endl;*/
    // Corrections
    RECOELE_ecalEnergy.push_back(cand->ecalEnergy()); 
    /*        RECOELE_isEcalEnergyCorrected[index] = cand->isEcalEnergyCorrected(); */
    /*        RECOELE_ecalEnergyError[index]       = cand->ecalEnergyError(); */
    /*        RECOELE_isMomentumCorrected[index]   = cand->isMomentumCorrected(); */
    /*        RECOELE_trackMomentumError[index]    = cand->trackMomentumError(); */
    /*        RECOELE_electronMomentumError[index] = cand->electronMomentumError(); */
    
    // Seed Collection
    if (useRECOformat) {
      edm::RefToBase<TrajectorySeed> seed = cand->gsfTrack()->extra()->seedRef();
      reco::ElectronSeedRef MyS = seed.castTo<reco::ElectronSeedRef>();
      ele_seedSubdet2.push_back(int(MyS->subDet2()));
      if(fabs(MyS->dPhi2()) < 100.) ele_seedDphi2.push_back(double(MyS->dPhi2()));
      if(fabs(MyS->dRz2()) < 100.)  ele_seedDrz2.push_back(double(MyS->dRz2()));
      
      ele_seedSubdet1.push_back(int(MyS->subDet1()));
      if(fabs(MyS->dPhi1()) < 100.) ele_seedDphi1.push_back(double(MyS->dPhi1()));
      if(fabs(MyS->dRz1()) < 100.)  ele_seedDrz1.push_back(double(MyS->dRz1()));
    }
    
    //
    edm::Handle<edm::View<reco::GsfElectron> > elecoll;
    iEvent.getByToken(mvaElectronTag_, elecoll);
    int indextris=0;
    for (edm::View<reco::GsfElectron>::const_iterator candid = elecoll->begin(); candid!=elecoll->end(); ++candid) {
      edm::Ref<edm::View<reco::GsfElectron> > eletrackrefa(elecoll,indextris);
      if ( cand->sigmaEtaEta()==candid->sigmaEtaEta() && cand->sigmaIetaIeta()==candid->sigmaIetaIeta() ){
	RECOELE_mvaTrigV0.push_back((*mvaTrigV0map)[eletrackrefa]);
	RECOELE_mvaNonTrigV0.push_back((*mvaNonTrigV0map)[eletrackrefa]);
	/*std::cout << "BDT MVA eleID flag = "
	  << RECOELE_mvaTrigV0[index] << " "
	  << RECOELE_mvaNonTrigV0[index] << " "
	  << std::endl;*/
      }
      indextris++;
    }
    // Conversion Finder
    ConvMapDist.push_back((*conversiondistmap)[eletrackref]);
    ConvMapDcot.push_back((*conversiondcotmap)[eletrackref]);
    /*std::cout << "Conversion finder = "
      << ConvMapDist[index] << " " 
      << ConvMapDcot[index] << " " 
      << std::endl;*/
      
    // Matching
    if (fillMCTruth==true){
      int i=0;
      for ( reco::CandidateCollection::const_iterator hIter=CollEle->begin(); hIter!= CollEle->end(); ++hIter ){
	//cout << "Reco Electron with pT= " << hIter->pt() << " " << RECOELE_PT[index] << " " << fabs(hIter->pt()-RECOELE_PT[index]) << " and mass="<< hIter->mass()<< endl;
	if (fabs(hIter->pt()-cand->p4().pt())<0.01){
	  i=hIter-(CollEle->begin());
	  CandidateRef Ref( CollEle, i );
	  edm::Ref<std::vector<reco::GenParticle> > genrefEle = (*GenParticlesMatchEle)[Ref];
	  if (!genrefEle.isNull()){
	    cout << "GenElectron with pT= " << genrefEle->p4().pt() << " and mass="<< genrefEle->p4().mass()<< endl;
	    RECOELE_MatchingMCTruth.push_back(true);
	    RECOELE_MatchingMCpT.push_back(genrefEle->p4().pt());
	    RECOELE_MatchingMCEta.push_back(genrefEle->p4().eta());
	    RECOELE_MatchingMCPhi.push_back(genrefEle->p4().phi());	    
	  } 
	  else {
	    cout << "There is no a reference to a genElectron" << endl;
	  }
	}   
      }    
    }
    index ++;
  }    
}
  

void HZZ4LeptonsRootTree::fillAdditionalRECO(const edm::Event& iEvent)
{
  RECO_ZMM_MASS.clear();
  RECO_ZEE_MASS.clear();
  RECO_DiLep_MASS.clear();
  RECO_ZMM_PT.clear();
  RECO_ZEE_PT.clear();
  RECO_DiLep_PT.clear();
  RECO_ZMM_ETA.clear();
  RECO_ZEE_ETA.clear();
  RECO_DiLep_ETA.clear();
  RECO_ZMM_PHI.clear();
  RECO_ZEE_PHI.clear();
  RECO_DiLep_PHI.clear();
  RECO_ZMMss_MASS.clear();
  RECO_ZEEss_MASS.clear();
  RECO_ZEM_MASS.clear();
  RECO_ZMMss_PT.clear();
  RECO_ZEEss_PT.clear();
  RECO_ZEM_PT.clear();
  RECO_ZMMss_ETA.clear();
  RECO_ZEEss_ETA.clear();
  RECO_ZEM_ETA.clear();
  RECO_ZMMss_PHI.clear();
  RECO_ZEEss_PHI.clear();
  RECO_ZEM_PHI.clear();
  RECO_MMMM_MASS.clear();
  RECO_MMMM_PT.clear();
  RECO_MMMM_ETA.clear();
  RECO_MMMM_PHI.clear();
  RECO_MMMM_MASS_REFIT.clear();
  RECO_EEEE_MASS.clear();
  RECO_EEEE_PT.clear();
  RECO_EEEE_ETA.clear();
  RECO_EEEE_PHI.clear();
  RECO_EEEE_MASS_REFIT.clear();
  RECO_EEMM_MASS.clear();
  RECO_EEMM_PT.clear();
  RECO_EEMM_ETA.clear();
  RECO_EEMM_PHI.clear();
  RECO_EEMM_MASS_REFIT.clear();
  RECO_LLL0_MASS.clear();
  RECO_LLL1_MASS.clear();
  RECO_LLL2_MASS.clear();
  RECO_LLL3_MASS.clear();
  RECO_LLL0_PT.clear();
  RECO_LLL1_PT.clear();
  RECO_LLL2_PT.clear();
  RECO_LLL3_PT.clear();
  RECO_LLLl0_MASS.clear();
  RECO_LLLl1_MASS.clear();
  RECO_LLLl0_PT.clear();
  RECO_LLLl1_PT.clear();
  RECO_LLLL0ss_MASS.clear();
  RECO_LLLL0ss_PT.clear();
  RECO_LLLL1ss_MASS.clear();
  RECO_LLLL1ss_PT.clear();
  RECO_LLLL2ss_MASS.clear();
  RECO_LLLL2ss_PT.clear();
  //RECOcollNameLLLLssos_MASS.clear();
  //RECOcollNameLLLLssos_PT.clear();
  RECO_LLLL_MASS.clear();
  RECO_LLLL_PT.clear();
  RECO_LLLL_ETA.clear();
  RECO_LLLL_PHI.clear();
  // Matching ZtoMuMu
  RECOzMuMu_MatchingMCTruth.clear();
  RECOzMuMu_MatchingMCpT.clear();
  RECOzMuMu_MatchingMCmass.clear();
  RECOzMuMu_MatchingMCEta.clear();
  RECOzMuMu_MatchingMCPhi.clear();
  // Matching ZtoEE
  RECOzEE_MatchingMCTruth.clear();
  RECOzEE_MatchingMCpT.clear();
  RECOzEE_MatchingMCmass.clear();
  RECOzEE_MatchingMCEta.clear();
  RECOzEE_MatchingMCPhi.clear();
  // Matching HiggsToMMMM
  RECOHzzMMMM_MatchingMCTruth.clear();
  RECOHzzMMMM_MatchingMCpT.clear();
  RECOHzzMMMM_MatchingMCmass.clear();
  RECOHzzMMMM_MatchingMCEta.clear();
  RECOHzzMMMM_MatchingMCPhi.clear();
  // Matching HiggsToEEEE
  RECOHzzEEEE_MatchingMCTruth.clear();
  RECOHzzEEEE_MatchingMCpT.clear();
  RECOHzzEEEE_MatchingMCmass.clear();
  RECOHzzEEEE_MatchingMCEta.clear();
  RECOHzzEEEE_MatchingMCPhi.clear();
  // Matching HiggsToEEMM
  RECOHzzEEMM_MatchingMCTruth.clear();
  RECOHzzEEMM_MatchingMCpT.clear();
  RECOHzzEEMM_MatchingMCmass.clear();
  RECOHzzEEMM_MatchingMCEta.clear();
  RECOHzzEEMM_MatchingMCPhi.clear();
  //tmp candidate collections
  /*reco::CandidateCollection *leptonscands2e2mu_;
    reco::CandidateCollection *leptonscands2e2murf_;
    reco::CandidateCollection *leptonscands4mu_;
    reco::CandidateCollection *leptonscands4murf_;
    reco::CandidateCollection *leptonscands4e_;
    reco::CandidateCollection *leptonscands4erf_;*/
  
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
  
  /*leptonscands2e2mu_= new (CandidateCollection);
  leptonscands2e2murf_= new (CandidateCollection);
  leptonscands4mu_= new (CandidateCollection);
  leptonscands4murf_= new (CandidateCollection);
  leptonscands4e_= new (CandidateCollection);
  leptonscands4erf_= new (CandidateCollection);*/
  
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
  //Matching ZtoMuMu:
  edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchZMuMu;
  iEvent.getByLabel("goodZtoMuMuMCMatch", GenParticlesMatchZMuMu);
  if (GenParticlesMatchZMuMu.isValid() ){
    cout << "The matched map collection has size= " <<  GenParticlesMatchZMuMu->size() << endl;
  }
  //Matching ZtoEE:
  edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchZEE;
  iEvent.getByLabel("goodZtoEEMCMatch", GenParticlesMatchZEE);
  if (GenParticlesMatchZEE.isValid() ){
    cout << "The matched map collection has size= " <<  GenParticlesMatchZEE->size() << endl;
  }
  // di-leptons OS
  leptonscands_Z0->clear();
  leptonscands_Z1->clear();
  for (unsigned int i=0; i<RECOcollNameZ.size(); i++) {  
    edm::Handle<edm::View<Candidate> > CandidatesZ;
    iEvent.getByLabel(RECOcollNameZ.at(i), CandidatesZ); 
    int kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesZ->begin();cand != CandidatesZ->end(); ++ cand ) { 
      if (i==0) RECO_ZMM_MASS.push_back(cand->p4().mass());
      if (i==0) RECO_ZMM_PT.push_back(cand->p4().pt());
      if (i==0) RECO_ZMM_ETA.push_back(cand->p4().eta());
      if (i==0) RECO_ZMM_PHI.push_back(cand->p4().phi());
      if (i==1) RECO_ZEE_MASS.push_back(cand->p4().mass());
      if (i==1) RECO_ZEE_PT.push_back(cand->p4().pt());
      if (i==1) RECO_ZEE_ETA.push_back(cand->p4().eta());
      if (i==1) RECO_ZEE_PHI.push_back(cand->p4().phi());
      cout << "di-lepton candidate of type=" << RECOcollNameZ.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl; 
	if (i==0) {
	  RECO_ZMM_PT.push_back(cand->daughter(j)->p4().pt());
	  RECO_ZMM_ETA.push_back(cand->daughter(j)->p4().eta());
	  RECO_ZMM_PHI.push_back(cand->daughter(j)->p4().phi());
	  leptonscands_Z0->push_back( cand->daughter(j)->clone());
	}
	if (i==1) {
	  RECO_ZEE_PT.push_back(cand->daughter(j)->p4().pt());
	  RECO_ZEE_ETA.push_back(cand->daughter(j)->p4().eta());
	  RECO_ZEE_PHI.push_back(cand->daughter(j)->p4().phi());
	  leptonscands_Z1->push_back(cand->daughter(j)->clone());
	}
      }
      if (i==0 && fillMCTruth==true){
	// Matching ZtoMuMu
	edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesZ,kk);
	edm::Ref<std::vector<reco::GenParticle> > genrefZMuMu = (*GenParticlesMatchZMuMu)[Ref]; 
	if (!genrefZMuMu.isNull()){
	  cout << "Gen Z with pT= " << genrefZMuMu->p4().pt() << " and mass="<< genrefZMuMu->p4().mass() << endl;	  
	  RECOzMuMu_MatchingMCTruth.push_back(true);
	  RECOzMuMu_MatchingMCpT.push_back(genrefZMuMu->p4().pt());
	  RECOzMuMu_MatchingMCmass.push_back(genrefZMuMu->p4().mass());
	  RECOzMuMu_MatchingMCEta.push_back(genrefZMuMu->p4().eta());
	  RECOzMuMu_MatchingMCPhi.push_back(genrefZMuMu->p4().phi());
	}
      }
      if (i==1 && fillMCTruth==true){
	// Matching ZtoEE
	edm::Ref<edm::View<reco::Candidate> > Ref(CandidatesZ,kk);
	edm::Ref<std::vector<reco::GenParticle> > genrefZEE = (*GenParticlesMatchZEE)[Ref]; 
	if (!genrefZEE.isNull()){
	  cout << "Gen Z with pT= " << genrefZEE->p4().pt() << " and mass="<< genrefZEE->p4().mass() << endl;	  
	  RECOzEE_MatchingMCTruth.push_back(true);
	  RECOzEE_MatchingMCpT.push_back(genrefZEE->p4().pt());
	  RECOzEE_MatchingMCmass.push_back(genrefZEE->p4().mass());
	  RECOzEE_MatchingMCEta.push_back(genrefZEE->p4().eta());
	  RECOzEE_MatchingMCPhi.push_back(genrefZEE->p4().phi());
	}
      }
      kk++;
    }
  }
  // di-leptons SS and cross-leptons
  leptonscands_Zss0->clear();
  leptonscands_Zss1->clear();
  leptonscands_Zcross->clear();
  for (unsigned int i=0; i<RECOcollNameZss.size(); i++) {  
    edm::Handle<edm::View<Candidate> > CandidatesZss;
    iEvent.getByLabel(RECOcollNameZss.at(i), CandidatesZss); 
    int kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesZss->begin();cand != CandidatesZss->end(); ++ cand ) { 
      if (i==0)  RECO_ZMMss_MASS.push_back(cand->p4().mass());
      if (i==0)  RECO_ZMMss_PT.push_back(cand->p4().pt());
      if (i==0)  RECO_ZMMss_ETA.push_back(cand->p4().eta());
      if (i==0)  RECO_ZMMss_PHI.push_back(cand->p4().phi());
      
      if (i==1)  RECO_ZEEss_MASS.push_back(cand->p4().mass());
      if (i==1)  RECO_ZEEss_PT.push_back(cand->p4().pt());
      if (i==1)  RECO_ZEEss_ETA.push_back(cand->p4().eta());
      if (i==1)  RECO_ZEEss_PHI.push_back(cand->p4().phi());
      
      if (i==2)  RECO_ZEM_MASS.push_back(cand->p4().mass());
      if (i==2)  RECO_ZEM_PT.push_back(cand->p4().pt());
      if (i==2)  RECO_ZEM_ETA.push_back(cand->p4().eta());
      if (i==2)  RECO_ZEM_PHI.push_back(cand->p4().phi());
      
      cout << "di-lepton candidate of type=" << RECOcollNameZss.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl; 
	if (i==0) {
	  RECO_ZMMss_PT.push_back(cand->daughter(j)->p4().pt());
	  RECO_ZMMss_ETA.push_back(cand->daughter(j)->p4().eta());
	  RECO_ZMMss_PHI.push_back(cand->daughter(j)->p4().phi());
	  leptonscands_Zss0->push_back( cand->daughter(j)->clone());
	}
	if (i==1) {
	  RECO_ZEEss_PT.push_back(cand->daughter(j)->p4().pt());
	  RECO_ZEEss_ETA.push_back(cand->daughter(j)->p4().eta());
	  RECO_ZEEss_PHI.push_back(cand->daughter(j)->p4().phi());
	  leptonscands_Zss1->push_back( cand->daughter(j)->clone());
	}
	if (i==2) {
	  RECO_ZEM_PT.push_back(cand->daughter(j)->p4().pt());
	  RECO_ZEM_ETA.push_back(cand->daughter(j)->p4().eta());
	  RECO_ZEM_PHI.push_back(cand->daughter(j)->p4().phi());
	  leptonscands_Zcross->push_back( cand->daughter(j)->clone());	  
	}
      }
      kk++;
    }
  }
  // di-leptons ALL
  leptonscands_DiLep->clear();
  edm::Handle<edm::View<Candidate> > CandidatesDiLep;
  iEvent.getByLabel(RECOcollNameDiLep, CandidatesDiLep); 
  int kkk=0;
  for( edm::View<Candidate>::const_iterator cand = CandidatesDiLep->begin();cand != CandidatesDiLep->end(); ++ cand ) { 
    RECO_DiLep_MASS.push_back(cand->p4().mass());
    RECO_DiLep_PT.push_back(cand->p4().pt());
    RECO_DiLep_ETA.push_back(cand->p4().eta());
    RECO_DiLep_PHI.push_back(cand->p4().phi());
    cout << "di-lepton candidate of type=" << RECOcollNameDiLep.label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
    for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
      cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl;
      RECO_DiLep_PT.push_back(cand->daughter(j)->p4().pt());
      RECO_DiLep_ETA.push_back(cand->daughter(j)->p4().eta());
      RECO_DiLep_PHI.push_back(cand->daughter(j)->p4().phi());
      leptonscands_DiLep->push_back(cand->daughter(j)->clone());
    }
    kkk++;
  }
  // MuMuMuMu
  //int i=1; 
  leptonscands_MMMM->clear();
  edm::Handle<edm::View<Candidate> > CandidatesMMMM;
  iEvent.getByLabel(RECOcollNameMMMM.at(0), CandidatesMMMM);
  edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchHMMMM;
  iEvent.getByLabel("goodHiggsTozzToMMMMMCMatch", GenParticlesMatchHMMMM);
  if (GenParticlesMatchHMMMM.isValid()){
    cout << "The matched map collection has size= " <<  GenParticlesMatchHMMMM->size() << endl;
  }
  int kk=0;
  for( edm::View<Candidate>::const_iterator cand = CandidatesMMMM->begin();cand != CandidatesMMMM->end(); ++ cand ) { 
    RECO_MMMM_MASS.push_back(cand->p4().mass());
    RECO_MMMM_PT.push_back(cand->p4().pt());
    RECO_MMMM_ETA.push_back(cand->p4().eta());
    RECO_MMMM_PHI.push_back(cand->p4().phi());
    int l=0;
    //cout << "index" << i-1 << endl;
    for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
      RECO_MMMM_MASS.push_back(cand->daughter(j)->p4().mass());
      RECO_MMMM_PT.push_back(cand->daughter(j)->p4().pt());
      RECO_MMMM_ETA.push_back(cand->daughter(j)->p4().eta());
      RECO_MMMM_PHI.push_back(cand->daughter(j)->p4().phi());
      //cout << "index" << i+j <<endl;
      for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	RECO_MMMM_MASS.push_back(cand->daughter(j)->daughter(k)->p4().mass());
	RECO_MMMM_PT.push_back(cand->daughter(j)->daughter(k)->p4().pt());
	RECO_MMMM_ETA.push_back(cand->daughter(j)->daughter(k)->p4().eta());
	RECO_MMMM_PHI.push_back(cand->daughter(j)->daughter(k)->p4().phi());
	leptonscands_MMMM->push_back(cand->daughter(j)->daughter(k)->clone());
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
	RECOHzzMMMM_MatchingMCTruth.push_back(true);
	RECOHzzMMMM_MatchingMCpT.push_back(genrefHzzMMMM->p4().pt());
	RECOHzzMMMM_MatchingMCmass.push_back(genrefHzzMMMM->p4().mass());
	RECOHzzMMMM_MatchingMCEta.push_back(genrefHzzMMMM->p4().eta());
	RECOHzzMMMM_MatchingMCPhi.push_back(genrefHzzMMMM->p4().phi());
      }
    }
    kk++;
  }
  // EEEE
  //i=1;
  leptonscands_EEEE->clear();
  edm::Handle<edm::View<Candidate> > CandidatesEEEE;
  iEvent.getByLabel(RECOcollNameEEEE.at(0), CandidatesEEEE);
  edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchHEEEE;
  iEvent.getByLabel("goodHiggsTozzToEEEEMCMatch", GenParticlesMatchHEEEE);
  if (GenParticlesMatchHEEEE.isValid()){
    cout << "The matched map collection has size= " <<  GenParticlesMatchHEEEE->size() << endl;
  }
  kk=0;
  for( edm::View<Candidate>::const_iterator cand = CandidatesEEEE->begin();cand != CandidatesEEEE->end(); ++ cand ) {
    RECO_EEEE_MASS.push_back(cand->p4().mass());
    RECO_EEEE_PT.push_back(cand->p4().pt());
    RECO_EEEE_ETA.push_back(cand->p4().eta());
    RECO_EEEE_PHI.push_back(cand->p4().phi());
    int l=0;
    //cout << "index" << i-1 << endl;
    for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
      RECO_EEEE_MASS.push_back(cand->daughter(j)->p4().mass());
      RECO_EEEE_PT.push_back(cand->daughter(j)->p4().pt());
      RECO_EEEE_ETA.push_back(cand->daughter(j)->p4().eta());
      RECO_EEEE_PHI.push_back(cand->daughter(j)->p4().phi());
      //cout << "index" << i+j <<endl;
      for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	RECO_EEEE_MASS.push_back(cand->daughter(j)->daughter(k)->p4().mass());
	RECO_EEEE_PT.push_back(cand->daughter(j)->daughter(k)->p4().pt());
	RECO_EEEE_ETA.push_back(cand->daughter(j)->daughter(k)->p4().eta());
	RECO_EEEE_PHI.push_back(cand->daughter(j)->daughter(k)->p4().phi());
	leptonscands_EEEE->push_back(cand->daughter(j)->daughter(k)->clone());
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
	RECOHzzEEEE_MatchingMCTruth.push_back(true);
	RECOHzzEEEE_MatchingMCpT.push_back(genrefHzzEEEE->p4().pt());
	RECOHzzEEEE_MatchingMCmass.push_back(genrefHzzEEEE->p4().mass());
	RECOHzzEEEE_MatchingMCEta.push_back(genrefHzzEEEE->p4().eta());
	RECOHzzEEEE_MatchingMCPhi.push_back(genrefHzzEEEE->p4().phi());
      }
    }
    kk++;
  }
  // EEMM
  //i=1;
  leptonscands_EEMM->clear();
  edm::Handle<edm::View<Candidate> > CandidatesEEMM;
  iEvent.getByLabel(RECOcollNameEEMM.at(0), CandidatesEEMM);
  edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchHEEMM;
  iEvent.getByLabel("goodHiggsTozzToEEMMMCMatch", GenParticlesMatchHEEMM);
  if (GenParticlesMatchHEEMM.isValid()){
    cout << "The matched map collection has size= " <<  GenParticlesMatchHEEMM->size() << endl;
  }
  kk=0;
  for( edm::View<Candidate>::const_iterator cand = CandidatesEEMM->begin();cand != CandidatesEEMM->end(); ++ cand ) {
    RECO_EEMM_MASS.push_back(cand->p4().mass());
    RECO_EEMM_PT.push_back(cand->p4().pt());
    RECO_EEMM_ETA.push_back(cand->p4().eta());
    RECO_EEMM_PHI.push_back(cand->p4().phi());
    int l=0;
    //cout << "index" << i-1 << endl;
    for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
      RECO_EEMM_MASS.push_back(cand->daughter(j)->p4().mass());
      RECO_EEMM_PT.push_back(cand->daughter(j)->p4().pt());
      RECO_EEMM_ETA.push_back(cand->daughter(j)->p4().eta());
      RECO_EEMM_PHI.push_back(cand->daughter(j)->p4().phi());
      //cout << "index" << i+j <<endl;
      for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	RECO_EEMM_MASS.push_back(cand->daughter(j)->daughter(k)->p4().mass());
	RECO_EEMM_PT.push_back(cand->daughter(j)->daughter(k)->p4().pt());
	RECO_EEMM_ETA.push_back(cand->daughter(j)->daughter(k)->p4().eta());
	RECO_EEMM_PHI.push_back(cand->daughter(j)->daughter(k)->p4().phi());
	leptonscands_EEMM->push_back(cand->daughter(j)->daughter(k)->clone());
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
	RECOHzzEEMM_MatchingMCTruth.push_back(true);
	RECOHzzEEMM_MatchingMCpT.push_back(genrefHzzEEMM->p4().pt());
	RECOHzzEEMM_MatchingMCmass.push_back(genrefHzzEEMM->p4().mass());
	RECOHzzEEMM_MatchingMCEta.push_back(genrefHzzEEMM->p4().eta());
	RECOHzzEEMM_MatchingMCPhi.push_back(genrefHzzEEMM->p4().phi());
      }
    }
    kk++;
  }
  // tri-leptons
  leptonscands_LLL0->clear();
  leptonscands_LLL1->clear();
  leptonscands_LLL2->clear();
  leptonscands_LLL3->clear();
  for (unsigned int i=0; i<RECOcollNameLLL.size(); i++) {  
    edm::Handle<edm::View<Candidate> > CandidatesLLL;
    iEvent.getByLabel(RECOcollNameLLL.at(i), CandidatesLLL); 
    int k=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesLLL->begin();cand != CandidatesLLL->end(); ++ cand ) { 
      if (i==0) RECO_LLL0_MASS.push_back(cand->p4().mass());
      if (i==0) RECO_LLL0_PT.push_back(cand->p4().pt());
      if (i==1) RECO_LLL1_MASS.push_back(cand->p4().mass());
      if (i==1) RECO_LLL1_PT.push_back(cand->p4().pt());
      if (i==2) RECO_LLL2_MASS.push_back(cand->p4().mass());
      if (i==2) RECO_LLL2_PT.push_back(cand->p4().pt());
      if (i==3) RECO_LLL3_MASS.push_back(cand->p4().mass());
      if (i==4) RECO_LLL3_PT.push_back(cand->p4().pt());
      cout << "Tri-lepton candidate of type=" << RECOcollNameLLL.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	if (i==0) {
	  RECO_LLL0_PT.push_back(cand->daughter(j)->p4().pt());
	  leptonscands_LLL0->push_back( cand->daughter(j)->clone());
	}
	if (i==1) {
	  RECO_LLL1_PT.push_back(cand->daughter(j)->p4().pt());
	  leptonscands_LLL1->push_back(cand->daughter(j)->clone());
	}
	if (i==2) {
	  RECO_LLL2_PT.push_back(cand->daughter(j)->p4().pt());
	  leptonscands_LLL2->push_back(cand->daughter(j)->clone()); 
	}
	if (i==3) {
	  RECO_LLL3_PT.push_back(cand->daughter(j)->p4().pt());
	  leptonscands_LLL3->push_back(cand->daughter(j)->clone());
	}
      }
      k++;
    }
  }
  // 4-leptons SS
  leptonscands_LLLLss0->clear();
  leptonscands_LLLLss1->clear();
  leptonscands_LLLLss2->clear();
  for (unsigned int i=0; i<RECOcollNameLLLLss.size(); i++) {  
    edm::Handle<edm::View<Candidate> > CandidatesLLLLss;
    iEvent.getByLabel(RECOcollNameLLLLss.at(i), CandidatesLLLLss);
    int kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesLLLLss->begin();cand != CandidatesLLLLss->end(); ++ cand ) { 
      if (i==0) RECO_LLLL0ss_MASS.push_back(cand->p4().mass());
      if (i==0) RECO_LLLL0ss_PT.push_back(cand->p4().pt());
      if (i==1) RECO_LLLL1ss_MASS.push_back(cand->p4().mass());
      if (i==1) RECO_LLLL1ss_PT.push_back(cand->p4().pt());
      if (i==2) RECO_LLLL2ss_MASS.push_back(cand->p4().mass());
      if (i==2) RECO_LLLL2ss_PT.push_back(cand->p4().pt());
      cout << "4-lepton ss candidate of type=" << RECOcollNameLLLLss.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  if (i==0) {
	    if (j==0) RECO_LLLL0ss_PT.push_back(cand->daughter(j)->daughter(k)->p4().pt());
	    if (j==1) RECO_LLLL0ss_PT.push_back(cand->daughter(j)->daughter(k)->p4().pt());
	    leptonscands_LLLLss0->push_back(cand->daughter(j)->daughter(k)->clone());
	  }
	  if (i==1) {
	    if (j==0) RECO_LLLL1ss_PT.push_back(cand->daughter(j)->daughter(k)->p4().pt());
	    if (j==1) RECO_LLLL1ss_PT.push_back(cand->daughter(j)->daughter(k)->p4().pt());
	    leptonscands_LLLLss1->push_back(cand->daughter(j)->daughter(k)->clone());
	  }
	  if (i==2) {
	    if (j==0) RECO_LLLL2ss_PT.push_back(cand->daughter(j)->daughter(k)->p4().pt());
	    if (j==1) RECO_LLLL2ss_PT.push_back(cand->daughter(j)->daughter(k)->p4().pt());
	    leptonscands_LLLLss2->push_back(cand->daughter(j)->daughter(k)->clone());
	  }
	}
      }
      kk++;
    }      
  }
  // 4-leptons 3l+l
  leptonscands_LLLl0->clear();
  leptonscands_LLLl1->clear(); 
  for (unsigned int i=0; i<RECOcollNameLLLl.size(); i++){ 
    edm::Handle<edm::View<Candidate> > CandidatesLLLl;  
    iEvent.getByLabel(RECOcollNameLLLl.at(i), CandidatesLLLl);  
    int kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesLLLl->begin();cand != CandidatesLLLl->end(); ++ cand ) { 
      if (i==0) RECO_LLLl0_MASS.push_back(cand->p4().mass());
      if (i==0) RECO_LLLl0_PT.push_back(cand->p4().pt());
      if (i==1) RECO_LLLl1_MASS.push_back(cand->p4().mass());
      if (i==1) RECO_LLLl1_PT.push_back(cand->p4().pt());
      cout << "4-lepton from 3l+l candidate of type=" << RECOcollNameLLLl.at(i).label()
	   << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	//cout << "j= " << j << endl;
	if ( i==0 && j==1) {
	  RECO_LLLl0_PT.push_back(cand->daughter(j)->p4().pt());
	  leptonscands_LLLl0->push_back(cand->daughter(j)->clone());
	}
	if ( i==1 && j==1) {
	  RECO_LLLl1_PT.push_back(cand->daughter(j)->p4().pt()); 
	  leptonscands_LLLl1->push_back(cand->daughter(j)->clone());
	}
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) { 
	  //cout << "k= " << k << endl;
	  if (i==0) {
	    if (j==0) RECO_LLLl0_PT.push_back(cand->daughter(j)->daughter(k)->p4().pt());
	    leptonscands_LLLl0->push_back(cand->daughter(j)->daughter(k)->clone());
	  }
	  if (i==1) {
	    if (j==0) RECO_LLLl1_PT.push_back(cand->daughter(j)->daughter(k)->p4().pt());
	    leptonscands_LLLl1->push_back(cand->daughter(j)->daughter(k)->clone());
	  }
	} 
      } 
      kk++;
    }
  }
  // LLLL merged (no flavour, no charge)
  //i=1;
  leptonscands_LLLL->clear();
  edm::Handle<edm::View<Candidate> > CandidatesLLLL;
  iEvent.getByLabel(RECOcollNameLLLL, CandidatesLLLL);
  kk=0;    
  for( edm::View<Candidate>::const_iterator cand = CandidatesLLLL->begin();cand != CandidatesLLLL->end(); ++ cand ) {
    cout << "4lepton (any flavour and charge) candidate of type=" << RECOcollNameLLLL.label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
    
    RECO_LLLL_MASS.push_back(cand->p4().mass());
    RECO_LLLL_PT.push_back(cand->p4().pt());
    RECO_LLLL_ETA.push_back(cand->p4().eta());
    RECO_LLLL_PHI.push_back(cand->p4().phi());
    int l=0;
    //cout << "index" << i-1 << endl;
    for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
      RECO_LLLL_MASS.push_back(cand->daughter(j)->p4().mass());
      RECO_LLLL_PT.push_back(cand->daughter(j)->p4().pt());
      RECO_LLLL_ETA.push_back(cand->daughter(j)->p4().eta());
      RECO_LLLL_PHI.push_back(cand->daughter(j)->p4().phi());
      //cout << "index" << i+j <<endl;
      for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	RECO_LLLL_MASS.push_back(cand->daughter(j)->daughter(k)->p4().mass());
	RECO_LLLL_PT.push_back(cand->daughter(j)->daughter(k)->p4().pt());
	RECO_LLLL_ETA.push_back(cand->daughter(j)->daughter(k)->p4().eta());
	RECO_LLLL_PHI.push_back(cand->daughter(j)->daughter(k)->p4().phi());
	leptonscands_LLLL->push_back(cand->daughter(j)->daughter(k)->clone());
	//cout << "index" << i+j+k+l+2 <<endl;
      }
      l++;
    }
    kk++;
  }
}


void HZZ4LeptonsRootTree::fillP3Covariance(const reco::PFCandidate &c, TMatrixDSym &bigCov, int offset) const 
{
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
    } 
  } 
}
  
void HZZ4LeptonsRootTree::fillConstraintVtx2e2mu(const edm::Event& iEvent)
{
  StdFitVertexX.clear();
  StdFitVertexY.clear();
  StdFitVertexZ.clear();
  StdFitVertexChi2r.clear();
  StdFitVertexProb.clear();
  StdFitVertexTrack_PT.clear();
  StdFitVertexTrack_ETA.clear();
  StdFitVertexTrack_PHI.clear();
  KinFitVertexX.clear();
  KinFitVertexY.clear();
  KinFitVertexZ.clear();
  KinFitVertexChi2r.clear();
  KinFitVertexProb.clear();
  edm::Handle<reco::VertexCollection> StandardFitVtx_;
  iEvent.getByLabel(StandardFitVertex, StandardFitVtx_);
  int jjj=0;
  for (vector<reco::Vertex>::const_iterator cand=StandardFitVtx_->begin(); cand!=StandardFitVtx_->end(); ++cand){
    StdFitVertexX.push_back(cand->position().x());
    StdFitVertexY.push_back(cand->position().y());
    StdFitVertexZ.push_back(cand->position().z());
    StdFitVertexChi2r.push_back(cand->chi2()/cand->ndof());
    StdFitVertexProb.push_back(TMath::Prob(cand->chi2(),cand->ndof()));
    /*cout << "Std Fit: " 
      <<  StdFitVertexX[jjj] << " " 
      <<  StdFitVertexY[jjj] << " " 
      <<  StdFitVertexZ[jjj] << " " 
      <<  StdFitVertexChi2r[jjj] << " " 
      <<  StdFitVertexProb[jjj] << endl;*/
    // Refitted tracks     
    bool hasRefittedTracks = cand->hasRefittedTracks();
    if(hasRefittedTracks){	
      vector<reco::Track> refit_tks= cand->refittedTracks(); 	
      for(unsigned int i=0; i< refit_tks.size(); i++){	 
	if (i > 3) break;
	//cout << "Track momentum Refit=" << refit_tks[i].pt() << endl;
	StdFitVertexTrack_PT.push_back(refit_tks[i].pt());
	StdFitVertexTrack_ETA.push_back(refit_tks[i].eta());
	StdFitVertexTrack_PHI.push_back(refit_tks[i].phi());
      }
    }
    jjj++;
  }
  edm::Handle<edm::View<Candidate> > CandidatesEEMM;
  iEvent.getByLabel(RECOcollNameEEMM.at(0), CandidatesEEMM);
  edm::Handle<edm::ValueMap<float> > refmassmap;
  iEvent.getByLabel(RefittedMass, refmassmap);
  int kk=0;
  for( edm::View<Candidate>::const_iterator cand = CandidatesEEMM->begin();cand != CandidatesEEMM->end(); ++ cand ) {
    edm::Ref<edm::View<Candidate> > Ref(CandidatesEEMM,kk);
    //cout << "Original 2e2mu mass is= " << cand->p4().mass() << " Refitted mass is= " << (*refmassmap)[Ref] << endl;
    RECO_EEMM_MASS_REFIT.push_back((*refmassmap)[Ref]);
    kk++;
  }
  edm::Handle<reco::VertexCollection> KinematicFitVtx_;
  iEvent.getByLabel(KinematicFitVertex, KinematicFitVtx_);
  jjj=0;
  for (vector<reco::Vertex>::const_iterator cand=KinematicFitVtx_->begin(); cand!=KinematicFitVtx_->end(); ++cand){
    KinFitVertexX.push_back(cand->position().x());
    KinFitVertexY.push_back(cand->position().y());
    KinFitVertexZ.push_back(cand->position().z());
    KinFitVertexChi2r.push_back(cand->chi2()/cand->ndof());
    KinFitVertexProb.push_back(TMath::Prob(cand->chi2(),cand->ndof()));
    /*cout << "Kin Fit: " 
      <<  KinFitVertexX[jjj] << " " 
      <<  KinFitVertexY[jjj] << " " 
      <<  KinFitVertexZ[jjj] << " " 
      <<  KinFitVertexChi2r[jjj] << " " 
      <<  KinFitVertexProb[jjj] << endl;*/
    jjj++;
  }
}

void HZZ4LeptonsRootTree::fillConstraintVtx4mu(const edm::Event& iEvent)
{
  StdFitVertexXMMMM.clear();
  StdFitVertexYMMMM.clear();
  StdFitVertexZMMMM.clear();
  StdFitVertexChi2rMMMM.clear();
  StdFitVertexProbMMMM.clear();
  StdFitVertexTrackMMMM_PT.clear();
  StdFitVertexTrackMMMM_ETA.clear();
  StdFitVertexTrackMMMM_PHI.clear();
  KinFitVertexXMMMM.clear();
  KinFitVertexYMMMM.clear();
  KinFitVertexZMMMM.clear();
  KinFitVertexChi2rMMMM.clear();
  KinFitVertexProbMMMM.clear();
  edm::Handle<reco::VertexCollection> StandardFitVtx_;
  iEvent.getByLabel(StandardFitVertexMMMM, StandardFitVtx_);
  int jjj=0;
  for (vector<reco::Vertex>::const_iterator cand=StandardFitVtx_->begin(); cand!=StandardFitVtx_->end(); ++cand){
    StdFitVertexXMMMM.push_back(cand->position().x());
    StdFitVertexYMMMM.push_back(cand->position().y());
    StdFitVertexZMMMM.push_back(cand->position().z());
    StdFitVertexChi2rMMMM.push_back(cand->chi2()/cand->ndof());
    StdFitVertexProbMMMM.push_back(TMath::Prob(cand->chi2(),cand->ndof()));
    /*cout << "Std Fit MMMM: " 
      <<  StdFitVertexXMMMM[jjj] << " " 
      <<  StdFitVertexYMMMM[jjj] << " " 
      <<  StdFitVertexZMMMM[jjj] << " " 
      <<  StdFitVertexChi2rMMMM[jjj] << " " 
      <<  StdFitVertexProbMMMM[jjj] << endl;*/
    // Refitted tracks     
    bool hasRefittedTracks = cand->hasRefittedTracks();
    if(hasRefittedTracks){	
      vector<reco::Track> refit_tks= cand->refittedTracks(); 	
      for(unsigned int i=0; i< refit_tks.size(); i++){
	if (i > 3) break;
	//cout << "Track momentum Refit=" << refit_tks[i].pt() << endl;
	StdFitVertexTrackMMMM_PT.push_back(refit_tks[i].pt()) ;
	StdFitVertexTrackMMMM_ETA.push_back(refit_tks[i].eta()) ;
	StdFitVertexTrackMMMM_PHI.push_back(refit_tks[i].phi());
      }
    }
    jjj++;
  }
  edm::Handle<edm::View<Candidate> > CandidatesMMMM;
  iEvent.getByLabel(RECOcollNameMMMM.at(0), CandidatesMMMM);
  edm::Handle<edm::ValueMap<float> > refmassmap;
  iEvent.getByLabel(RefittedMassMMMM, refmassmap);
  int kk=0;
  for( edm::View<Candidate>::const_iterator cand = CandidatesMMMM->begin();cand != CandidatesMMMM->end(); ++ cand ) {
    edm::Ref<edm::View<Candidate> > Ref(CandidatesMMMM,kk);
    //cout << "Original 4mu mass is= " << cand->p4().mass() << " Refitted mass is= " << (*refmassmap)[Ref] << endl;
    RECO_MMMM_MASS_REFIT.push_back((*refmassmap)[Ref]);
    kk++;
  }
  edm::Handle<reco::VertexCollection> KinematicFitVtx_;
  iEvent.getByLabel(KinematicFitVertexMMMM, KinematicFitVtx_);
  jjj=0;
  for (vector<reco::Vertex>::const_iterator cand=KinematicFitVtx_->begin(); cand!=KinematicFitVtx_->end(); ++cand){
    KinFitVertexXMMMM.push_back(cand->position().x());
    KinFitVertexYMMMM.push_back(cand->position().y());
    KinFitVertexZMMMM.push_back(cand->position().z());
    KinFitVertexChi2rMMMM.push_back(cand->chi2()/cand->ndof());
    KinFitVertexProbMMMM.push_back(TMath::Prob(cand->chi2(),cand->ndof()));
    /*cout << "Kin Fit MMMM: " 
      <<  KinFitVertexXMMMM[jjj] << " " 
      <<  KinFitVertexYMMMM[jjj] << " " 
      <<  KinFitVertexZMMMM[jjj] << " " 
      <<  KinFitVertexChi2rMMMM[jjj] << " " 
      <<  KinFitVertexProbMMMM[jjj] << endl;*/
    jjj++;
  }
}

void HZZ4LeptonsRootTree::fillConstraintVtx4e(const edm::Event& iEvent)
{
  StdFitVertexXEEEE.clear();
  StdFitVertexYEEEE.clear();
  StdFitVertexZEEEE.clear();
  StdFitVertexChi2rEEEE.clear();
  StdFitVertexProbEEEE.clear();
  StdFitVertexTrackEEEE_PT.clear();
  StdFitVertexTrackEEEE_ETA.clear();
  StdFitVertexTrackEEEE_PHI.clear();
  KinFitVertexXEEEE.clear();
  KinFitVertexYEEEE.clear();
  KinFitVertexZEEEE.clear();
  KinFitVertexChi2rEEEE.clear();
  KinFitVertexProbEEEE.clear();
  edm::Handle<reco::VertexCollection> StandardFitVtx_;
  iEvent.getByLabel(StandardFitVertexEEEE, StandardFitVtx_);
  int jjj=0;
  for (vector<reco::Vertex>::const_iterator cand=StandardFitVtx_->begin(); cand!=StandardFitVtx_->end(); ++cand){
    StdFitVertexXEEEE.push_back(cand->position().x());
    StdFitVertexYEEEE.push_back(cand->position().y());
    StdFitVertexZEEEE.push_back(cand->position().z());
    StdFitVertexChi2rEEEE.push_back(cand->chi2()/cand->ndof());
    StdFitVertexProbEEEE.push_back(TMath::Prob(cand->chi2(),cand->ndof()));
    /*cout << "Std Fit EEEE: " 
      <<  StdFitVertexXEEEE[jjj] << " " 
      <<  StdFitVertexYEEEE[jjj] << " " 
      <<  StdFitVertexZEEEE[jjj] << " " 
      <<  StdFitVertexChi2rEEEE[jjj] << " " 
      <<  StdFitVertexProbEEEE[jjj] << endl;*/
    // Refitted tracks     
    bool hasRefittedTracks = cand->hasRefittedTracks();
    if(hasRefittedTracks){	
      vector<reco::Track> refit_tks= cand->refittedTracks(); 	
      for(unsigned int i=0; i< refit_tks.size(); i++){
	if (i > 3) break;
	//cout << "Track momentum Refit=" << refit_tks[i].pt() << endl;
	StdFitVertexTrackEEEE_PT.push_back(refit_tks[i].pt()) ;
	StdFitVertexTrackEEEE_ETA.push_back(refit_tks[i].eta()) ;
	StdFitVertexTrackEEEE_PHI.push_back(refit_tks[i].phi());
      }
    }
    jjj++;      
  }
  edm::Handle<edm::View<Candidate> > CandidatesEEEE;
  iEvent.getByLabel(RECOcollNameEEEE.at(0), CandidatesEEEE);
  edm::Handle<edm::ValueMap<float> > refmassmap;
  iEvent.getByLabel(RefittedMassEEEE, refmassmap);
  int kk=0;
  for( edm::View<Candidate>::const_iterator cand = CandidatesEEEE->begin();cand != CandidatesEEEE->end(); ++ cand ) {
    edm::Ref<edm::View<Candidate> > Ref(CandidatesEEEE,kk);
    cout << "Original 4e mass is= " << cand->p4().mass() << " Refitted mass is= " << (*refmassmap)[Ref] << endl;
    RECO_EEEE_MASS_REFIT.push_back((*refmassmap)[Ref]);
    kk++;
  }
  edm::Handle<reco::VertexCollection> KinematicFitVtx_;
  iEvent.getByLabel(KinematicFitVertexEEEE, KinematicFitVtx_);
  jjj=0;
  for (vector<reco::Vertex>::const_iterator cand=KinematicFitVtx_->begin(); cand!=KinematicFitVtx_->end(); ++cand){
    KinFitVertexXEEEE.push_back(cand->position().x());
    KinFitVertexYEEEE.push_back(cand->position().y());
    KinFitVertexZEEEE.push_back(cand->position().z());
    KinFitVertexChi2rEEEE.push_back(cand->chi2()/cand->ndof());
    KinFitVertexProbEEEE.push_back(TMath::Prob(cand->chi2(),cand->ndof()));
    /*cout << "Kin Fit EEEE: " 
      <<  KinFitVertexXEEEE[jjj] << " " 
      <<  KinFitVertexYEEEE[jjj] << " " 
      <<  KinFitVertexZEEEE[jjj] << " " 
      <<  KinFitVertexChi2rEEEE[jjj] << " " 
      <<  KinFitVertexProbEEEE[jjj] << endl;*/
    jjj++;
  }
}

void HZZ4LeptonsRootTree::fillConstraintVtxDiLeptons(const edm::Event& iEvent)
{
  StdFitVertexChi2rDiLep.clear();
  StdFitVertexProbDiLep.clear();
  edm::Handle<reco::VertexCollection> StandardFitVtxDiLep_;
  iEvent.getByLabel(StandardFitVertexDiLep, StandardFitVtxDiLep_);
  int jjj=0;
  for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxDiLep_->begin(); cand!=StandardFitVtxDiLep_->end(); ++cand){
    StdFitVertexChi2rDiLep.push_back(cand->chi2()/cand->ndof());
    StdFitVertexProbDiLep.push_back(TMath::Prob(cand->chi2(),cand->ndof()));
    /*cout << "Std Fit DiLeptons: " 
      <<  StdFitVertexChi2rDiLep[jjj] << " " 
      <<  StdFitVertexProbDiLep[jjj] << endl;*/
    jjj++;
  }
}

void HZZ4LeptonsRootTree::fillConstraintVtxTriLeptons(const edm::Event& iEvent)
{
  StdFitVertexChi2rMMM.clear();
  StdFitVertexProbMMM.clear();
  StdFitVertexChi2rMME.clear();
  StdFitVertexProbMME.clear();
  StdFitVertexChi2rEEE.clear();
  StdFitVertexProbEEE.clear();
  StdFitVertexChi2rMEE.clear();
  StdFitVertexProbMEE.clear();
  edm::Handle<reco::VertexCollection> StandardFitVtxMMM_;
  iEvent.getByLabel(StandardFitVertexMMM, StandardFitVtxMMM_);
  int jjj=0;
  for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxMMM_->begin(); cand!=StandardFitVtxMMM_->end(); ++cand){
    StdFitVertexChi2rMMM.push_back(cand->chi2()/cand->ndof());
    StdFitVertexProbMMM.push_back(TMath::Prob(cand->chi2(),cand->ndof()));
    /*cout << "Std Fit MMM: " 
      <<  StdFitVertexChi2rMMM[jjj] << " " 
      <<  StdFitVertexProbMMM[jjj] << endl;*/
    jjj++;
  }
  edm::Handle<reco::VertexCollection> StandardFitVtxMME_;
  iEvent.getByLabel(StandardFitVertexMME, StandardFitVtxMME_);
  jjj=0;
  for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxMME_->begin(); cand!=StandardFitVtxMME_->end(); ++cand){
    StdFitVertexChi2rMME.push_back(cand->chi2()/cand->ndof());
    StdFitVertexProbMME.push_back(TMath::Prob(cand->chi2(),cand->ndof()));
    /*cout << "Std Fit MME: " 
      <<  StdFitVertexChi2rMME[jjj] << " " 
      <<  StdFitVertexProbMME[jjj] << endl;*/
    jjj++;
  }
  edm::Handle<reco::VertexCollection> StandardFitVtxEEE_;
  iEvent.getByLabel(StandardFitVertexEEE, StandardFitVtxEEE_);
  jjj=0;
  for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxEEE_->begin(); cand!=StandardFitVtxEEE_->end(); ++cand){
    StdFitVertexChi2rEEE.push_back(cand->chi2()/cand->ndof());
    StdFitVertexProbEEE.push_back(TMath::Prob(cand->chi2(),cand->ndof()));
    /*cout << "Std Fit EEE: " 
      <<  StdFitVertexChi2rEEE[jjj] << " " 
      <<  StdFitVertexProbEEE[jjj] << endl;*/
    jjj++;
  }
  edm::Handle<reco::VertexCollection> StandardFitVtxMEE_;
  iEvent.getByLabel(StandardFitVertexMEE, StandardFitVtxMEE_);
  jjj=0;
  for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxMEE_->begin(); cand!=StandardFitVtxMEE_->end(); ++cand){
    StdFitVertexChi2rMEE.push_back(cand->chi2()/cand->ndof());
    StdFitVertexProbMEE.push_back(TMath::Prob(cand->chi2(),cand->ndof()));
    /*cout << "Std Fit MEE: " 
      <<  StdFitVertexChi2rMEE[jjj] << " " 
      <<  StdFitVertexProbMEE[jjj] << endl;*/
    jjj++;
  }
}

void HZZ4LeptonsRootTree::fillGD2e2mu(const edm::Event& iEvent)
{
  ftsigma.clear();
  gdX.clear();
  gdY.clear();
  gdZ.clear();
  ftsigmalag.clear();
  gdlagX.clear();
  gdlagY.clear();
  gdlagZ.clear();
  gdlagProb.clear();
  gdlagNdof.clear();
  edm::Handle<vector<double> > GeomD;
  iEvent.getByLabel(ftsigma_Vert, GeomD);
  int jjj=0;
  for (vector<double>::const_iterator cand=GeomD->begin(); cand!=GeomD->end(); ++cand){
    ftsigma.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > GeomDlag;
  iEvent.getByLabel(ftsigmalag_Vert, GeomDlag);
  jjj=0;
  for (vector<double>::const_iterator cand=GeomDlag->begin(); cand!=GeomDlag->end(); ++cand){
    ftsigmalag.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdX_;
  iEvent.getByLabel(gdX_Vert, gdX_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdX_->begin(); cand!=gdX_->end(); ++cand){
    gdX.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagX_;
  iEvent.getByLabel(gdlagX_Vert, gdlagX_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagX_->begin(); cand!=gdlagX_->end(); ++cand){
    gdlagX.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdY_;
  iEvent.getByLabel(gdY_Vert, gdY_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdY_->begin(); cand!=gdY_->end(); ++cand){
    gdY.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagY_;
  iEvent.getByLabel(gdlagY_Vert, gdlagY_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagY_->begin(); cand!=gdlagY_->end(); ++cand){
    gdlagY.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdZ_;
  iEvent.getByLabel(gdZ_Vert, gdZ_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdZ_->begin(); cand!=gdZ_->end(); ++cand){
    gdZ.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagZ_;
  iEvent.getByLabel(gdlagZ_Vert, gdlagZ_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagZ_->begin(); cand!=gdlagZ_->end(); ++cand){
    gdlagZ.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagProb_;
  iEvent.getByLabel(gdlagProb_Vert, gdlagProb_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagProb_->begin(); cand!=gdlagProb_->end(); ++cand){
    gdlagProb.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagNdof_;
  iEvent.getByLabel(gdlagNdof_Vert, gdlagNdof_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagNdof_->begin(); cand!=gdlagNdof_->end(); ++cand){
    gdlagNdof.push_back((*cand));
    jjj++;
  }
}

void HZZ4LeptonsRootTree::fillGD4mu(const edm::Event& iEvent)
{
  ftsigmaMMMM.clear();
  gdXMMMM.clear();
  gdYMMMM.clear();
  gdZMMMM.clear();
  ftsigmalagMMMM.clear();
  gdlagXMMMM.clear();
  gdlagYMMMM.clear();
  gdlagZMMMM.clear();
  gdlagProbMMMM.clear();
  gdlagNdofMMMM.clear();
  edm::Handle<vector<double> > GeomD_MMMM;
  iEvent.getByLabel(ftsigma_VertMMMM, GeomD_MMMM);
  int jjj=0;
  for (vector<double>::const_iterator cand=GeomD_MMMM->begin(); cand!=GeomD_MMMM->end(); ++cand){
    ftsigmaMMMM.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > GeomDlag_MMMM;
  iEvent.getByLabel(ftsigmalag_VertMMMM, GeomDlag_MMMM);
  jjj=0;
  for (vector<double>::const_iterator cand=GeomDlag_MMMM->begin(); cand!=GeomDlag_MMMM->end(); ++cand){
    ftsigmalagMMMM.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdXMMMM_;
  iEvent.getByLabel(gdX_VertMMMM, gdXMMMM_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdXMMMM_->begin(); cand!=gdXMMMM_->end(); ++cand){
    gdXMMMM.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagXMMMM_;
  iEvent.getByLabel(gdlagX_VertMMMM, gdlagXMMMM_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagXMMMM_->begin(); cand!=gdlagXMMMM_->end(); ++cand){
    gdlagXMMMM.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdYMMMM_;
  iEvent.getByLabel(gdY_VertMMMM, gdYMMMM_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdYMMMM_->begin(); cand!=gdYMMMM_->end(); ++cand){
    gdYMMMM.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagYMMMM_;
  iEvent.getByLabel(gdlagY_VertMMMM, gdlagYMMMM_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagYMMMM_->begin(); cand!=gdlagYMMMM_->end(); ++cand){
    gdlagYMMMM.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdZMMMM_;
  iEvent.getByLabel(gdZ_VertMMMM, gdZMMMM_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdZMMMM_->begin(); cand!=gdZMMMM_->end(); ++cand){
    gdZMMMM.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagZMMMM_;
  iEvent.getByLabel(gdlagZ_VertMMMM, gdlagZMMMM_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagZMMMM_->begin(); cand!=gdlagZMMMM_->end(); ++cand){
    gdlagZMMMM.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagProbMMMM_;
  iEvent.getByLabel(gdlagProb_VertMMMM, gdlagProbMMMM_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagProbMMMM_->begin(); cand!=gdlagProbMMMM_->end(); ++cand){
    gdlagProbMMMM.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagNdofMMMM_;
  iEvent.getByLabel(gdlagNdof_VertMMMM, gdlagNdofMMMM_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagNdofMMMM_->begin(); cand!=gdlagNdofMMMM_->end(); ++cand){
    gdlagNdofMMMM.push_back((*cand));
    jjj++;
  }
}

void HZZ4LeptonsRootTree::fillGD4e(const edm::Event& iEvent)
{
  ftsigmaEEEE.clear();
  gdXEEEE.clear();
  gdYEEEE.clear();
  gdZEEEE.clear();
  ftsigmalagEEEE.clear();
  gdlagXEEEE.clear();
  gdlagYEEEE.clear();
  gdlagZEEEE.clear();
  gdlagProbEEEE.clear();
  gdlagNdofEEEE.clear();
  edm::Handle<vector<double> > GeomD_EEEE;
  iEvent.getByLabel(ftsigma_VertEEEE, GeomD_EEEE);
  int jjj=0;
  for (vector<double>::const_iterator cand=GeomD_EEEE->begin(); cand!=GeomD_EEEE->end(); ++cand){
    ftsigmaEEEE.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > GeomDlag_EEEE;
  iEvent.getByLabel(ftsigmalag_VertEEEE, GeomDlag_EEEE);
  jjj=0;
  for (vector<double>::const_iterator cand=GeomDlag_EEEE->begin(); cand!=GeomDlag_EEEE->end(); ++cand){
    ftsigmalagEEEE.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdXEEEE_;
  iEvent.getByLabel(gdX_VertEEEE, gdXEEEE_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdXEEEE_->begin(); cand!=gdXEEEE_->end(); ++cand){
    gdXEEEE.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagXEEEE_;
  iEvent.getByLabel(gdlagX_VertEEEE, gdlagXEEEE_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagXEEEE_->begin(); cand!=gdlagXEEEE_->end(); ++cand){
    gdlagXEEEE.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdYEEEE_;
  iEvent.getByLabel(gdY_VertEEEE, gdYEEEE_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdYEEEE_->begin(); cand!=gdYEEEE_->end(); ++cand){
    gdYEEEE.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagYEEEE_;
  iEvent.getByLabel(gdlagY_VertEEEE, gdlagYEEEE_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagYEEEE_->begin(); cand!=gdlagYEEEE_->end(); ++cand){
    gdlagYEEEE.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdZEEEE_;
  iEvent.getByLabel(gdZ_VertEEEE, gdZEEEE_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdZEEEE_->begin(); cand!=gdZEEEE_->end(); ++cand){
    gdZEEEE.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagZEEEE_;
  iEvent.getByLabel(gdlagZ_VertEEEE, gdlagZEEEE_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagZEEEE_->begin(); cand!=gdlagZEEEE_->end(); ++cand){
    gdlagZEEEE.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagProbEEEE_;
  iEvent.getByLabel(gdlagProb_VertEEEE, gdlagProbEEEE_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagProbEEEE_->begin(); cand!=gdlagProbEEEE_->end(); ++cand){
    gdlagProbEEEE.push_back((*cand));
    jjj++;
  }
  edm::Handle<vector<double> > gdlagNdofEEEE_;
  iEvent.getByLabel(gdlagNdof_VertEEEE, gdlagNdofEEEE_);
  jjj=0;
  for (vector<double>::const_iterator cand=gdlagNdofEEEE_->begin(); cand!=gdlagNdofEEEE_->end(); ++cand){
    gdlagNdofEEEE.push_back((*cand));
    jjj++;
  }
}

  
//Vertices
void HZZ4LeptonsRootTree::fillVertices(const edm::Event& iEvent)
{
  RECO_NVTX=0;
  RECO_VERTEX_x.clear();
  RECO_VERTEX_y.clear();
  RECO_VERTEX_z.clear();
  RECO_VERTEX_ndof.clear();
  RECO_VERTEX_chi2.clear();
  RECO_VERTEX_ntracks.clear();
  RECO_VERTEXPROB.clear();
  RECO_VERTEX_TRACK_PT.clear();
  RECO_VERTEX_isValid.clear();

  RECO_NVTX=PV.size();
  //cout << "Number of Vertices in the event= " << RECO_NVTX << endl;
  int index_vertex = 0;
  for (VertexCollection::const_iterator i=PV.begin(); i!=PV.end();i++) {
    //RECO_NVTX.push_back(index_vertex);
    RECO_VERTEX_x.push_back(i->x());
    RECO_VERTEX_y.push_back(i->y());
    RECO_VERTEX_z.push_back(i->z());
    RECO_VERTEX_ndof.push_back(i->ndof());
    RECO_VERTEX_chi2.push_back(i->chi2());
    RECO_VERTEX_ntracks.push_back(i->tracksSize());
    RECO_VERTEXPROB.push_back(ChiSquaredProbability(i->chi2(),i->ndof()));
    RECO_VERTEX_isValid.push_back(i->isValid());
    cout << "Vertex made by " << i->tracksSize() << " tracks with chi2="<< i->chi2() << " and ndof=" << i->ndof() << " and prob=" << RECO_VERTEXPROB[index_vertex] << endl;
    int indice=0;
    for(std::vector<reco::TrackBaseRef>::const_iterator iter = i->tracks_begin();
	iter != i->tracks_end(); iter++) {
      cout << "pT of tracks building the vertex= " << (**iter).pt() << endl; 
      RECO_VERTEX_TRACK_PT.push_back((**iter).pt());
      indice++;
    }
    index_vertex++;
  } // loop on vertices
}

//Tracks
void HZZ4LeptonsRootTree::fillTracks(const edm::Event& iEvent)
{
  RECO_NTRACK=0;
  RECO_TRACK_PT.clear();
  RECO_TRACK_ETA.clear();
  RECO_TRACK_PHI.clear();
  RECO_TRACK_CHI2.clear();
  RECO_TRACK_CHI2RED.clear();
  RECO_TRACK_CHI2PROB.clear();
  RECO_TRACK_NHITS.clear();
  RECO_TRACK_DXY.clear();
  RECO_TRACK_DXYERR.clear();
  RECO_TRACK_DZ.clear();
  RECO_TRACK_DZERR.clear();
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(tracksTag_, tracks);
  
  RECO_NTRACK=tracks->size();
  //
  cout << "Number of Tracks in the event= " << RECO_NTRACK << endl;
  int countk=0;
  for ( TrackCollection::const_iterator i=tracks->begin(); i!=tracks->end(); i++) { 
    //RECO_NTRACK=countk;
    RECO_TRACK_PT.push_back(i->pt());
    RECO_TRACK_ETA.push_back(i->eta());
    RECO_TRACK_PHI.push_back(i->phi());
    RECO_TRACK_CHI2.push_back(i->chi2());
    RECO_TRACK_CHI2RED.push_back(i->normalizedChi2());
    //RECO_TRACK_CHI2PROB.push_back(TMath::Prob(i->chi2(),i->ndof()));
    RECO_TRACK_CHI2PROB.push_back(ChiSquaredProbability(i->chi2(),i->ndof()));
    RECO_TRACK_NHITS.push_back(i->numberOfValidHits());
    RECO_TRACK_DXY.push_back(i->dxy());
    RECO_TRACK_DXYERR.push_back(i->dxyError());
    RECO_TRACK_DZ.push_back(i->dz());
    RECO_TRACK_DZERR.push_back(i->dzError());
    countk++;
  }
}

