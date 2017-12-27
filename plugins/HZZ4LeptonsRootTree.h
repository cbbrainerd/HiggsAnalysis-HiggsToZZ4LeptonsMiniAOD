/** \class  HZZ4LeptonsRootTree
 *
 *  Root Tree for H->ZZ->4l analysis.
 *
 *  Author: N. De Filippis - Politecnico and INFN Bari
 *
 *  Modified (from array to vector format) by Sherif Elgammal / Nicola De Filippis / Giorgia Miniello
 *  11/3/2016    
 *  CMSSW_7_6_X     
 */

#ifndef MyCodeArea_Analyzer_HZZ4LeptonsRootTree_h
#define MyCodeArea_Analyzer_HZZ4LeptonsRootTree_h
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include <FWCore/Framework/interface/ESHandle.h>
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//====================== Get PixelMatchGsfElectronCollection =========
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
//====================== Get PixelMatchGsfElectronCollection =========
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
//#include "DQMOffline/JetMET/interface/CaloMETAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/METReco/interface/SpecificCaloMETData.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CorrMETData.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
//==================== start a part for Jet ================================
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
//==================== end a part for Jet ================================
//==================== start a part for photons ==========================
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
//#include "tmp/CaloGeometryBuilder/interface/CaloGeometryBuilder.h"
#include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Common/interface/RefCore.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/TypeID.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
 
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "CommonTools/TriggerUtils/interface/GenericTriggerEventFlag.h"
 
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h" 
#include "DataFormats/MuonReco/interface/MuonEnergy.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

// Trigger
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "CommonTools/ParticleFlow/interface/PFMETAlgo.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//====================== end a part for photons ==========================
#include <vector>
#include <set>
#include <stdio.h>
#include "TProfile.h"
#include "TProfile2D.h"
//#include <iostream.h>
//#include <ostream.h>
//#include <fstream.h>
//#include <vector.h>
#include "TFile.h"
#include <math.h>
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "Math/LorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"


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


//class TFile;

class MultiTrajectoryStateMode ;
class EgammaTowerIsolation ;


// Class to create TTree variables
#include <TFile.h> 
#include <TTree.h> 

// Namespaces
using namespace reco;
using namespace std;
using namespace pat;

//
// class declaration
//


class HZZ4LeptonsRootTree : public edm::EDAnalyzer {
   public:
      explicit HZZ4LeptonsRootTree( const edm::ParameterSet& );  //constructor
      ~HZZ4LeptonsRootTree();  //distructor
      virtual void analyze( const edm::Event&, const edm::EventSetup& );
      virtual void beginJob();
      virtual void endJob();
      void IntialValues();
      void EventsReWeighting(const edm::Event& evt);
      float delR(float eta1,float phi1,float eta2,float phi2);
      void triggermatching(const edm::Event& iEvent);
      bool IsMuMatchedToHLTMu(const reco::Muon &mu,std::vector<reco::Particle> HLTMu,std::vector<string> HLTMuNames,double DR,double DPtRel);
      bool IsEleMatchedToHLTEle(const reco::GsfElectron &ele,std::vector<reco::Particle> HLTEle,std::vector<string> HLTEleNames,double DR,double DPtRel);
      const std::string getParticleName(int id);
      bool match(double mass, double pt, int charge,const reco::CandidateCollection *c1Coll);
      bool matchParticle(double mass, double pt, int charge, const reco::Candidate *c1);
      struct SortCandByDecreasingPt {
	bool operator()( const reco::Candidate &c1, const reco::Candidate &c2) const {
	  return c1.pt() > c2.pt();
	}
      };
      void SetMCValues(const reco::Candidate& cand, int nMC);
      void fillgenparticles(const edm::Event& iEvent, const edm::EventSetup &es);
      void fillmc(const edm::Event& iEvent);
      void fillBTagging(const edm::Event& iEvent);
      void fillMET(const edm::Event& evt);
      void fillPU(const edm::Event& iEvent);
      void fillRho(const edm::Event& evt);
      void fillPhotons(const edm::Event& evt);
      void fillMuons(const edm::Event& iEvent,const edm::EventSetup& iSetup);
      void fillElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      void EventsMCReWeighting(const edm::Event& iEvent);
      void fillHLTFired(const edm::Event& iEvent);
      void filljets(const edm::Event& iEvent);
      void fillAdditionalRECO(const edm::Event& iEvent);
      void fillP3Covariance(const reco::PFCandidate &c, TMatrixDSym &bigCov, int offset) const;
      void fillConstraintVtx2e2mu(const edm::Event& iEvent);
      void fillConstraintVtx4e(const edm::Event& iEvent);
      void fillConstraintVtx4mu(const edm::Event& iEvent);
      void fillConstraintVtxDiLeptons(const edm::Event& iEvent);
      void fillConstraintVtxTriLeptons(const edm::Event& iEvent);
      void fillGD2e2mu(const edm::Event& iEvent);
      void fillGD4e(const edm::Event& iEvent);
      void fillGD4mu(const edm::Event& iEvent);
      void fillVertices(const edm::Event& iEvent);
      void fillTracks(const edm::Event& iEvent);
private:

      TTree* mytree;
      //====================================================
      //
      //    Create variable for Nb of event,run,lumi
      //
      //====================================================
      TH1F* h1_DaughterID_;

      GlobalPoint eventVertexPosition_;
      int Run;
      int Event;
      int lumi;
      int bunch;
      int event_runNo;
      int event_evtNo;
      int event_lumi;
      int event_bunch;
      float EffectiveArea;
      //==================================================
      //
      //     Create Branches for Vertices variables
      //
      //==================================================
      int RECO_NVTX;
      std::vector<float> RECO_VERTEX_x;
      std::vector<float> RECO_VERTEX_y;
      std::vector<float> RECO_VERTEX_z;
      std::vector<float> RECO_VERTEX_ndof;
      std::vector<float> RECO_VERTEX_chi2;
      std::vector<float> RECO_VERTEX_ntracks;
      std::vector<float> RECO_VERTEXPROB;
      std::vector<float> RECO_VERTEX_TRACK_PT;
      std::vector<bool> RECO_VERTEX_isValid;
      //=============================================================
      //                   
      //           Create Branches for Tracks variables
      //
      //=============================================================
      int RECO_NTRACK;
      std::vector<float> RECO_TRACK_PT;
      std::vector<float> RECO_TRACK_ETA;
      std::vector<float> RECO_TRACK_PHI;
      std::vector<float> RECO_TRACK_CHI2;
      std::vector<float> RECO_TRACK_CHI2RED;
      std::vector<float> RECO_TRACK_CHI2PROB;
      std::vector<int> RECO_TRACK_NHITS;
      std::vector<float> RECO_TRACK_DXY;
      std::vector<float> RECO_TRACK_DXYERR;
      std::vector<float> RECO_TRACK_DZ;
      std::vector<float> RECO_TRACK_DZERR;
      //==================================================
      //
      //           Create vectors for Muons variables
      //
      //==================================================
      int RECO_NMU;
      std::vector<float> RECOMU_isPFMu;
      std::vector<float> RECOMU_isGlobalMu;
      std::vector<float> RECOMU_isStandAloneMu;
      std::vector<float> RECOMU_isTrackerMu;
      std::vector<float> RECOMU_isCaloMu;
      // Kinematic of muon
      std::vector<float> RECOMU_E;
      std::vector<float> RECOMU_PT;
      std::vector<float> RECOMU_P;
      std::vector<float> RECOMU_ETA;
      std::vector<float> RECOMU_THETA;
      std::vector<float> RECOMU_PHI;
      std::vector<float> RECOMU_MASS;
      std::vector<int> RECOMU_CHARGE;
      // Covariance matrix
      std::vector<float> RECOMU_COV;
      // Isolation
      std::vector<float> RECOMU_TRACKISO;
      std::vector<float> RECOMU_TRACKISO_SUMPT;
      //temporary solution for reducedrechit problem
      std::vector<float> RECOMU_ECALISO;
      std::vector<float> RECOMU_HCALISO;
      std::vector<float> RECOMU_X;
      std::vector<float> RECOMU_PFchHad;
      std::vector<float> RECOMU_PFneuHad;
      std::vector<float> RECOMU_PFphoton;
      std::vector<float> RECOMU_PFPUchAllPart;
      // Vertexing
      std::vector<float> RECOMU_SIP;
      std::vector<float> RECOMU_IP;
      std::vector<float> RECOMU_IPERROR;
      std::vector<float> RECOMU_SIP_KF;
      std::vector<float> RECOMU_IP_KF;
      std::vector<float> RECOMU_IPERROR_KF;
      std::vector<float> RECOMU_PFX_dB;
      std::vector<float> RECOMU_PFX_rho;
      std::vector<float> RECOMU_STIP;
      std::vector<float> RECOMU_SLIP;
      std::vector<float> RECOMU_TIP;
      std::vector<float> RECOMU_LIP;
      std::vector<float> RECOMU_TIPERROR;
      std::vector<float> RECOMU_LIPERROR;
      // Other properties
      std::vector<float> RECOMU_numberOfMatches;
      std::vector<float> RECOMU_numberOfMatchedStations;
      std::vector<float> RECOMU_caloCompatibility;
      std::vector<float> RECOMU_segmentCompatibility;
      std::vector<float> RECOMU_glbmuPromptTight;
      // Track properties
      std::vector<float> RECOMU_mubesttrkType;
      std::vector<float> RECOMU_mubesttrkDxy;
      std::vector<float> RECOMU_mubesttrkDxyB;
      std::vector<float> RECOMU_mubesttrkDxyError;
      std::vector<float> RECOMU_mubesttrkDz;
      std::vector<float> RECOMU_mubesttrkDzB;
      std::vector<float> RECOMU_mubesttrkDzError;
      std::vector<float> RECOMU_mubesttrkPTError;

      std::vector<float> RECOMU_mutrkPT;
      std::vector<float> RECOMU_mutrkPTError;
      //std::vector<float> RECOMU_mutrkPTError;
      std::vector<float> RECOMU_mutrkDxy;
      std::vector<float> RECOMU_mutrkDxyError;
      std::vector<float> RECOMU_mutrkDxyB;
      std::vector<float> RECOMU_mutrkDz;
      std::vector<float> RECOMU_mutrkDzError;
      std::vector<float> RECOMU_mutrkDzB;
      std::vector<float> RECOMU_mutrkChi2PerNdof;
      std::vector<int> RECOMU_mutrkCharge;
      std::vector<int> RECOMU_mutrkNHits;
      std::vector<int> RECOMU_mutrkNPixHits;
      std::vector<int> RECOMU_mutrkNStripHits;
      std::vector<int> RECOMU_mutrkNMuonHits;
      std::vector<int> RECOMU_mutrktrackerLayersWithMeasurement;
      
      std::vector<float> RECOMU_muInnertrkDxy;
      std::vector<float> RECOMU_muInnertrkDxyError;
      std::vector<float> RECOMU_muInnertrkDxyB;
      std::vector<float> RECOMU_muInnertrkDz;
      std::vector<float> RECOMU_muInnertrkDzError;
      std::vector<float> RECOMU_muInnertrkDzB;
      std::vector<float> RECOMU_muInnertrkChi2PerNdof;
      std::vector<float> RECOMU_muInnertrktrackerLayersWithMeasurement;
      std::vector<float> RECOMU_muInnertrkPT;
      //std::vector<float> RECOMU_muInnertrkPTError;
      std::vector<float> RECOMU_muInnertrkPTError;
      std::vector<int> RECOMU_muInnertrkCharge;
      std::vector<int> RECOMU_muInnertrkNHits;
      std::vector<int> RECOMU_muInnertrkNPixHits;
      std::vector<int> RECOMU_muInnertrkNStripHits;
      
      // Tracker muon properties
      // std::vector<float> RECOMU_trkmuArbitration;
      std::vector<float> RECOMU_trkmuArbitration;
      std::vector<float> RECOMU_trkmu2DCompatibilityLoose;
      std::vector<float> RECOMU_trkmu2DCompatibilityTight;
      std::vector<float> RECOMU_trkmuOneStationLoose;
      std::vector<float> RECOMU_trkmuOneStationTight;
      std::vector<float> RECOMU_trkmuLastStationLoose;
      std::vector<float> RECOMU_trkmuLastStationTight;
      std::vector<float> RECOMU_trkmuOneStationAngLoose;
      std::vector<float> RECOMU_trkmuOneStationAngTight;
      std::vector<float> RECOMU_trkmuLastStationAngLoose;
      std::vector<float> RECOMU_trkmuLastStationAngTight;
      std::vector<float> RECOMU_trkmuLastStationOptimizedLowPtLoose;
      std::vector<float> RECOMU_trkmuLastStationOptimizedLowPtTight;

      std::vector<int> RECOMU_MatchingMCTruth;
      std::vector<float> RECOMU_MatchingMCpT;
      std::vector<float> RECOMU_MatchingMCEta;
      std::vector<float> RECOMU_MatchingMCPhi;

      //===================================================  
      //
      //    Create vectors for gen particles variables
      //
      //===================================================
      int kk,ii,l;
      std::vector<float> MC_LEPT_PT;
      std::vector<float> MC_LEPT_ETA;
      std::vector<float> MC_LEPT_PHI;
      std::vector<float> MC_LEPT_THETA;
      std::vector<int> MC_LEPT_PDGID;
      std::vector<float> MC_Z_PT;
      std::vector<float> MC_Z_ETA;
      std::vector<float> MC_Z_PHI;
      std::vector<float> MC_Z_THETA;
      std::vector<float> MC_Z_MASS;
      std::vector<int> MC_Z_PDGID;
      std::vector<float> MC_fourl_MASS;
      std::vector<float> MC_fourl_PT;
      std::vector<int> MC_fourl_PDGID;
      std::vector<float> MC_ZZ_MASS;
      std::vector<float> MC_ZZ_PT;
      std::vector<float> MC_ZZ_ETA;
      std::vector<float> MC_ZZ_PHI;
      std::vector<float> MC_ZZ_THETA;
      std::vector<int> MC_ZZ_PDGID;
      std::vector<float> MC_E;
      std::vector<float> MC_PT;
      std::vector<float> MC_ETA;
      std::vector<float> MC_THETA;
      std::vector<float> MC_PHI;
      std::vector<float> MC_MASS;
      std::vector<int> MC_PDGID;
      //====================================================
      //
      //   Create vectors for Pimary Vertice variables
      //
      //====================================================
      int value3_;
      std::vector<int> nbPv;
      std::vector<int> Nbdof;
      std::vector<float> PositionX;
      std::vector<float> PositionY;
      std::vector<float> PositionZ;
      std::vector<float> PositionRho;
      //====================================================
      //
      //   Create vectors for Jets variables
      //
      //====================================================
      double RHO,RHO_ele,RHO_mu;
      int RECO_PFJET_N;
      std::vector<int> RECO_PFJET_CHARGE;
      std::vector<float> RECO_PFJET_ET;
      std::vector<float> RECO_PFJET_PT;
      std::vector<float> RECO_PFJET_ETA;
      std::vector<float> RECO_PFJET_PHI;
      std::vector<int> RECO_PFJET_PUID;
      std::vector<float> RECO_PFJET_PUID_MVA;
      
      //=============================================================
      //
      //  Create vectors for Photons Tree
      //
      //=============================================================
      int iphot;
      int RECO_NPHOT;
      std::vector<float> RECOPHOT_PT;
      std::vector<float> RECOPHOT_ETA;
      std::vector<float> RECOPHOT_PHI;
      std::vector<float> RECOPHOT_THETA;
      int RECO_NPFPHOT;
      std::vector<float> RECOPFPHOT_PT;
      std::vector<float> RECOPFPHOT_PTError;
      std::vector<float> RECOPFPHOT_ETA;
      std::vector<float> RECOPFPHOT_PHI;
      std::vector<float> RECOPFPHOT_THETA;
      std::vector<float> RECOPFPHOT_PFchAllPart;
      std::vector<float> RECOPFPHOT_PFchHad;
      std::vector<float> RECOPFPHOT_PFneuHad;
      std::vector<float> RECOPFPHOT_PFphoton;
      std::vector<float> RECOPFPHOT_PFPUchAllPart;
      std::vector<float> RECOPFPHOT_PFX_rho;

      std::vector<int> RECOPHOT_MatchingMCTruth;
      std::vector<float> RECOPHOT_MatchingMCpT;
      std::vector<float> RECOPHOT_MatchingMCEta;
      std::vector<float> RECOPHOT_MatchingMCPhi;
      //=============================================================
      //
      //            Method for BTagging Tree
      //
      //=============================================================
      int iBtag;
      std::vector<int> tCHighEff_nb;
      std::vector<float> tCHighEff_BTagJet_PT;
      std::vector<float> tCHighEff_BTagJet_ETA;
      std::vector<float> tCHighEff_BTagJet_PHI;
      std::vector<float> tCHighEff_BTagJet_DISCR;
      std::vector<int> tCHighPur_nb;
      std::vector<float> tCHighPur_BTagJet_PT;
      std::vector<float> tCHighPur_BTagJet_ETA;
      std::vector<float> tCHighPur_BTagJet_PHI;
      std::vector<float> tCHighPur_BTagJet_DISCR;
      std::vector<int> cSV_nb;
      std::vector<float> cSV_BTagJet_PT;
      std::vector<float> cSV_BTagJet_ETA;
      std::vector<float> cSV_BTagJet_PHI;
      std::vector<float> cSV_BTagJet_DISCR;
      std::vector<float> cSV_BTagJet_ET;
      //=============================================================
      //
      //           Create Branchs for Muons match HLT variables
      //
      //=============================================================
      char HLTPathsFired[20000];
      int RECO_nMuHLTMatch,RECO_nEleHLTMatch;
      std::vector<float> RECOMU_PT_MuHLTMatch,RECOMU_ETA_MuHLTMatch,RECOMU_PHI_MuHLTMatch;
      std::vector<float> RECOELE_PT_EleHLTMatch,RECOELE_ETA_EleHLTMatch,RECOELE_PHI_EleHLTMatch;

      //=============================================================
      //                   
      //           Create Branchs for PF MET
      //
      //=============================================================
      float PFMet_pt;
      float PFMet_eta;
      float PFMet_phi;
      float PFMet_en;
      float PFMet_px;
      float PFMet_py;
      float PFMet_pz;
      float PFMet_sumEt; 
      double METSign;
      //==================================================
      //
      //     Create vectors for Electrons variables
      //
      //==================================================
      int RECO_NELE;
      std::vector<int> RECOELE_isEcalDriven;
      std::vector<int> RECOELE_isTrackerDriven;
      std::vector<int> RECOELE_CHARGE;
      std::vector<float> RECOELE_E;
      std::vector<float> RECOELE_PT;
      std::vector<float> RECOELE_P;
      std::vector<float> RECOELE_ETA;
      std::vector<float> RECOELE_THETA;
      std::vector<float> RECOELE_PHI;
      std::vector<float> RECOELE_MASS;
      std::vector<float> ele_sclRawE;
      std::vector<float> RECOELE_scl_E;
      std::vector<float> RECOELE_scl_Et;
      std::vector<float> RECOELE_scl_Eta;
      std::vector<float> RECOELE_scl_Phi;
      std::vector<float> ele_sclX;
      std::vector<float> ele_sclY;
      std::vector<float> ele_sclZ;
      std::vector<float> RECOELE_PTError;
      std::vector<float> RECOELE_COV;
      std::vector<float> RECOELE_EGMECALISO;
      std::vector<float> RECOELE_EGMHCALISO;
      std::vector<float> RECOELE_EGMX;
      // PF isolation
      std::vector<float> RECOELE_PFchAllPart;
      std::vector<float> RECOELE_PFchHad;
      std::vector<float> RECOELE_PFneuHad;
      std::vector<float> RECOELE_PFphoton;
      std::vector<float> RECOELE_PFPUchAllPart;
      std::vector<float> RECOELE_PFX_dB;
      std::vector<float> RECOELE_PFX_rho;
      // Vertexing DA
      std::vector<float> RECOELE_SIP;
      std::vector<float> RECOELE_IP;
      std::vector<float> RECOELE_IPERROR;
      // KF
      std::vector<float> RECOELE_SIP_KF;
      std::vector<float> RECOELE_IP_KF;
      std::vector<float> RECOELE_IPERROR_KF;
      std::vector<float> RECOELE_STIP;
      std::vector<float> RECOELE_SLIP;
      std::vector<float> RECOELE_TIP;
      std::vector<float> RECOELE_LIP;
      std::vector<float> RECOELE_TIPERROR;
      std::vector<float> RECOELE_LIPERROR;
      // GsfTrack
      std::vector<float> RECOELE_gsftrack_NPixHits;
      std::vector<float> RECOELE_gsftrack_NStripHits;
      std::vector<float> RECOELE_gsftrack_chi2;
      std::vector<float> RECOELE_gsftrack_dxyB;
      std::vector<float> RECOELE_gsftrack_dxy;
      std::vector<float> RECOELE_gsftrack_dxyError;
      std::vector<float> RECOELE_gsftrack_dzB;
      std::vector<float> RECOELE_gsftrack_dz;
      std::vector<float> RECOELE_gsftrack_dzError;
      //Conversion variables
      std::vector<int> RECOELE_gsftrack_losthits;
      std::vector<int> RECOELE_gsftrack_validhits;
      std::vector<int> RECOELE_gsftrack_expected_inner_hits;
      // Track-Cluster matching attributes
      std::vector<float> RECOELE_ep;
      std::vector<float> RECOELE_eSeedp;
      std::vector<float> RECOELE_eSeedpout;
      std::vector<float> RECOELE_eElepout;
      std::vector<float> RECOELE_deltaEtaIn;
      std::vector<float> RECOELE_deltaEtaSeed;
      std::vector<float> RECOELE_deltaEtaEle;
      std::vector<float> RECOELE_deltaPhiIn;
      std::vector<float> RECOELE_deltaPhiSeed;
      std::vector<float> RECOELE_deltaPhiEle;
      // Fiducial flags
      std::vector<int> RECOELE_isbarrel;
      std::vector<int> RECOELE_isendcap;
      std::vector<int> RECOELE_isGap;
      std::vector<int> RECOELE_isEBetaGap;
      std::vector<int> RECOELE_isEBphiGap;
      std::vector<int> RECOELE_isEEdeeGap;
      std::vector<int> RECOELE_isEEringGap;
      // Shower shape
      std::vector<float> RECOELE_sigmaIetaIeta;
      std::vector<float> RECOELE_sigmaEtaEta;
      std::vector<float> RECOELE_e15;
      std::vector<float> RECOELE_e25max;
      std::vector<float> RECOELE_e55;
      std::vector<float> RECOELE_he;
      //std::vector<float> RECOELE_r9;
      // Brem & Classifaction
      std::vector<float> RECOELE_fbrem;
      std::vector<float> RECOELE_nbrems;
      std::vector<float> RECOELE_Class;
      std::vector<float> RECOELE_fbrem_mean;
      std::vector<float> RECOELE_fbrem_mode;
      // Corrections
      std::vector<float> RECOELE_ecalEnergy;
      // Seed Collection
      std::vector<int> ele_seedSubdet2;
      std::vector<double> ele_seedDphi2;
      std::vector<double> ele_seedDrz2;
      std::vector<int> ele_seedSubdet1;
      std::vector<double> ele_seedDphi1;
      std::vector<double> ele_seedDrz1;
      std::vector<float> RECOELE_mvaTrigV0;
      std::vector<float> RECOELE_mvaNonTrigV0;
      // Conversion Finder
      std::vector<float> ConvMapDist;
      std::vector<float> ConvMapDcot;
      // Matching
      std::vector<int> RECOELE_MatchingMCTruth;
      std::vector<float> RECOELE_MatchingMCpT;
      std::vector<float> RECOELE_MatchingMCEta;
      std::vector<float> RECOELE_MatchingMCPhi;
      //=============================================================
      //
      //  Create Branches for reco leptons Tree
      //
      //=============================================================
      // RECO additional block for reconstructed higgs, Z and their daughters
      std::vector<float> RECO_ZMM_MASS;
      std::vector<float> RECO_ZEE_MASS;
      std::vector<float> RECO_DiLep_MASS;
      std::vector<float> RECO_ZMM_PT;
      std::vector<float> RECO_ZEE_PT;
      std::vector<float> RECO_DiLep_PT;
      std::vector<float> RECO_ZMM_ETA;
      std::vector<float> RECO_ZEE_ETA;
      std::vector<float> RECO_DiLep_ETA;
      std::vector<float> RECO_ZMM_PHI;
      std::vector<float> RECO_ZEE_PHI;
      std::vector<float> RECO_DiLep_PHI;
      std::vector<float> RECO_ZMMss_MASS;
      std::vector<float> RECO_ZEEss_MASS;
      std::vector<float> RECO_ZEM_MASS;
      std::vector<float> RECO_ZMMss_PT;
      std::vector<float> RECO_ZEEss_PT;
      std::vector<float> RECO_ZEM_PT;
      std::vector<float> RECO_ZMMss_ETA;
      std::vector<float> RECO_ZEEss_ETA;
      std::vector<float> RECO_ZEM_ETA;
      std::vector<float> RECO_ZMMss_PHI;
      std::vector<float> RECO_ZEEss_PHI;
      std::vector<float> RECO_ZEM_PHI;
      std::vector<float> RECO_MMMM_MASS;
      std::vector<float> RECO_MMMM_PT;
      std::vector<float> RECO_MMMM_ETA;
      std::vector<float> RECO_MMMM_PHI;
      std::vector<float> RECO_MMMM_MASS_REFIT;
      std::vector<float> RECO_EEEE_MASS;
      std::vector<float> RECO_EEEE_PT;
      std::vector<float> RECO_EEEE_ETA;
      std::vector<float> RECO_EEEE_PHI;
      std::vector<float> RECO_EEEE_MASS_REFIT;
      std::vector<float> RECO_EEMM_MASS;
      std::vector<float> RECO_EEMM_PT;
      std::vector<float> RECO_EEMM_ETA;
      std::vector<float> RECO_EEMM_PHI;
      std::vector<float> RECO_EEMM_MASS_REFIT;
      std::vector<float> RECO_LLL0_MASS;
      std::vector<float> RECO_LLL1_MASS;
      std::vector<float> RECO_LLL2_MASS;
      std::vector<float> RECO_LLL3_MASS;
      std::vector<float> RECO_LLL0_PT;
      std::vector<float> RECO_LLL1_PT;
      std::vector<float> RECO_LLL2_PT;
      std::vector<float> RECO_LLL3_PT;
      std::vector<float> RECO_LLLl0_MASS;
      std::vector<float> RECO_LLLl1_MASS;
      std::vector<float> RECO_LLLl0_PT;
      std::vector<float> RECO_LLLl1_PT;
      std::vector<float> RECO_LLLL0ss_MASS;
      std::vector<float> RECO_LLLL0ss_PT;
      std::vector<float> RECO_LLLL1ss_MASS;
      std::vector<float> RECO_LLLL1ss_PT;
      std::vector<float> RECO_LLLL2ss_MASS;
      std::vector<float> RECO_LLLL2ss_PT;
      //std::vector<float> RECOcollNameLLLLssos_MASS;
      //std::vector<float> RECOcollNameLLLLssos_PT;
      std::vector<float> RECO_LLLL_MASS;
      std::vector<float> RECO_LLLL_PT;
      std::vector<float> RECO_LLLL_ETA;
      std::vector<float> RECO_LLLL_PHI;
      // Matching ZtoMuMu
      std::vector<bool> RECOzMuMu_MatchingMCTruth;
      std::vector<float> RECOzMuMu_MatchingMCpT;
      std::vector<float> RECOzMuMu_MatchingMCmass;
      std::vector<float> RECOzMuMu_MatchingMCEta;
      std::vector<float> RECOzMuMu_MatchingMCPhi;
      // Matching ZtoEE
      std::vector<bool> RECOzEE_MatchingMCTruth;
      std::vector<float> RECOzEE_MatchingMCpT;
      std::vector<float> RECOzEE_MatchingMCmass;
      std::vector<float> RECOzEE_MatchingMCEta;
      std::vector<float> RECOzEE_MatchingMCPhi;
      // Matching HiggsToMMMM
      std::vector<bool> RECOHzzMMMM_MatchingMCTruth;
      std::vector<float> RECOHzzMMMM_MatchingMCpT;
      std::vector<float> RECOHzzMMMM_MatchingMCmass;
      std::vector<float> RECOHzzMMMM_MatchingMCEta;
      std::vector<float> RECOHzzMMMM_MatchingMCPhi;
      // Matching HiggsToEEEE
      std::vector<bool> RECOHzzEEEE_MatchingMCTruth;
      std::vector<float> RECOHzzEEEE_MatchingMCpT;
      std::vector<float> RECOHzzEEEE_MatchingMCmass;
      std::vector<float> RECOHzzEEEE_MatchingMCEta;
      std::vector<float> RECOHzzEEEE_MatchingMCPhi;
      // Matching HiggsToEEMM
      std::vector<bool> RECOHzzEEMM_MatchingMCTruth;
      std::vector<float> RECOHzzEEMM_MatchingMCpT;
      std::vector<float> RECOHzzEEMM_MatchingMCmass;
      std::vector<float> RECOHzzEEMM_MatchingMCEta;
      std::vector<float> RECOHzzEEMM_MatchingMCPhi;
      //=============================================================
      //                   
      //           Create Branchs for PileUp tree
      //
      //=============================================================
      //int num_PU_vertices;
      //int PU_BunchCrossing;
      // Beam Spot
      BeamSpot bs;
      double BeamSpot_X,BeamSpot_Y,BeamSpot_Z;
      //=============================================================
      //                   
      //           Create Branch for Rho
      //
      //=============================================================
      std::vector<float> Rho2;
      //=============================================================
      //                   
      //           Create Branch for events reweighting
      //
      //=============================================================
      std::vector<float> MC_weighting;
      //=============================================================
      //                   
      //           Create Branchs for ConstraintVtx2e2mu tree
      //
      //=============================================================
      std::vector<double> StdFitVertexX;
      std::vector<double> StdFitVertexY;
      std::vector<double> StdFitVertexZ;
      std::vector<double> StdFitVertexChi2r;
      std::vector<double> StdFitVertexProb;
      std::vector<double> StdFitVertexTrack_PT;
      std::vector<double> StdFitVertexTrack_ETA;
      std::vector<double> StdFitVertexTrack_PHI;
      std::vector<double> KinFitVertexX;
      std::vector<double> KinFitVertexY;
      std::vector<double> KinFitVertexZ;
      std::vector<double> KinFitVertexChi2r;
      std::vector<double> KinFitVertexProb;
      //=============================================================
      //                   
      //           Create Branchs for fillConstraintVtx4mu tree
      //
      //=============================================================
      std::vector<double> StdFitVertexXMMMM;
      std::vector<double> StdFitVertexYMMMM;
      std::vector<double> StdFitVertexZMMMM;
      std::vector<double> StdFitVertexChi2rMMMM;
      std::vector<double> StdFitVertexProbMMMM;
      std::vector<double> StdFitVertexTrackMMMM_PT;
      std::vector<double> StdFitVertexTrackMMMM_ETA;
      std::vector<double> StdFitVertexTrackMMMM_PHI;
      std::vector<double> KinFitVertexXMMMM;
      std::vector<double> KinFitVertexYMMMM;
      std::vector<double> KinFitVertexZMMMM;
      std::vector<double> KinFitVertexChi2rMMMM;
      std::vector<double> KinFitVertexProbMMMM;
      //=============================================================
      //                   
      //           Create Branchs for fillConstraintVtx4e tree
      //
      //=============================================================
      std::vector<double> StdFitVertexXEEEE;
      std::vector<double> StdFitVertexYEEEE;
      std::vector<double> StdFitVertexZEEEE;
      std::vector<double> StdFitVertexChi2rEEEE;
      std::vector<double> StdFitVertexProbEEEE;
      std::vector<double> StdFitVertexTrackEEEE_PT;
      std::vector<double> StdFitVertexTrackEEEE_ETA;
      std::vector<double> StdFitVertexTrackEEEE_PHI;
      std::vector<double> KinFitVertexXEEEE;
      std::vector<double> KinFitVertexYEEEE;
      std::vector<double> KinFitVertexZEEEE;
      std::vector<double> KinFitVertexChi2rEEEE;
      std::vector<double> KinFitVertexProbEEEE;
      //=============================================================
      //                   
      //    Create Branchs for fillConstraintVtxDiLeptons tree
      //
      //=============================================================
      std::vector<double> StdFitVertexChi2rDiLep;
      std::vector<double> StdFitVertexProbDiLep;
      //=============================================================
      //                   
      //    Create Branchs for fillConstraintVtxTripLeptons tree
      //
      //=============================================================
      std::vector<double> StdFitVertexChi2rMMM;
      std::vector<double> StdFitVertexProbMMM;
      std::vector<double> StdFitVertexChi2rMME;
      std::vector<double> StdFitVertexProbMME;
      std::vector<double> StdFitVertexChi2rEEE;
      std::vector<double> StdFitVertexProbEEE;
      std::vector<double> StdFitVertexChi2rMEE;
      std::vector<double> StdFitVertexProbMEE;
      //=============================================================
      //                   
      //           Create Branchs for Geom. Discri. tree
      //
      //=============================================================
      std::vector<double> ftsigma;
      std::vector<double> gdX;
      std::vector<double> gdY;
      std::vector<double> gdZ;
      std::vector<double> ftsigmalag;
      std::vector<double> gdlagX;
      std::vector<double> gdlagY;
      std::vector<double> gdlagZ;
      std::vector<double> gdlagProb;
      std::vector<double> gdlagNdof;
      std::vector<double> ftsigmaMMMM;
      std::vector<double> gdXMMMM;
      std::vector<double> gdYMMMM;
      std::vector<double> gdZMMMM;
      std::vector<double> ftsigmalagMMMM;
      std::vector<double> gdlagXMMMM;
      std::vector<double> gdlagYMMMM;
      std::vector<double> gdlagZMMMM;
      std::vector<double> gdlagProbMMMM;
      std::vector<double> gdlagNdofMMMM;
      std::vector<double> ftsigmaEEEE;
      std::vector<double> gdXEEEE;
      std::vector<double> gdYEEEE;
      std::vector<double> gdZEEEE;
      std::vector<double> ftsigmalagEEEE;
      std::vector<double> gdlagXEEEE;
      std::vector<double> gdlagYEEEE;
      std::vector<double> gdlagZEEEE;
      std::vector<double> gdlagProbEEEE;
      std::vector<double> gdlagNdofEEEE;
      //=============================================================
      //                   
      //           Create Branchs for  RECORF block 2e2mu tree
      //
      //=============================================================
      std::vector<double> RECORF_2e2mu_cosTheta1_spin;
      std::vector<double> RECORF_2e2mu_cosTheta2_spin;
      std::vector<double> RECORF_2e2mu_cosThetaStar_spin;
      std::vector<double> RECORF_2e2mu_Phi_spin;
      std::vector<double> RECORF_2e2mu_Phi1_spin;
      std::vector<double> RECORF_2e2mu_Phi2_spin;
      std::vector<double> RECORF_2e2mu_phi1RF_spin;
      std::vector<double> RECORF_2e2mu_phi2RF_spin;
      std::vector<double> RECORF_2e2mu_MELA;
      
      std::vector<double> RECORF_4e_cosTheta1_spin;
      std::vector<double> RECORF_4e_cosTheta2_spin;
      std::vector<double> RECORF_4e_cosThetaStar_spin;
      std::vector<double> RECORF_4e_Phi_spin;
      std::vector<double> RECORF_4e_Phi1_spin;
      std::vector<double> RECORF_4e_Phi2_spin;
      std::vector<double> RECORF_4e_phi1RF_spin;
      std::vector<double> RECORF_4e_phi2RF_spin;
      std::vector<double> RECORF_4e_MELA;
      
      std::vector<double> RECORF_4mu_cosTheta1_spin;
      std::vector<double> RECORF_4mu_cosTheta2_spin;
      std::vector<double> RECORF_4mu_cosThetaStar_spin;
      std::vector<double> RECORF_4mu_Phi_spin;
      std::vector<double> RECORF_4mu_Phi1_spin;
      std::vector<double> RECORF_4mu_Phi2_spin;
      std::vector<double> RECORF_4mu_phi1RF_spin;
      std::vector<double> RECORF_4mu_phi2RF_spin;
      std::vector<double> RECORF_4mu_MELA;

      std::vector<double> MCRF_cosTheta1_spin;
      std::vector<double> MCRF_cosTheta2_spin;
      std::vector<double> MCRF_cosThetaStar_spin;
      std::vector<double> MCRF_Phi_spin;
      std::vector<double> MCRF_Phi1_spin;
      std::vector<double> MCRF_Phi2_spin;
      std::vector<double> MCRF_phi1RF_spin;
      std::vector<double> MCRF_phi2RF_spin;
      std::vector<double> MCRF_MELA;
      //=============================================================
      // RECO MET
      float genmet;
      float calomet;
      float LxyPhoConv, FitProbPhoConv;
      unsigned int NHitsBeforeVtxPhoConv;
      double lxyMin_;
      double probMin_;
      int nHitsBeforeVtxMax_;
      bool allowCkfMatch_;
      double maxAbsZ_;
      double maxd0_;
      int minNdof_;
      int NbGoodPv_;
      // counters
      //int irun, ievt,ils,nevt;
      //float Avginstlumi;
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
      
      // Generator
      edm::EDGetTokenT<GenEventInfoProduct> generator_;

      // HTL
      bool fillHLTinfo;
      edm::InputTag HLTInfoFired;
      std::string HLTAnalysisinst;
      std::vector<edm::InputTag> flagHLTnames; 
      
      edm::EDGetTokenT<trigger::TriggerEvent> triggerEvent;
      edm::EDGetTokenT<edm::Association<std::vector<pat::TriggerObjectStandAlone> > > triggerMatchObject;
      edm::InputTag triggerMatchObjectEle;
      std::string triggerFilter,triggerEleFilter,triggerHLTcollection;
      
      // SkimEarlyData
      std::string SkimEarlyDataAnalysisinst;
      std::vector<edm::InputTag> flagSkimEarlyDatanames; 
      
      // MC truth
      bool fillMCTruth;
      edm::EDGetTokenT<edm::View<reco::Candidate> > MCcollName;
      edm::EDGetTokenT<vector<reco::GenParticle> > genParticles_;
      edm::EDGetTokenT<edm::View<reco::Candidate> > fourgenleptons_;
      edm::EDGetTokenT<edm::View<reco::Candidate> > digenZ_;
	  
      // RECO
      bool useAdditionalRECO;
      bool use2011EA;
      std::vector<edm::InputTag> RECOcollNameBest2e2mu,RECOcollNameBest4e,RECOcollNameBest4mu,
	RECOcollNameBestRestFrame2e2mu,RECOcollNameBestRestFrame4mu,RECOcollNameBestRestFrame4e;
      std::vector<edm::InputTag> RECOcollNameZ,RECOcollNameZss,
	RECOcollNameMMMM,RECOcollNameEEEE,RECOcollNameEEMM,
	RECOcollNameLLLLss,RECOcollNameLLL,RECOcollNameLLLl,RECOcollNameLLLLssos;
      edm::InputTag RECOcollNameLLLL,RECOcollNameDiLep;
      
      // electron and muon tags
      bool useBestCandidate;
      edm::InputTag BestCandidatesLeptonsTag_;
      edm::InputTag clusterCollectionTag_,gsftrackCollection_;
      edm::EDGetTokenT<edm::View<reco::Muon> > muonPFTag_;
      edm::EDGetTokenT<edm::View<reco::Muon> > muonTag_;
      edm::EDGetTokenT<edm::ValueMap<float> > muonCorrPtErrorMapTag_;
      edm::EDGetTokenT<edm::View<reco::GsfElectron> > electronEgmTag_;
      edm::InputTag electronMapTag_;
      edm::InputTag electronEgmTkMapTag_;
      edm::InputTag electronEgmEcalMapTag_;
      edm::InputTag electronEgmHcalMapTag_;
      edm::EDGetTokenT<edm::View<reco::GsfElectron> > mvaElectronTag_;
      edm::EDGetTokenT<edm::ValueMap<float> > mvaTrigV0MapTag_,mvaNonTrigV0MapTag_;
      
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
      edm::EDGetTokenT<vector<reco::Track> > tracksTag_;
      edm::EDGetTokenT<edm::View<reco::PFCandidate> > pfphotonsTag_;
      edm::EDGetTokenT<edm::View<reco::Photon> > photonsTag_;
      vector<reco::Vertex> PV;
      edm::EDGetTokenT<std::vector<reco::Vertex> >verticesTag_;
      edm::EDGetTokenT<double> rhojetsTag_;
      edm::EDGetTokenT<std::vector<reco::PFJet> > jetsTag_, jetsDataTag_, jetsMVATag_;
      edm::EDGetTokenT<edm::ValueMap<float> >  PuJetMvaMCfullDiscr_,PuJetMvaMCfullId_;
      edm::EDGetTokenT<edm::ValueMap<float> >  PuJetMvaDatafullDiscr_,PuJetMvaDatafullId_;
      // MET
      edm::InputTag trackermetTag_;
      edm::EDGetTokenT<vector<reco::PFMET> >  pfmetTag_;
      edm::EDGetTokenT<vector<reco::GenMET> > genmetTag_;
      edm::InputTag calometTag_,calometoptTag_,calometoptnohfTag_,calometoptnohfhoTag_;
      edm::InputTag calometopthoTag_,calometnohfTag_,calometnohfhoTag_,calomethoTag_;
      bool useAdditionalMET_;
      edm::InputTag htmetic5Tag_,htmetkt4Tag_,htmetkt6Tag_,htmetsc5Tag_,htmetsc7Tag_;
      edm::InputTag jescormetic5Tag_,jescormetkt4Tag_,jescormetkt6Tag_,jescormetsc5Tag_,jescormetsc7Tag_;
      edm::InputTag cormetMuTag_;
        
      // Conversion
      edm::EDGetTokenT<edm::ValueMap<float> > ConvMapDistTag_,ConvMapDcotTag_;

      // Matching
      edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > > goodElectronMCMatch_;
      edm::EDGetTokenT<reco::CandidateCollection> myElectrons_;
      edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > > goodMuonMCMatch_;
      edm::EDGetTokenT<reco::CandidateCollection> myMuons_;
      edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > > goodGammaMCMatch_;
      edm::EDGetTokenT<reco::CandidateCollection> myGammas_;

      // Beam Spot
      edm::EDGetTokenT<reco::BeamSpot> offlineBeamSpot_;
      
      // bTagging
      edm::EDGetTokenT<edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,std::vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference> >  
	tCHighEff_bTag_, tCHighPur_bTag_,cSV_bTag_;

      
      // counters
      int irun, ievt,ils,nevt;
      float Avginstlumi;
      
      
      // MC info
      edm::ESHandle<ParticleDataTable>  pdt_;

      std::string outputFile_; // output file
  
    
      
      /// number of HLT trigger paths requested in configuration
      unsigned int n_;
      bool firstevent_;

      // PG and FRC 06-07-11 try to reduce printout!
      bool debug;
      
      /// list of required HLT triggers by HLT name
      //std::vector<std::string > HLTPathsByName_;
      /// list of required HLT triggers by HLT index
      //std::vector<unsigned int> HLTPathsByIndex_;

      //std::vector<std::string> triggerFilterMuon;
      //std::vector<std::string> triggerFilterEle;

      
      //===============================
      // root file to store histograms
      TFile*  rootFile_;
      // min and max of energy histograms
      
      
      

};


#endif



