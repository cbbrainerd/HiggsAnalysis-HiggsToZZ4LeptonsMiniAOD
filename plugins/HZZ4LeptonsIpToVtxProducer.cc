
/**\class HZZ4LeptonsIpToVtxProducer (from HZZ2e2muIpToVtxProducer)
 *
 *
 * Original Author:  Alexis Pompili - Univ. of Bari & INFN
 *    modifed by     Nicola De Filippis - Polytechnic. of Bari & INFN               
 *
 * Computes for each lepton the IP significance w.r.t to PrimaryVertex
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsIpToVtxProducer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

// Tracker tracks 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GeometryVector/interface/Basic3DVector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"


// Transient tracks
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

// Vertex utilities
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// Beam Spot
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// Geometry
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Common/interface/AssociationVector.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/RefToBase.h"

// Other Tracking Tools
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/GsfTools/interface/GSUtilities.h"
#include "TrackingTools/GsfTools/interface/GsfPropagatorAdapter.h"
#include "TrackingTools/GsfTools/interface/GaussianSumUtilities1D.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianState1D.h"
//
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

// Other specific
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "TrackingTools/IPTools/interface/IPTools.h"


// Maps
#include "DataFormats/Common/interface/ValueMap.h"

//-- system include files
#include <memory>
#include <vector>
#include <string>

//-- namespaces
using namespace edm;
using namespace std;
using namespace reco;
using namespace IPTools;

bool cmpf( float a, float b ) {
  return fabs(a) > fabs(b);
}

//
//-- constructor
//
HZZ4LeptonsIpToVtxProducer::HZZ4LeptonsIpToVtxProducer(const edm::ParameterSet& pset) {

  // Decay Channel
  decaychannel                                                    = pset.getParameter<std::string>("decaychannel");
  if (decaychannel=="2e2mu" || decaychannel=="4mu")  muonTag_     = consumes<edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("MuonsLabel"));
  if (decaychannel=="2e2mu" || decaychannel=="4e")   electronTag_ = consumes<edm::View<pat::Electron> >(pset.getParameter<edm::InputTag>("ElectronsLabel"));
  vertexTag_                                                      = consumes<std::vector<reco::Vertex> >(pset.getParameter<edm::InputTag>("VertexLabel"));
  offlineBeamSpot_                                                = consumes<reco::BeamSpot>(pset.getParameter<edm::InputTag>("BeamSpotLabel"));  
  useBeamSpot_                                                    = pset.getParameter<bool>("useBeamSpot");

  // PG and FRC 06-07-11
  debug	=	pset.getUntrackedParameter<bool> ("debug");
  
  //
  string alias;
  string iName = "IpToVtx";
  //
  if (decaychannel=="2e2mu" || decaychannel=="4e")  {
    produces<edm::ValueMap<float> >( alias = "VertexEleMap");  
    produces<edm::ValueMap<float> >( alias = "VertexValueEleMap");  
    produces<edm::ValueMap<float> >( alias = "VertexErrorEleMap"); 
  }
  if (decaychannel=="2e2mu" || decaychannel=="4mu"){ 
    produces<edm::ValueMap<float> >( alias = "VertexMuMap");
    produces<edm::ValueMap<float> >( alias = "VertexValueMuMap");
    produces<edm::ValueMap<float> >( alias = "VertexErrorMuMap");
  }
}

//
//-- destructor
//
HZZ4LeptonsIpToVtxProducer::~HZZ4LeptonsIpToVtxProducer() {
}

void HZZ4LeptonsIpToVtxProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<edm::ValueMap<float> >     VertexMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler     fillerMu(*VertexMuMap);
  auto_ptr<edm::ValueMap<float> > VertexEleMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerEle(*VertexEleMap);
  
  auto_ptr<edm::ValueMap<float> >     VertexValueMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler     fillerValueMu(*VertexValueMuMap);
  auto_ptr<edm::ValueMap<float> > VertexValueEleMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerValueEle(*VertexValueEleMap);
  
  auto_ptr<edm::ValueMap<float> >     VertexErrorMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler     fillerErrorMu(*VertexErrorMuMap);
  auto_ptr<edm::ValueMap<float> > VertexErrorEleMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerErrorEle(*VertexErrorEleMap);


  // Get the B-field
  ESHandle<MagneticField> B;
  iSetup.get<IdealMagneticFieldRecord>().get( B );

  // Get the geometry
  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  // get the track builder
  ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);

  // Beamspot 
  Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(offlineBeamSpot_,recoBeamSpotHandle);
  const BeamSpot bs = *recoBeamSpotHandle;

  if(debug) cout << "BeamSpot position: X,Y,Z=" << bs.position().x() << " " << bs.position().y() << " " << bs.position().z() << " " << bs.x0() << " " << bs.y0() << " " << bs.z0() << endl;
  GlobalPoint BSVertex(bs.position().x(),bs.position().y(),bs.position().z());
  Basic3DVector<double> BSVertexErr(bs.x0Error(),bs.y0Error(),bs.z0Error());

  reco::Vertex::Point BSp(bs.position().x(),bs.position().y(),bs.position().z());
  reco::Vertex::Error BSe;
  
  BSe(0,0) = bs.x0Error()*bs.x0Error();
  BSe(1,1) = bs.y0Error()*bs.y0Error();
  BSe(2,2) = bs.z0Error()*bs.z0Error();
  reco::Vertex BSprimVertex = reco::Vertex(BSp,BSe,1,1,1);


  // Get primary vertex collection
  Handle<reco::VertexCollection> recoPVCollection;
  iEvent.getByToken(vertexTag_,recoPVCollection);
  //
  reco::Vertex primVertex;
  GlobalPoint pVertex(0., 0., 0.);
  Basic3DVector<double> pVertexErr(0.,0.,0.);
  
  bool pvfound = (recoPVCollection->size() != 0);

  if (pvfound){
    for(reco::VertexCollection::const_iterator it=recoPVCollection->begin() ; it!=recoPVCollection->end() ; ++it){
      //if (it==recoPVCollection->begin()){
      if(!it->isFake() && it->ndof() > 4 && fabs(it->position().z()) <= 24 && fabs(it->position().rho()) <= 2){
	primVertex=*it;
	pVertex = GlobalPoint(it->position().x(), it->position().y(), it->position().z());
	pVertexErr = Basic3DVector<double>(it->xError(),it->yError(),it->zError());
	cout << "P vertex position used for computing SIP for electron and muons is (x,y,z)= " << 
	  pVertex.x() << " " <<
	  pVertex.y() << " " <<
	  pVertex.z() << endl;
	break;
      }
    }     
  }

  //
  /////////////////////////////////////////////
  //
  //--track refs
  TrackRef mutrack;
  GsfTrackRef eletrack;
  //
  //--transient tracks
  reco::TransientTrack mutranstrack;
  reco::TransientTrack eletranstrack;
  //--signif. vars
  float muSignificance3D;
  float eleSignificance3D;
  float  muValue3D,muError3D;
  float eleValue3D,eleError3D;

  //
  // Muons first
  //============ 
  //
 
  if (decaychannel=="2e2mu" || decaychannel=="4mu"){
    edm::Handle<edm::View<pat::Muon> > muCandidates;
    iEvent.getByToken(muonTag_,muCandidates);

    size_t nmu = muCandidates->size();
    std::vector<float> sigmu(nmu),valuemu(nmu),errormu(nmu);
    unsigned int indexmu=0;	
    for (edm::View<pat::Muon>::const_iterator muCand = muCandidates->begin(); muCand != muCandidates->end(); ++muCand){
      edm::Ref<edm::View<pat::Muon> > muonBaseRef(muCandidates,indexmu); 
      reco::TrackRef innerTrack = muonBaseRef->innerTrack();
      reco::TrackRef globalTrack= muonBaseRef->globalTrack();
      reco::TrackRef bestTrack  = muonBaseRef->muonBestTrack();
      reco::TrackRef mutrack = innerTrack;
      // Make sure the collection it points to is there
      if ( bestTrack.isNonnull() && bestTrack.isAvailable() ) mutrack = bestTrack;
      
      if (! (mutrack.isNonnull() && mutrack.isAvailable() )) {	
	if(debug) cout <<"tracker trackref is null since muon is STA" << endl;
	mutrack=muCand->get<TrackRef,reco::StandAloneMuonTag>();
      }  
      //      
      mutranstrack = trackBuilder->build( mutrack ) ;
      //  

      // TrajectoryStateOnSurface muTSOS;
      // if (useBeamSpot_==true){ 
      // 	muTSOS = IPTools::transverseExtrapolate(mutranstrack.impactPointState(), BSVertex, mutranstrack.field());	
      // } 
      // else {
      // 	muTSOS = IPTools::transverseExtrapolate(mutranstrack.impactPointState(), pVertex, mutranstrack.field());
      // }

      //
      //if (muTSOS.isValid()){
	std::pair<bool,Measurement1D> muIPpair;
	
	if (useBeamSpot_==true){ 		  
	  //muIPpair = IPTools::absoluteImpactParameter3D(mutranstrack, BSprimVertex);
	  muIPpair = IPTools::signedImpactParameter3D(mutranstrack,GlobalVector(mutrack->px(),mutrack->py(),mutrack->pz()),BSprimVertex);
	} 
	else {	 
	  //muIPpair = IPTools::absoluteImpactParameter3D(mutranstrack, primVertex);
	  muIPpair = IPTools::signedImpactParameter3D(mutranstrack,GlobalVector(mutrack->px(),mutrack->py(),mutrack->pz()),primVertex);	  
	}
	//
	if (muIPpair.first){
	  muSignificance3D = muIPpair.second.significance();
	  muValue3D = muIPpair.second.value();
	  muError3D = muIPpair.second.error();
	  if(debug) cout << "3DIP Significance for muons= " << muSignificance3D << endl;	
	  sigmu[indexmu]=float(muSignificance3D);
	  valuemu[indexmu]=float(muValue3D);
	  errormu[indexmu]=float(muError3D);
	} 	       
	//}
      //
      indexmu++;
    } //-- muon loop closed
    //

    // filling map
    fillerMu.insert(muCandidates, sigmu.begin(), sigmu.end());
    fillerMu.fill();

    fillerValueMu.insert(muCandidates, valuemu.begin(), valuemu.end());
    fillerValueMu.fill();
    fillerErrorMu.insert(muCandidates, errormu.begin(), errormu.end());
    fillerErrorMu.fill();

  }
  
  //
  // Electrons now
  //  //
  if (decaychannel=="2e2mu" || decaychannel=="4e"){    
    
    Handle<edm::View<pat::Electron> > eleCandidates; 
    iEvent.getByToken(electronTag_,eleCandidates);

    size_t nele = eleCandidates->size();
    std::vector<float> sigele(nele),valueele(nele),errorele(nele);
    unsigned int indexele=0;
    for (edm::View<pat::Electron>::const_iterator eleCand = eleCandidates->begin(); eleCand != eleCandidates->end(); ++eleCand){

      eletrack = eleCand->get<GsfTrackRef>();

      if ( eletrack.isNonnull() && eletrack.isAvailable() ) 
	eletranstrack = trackBuilder->build( eletrack ) ;

      // TrajectoryStateOnSurface eleTSOS;

      // if (useBeamSpot_==true){ 
      // 	eleTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), BSVertex, eletranstrack.field());
      // } 
      // else {
      // 	eleTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), pVertex, eletranstrack.field());
      // }


      //if (eleTSOS.isValid()){
	std::pair<bool,Measurement1D> eleIPpair;
	if (useBeamSpot_==true){ 
	  //eleIPpair = IPTools::signedImpactParameter3D(eletranstrack, eleTSOS.globalDirection(), BSprimVertex);
	  eleIPpair = IPTools::signedImpactParameter3D(eletranstrack,GlobalVector(eletrack->px(),eletrack->py(),eletrack->pz()),BSprimVertex);
	}
	else {
	  //eleIPpair = IPTools::signedImpactParameter3D(eletranstrack, eleTSOS.globalDirection(), primVertex);
	  eleIPpair = IPTools::signedImpactParameter3D(eletranstrack,GlobalVector(eletrack->px(),eletrack->py(),eletrack->pz()),primVertex);
	}
	//
	if (eleIPpair.first){
	  eleSignificance3D = eleIPpair.second.significance();
	  eleValue3D = eleIPpair.second.value();
	  eleError3D = eleIPpair.second.error();
	  if(debug) cout << "3DIP Significance for electrons= " <<  eleSignificance3D << endl;
	  sigele[indexele]=float(eleSignificance3D);
	  valueele[indexele]=float(eleValue3D);
	  errorele[indexele]=float(eleError3D);
	} 	
	//}            
      indexele++;
    } //-- ele loop closed
      
    fillerEle.insert(eleCandidates, sigele.begin(), sigele.end());
    fillerEle.fill();
    fillerValueEle.insert(eleCandidates, valueele.begin(), valueele.end());
    fillerValueEle.fill();
    fillerErrorEle.insert(eleCandidates, errorele.begin(), errorele.end());
    fillerErrorEle.fill();
  }

  if (decaychannel=="2e2mu" || decaychannel=="4mu"){
    // iEvent.put( VertexMuMap, "VertexMuMap" );
    // iEvent.put( VertexValueMuMap, "VertexValueMuMap" );
    // iEvent.put( VertexErrorMuMap, "VertexErrorMuMap" );

    iEvent.put(std::make_unique<edm::ValueMap<float>>(*VertexMuMap), "VertexMuMap" );
    iEvent.put(std::make_unique<edm::ValueMap<float>>(*VertexValueMuMap), "VertexValueMuMap" );
    iEvent.put(std::make_unique<edm::ValueMap<float>>(*VertexErrorMuMap), "VertexErrorMuMap" );

  }
  if (decaychannel=="2e2mu" || decaychannel=="4e"){
    // iEvent.put( VertexEleMap, "VertexEleMap" );
    // iEvent.put( VertexValueEleMap, "VertexValueEleMap" );
    // iEvent.put( VertexErrorEleMap, "VertexErrorEleMap" );

    iEvent.put(std::make_unique<edm::ValueMap<float>>(*VertexEleMap), "VertexEleMap" );
    iEvent.put(std::make_unique<edm::ValueMap<float>>(*VertexValueEleMap), "VertexValueEleMap" );
    iEvent.put(std::make_unique<edm::ValueMap<float>>(*VertexErrorEleMap), "VertexErrorEleMap" );
  }
  
  //
  
}
   
