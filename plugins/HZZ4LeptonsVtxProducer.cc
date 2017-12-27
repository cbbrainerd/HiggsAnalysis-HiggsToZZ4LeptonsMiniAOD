/**\class HZZ4LeptonsVtxProducer
 *
 *
 * Original Author:  Alexis Pompili - Univ. of Bari & INFN           
 *
 * Computes for each lepton the IP significance w.r.t to PrimaryVertex and other VTX stuff
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsVtxProducer.h"

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

// JET STUFF
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SoftLeptonTagInfo.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
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

bool cmpf_new( float a, float b ) {
  return fabs(a) > fabs(b);
}

//
//-- constructor
//
HZZ4LeptonsVtxProducer::HZZ4LeptonsVtxProducer(const edm::ParameterSet& pset) {

  // Decay Channel
  decaychannel                                                    = pset.getParameter<std::string>("decaychannel");
  if (decaychannel=="2e2mu" || decaychannel=="4mu")  muonTag_     = pset.getParameter<edm::InputTag>("MuonsLabel");
  if (decaychannel=="2e2mu" || decaychannel=="4e")   electronTag_ = pset.getParameter<edm::InputTag>("ElectronsLabel");
  //  
  //bTagAlgo_                                                       = pset.getParameter<>;
  //
  string alias;
  string iName = "Vtx";
  //
  if (decaychannel=="2e2mu" || decaychannel=="4e")  {
    //
    produces<edm::ValueMap<float> >( alias = "Sip3DPVEleMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip3DPVCheckEleMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip3DPVwbsEleMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip3DBSEleMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip3DcloPVEleMap");
    produces<edm::ValueMap<float> >( alias = "Sip3DcloPVwbsEleMap");
    //
    produces<edm::ValueMap<float> >( alias = "Sip2DPVEleMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip2DPVwbsEleMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip2DBSEleMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip2DcloPVEleMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip2DcloPVwbsEleMap"); 
    //
    produces<edm::ValueMap<float> >( alias = "SipXYPVCheckNoIPToolsEleMap");
    produces<edm::ValueMap<float> >( alias = "SipXYBSCheckNoIPToolsEleMap");
    produces<edm::ValueMap<float> >( alias = "SipZPVCheckNoIPToolsEleMap");
    produces<edm::ValueMap<float> >( alias = "SipZBSCheckNoIPToolsEleMap");
  }
  //
  if (decaychannel=="2e2mu" || decaychannel=="4mu"){ 
    //
    produces<edm::ValueMap<float> >( alias = "Sip3DPVMuMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip3DPVCheckMuMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip3DPVwbsMuMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip3DBSMuMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip3DcloPVMuMap");
    produces<edm::ValueMap<float> >( alias = "Sip3DcloPVwbsMuMap");
    //
    produces<edm::ValueMap<float> >( alias = "Sip2DPVMuMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip2DPVwbsMuMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip2DBSMuMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip2DcloPVMuMap"); 
    produces<edm::ValueMap<float> >( alias = "Sip2DcloPVwbsMuMap"); 
    //
    produces<edm::ValueMap<float> >( alias = "SipXYPVCheckNoIPToolsMuMap");
    produces<edm::ValueMap<float> >( alias = "SipXYBSCheckNoIPToolsMuMap");
    produces<edm::ValueMap<float> >( alias = "SipZPVCheckNoIPToolsMuMap");
    produces<edm::ValueMap<float> >( alias = "SipZBSCheckNoIPToolsMuMap");
  }
}

//
//-- destructor
//
HZZ4LeptonsVtxProducer::~HZZ4LeptonsVtxProducer() {
}

////////////////////////////////////////////////////////////////////////////////////

void HZZ4LeptonsVtxProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //
  // muons
  //
  auto_ptr<edm::ValueMap<float> > sip3DPVMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSip3DPV(*sip3DPVMuMap);
  auto_ptr<edm::ValueMap<float> > sip3DPVCheckMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSip3DPVCheck(*sip3DPVCheckMuMap);
  auto_ptr<edm::ValueMap<float> > sip3DPVwbsMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSip3DPVwbs(*sip3DPVwbsMuMap);
  auto_ptr<edm::ValueMap<float> > sip3DBSMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSip3DBS(*sip3DBSMuMap);
  auto_ptr<edm::ValueMap<float> > sip3DcloPVMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSip3DcloPV(*sip3DcloPVMuMap);
  auto_ptr<edm::ValueMap<float> > sip3DcloPVwbsMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSip3DcloPVwbs(*sip3DcloPVwbsMuMap);
  //
  auto_ptr<edm::ValueMap<float> > sip2DPVMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSip2DPV(*sip2DPVMuMap);
  auto_ptr<edm::ValueMap<float> > sip2DPVwbsMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSip2DPVwbs(*sip2DPVwbsMuMap);
  auto_ptr<edm::ValueMap<float> > sip2DBSMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSip2DBS(*sip2DBSMuMap);
  auto_ptr<edm::ValueMap<float> > sip2DcloPVMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSip2DcloPV(*sip2DcloPVMuMap);
  auto_ptr<edm::ValueMap<float> > sip2DcloPVwbsMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSip2DcloPVwbs(*sip2DcloPVwbsMuMap);
  //
  auto_ptr<edm::ValueMap<float> > sipXYPVCheckNoIPToolsMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSipXYPVCheckNoIPTools(*sipXYPVCheckNoIPToolsMuMap);
  auto_ptr<edm::ValueMap<float> > sipXYBSCheckNoIPToolsMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSipXYBSCheckNoIPTools(*sipXYBSCheckNoIPToolsMuMap);
  auto_ptr<edm::ValueMap<float> > sipZPVCheckNoIPToolsMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSipZPVCheckNoIPTools(*sipZPVCheckNoIPToolsMuMap);
  auto_ptr<edm::ValueMap<float> > sipZBSCheckNoIPToolsMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuSipZBSCheckNoIPTools(*sipZBSCheckNoIPToolsMuMap);
  //
  ///// electrons
  //
  auto_ptr<edm::ValueMap<float> > sip3DPVEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSip3DPV(*sip3DPVEleMap);
  auto_ptr<edm::ValueMap<float> > sip3DPVCheckEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSip3DPVCheck(*sip3DPVCheckEleMap);
  auto_ptr<edm::ValueMap<float> > sip3DPVwbsEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSip3DPVwbs(*sip3DPVwbsEleMap);
  auto_ptr<edm::ValueMap<float> > sip3DBSEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSip3DBS(*sip3DBSEleMap);
  auto_ptr<edm::ValueMap<float> > sip3DcloPVEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSip3DcloPV(*sip3DcloPVEleMap);
  auto_ptr<edm::ValueMap<float> > sip3DcloPVwbsEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSip3DcloPVwbs(*sip3DcloPVwbsEleMap);
  //
  auto_ptr<edm::ValueMap<float> > sip2DPVEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSip2DPV(*sip2DPVEleMap);
  auto_ptr<edm::ValueMap<float> > sip2DPVwbsEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSip2DPVwbs(*sip2DPVwbsEleMap);
  auto_ptr<edm::ValueMap<float> > sip2DBSEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSip2DBS(*sip2DBSEleMap);
  auto_ptr<edm::ValueMap<float> > sip2DcloPVEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSip2DcloPV(*sip2DcloPVEleMap);
  auto_ptr<edm::ValueMap<float> > sip2DcloPVwbsEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSip2DcloPVwbs(*sip2DcloPVwbsEleMap);
  //
  auto_ptr<edm::ValueMap<float> > sipXYPVCheckNoIPToolsEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSipXYPVCheckNoIPTools(*sipXYPVCheckNoIPToolsEleMap);
  auto_ptr<edm::ValueMap<float> > sipXYBSCheckNoIPToolsEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSipXYBSCheckNoIPTools(*sipXYBSCheckNoIPToolsEleMap);
  auto_ptr<edm::ValueMap<float> > sipZPVCheckNoIPToolsEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSipZPVCheckNoIPTools(*sipZPVCheckNoIPToolsEleMap);
  auto_ptr<edm::ValueMap<float> > sipZBSCheckNoIPToolsEleMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler fillerEleSipZBSCheckNoIPTools(*sipZBSCheckNoIPToolsEleMap);
  //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Get the B-field
  ESHandle<MagneticField> B;
  iSetup.get<IdealMagneticFieldRecord>().get( B );

  // Get the geometry
  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  // get the track builder
  ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);
  // this is needed by the IPTools methods 
  //
  ////////////////////////////////////////////////////////////////
  //
  // OFFLINE BEST PRIMARY VERTEX ("offlinePrimaryVertices")
  //
  Handle<reco::VertexCollection> recoPVCollection;
  //iEvent.getByLabel(vertexTag_.label(),recoPVCollection);
  iEvent.getByLabel("offlineSlimmedPrimaryVertices",recoPVCollection);
  //
  reco::Vertex primVertex;
  bool isPVfound = (recoPVCollection->size() != 0);
  bool isPVgood = false;
  //
  if ( isPVfound )
    {
      for (size_t idx = 0; idx != recoPVCollection->size(); idx++)
	 {
	   const reco::Vertex &vtx = (*recoPVCollection)[idx];
	   if( vtx.ndof()>4 && fabs(vtx.z())<=25 && vtx.position().Rho()<=2 ) isPVgood = true;
	 }
    }
  //
  //if (isPVfound && isPVgood)
  if (isPVfound)
    {
     PrimaryVertexSorter pvs;
     vector<reco::Vertex> sortedList = pvs.sortedList(*(recoPVCollection.product()) );
     primVertex = (sortedList.front());
    } else {
            //creating a dummy PV
            reco::Vertex::Point p(0,0,0);  
     	    reco::Vertex::Error e;
     	    e(0,0) = 0.0015*0.0015;
     	    e(1,1) = 0.0015*0.0015;
     	    e(2,2) = 15.*15.;
     	    primVertex = reco::Vertex(p,e,1,1,1);
    }
  //
  GlobalPoint pVertex(primVertex.position().x(),primVertex.position().y(),primVertex.position().z());
  //Basic3DVector<double> pVertexErr(primVertex.xError(),primVertex.yError(),primVertex.zError());
  //
  cout << "PV FOUND = " << isPVfound << " and PV GOOD = " << isPVgood << endl;
  cout << "Best/front PV has Chi2 = " << primVertex.chi2() << " and NDof = " << primVertex.ndof() << endl;
  cout << "Best/front PV has Xpos = " << primVertex.position().x() << " with Xerror = " << primVertex.xError() << endl;
  cout << "Best/front PV has Ypos = " << primVertex.position().y() << " with Yerror = " << primVertex.yError() << endl;
  cout << "Best/front PV has Zpos = " << primVertex.position().z() << " with Zerror = " << primVertex.zError() << endl;
  //
  ///////////////////////////////////////////////////////////////////////////////////
  //
  // OFFLINE BEST PRIMARY VERTEX WITH BS CONSTRAINT ("offlinePrimaryVerticesWithBS")
  //
  Handle<reco::VertexCollection> recoPVCollectionWithBS;
  //iEvent.getByLabel(vertexTag_.label(),recoPVCollection);
  iEvent.getByLabel("offlinePrimaryVerticesWithBS",recoPVCollectionWithBS);
  //
  reco::Vertex primVertexWithBS;
  bool isPVwBSfound = (recoPVCollectionWithBS->size() != 0);
  bool isPVwBSgood = false;
  //
  if ( isPVwBSfound )
    {
      for (size_t jdx = 0; jdx != recoPVCollectionWithBS->size(); jdx++)
	 {
	   const reco::Vertex &vtxbs = (*recoPVCollectionWithBS)[jdx];
	   if( vtxbs.ndof()>4 && fabs(vtxbs.z())<=25 && vtxbs.position().Rho()<=2 ) isPVwBSgood = true;
	 }
    }
  //
  //if (isPVwBSfound && isPVwBSgood)
  if (isPVwBSfound)
    {
     PrimaryVertexSorter pvsbs;
     vector<reco::Vertex> sortedListWithBS = pvsbs.sortedList(*(recoPVCollectionWithBS.product()) );
     primVertexWithBS = (sortedListWithBS.front());
    } else {
            //creating a dummy PV
            reco::Vertex::Point p(0,0,0);  
     	    reco::Vertex::Error e;
     	    e(0,0) = 0.0015*0.0015;
     	    e(1,1) = 0.0015*0.0015;
     	    e(2,2) = 15.*15.;
     	    primVertexWithBS = reco::Vertex(p,e,1,1,1);
    }
  //
  //GlobalPoint pVertexWithBS(primVertexWithBS.position().x(),primVertexWithBS.position().y(),primVertexWithBS.position().z());
  //Basic3DVector<double> pVertexWithBSErr(primVertexWithBS.xError(),primVertexWithBS.yError(),primVertexWithBS.zError());  
  //
  cout << "PVBS FOUND = " << isPVwBSfound << " and PV GOOD = " << isPVwBSgood << endl;
  cout << "Best/front PVBS has Chi2 = " << primVertexWithBS.chi2() << " and NDof = " << primVertexWithBS.ndof() << endl;
  cout << "Best/front PVBS has Xpos = " << primVertexWithBS.position().x() << " with Xerror = " << primVertexWithBS.xError() << endl;
  cout << "Best/front PVBS has Ypos = " << primVertexWithBS.position().y() << " with Yerror = " << primVertexWithBS.yError() << endl;
  cout << "Best/front PVBS has Zpos = " << primVertexWithBS.position().z() << " with Zerror = " << primVertexWithBS.zError() << endl;
  //
  //////////////////////////////////////////////////////////////////////////////////////////
  //
  // BEAMSPOT
  // 
  Handle<reco::BeamSpot> recoBeamSpotHandle;
  //iEvent.getByType(recoBeamSpotHandle);
  iEvent.getByLabel("offlineBeamSpot",recoBeamSpotHandle);
  //
  bool isBSvalid = false;
  //
  reco::Vertex BSribbon;
  //
  if (recoBeamSpotHandle.isValid() )
    {
     isBSvalid = true;
     const reco::BeamSpot bs = *recoBeamSpotHandle;
     BSribbon = reco::Vertex(bs.position(),bs.covariance3D());
    } else {
            cout << "No beam spot available from EventSetup !" << endl;
            // creating a dummy beamspot 
            reco::Vertex::Point p(0,0,0);  
     	    reco::Vertex::Error e;
     	    e(0,0) = 0.005*0.005;
     	    e(1,1) = 0.005*0.005;
     	    e(2,2) = 15.*15.;
            BSribbon = reco::Vertex(p,e,1,1,1);
  }
  //
  GlobalPoint BSVertex(BSribbon.position().x(),BSribbon.position().y(),BSribbon.position().z());
  //Basic3DVector<double> BSVertexErr(bs.x0Error(),bs.y0Error(),bs.z0Error());
  //
  cout << "Valid BS FOUND = " << isBSvalid << endl;
  cout << "BS has Xpos = " << BSribbon.position().x() << " with Xerror = " << BSribbon.xError() << endl;
  cout << "BS has Ypos = " << BSribbon.position().y() << " with Yerror = " << BSribbon.yError() << endl;
  cout << "BS has Zpos = " << BSribbon.position().z() << " with Zerror = " << BSribbon.zError() << endl;
  //
  //////////////////////////////////////////////////////////////////////////////////////////
  //
  ///////////////////////////////////////////// JETS /////////////////////////////////
  //
  // get btag info
  //
  //  edm::Handle<reco::JetTagCollection> bTagHandle;
  //event.getByLabel(bTagAlgo_,bTagHandle);
  //if(not iEvent.getByLabel(bTagLabel_,bTagHandle)){
  //                                                  std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
  //                                                  <<bTagLabel_<<std::endl; exit(0);
  //                                                }
  //const reco::JetTagCollection & bTagColl= *(bTagHandle.product());


  //
  ///////////////////////////////////////////////////////////////////////////////////////
  //
  //--track refs
  //TrackRef mutrack;
  //GsfTrackRef eletrack;
  //
  //--transient tracks
  //reco::TransientTrack mutranstrack;
  //reco::TransientTrack eletranstrack;
  //
  // Muons
  //======= 
  //
  if (decaychannel=="2e2mu" || decaychannel=="4mu")
    {
     edm::Handle<edm::View<Muon> > muCandidates;
     iEvent.getByLabel(muonTag_.label(),muCandidates);
     //
     size_t nmu = muCandidates->size();
     //
     std::vector<float> sip3Dmu(nmu);
     std::vector<float> sip3DmuCheck(nmu);
     std::vector<float> sip3Dmu_wbs(nmu);
     std::vector<float> sip3Dmu_bs(nmu);
     std::vector<float> sip3Dmu_mindz(nmu);
     std::vector<float> sip3Dmu_wbs_mindz(nmu);
     //
     std::vector<float> sip2Dmu(nmu);
     std::vector<float> sip2Dmu_wbs(nmu);
     std::vector<float> sip2Dmu_bs(nmu);
     std::vector<float> sip2Dmu_mindz(nmu);
     std::vector<float> sip2Dmu_wbs_mindz(nmu);
     //
     std::vector<float> sipXYmuCheckNoIPTools(nmu);
     std::vector<float> sipXYBSmuCheckNoIPTools(nmu);
     std::vector<float> sipZmuCheckNoIPTools(nmu);
     std::vector<float> sipZBSmuCheckNoIPTools(nmu);
     //
     unsigned int indexmu=0;
     //	
     for (edm::View<reco::Muon>::const_iterator muCand = muCandidates->begin(); muCand != muCandidates->end(); ++muCand)
        {
	 TrackRef mutrack = muCand->get<TrackRef>(); 
         reco::TransientTrack muTransTrack = trackBuilder->build( mutrack ) ;
	 // protection for STA muons follows:
	 if (mutrack.isNull()){
	                       cout <<"tracker trackref is null since muon is STA" << endl;
	                       mutrack=muCand->get<TrackRef,reco::StandAloneMuonTag>();
	                      }  
	 //
	 const GlobalVector pmuvec = GlobalVector(mutrack->px(),mutrack->py(),mutrack->pz());
	 //
	 //-----------------checks--------------------------------
	 std::pair<bool,Measurement1D> muIP3DmindzCheck;
	 float muSignif3DmindzCheck = -9999;
	 //if (isPVfound && isPVgood)
	 if (isPVfound)
	   {
	     std::pair<bool,Measurement1D> muIP3DmindzCheck = IPTools::signedImpactParameter3D(muTransTrack,pmuvec,primVertex);
	     muSignif3DmindzCheck = muIP3DmindzCheck.second.significance();
	   }
	 sip3DmuCheck[indexmu] = muSignif3DmindzCheck; 
	 //
	 //---------------other check / PV ------------------------------
	 float mudxySignif =-9999;
	 float mudzSignif = -9999;
	 if (isPVfound) 
	   {
	     mudxySignif = mutrack->dxy(primVertex.position())/(mutrack->dxyError());
	     mudzSignif = mutrack->dz(primVertex.position())/(mutrack->dzError()); 
	   }
	 sipXYmuCheckNoIPTools[indexmu] = float(mudxySignif);
	 sipZmuCheckNoIPTools[indexmu] = float(mudzSignif);
	 //
	 //---------------other check / BS ------------------------------
	 float mudxySignifBS =-9999;
	 float mudzSignifBS = -9999;
	 if (isBSvalid) 
	   {
	     mudxySignifBS = mutrack->dxy(BSribbon.position())/(mutrack->dxyError());
	     mudzSignifBS = mutrack->dz(BSribbon.position())/(mutrack->dzError()); 
	   }
	 sipXYBSmuCheckNoIPTools[indexmu] = float(mudxySignifBS);
	 sipZBSmuCheckNoIPTools[indexmu] = float(mudzSignifBS);
	 //
	 //---------------------------------------------------------
	 //
	 // Looking for the associated closest PV :
	 //----------------------------------------
	 //
	 double dzMin = -9999.;
	 float muSignif3Dmindz = -9999;
	 float muSignif2Dmindz = -9999;
	 //
	 //if (isPVfound && isPVgood)
	 if (isPVfound)
	   {
            for (size_t indx = 0; indx != recoPVCollection->size(); indx++)
    	       {
	        const reco::Vertex &tempVtx = (*recoPVCollection)[indx];
		//
	        if ( fabs(mutrack->dz(tempVtx.position())) < fabs(dzMin) )
	          {
	           dzMin = mutrack->dz(tempVtx.position());
		   cout << "dzMin = " << dzMin << endl;
		   muSignif3Dmindz = float( getMuDauIP3Dsignif(mutrack,trackBuilder,tempVtx) );
		   muSignif2Dmindz = float( getMuDauIP2Dsignif(mutrack,trackBuilder,tempVtx) );
	          }
		//
	       }
	    //
	   }
	 //
	 double dzMinWBS = -9999.;
	 float muSignif3DmindzWithBS = -9999;
	 float muSignif2DmindzWithBS = -9999;
	 //
	 //if (isPVwBSfound && isPVwBSgood)
	 if (isPVwBSfound )
	   {
            for (size_t jndx = 0; jndx != recoPVCollectionWithBS->size(); jndx++)
    	       {
	        const reco::Vertex &tempVtxWithBS = (*recoPVCollectionWithBS)[jndx];
		//
	        if ( fabs(mutrack->dz(tempVtxWithBS.position())) < fabs(dzMinWBS))
	          {
	           dzMinWBS = mutrack->dz(tempVtxWithBS.position());
		   muSignif3DmindzWithBS = float( getMuDauIP3Dsignif(mutrack,trackBuilder,tempVtxWithBS) );
		   muSignif2DmindzWithBS = float( getMuDauIP2Dsignif(mutrack,trackBuilder,tempVtxWithBS) );
	          }
		//
	       }
	    //
	   }
	 //
	 //////////////////////////////////////////
	 //
	 // Best PV:
	 sip3Dmu[indexmu] = float( getMuDauIP3Dsignif(mutrack,trackBuilder,primVertex) );
	 sip2Dmu[indexmu] = float( getMuDauIP2Dsignif(mutrack,trackBuilder,primVertex) );
	 //
	 // Best PV with BS constraint:
	 sip3Dmu_wbs[indexmu] = float( getMuDauIP3Dsignif(mutrack,trackBuilder,primVertexWithBS) );
	 sip2Dmu_wbs[indexmu] = float( getMuDauIP2Dsignif(mutrack,trackBuilder,primVertexWithBS) );
	 //
	 // BeamSpot:
         sip3Dmu_bs[indexmu] = float( getMuDauIP3Dsignif(mutrack,trackBuilder,BSribbon) );
	 sip2Dmu_bs[indexmu] = float( getMuDauIP2Dsignif(mutrack,trackBuilder,BSribbon) );
	 //
	 // Closest PV:
         sip3Dmu_mindz[indexmu] = muSignif3Dmindz;
         sip2Dmu_mindz[indexmu] = muSignif2Dmindz;
         //
	 // Closest PV with BS constraint
	 sip3Dmu_wbs_mindz[indexmu] = muSignif3DmindzWithBS;
	 sip2Dmu_wbs_mindz[indexmu] = muSignif2DmindzWithBS;
         //
	 cout << "MUON SIP3D (bestPV, bestPVwBS, BS, closestPV, closestPVwBS) = " << sip3Dmu[indexmu] << 
	          ", " << sip3Dmu_wbs[indexmu] << ", "<< sip3Dmu_bs[indexmu] << 
                  ", " << sip3Dmu_mindz[indexmu] << ", " << sip3Dmu_wbs_mindz[indexmu] << endl;
	 //
	 cout << "MUON SIP2D (bestPV, bestPVwBS, BS, closestPV, closestPVwBS) = " << sip2Dmu[indexmu] << 
	          ", " << sip2Dmu_wbs[indexmu] << ", "<< sip2Dmu_bs[indexmu] << 
                  ", " << sip2Dmu_mindz[indexmu] << ", " << sip2Dmu_wbs_mindz[indexmu] << endl;
	 //
	 cout << "CHECK ALGO: MUON SIP3D (bestPV,bestPV-check) = " << sip3Dmu[indexmu] << ", " << sip3DmuCheck[indexmu] << endl;
	 cout << "CHECK ALGO: MUON SIP2D (bestPV,bestPV-checkNOIPTools,BS-checkNOIPTools) = " << sip2Dmu[indexmu] << 
	         ", " << sipXYmuCheckNoIPTools[indexmu] << ", " << sipXYBSmuCheckNoIPTools[indexmu] << endl;
	 //
         indexmu++;
	 //
        }
      //
      // inserting & filling maps
      //
      fillerMuSip3DPV.insert(muCandidates, sip3Dmu.begin(), sip3Dmu.end());
      fillerMuSip3DPVCheck.insert(muCandidates, sip3DmuCheck.begin(), sip3DmuCheck.end());
      fillerMuSip3DPVwbs.insert(muCandidates, sip3Dmu_wbs.begin(), sip3Dmu_wbs.end());
      fillerMuSip3DBS.insert(muCandidates, sip3Dmu_bs.begin(), sip3Dmu_bs.end());
      fillerMuSip3DcloPV.insert(muCandidates, sip3Dmu_mindz.begin(), sip3Dmu_mindz.end());
      fillerMuSip3DcloPVwbs.insert(muCandidates, sip3Dmu_wbs_mindz.begin(), sip3Dmu_wbs_mindz.end());
      //
      fillerMuSip2DPV.insert(muCandidates, sip2Dmu.begin(), sip2Dmu.end());
      fillerMuSip2DPVwbs.insert(muCandidates, sip2Dmu_wbs.begin(), sip2Dmu_wbs.end());
      fillerMuSip2DBS.insert(muCandidates, sip2Dmu_bs.begin(), sip2Dmu_bs.end());
      fillerMuSip2DcloPV.insert(muCandidates, sip2Dmu_mindz.begin(), sip2Dmu_mindz.end());
      fillerMuSip2DcloPVwbs.insert(muCandidates, sip2Dmu_wbs_mindz.begin(), sip2Dmu_wbs_mindz.end());
      //
      fillerMuSipXYPVCheckNoIPTools.insert(muCandidates, sipXYmuCheckNoIPTools.begin(), sipXYmuCheckNoIPTools.end());    
      fillerMuSipXYBSCheckNoIPTools.insert(muCandidates, sipXYBSmuCheckNoIPTools.begin(), sipXYBSmuCheckNoIPTools.end());      
      fillerMuSipZPVCheckNoIPTools.insert(muCandidates, sipZmuCheckNoIPTools.begin(), sipZmuCheckNoIPTools.end());       
      fillerMuSipZBSCheckNoIPTools.insert(muCandidates, sipZBSmuCheckNoIPTools.begin(), sipZBSmuCheckNoIPTools.end());      
      //
      fillerMuSip3DPV.fill();
      fillerMuSip3DPVCheck.fill();
      fillerMuSip3DPVwbs.fill();
      fillerMuSip3DBS.fill();
      fillerMuSip3DcloPV.fill();
      fillerMuSip3DcloPVwbs.fill();
      //
      fillerMuSip2DPV.fill();
      fillerMuSip2DPVwbs.fill();
      fillerMuSip2DBS.fill();
      fillerMuSip2DcloPV.fill();
      fillerMuSip2DcloPVwbs.fill();
      //
      fillerMuSipXYPVCheckNoIPTools.fill();
      fillerMuSipXYBSCheckNoIPTools.fill();
      fillerMuSipZPVCheckNoIPTools.fill();
      fillerMuSipZBSCheckNoIPTools.fill();
      //
    }
  //
  // Electrons 
  //===========
  //
  if (decaychannel=="2e2mu" || decaychannel=="4e")
    {    
     Handle<edm::View<GsfElectron> > eleCandidates; 
     iEvent.getByLabel(electronTag_.label(),eleCandidates);
     //
     size_t nele = eleCandidates->size();
     //
     //std::vector<float> sigele(nele);
     std::vector<float> sip3Dele(nele);
     std::vector<float> sip3DeleCheck(nele);
     std::vector<float> sip3Dele_wbs(nele);
     std::vector<float> sip3Dele_bs(nele);
     std::vector<float> sip3Dele_mindz(nele);
     std::vector<float> sip3Dele_wbs_mindz(nele);
     //
     std::vector<float> sip2Dele(nele);
     std::vector<float> sip2Dele_wbs(nele);
     std::vector<float> sip2Dele_bs(nele);
     std::vector<float> sip2Dele_mindz(nele);
     std::vector<float> sip2Dele_wbs_mindz(nele);
     //
     std::vector<float> sipXYeleCheckNoIPTools(nele);
     std::vector<float> sipXYBSeleCheckNoIPTools(nele);
     std::vector<float> sipZeleCheckNoIPTools(nele);
     std::vector<float> sipZBSeleCheckNoIPTools(nele);
     //
     unsigned int indexele=0;
     //
     for (edm::View<reco::GsfElectron>::const_iterator eleCand = eleCandidates->begin(); eleCand != eleCandidates->end(); ++eleCand)
        {
         GsfTrackRef eletrack = eleCand->get<GsfTrackRef>();
	 reco::TransientTrack eleTransTrack = trackBuilder->build( eletrack ) ;
	 //	 
	 const GlobalVector pelevec = GlobalVector(eletrack->px(),eletrack->py(),eletrack->pz());
	 //
	 //-----------------checks--------------------------------
	 std::pair<bool,Measurement1D> eleIP3DmindzCheck;
	 float eleSignif3DmindzCheck = -9999;
	 //if (isPVfound && isPVgood)
	 if (isPVfound)
	   {
	     std::pair<bool,Measurement1D> eleIP3DmindzCheck = IPTools::signedImpactParameter3D(eleTransTrack,pelevec,primVertex);
	     eleSignif3DmindzCheck = eleIP3DmindzCheck.second.significance();
	   }
	 sip3DeleCheck[indexele] = float(eleSignif3DmindzCheck); 
	 //
	 //---------------other check / PV ------------------------------
	 float eledxySignif =-9999;
	 float eledzSignif = -9999;
	 if (isPVfound) 
	   {
	     eledxySignif = eletrack->dxy(primVertex.position())/(eletrack->dxyError());
	     eledzSignif = eletrack->dz(primVertex.position())/(eletrack->dzError()); 
	   }
	 sipXYeleCheckNoIPTools[indexele] = float(eledxySignif);
	 sipZeleCheckNoIPTools[indexele] = float(eledzSignif);
	 //
	 //---------------other check / BS ------------------------------
	 float eledxySignifBS =-9999;
	 float eledzSignifBS = -9999;
	 if (isBSvalid) 
	   {
	     eledxySignifBS = eletrack->dxy(BSribbon.position())/(eletrack->dxyError());
	     eledzSignifBS = eletrack->dz(BSribbon.position())/(eletrack->dzError()); 
	   }
	 sipXYBSeleCheckNoIPTools[indexele] = float(eledxySignifBS);
	 sipZBSeleCheckNoIPTools[indexele] = float(eledzSignifBS);
	 //
	 //
	 //---------------------------------------------------------
	 //
	 // Looking for the associated closest PV :
	 //----------------------------------------
	 //
	 double dzMin1 = -9999.;
	 float eleSignif3Dmindz = -9999;
	 float eleSignif2Dmindz = -9999;
	 //
	 //if (isPVfound && isPVgood)
	 if (isPVfound)
	   {
            for (size_t kndx = 0; kndx != recoPVCollection->size(); kndx++)
    	       {
	        const reco::Vertex &tempVtx = (*recoPVCollection)[kndx];
		//
	        if ( fabs(eletrack->dz(tempVtx.position())) < fabs(dzMin1) )
	          {
	           dzMin1 = eletrack->dz(tempVtx.position());
		   eleSignif3Dmindz = float( getEleDauIP3Dsignif(eletrack,trackBuilder,tempVtx) );
		   eleSignif2Dmindz = float( getEleDauIP2Dsignif(eletrack,trackBuilder,tempVtx) );
	          }
		//
	       }
	   }
	 //
	 cout << "eleSignif3Dmindz = " << eleSignif3Dmindz << endl;
	 //
	 double dzMinWBS1 = -9999.;
	 float eleSignif3DmindzWithBS = -9999;
	 float eleSignif2DmindzWithBS = -9999;
	 //
	 //if (isPVwBSfound && isPVwBSgood)
	 if (isPVwBSfound )
	   {
            for (size_t lndx = 0; lndx != recoPVCollectionWithBS->size(); lndx++)
    	       {
	        const reco::Vertex &tempVtxWithBS = (*recoPVCollectionWithBS)[lndx];
		//
	        if ( fabs(eletrack->dz(tempVtxWithBS.position())) < fabs(dzMinWBS1))
	          {
	           dzMinWBS1 = eletrack->dz(tempVtxWithBS.position());
		   eleSignif3DmindzWithBS = float( getEleDauIP3Dsignif(eletrack,trackBuilder,tempVtxWithBS) );
		   eleSignif2DmindzWithBS = float( getEleDauIP2Dsignif(eletrack,trackBuilder,tempVtxWithBS) );
	          }
		//
	       }
	    //
	   }
	 //
	 //////////////////////////////////////////
	 //
 	 // Best PV:
	 sip3Dele[indexele] = float( getEleDauIP3Dsignif(eletrack,trackBuilder,primVertex) );
	 sip2Dele[indexele] = float( getEleDauIP2Dsignif(eletrack,trackBuilder,primVertex) );
	 //
	 // Best PV with BS constraint:
	 sip3Dele_wbs[indexele] = float( getEleDauIP3Dsignif(eletrack,trackBuilder,primVertexWithBS) );
	 sip2Dele_wbs[indexele] = float( getEleDauIP2Dsignif(eletrack,trackBuilder,primVertexWithBS) );
	 //
	 // BeamSpot:
         sip3Dele_bs[indexele] = float( getEleDauIP3Dsignif(eletrack,trackBuilder,BSribbon) );
	 sip2Dele_bs[indexele] = float( getEleDauIP2Dsignif(eletrack,trackBuilder,BSribbon) );
	 //
	 // Closest PV:
         sip3Dele_mindz[indexele] = eleSignif3Dmindz;
         sip2Dele_mindz[indexele] = eleSignif2Dmindz;
         //
	 // Closest PV with BS constraint
	 sip3Dele_wbs_mindz[indexele] = eleSignif3DmindzWithBS;
	 sip2Dele_wbs_mindz[indexele] = eleSignif2DmindzWithBS;
         //
	 cout << "ELECTRON SIP3D (bestPV, bestPVwBS, BS, closestPV, closestPVwBS) = " << sip3Dele[indexele] << 
	          ", " << sip3Dele_wbs[indexele] << ", "<< sip3Dele_bs[indexele] << 
                  ", " << sip3Dele_mindz[indexele] << ", " << sip3Dele_wbs_mindz[indexele] << endl;
	 //
	 cout << "ELECTRON SIP2D (bestPV, bestPVwBS, BS, closestPV, closestPVwBS) = " << sip2Dele[indexele] << 
	          ", " << sip2Dele_wbs[indexele] << ", "<< sip2Dele_bs[indexele] << 
                  ", " << sip2Dele_mindz[indexele] << ", " << sip2Dele_wbs_mindz[indexele] << endl;
	 //
	 cout << "CHECK ALGO: ELECTRON SIP3D (bestPV,bestPV-check) = " << sip3Dele[indexele] << ", " << sip3DeleCheck[indexele] << endl;
	 cout << "CHECK ALGO: ELECTRON SIP2D (bestPV,bestPV-checkNOIPTools,BS-checkNOIPTools) = " << sip2Dele[indexele] << 
	         ", " << sipXYeleCheckNoIPTools[indexele] << ", " << sipXYBSeleCheckNoIPTools[indexele] << endl;
	 //   
         indexele++;
         //
        } 
     //
     // inserting & filling maps
     //
     fillerEleSip3DPV.insert(eleCandidates,sip3Dele.begin(), sip3Dele.end());
     fillerEleSip3DPVCheck.insert(eleCandidates,sip3DeleCheck.begin(), sip3DeleCheck.end());
     fillerEleSip3DPVwbs.insert(eleCandidates,sip3Dele_wbs.begin(), sip3Dele_wbs.end());
     fillerEleSip3DBS.insert(eleCandidates,sip3Dele_bs.begin(), sip3Dele_bs.end());
     fillerEleSip3DcloPV.insert(eleCandidates,sip3Dele_mindz.begin(), sip3Dele_mindz.end());
     fillerEleSip3DcloPVwbs.insert(eleCandidates,sip3Dele_wbs_mindz.begin(), sip3Dele_wbs_mindz.end());
     //
     fillerEleSip2DPV.insert(eleCandidates,sip2Dele.begin(), sip2Dele.end());
     fillerEleSip2DPVwbs.insert(eleCandidates,sip2Dele_wbs.begin(), sip2Dele_wbs.end());
     fillerEleSip2DBS.insert(eleCandidates,sip2Dele_bs.begin(), sip2Dele_bs.end());
     fillerEleSip2DcloPV.insert(eleCandidates,sip2Dele_mindz.begin(), sip2Dele_mindz.end());
     fillerEleSip2DcloPVwbs.insert(eleCandidates,sip2Dele_wbs_mindz.begin(), sip2Dele_wbs_mindz.end());
     //
     fillerEleSipXYPVCheckNoIPTools.insert(eleCandidates,sipXYeleCheckNoIPTools.begin(), sipXYeleCheckNoIPTools.end());
     fillerEleSipXYBSCheckNoIPTools.insert(eleCandidates,sipXYBSeleCheckNoIPTools.begin(), sipXYBSeleCheckNoIPTools.end());
     fillerEleSipZPVCheckNoIPTools.insert(eleCandidates,sipZeleCheckNoIPTools.begin(), sipZeleCheckNoIPTools.end());
     fillerEleSipZBSCheckNoIPTools.insert(eleCandidates,sipZBSeleCheckNoIPTools.begin(), sipZBSeleCheckNoIPTools.end());
     //
     fillerEleSip3DPV.fill();
     fillerEleSip3DPVCheck.fill();
     fillerEleSip3DPVwbs.fill();
     fillerEleSip3DBS.fill();
     fillerEleSip3DcloPV.fill();
     fillerEleSip3DcloPVwbs.fill();
     //
     fillerEleSip2DPV.fill();
     fillerEleSip2DPVwbs.fill();
     fillerEleSip2DBS.fill();
     fillerEleSip2DcloPV.fill();
     fillerEleSip2DcloPVwbs.fill();
     //
     fillerEleSipXYPVCheckNoIPTools.fill();
     fillerEleSipXYBSCheckNoIPTools.fill();
     fillerEleSipZPVCheckNoIPTools.fill();
     fillerEleSipZBSCheckNoIPTools.fill();
     //
  }
  //
  //=========
  //
  if (decaychannel=="2e2mu" || decaychannel=="4mu")
    {
      iEvent.put( sip3DPVMuMap, "Sip3DPVMuMap" );              // mu sip-3d-iptools w.r.t. PV
      iEvent.put( sip3DPVCheckMuMap, "Sip3DPVCheckMuMap" );    // mu sip-3d-iptools-check-method w.r.t. PV
      iEvent.put( sip3DPVwbsMuMap, "Sip3DPVwbsMuMap" );        // mu sip-3d-iptools w.r.t. PV with BS constraint
      iEvent.put( sip3DBSMuMap, "Sip3DBSMuMap" );              // mu sip-3d-iptools w.r.t. BS
      iEvent.put( sip3DcloPVMuMap, "Sip3DcloPVMuMap" );        // mu sip-3d-iptools w.r.t. closest-in-Z PV
      iEvent.put( sip3DcloPVwbsMuMap, "Sip3DcloPVwbsMuMap" );  // mu sip-3d-iptools w.r.t. closest-in-Z PV with BS constraint
      //
      iEvent.put( sip2DPVMuMap, "Sip2DPVMuMap" );              // mu sip-2d-iptools w.r.t. PV
      iEvent.put( sip2DPVwbsMuMap, "Sip2DPVwbsMuMap" );        // mu sip-2d-iptools w.r.t. PV with BS constraint
      iEvent.put( sip2DBSMuMap, "Sip2DBSMuMap" );              // mu sip-2d-iptools w.r.t. BS
      iEvent.put( sip2DcloPVMuMap, "Sip2DcloPVMuMap" );        // mu sip-2d-iptools w.r.t. closest-in-Z PV
      iEvent.put( sip2DcloPVwbsMuMap, "Sip2DcloPVwbsMuMap" );  // mu sip-2d-iptools w.r.t. closest-in-Z PV with BS constraint
      //
      iEvent.put( sipXYPVCheckNoIPToolsMuMap, "SipXYPVCheckNoIPToolsMuMap" );  // mu sip-XY-NOiptools w.r.t. PV
      iEvent.put( sipXYBSCheckNoIPToolsMuMap, "SipXYBSCheckNoIPToolsMuMap" );  // mu sip-XY-NOiptools w.r.t. BS
      iEvent.put( sipZPVCheckNoIPToolsMuMap, "SipZPVCheckNoIPToolsMuMap" );    // mu sip-Z-NOiptools w.r.t. PV
      iEvent.put( sipZBSCheckNoIPToolsMuMap, "SipZBSCheckNoIPToolsMuMap" );    // mu sip-Z-NOiptools w.r.t. BS
      //
    }
  if (decaychannel=="2e2mu" || decaychannel=="4e")
    {
      iEvent.put( sip3DPVEleMap, "Sip3DPVEleMap" );             // ele sip-3d-iptools w.r.t. PV
      iEvent.put( sip3DPVCheckEleMap, "Sip3DPVCheckEleMap" );   // ele sip-3d-iptools-check-method w.r.t. PV
      iEvent.put( sip3DPVwbsEleMap,"Sip3DPVwbsEleMap" );        // ele sip-3d-iptools w.r.t. PV with BS constraint
      iEvent.put( sip3DBSEleMap, "Sip3DBSEleMap" );             // ele sip-3d-iptools w.r.t. BS
      iEvent.put( sip3DcloPVEleMap, "Sip3DcloPVEleMap" );       // ele sip-3d-iptools w.r.t. closest-in-Z PV
      iEvent.put( sip3DcloPVwbsEleMap, "Sip3DcloPVwbsEleMap" ); // ele sip-3d-iptools w.r.t. closest-in-Z PV with BS constraint
      //
      iEvent.put( sip2DPVEleMap, "Sip2DPVEleMap" );             // ele sip-2d-iptools w.r.t. PV
      iEvent.put( sip2DPVwbsEleMap, "Sip2DPVwbsEleMap" );       // ele sip-2d-iptools w.r.t. PV with BS constraint
      iEvent.put( sip2DBSEleMap, "Sip2DBSEleMap" );             // ele sip-2d-iptools w.r.t. BS
      iEvent.put( sip2DcloPVEleMap, "Sip2DcloPVEleMap" );       // ele sip-2d-iptools w.r.t. closest-in-Z PV
      iEvent.put( sip2DcloPVwbsEleMap, "Sip2DcloPVwbsEleMap" ); // ele sip-2d-iptools w.r.t. closest-in-Z PV with BS constraint
      //
      iEvent.put( sipXYPVCheckNoIPToolsEleMap, "SipXYPVCheckNoIPToolsEleMap" ); // ele sip-XY-NOiptools w.r.t. PV
      iEvent.put( sipXYBSCheckNoIPToolsEleMap, "SipXYBSCheckNoIPToolsEleMap" ); // ele sip-XY-NOiptools w.r.t. BS
      iEvent.put( sipZPVCheckNoIPToolsEleMap, "SipZPVCheckNoIPToolsEleMap" );   // ele sip-Z-NOiptools w.r.t. PV
      iEvent.put( sipZBSCheckNoIPToolsEleMap, "SipZBSCheckNoIPToolsEleMap" );   // ele sip-Z-NOiptools w.r.t. BS
      //
    }
  //
}

///////////////////////////////////////////////////////////////////////////
//------------- helper functions:
///////////////////////////////////////////////////////////////////////////

// MUONS
////////

// the following function, given a muon, provides its signed IP3D significance

float HZZ4LeptonsVtxProducer::getMuDauIP3Dsignif(reco::TrackRef muDauTrack,
                                                 const ESHandle<TransientTrackBuilder> TrackBuilder,
                                                 reco::Vertex primaryVertex){
  //
  GlobalPoint pVertex(primaryVertex.position().x(),primaryVertex.position().y(),primaryVertex.position().z());
  // 
  reco::TransientTrack muDauTransTrack = TrackBuilder->build( muDauTrack ) ;
  //
  TrajectoryStateOnSurface muDauTSOS = muDauTransTrack.stateOnSurface(pVertex);
  //
  float muDauSignif3D;
  //
  if (!muDauTSOS.isValid())
    {
     muDauSignif3D = -9999;
     //
    } else {
            std::pair<bool,Measurement1D> muDauIPpair;
            muDauIPpair = IPTools::signedImpactParameter3D(muDauTransTrack, muDauTSOS.globalDirection(), primaryVertex);
            if (muDauIPpair.first)
              {
               // taking absolute value of significance: muDauSignif3D = fabs( muDauIPpair.second.significance() );
	       // no: get the sign as well:
	       muDauSignif3D = muDauIPpair.second.significance();
	       //
              } else {
                      muDauSignif3D = -999;
              }
    }
  //
  return muDauSignif3D;
}

/////////////

// the following function, given a muon, provides its signed IP2D significance

float HZZ4LeptonsVtxProducer::getMuDauIP2Dsignif(reco::TrackRef muDauTrack,
                                                 const ESHandle<TransientTrackBuilder> TrackBuilder,
                                                 reco::Vertex primaryVertex){
  //
  GlobalPoint pVertex(primaryVertex.position().x(),primaryVertex.position().y(),primaryVertex.position().z());
  // 
  reco::TransientTrack muDauTransTrack = TrackBuilder->build( muDauTrack ) ;
  //
  TrajectoryStateOnSurface muDauTSOS = muDauTransTrack.stateOnSurface(pVertex);
  //
  float muDauSignif2D;
  //
  if (!muDauTSOS.isValid())
    {
     muDauSignif2D = -9999;
     //
    } else {
            std::pair<bool,Measurement1D> muDauIPpair;
            muDauIPpair = IPTools::signedTransverseImpactParameter(muDauTransTrack, muDauTSOS.globalDirection(), primaryVertex);
            if (muDauIPpair.first)
              {
	       // getting the sign as well:
	       muDauSignif2D = muDauIPpair.second.significance();
	       //
              } else {
                      muDauSignif2D = -999;
              }
    }
  //
  return muDauSignif2D;
}

//

// ELECTRONS
////////////////////////////
//
float HZZ4LeptonsVtxProducer::getEleDauIP3Dsignif(reco::GsfTrackRef eleDauTrack,
                                                  const ESHandle<TransientTrackBuilder> TrackBuilder,
                                                  reco::Vertex primaryVertex){
  //
  GlobalPoint pVertex(primaryVertex.position().x(),primaryVertex.position().y(),primaryVertex.position().z());
  // 
  reco::TransientTrack eleDauTransTrack = TrackBuilder->build(eleDauTrack);
  TrajectoryStateOnSurface eleDauTSOS = eleDauTransTrack.stateOnSurface(pVertex);
  float eleDauSignif3D;
  //
  if (!eleDauTSOS.isValid())
    {
     eleDauSignif3D = -9999;
     //
    } else {
            std::pair<bool,Measurement1D> eleDauIPpair;
            eleDauIPpair = IPTools::signedImpactParameter3D(eleDauTransTrack, eleDauTSOS.globalDirection(), primaryVertex);
            if (eleDauIPpair.first)
              {
               // taking absolute value of significance: eleDauSignif3D = fabs( eleDauIPpair.second.significance() );
	       // no: get the sign as well:
	       eleDauSignif3D = eleDauIPpair.second.significance();
	       //
              } else {
                      eleDauSignif3D = -999;
              }
    }
  //
  return eleDauSignif3D;
}

//////////////////// the following function, given an electron, provides its signed IP2D significance

float HZZ4LeptonsVtxProducer::getEleDauIP2Dsignif(reco::GsfTrackRef eleDauTrack,
                                                  const ESHandle<TransientTrackBuilder> TrackBuilder,
                                                  reco::Vertex primaryVertex){
  //
  GlobalPoint pVertex(primaryVertex.position().x(),primaryVertex.position().y(),primaryVertex.position().z());
  // 
  reco::TransientTrack eleDauTransTrack = TrackBuilder->build(eleDauTrack);
  TrajectoryStateOnSurface eleDauTSOS = eleDauTransTrack.stateOnSurface(pVertex);
  //
  float eleDauSignif2D;
  //
  if (!eleDauTSOS.isValid())
    {
     eleDauSignif2D = -9999;
     //
    } else {
            std::pair<bool,Measurement1D> eleDauIPpair;
            eleDauIPpair = IPTools::signedTransverseImpactParameter(eleDauTransTrack, eleDauTSOS.globalDirection(), primaryVertex);
            if (eleDauIPpair.first)
              {
	       // getting the sign as well:
	       eleDauSignif2D = eleDauIPpair.second.significance();
	       //
              } else {
                      eleDauSignif2D = -999;
              }
    }
  //
  return eleDauSignif2D;
}

////////////////////////////


