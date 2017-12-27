/**\class HZZ4LeptonsTipLipToVtxProducer 
 *
 * Original Author:  Alexis Pompili - Univ. of Bari & INFN
 *       modified by N. De Filippis - Politecnico di Bari & INFN          
 * (from HiggsToZZ2e2m/src/CmsCandidateFiller.cc and some core code borrowed from C.Charlot-S.Baffioni)
 *
 * Computes for each lepton the Transverse & Longitudinal IP w.r.t to PrimaryVertex (Tip & Lip)
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsTipLipToVtxProducer.h"


// Tracker tracks 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GeometryVector/interface/Basic3DVector.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// Transient tracks
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

// Vertex utilities
#include "DataFormats/VertexReco/interface/Vertex.h"

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
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/GsfTools/interface/GsfPropagatorAdapter.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

// Other specific
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"

// Maps
#include "DataFormats/Common/interface/ValueMap.h"
#include "TrackingTools/IPTools/interface/IPTools.h"


//-- system include files
#include <memory>
#include <vector>
#include <string>
//
//-- namespaces
using namespace edm;
using namespace std;
using namespace reco;

//
//-- constructor
//
HZZ4LeptonsTipLipToVtxProducer::HZZ4LeptonsTipLipToVtxProducer(const edm::ParameterSet& pset) {
  
  // Decay Channel
  decaychannel                                                    = pset.getParameter<std::string>("decaychannel");
  if (decaychannel=="2e2mu" || decaychannel=="4mu")  muonTag_     = consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("MuonsLabel"));
  if (decaychannel=="2e2mu" || decaychannel=="4e")   electronTag_ = consumes<edm::View<reco::GsfElectron> >(pset.getParameter<edm::InputTag>("ElectronsLabel"));
  vertexTag_                                                      = consumes<std::vector<reco::Vertex> >(pset.getParameter<edm::InputTag>("VertexLabel"));
  offlineBeamSpot_                                                = consumes<reco::BeamSpot>(pset.getParameter<edm::InputTag>("BeamSpotLabel"));
  useBeamSpot_                                                    = pset.getParameter<bool>("useBeamSpot");
  debug	=	pset.getUntrackedParameter<bool> ("debug", false);

  //
  string alias;
  string iName = "TipLipToVtx";

  //
  if (decaychannel=="2e2mu" || decaychannel=="4e")  {
    produces<edm::ValueMap<float> >( "TipEleMap");
    produces<edm::ValueMap<float> >( "LipEleMap");

    produces<edm::ValueMap<float> >( "TipValueEleMap");
    produces<edm::ValueMap<float> >( "LipValueEleMap");
    produces<edm::ValueMap<float> >( "TipErrorEleMap");
    produces<edm::ValueMap<float> >( "LipErrorEleMap");

  }
  

  if (decaychannel=="2e2mu" || decaychannel=="4mu"){
    produces<edm::ValueMap<float> >( "TipMuMap");
    produces<edm::ValueMap<float> >( "LipMuMap");

    produces<edm::ValueMap<float> >( "TipValueMuMap");
    produces<edm::ValueMap<float> >( "LipValueMuMap");
    produces<edm::ValueMap<float> >( "TipErrorMuMap");
    produces<edm::ValueMap<float> >( "LipErrorMuMap");

  }

}

//-- destructor
//
HZZ4LeptonsTipLipToVtxProducer::~HZZ4LeptonsTipLipToVtxProducer() {
}

void HZZ4LeptonsTipLipToVtxProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //
  auto_ptr<edm::ValueMap<float> >     twodTipMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuTip(*twodTipMuMap);
  auto_ptr<edm::ValueMap<float> >     twodLipMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuLip(*twodLipMuMap);
  
  auto_ptr<edm::ValueMap<float> >     twodTipValueMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuTipValue(*twodTipValueMuMap);
  auto_ptr<edm::ValueMap<float> >     twodLipValueMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuLipValue(*twodLipValueMuMap);

  auto_ptr<edm::ValueMap<float> >     twodTipErrorMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuTipError(*twodTipErrorMuMap);
  auto_ptr<edm::ValueMap<float> >     twodLipErrorMuMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerMuLipError(*twodLipErrorMuMap);
  
  
  auto_ptr<edm::ValueMap<float> > twodTipEleMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerEleTip(*twodTipEleMap);
  auto_ptr<edm::ValueMap<float> > twodLipEleMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerEleLip(*twodLipEleMap);
  
  auto_ptr<edm::ValueMap<float> >     twodTipValueEleMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerEleTipValue(*twodTipValueEleMap);
  auto_ptr<edm::ValueMap<float> >     twodLipValueEleMap(new edm::ValueMap<float> ());
  edm::ValueMap<float> ::Filler fillerEleLipValue(*twodLipValueEleMap);
    
  auto_ptr<edm::ValueMap<float> >     twodTipErrorEleMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerEleTipError(*twodTipErrorEleMap);
  auto_ptr<edm::ValueMap<float> >     twodLipErrorEleMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerEleLipError(*twodLipErrorEleMap);

  //
  //////////////////////////////////
  //
  // Get the B-field
  ESHandle<MagneticField> B;
  iSetup.get<IdealMagneticFieldRecord>().get( B );
  
  // get the track builder
  ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);


  // Beamspot 
  Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(offlineBeamSpot_,recoBeamSpotHandle);
  const BeamSpot bs = *recoBeamSpotHandle;

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
  bool pvfound = (recoPVCollection->size() != 0);
  if (pvfound)
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
  Basic3DVector<double> pVertexErr(primVertex.xError(),primVertex.yError(),primVertex.zError());
  //
  /////////////////////////////////
  //
  //--------------- for propagation
  //
  // electrons:
  //===========
  const MultiTrajectoryStateTransform *theMtsTransform = new MultiTrajectoryStateTransform;
  const GsfPropagatorAdapter *theGeomPropBw = new GsfPropagatorAdapter(AnalyticalPropagator(B.product(),oppositeToMomentum)); 
  //
  // muons:
  //=======
  const TrajectoryStateTransform *theTransform = new TrajectoryStateTransform;
  Propagator *thePropBw = new AnalyticalPropagator(B.product(),oppositeToMomentum);                    
  // 
//   ESHandle<TrackerGeometry> trackerHandle;
//   iSetup.get<TrackerDigiGeometryRecord>().get(trackerHandle); 
  //
  float muTip,muLip,muSTip,muSLip,muTipSignif,muLipSignif;
  float eleTip,eleLip,eleSTip,eleSLip,eleTipSignif,eleLipSignif;
  //
  //
  //--track refs
  TrackRef mutrack;
  GsfTrackRef eletrack;
  reco::TransientTrack mutranstrack;
  reco::TransientTrack eletranstrack;

  //
  // Muons
  //====== 
  //
  if (decaychannel=="2e2mu" || decaychannel=="4mu"){
    
    edm::Handle<edm::View<Muon> > muCandidates;
    iEvent.getByToken(muonTag_,muCandidates);
    
    size_t nmu = muCandidates->size();
    std::vector<float> sigtipmu(nmu),siglipmu(nmu);
    std::vector<float> valuetipmu(nmu),valuelipmu(nmu);
    std::vector<float> errortipmu(nmu),errorlipmu(nmu);
    unsigned int indexmu=0;
    
    for (edm::View<reco::Muon>::const_iterator muCand = muCandidates->begin(); muCand != muCandidates->end(); ++muCand){
      //
      //mutrack = muCand->combinedMuon();
      mutrack = muCand->get<TrackRef>();
      if (mutrack.isNull()){
	if(debug) cout <<"tracker trackref is null since muon is STA" << endl;
	mutrack=muCand->get<TrackRef,reco::StandAloneMuonTag>();
      }  
      //
      ///////////////////////////
      //-- the following doesn't solve the problem for STA muon because 
      //   of "Invalid DetID: no GeomDet associated" when building TSOS
      //if (mutrack.isNull())
      //  {
      //   mutrack=muCand->get<TrackRef,reco::StandAloneMuonTag>();
      //  }
      ///////////////////////////
      //if(mutrack.isNull()){
      //mutrack = muCand->track();
      //if (mutrack.id().isValid()) {
      ////////////////////////////////	
      //

      mutranstrack = trackBuilder->build( mutrack );

      TrajectoryStateOnSurface innerMuTSOS;

      if (useBeamSpot_==true){ 
	//innerMuTSOS = mutranstrack.stateOnSurface(BSVertex);
	innerMuTSOS = IPTools::transverseExtrapolate(mutranstrack.impactPointState(), BSVertex, mutranstrack.field());
      } 
      else {
	//innerMuTSOS = mutranstrack.stateOnSurface(pVertex);
	innerMuTSOS = IPTools::transverseExtrapolate(mutranstrack.impactPointState(), pVertex, mutranstrack.field());
      } 

      if (innerMuTSOS.isValid() && !mutrack.isNull() ){    
	//
	//-- get propagated the inner TSOS to the PV:
	TrajectoryStateOnSurface vtxMuTSOS;
	if (useBeamSpot_==true){ 
	  vtxMuTSOS = TransverseImpactPointExtrapolator(*thePropBw).extrapolate(innerMuTSOS,BSVertex);
	} 
	else {
	  vtxMuTSOS = TransverseImpactPointExtrapolator(*thePropBw).extrapolate(innerMuTSOS,pVertex);
	}
	
	//		     
	if (!vtxMuTSOS.isValid()){		 
	  vtxMuTSOS = innerMuTSOS; //-protection for eventual failing extrapolation
	}
	//
	//
	//-- get the distances (transverse & longitudinal) between extrapolated position and PV position
	GlobalPoint mimpP = vtxMuTSOS.globalPosition();
	GlobalVector mdistV; 
	GlobalVector direction=vtxMuTSOS.globalDirection();

	if (useBeamSpot_==true){ 
	  mdistV = mimpP - BSVertex; 
	} 
	else {
	  mdistV = mimpP - pVertex; 
	}

	GlobalVector transverseIP(mdistV.x(),mdistV.y(),0.); 
	double ps = transverseIP.dot(direction);
	muTip = mdistV.perp()*((ps!=0)?ps/abs(ps):1.); // signed by definition
	muLip = mdistV.z();    // signed by definition
	//
	// compute full error calculation:
	//
	// - diagonal terms first:
	AlgebraicSymMatrix33 mvtxerrM; 
	if (useBeamSpot_==true){ 
	  mvtxerrM(0,0) = BSVertexErr.x()*BSVertexErr.x(); 
	  mvtxerrM(1,1) = BSVertexErr.y()*BSVertexErr.y();
	  mvtxerrM(2,2) = BSVertexErr.z()*BSVertexErr.z();
	} 
	else {
	  mvtxerrM(0,0) = pVertexErr.x()*pVertexErr.x(); 
	  mvtxerrM(1,1) = pVertexErr.y()*pVertexErr.y();
	  mvtxerrM(2,2) = pVertexErr.z()*pVertexErr.z();
	}
	//
	// - off-diagonal terms:
	AlgebraicSymMatrix33 merrorM = mvtxerrM + vtxMuTSOS.cartesianError().matrix().Sub<AlgebraicSymMatrix33>(0,0);
	AlgebraicVector2 mjacobianTranV;	
	AlgebraicVector1 mjacobianLongV;
	mjacobianTranV[0] = mdistV.x()/mdistV.perp();	
	mjacobianTranV[1] = mdistV.y()/mdistV.perp();
	mjacobianLongV[0] = 1.;	
	//			   
	//- errors:
	muSTip = sqrt(ROOT::Math::Similarity(merrorM.Sub<AlgebraicSymMatrix22>(0,0),mjacobianTranV));
	muSLip = sqrt(ROOT::Math::Similarity(merrorM.Sub<AlgebraicSymMatrix11>(2,2),mjacobianLongV));

	//
	muTipSignif=muTip/muSTip;
	muLipSignif=muLip/muSLip;
	sigtipmu[indexmu]=float(muTipSignif);
	siglipmu[indexmu]=float(muLipSignif);
	valuetipmu[indexmu]=float(muTip);
	valuelipmu[indexmu]=float(muLip);
	errortipmu[indexmu]=float(muSTip);
	errorlipmu[indexmu]=float(muSLip);
	if(debug) cout << "TipSignificance for muons= " << muTipSignif << endl;
      }
      //
      indexmu++;
    } //-- muon loop closed
    //
    // filling map
    fillerMuTip.insert(muCandidates, sigtipmu.begin(), sigtipmu.end());
    fillerMuLip.insert(muCandidates, siglipmu.begin(), siglipmu.end());
    fillerMuTip.fill();
    fillerMuLip.fill();

    fillerMuTipValue.insert(muCandidates, valuetipmu.begin(), valuetipmu.end());
    fillerMuLipValue.insert(muCandidates, valuelipmu.begin(), valuelipmu.end());
    fillerMuTipValue.fill();
    fillerMuLipValue.fill();

    fillerMuTipError.insert(muCandidates, errortipmu.begin(), errortipmu.end());
    fillerMuLipError.insert(muCandidates, errorlipmu.begin(), errorlipmu.end());
    fillerMuTipError.fill();
    fillerMuLipError.fill();

  } 
  
  //
  if (decaychannel=="2e2mu" || decaychannel=="4e"){

    Handle<edm::View<GsfElectron> > eleCandidates;
    iEvent.getByToken(electronTag_,eleCandidates);

    size_t nele = eleCandidates->size();
    std::vector<float> sigtipele(nele),siglipele(nele);
    std::vector<float> valuetipele(nele),valuelipele(nele);
    std::vector<float> errortipele(nele),errorlipele(nele);
    unsigned int indexele=0;
    
    for (edm::View<reco::GsfElectron>::const_iterator eleCand = eleCandidates->begin(); eleCand != eleCandidates->end(); ++eleCand){

      eletrack = eleCand->get<GsfTrackRef>();
      eletranstrack = trackBuilder->build( eletrack ) ;

      TrajectoryStateOnSurface innerTSOS;
      if (useBeamSpot_==true){ 
	//innerTSOS = eletranstrack.stateOnSurface(BSVertex);
	innerTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), BSVertex, eletranstrack.field());
      } 
      else {
	//innerTSOS = eletranstrack.stateOnSurface(pVertex);
	innerTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), pVertex, eletranstrack.field());
      }

      //
      if (innerTSOS.isValid()){
	//-- get propagated the inner TSOS to the PV:
	TrajectoryStateOnSurface vtxTSOS;
	
	if (useBeamSpot_==true){ 
	  vtxTSOS = TransverseImpactPointExtrapolator(*theGeomPropBw).extrapolate(innerTSOS,BSVertex);
	} 
	else {
	  vtxTSOS = TransverseImpactPointExtrapolator(*theGeomPropBw).extrapolate(innerTSOS,pVertex);
	}

	//		     
	if (!vtxTSOS.isValid()){		 
	  vtxTSOS = innerTSOS; //-protection for eventual failing extrapolation
	}       
	//
	//-- get the distances (transverse & longitudinal) between extrapolated position and PV position
	GlobalPoint impP = vtxTSOS.globalPosition();
	GlobalVector distV;
	GlobalVector direction=vtxTSOS.globalDirection();
      
	if (useBeamSpot_==true){ 
	  distV = impP - BSVertex; 	
	} 
	else {
	  distV = impP - pVertex; 
	}

	GlobalVector transverseIPele(distV.x(),distV.y(),0.); 
	double psele = transverseIPele.dot(direction);
	eleTip = distV.perp()*((psele!=0)?psele/abs(psele):1.); // signed by definition
	eleLip = distV.z();    // signed by definition
	//
	// compute full error calculation:
	//
	// - diagonal terms first:
	AlgebraicSymMatrix33 vtxerrM; 

	if (useBeamSpot_==true){ 
	  vtxerrM(0,0) = BSVertexErr.x()*BSVertexErr.x(); 
	  vtxerrM(1,1) = BSVertexErr.y()*BSVertexErr.y();
	  vtxerrM(2,2) = BSVertexErr.z()*BSVertexErr.z();
	} 
	else {
	  vtxerrM(0,0) = pVertexErr.x()*pVertexErr.x(); 
	  vtxerrM(1,1) = pVertexErr.y()*pVertexErr.y();
	  vtxerrM(2,2) = pVertexErr.z()*pVertexErr.z();
	}
		
	//
	// - off-diagonal terms:
	AlgebraicSymMatrix33 errorM = vtxerrM + vtxTSOS.cartesianError().matrix().Sub<AlgebraicSymMatrix33>(0,0);
	AlgebraicVector2 jacobianTranV;	
	AlgebraicVector1 jacobianLongV;
	jacobianTranV[0] = distV.x()/distV.perp();	
	jacobianTranV[1] = distV.y()/distV.perp();
	jacobianLongV[0] = 1.;	
	//			   
	//- errors:
	eleSTip = sqrt(ROOT::Math::Similarity(errorM.Sub<AlgebraicSymMatrix22>(0,0),jacobianTranV));
	eleSLip = sqrt(ROOT::Math::Similarity(errorM.Sub<AlgebraicSymMatrix11>(2,2),jacobianLongV));
	//
	eleTipSignif=eleTip/eleSTip;
	eleLipSignif=eleLip/eleSLip;
	sigtipele[indexele]=float(eleTipSignif);
	siglipele[indexele]=float(eleLipSignif);
	valuetipele[indexele]=float(eleTip);
	valuelipele[indexele]=float(eleLip);
	errortipele[indexele]=float(eleSTip);
	errorlipele[indexele]=float(eleSLip);

	if(debug) cout << "TipSignificance for electrons= " << eleTipSignif << endl;
      }
      //
      indexele++;
    } //-- ele loop closed
    //
    // filling map
    fillerEleTip.insert(eleCandidates, sigtipele.begin(), sigtipele.end());
    fillerEleLip.insert(eleCandidates, siglipele.begin(), siglipele.end());
    fillerEleTip.fill();
    fillerEleLip.fill();

    fillerEleTipValue.insert(eleCandidates, valuetipele.begin(), valuetipele.end());
    fillerEleLipValue.insert(eleCandidates, valuelipele.begin(), valuelipele.end());
    fillerEleTipValue.fill();
    fillerEleLipValue.fill();

    fillerEleTipError.insert(eleCandidates, errortipele.begin(), errortipele.end());
    fillerEleLipError.insert(eleCandidates, errorlipele.begin(), errorlipele.end());
    fillerEleTipError.fill();
    fillerEleLipError.fill();


 }
  
  //
  //////////////////////////////////////////////////////////
  //  
  //
  if (decaychannel=="2e2mu" || decaychannel=="4e"){ 
    iEvent.put( twodTipEleMap, "TipEleMap");
    iEvent.put( twodLipEleMap, "LipEleMap");

    iEvent.put( twodTipValueEleMap, "TipValueEleMap");
    iEvent.put( twodLipValueEleMap, "LipValueEleMap");
    
    iEvent.put( twodTipErrorEleMap, "TipErrorEleMap");
    iEvent.put( twodLipErrorEleMap, "LipErrorEleMap");

  }
  //
  if (decaychannel=="2e2mu" || decaychannel=="4mu"){
    iEvent.put( twodTipMuMap, "TipMuMap");
    iEvent.put( twodLipMuMap, "LipMuMap");

    iEvent.put( twodTipValueMuMap, "TipValueMuMap");
    iEvent.put( twodLipValueMuMap, "LipValueMuMap");
    
    iEvent.put( twodTipErrorMuMap, "TipErrorMuMap");
    iEvent.put( twodLipErrorMuMap, "LipErrorMuMap");
  }
  //
  delete theGeomPropBw;
  delete theMtsTransform; 
  //
  delete thePropBw;
  delete theTransform;
  //
}


      
