#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include <DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include <DataFormats/GeometryVector/interface/Point3DBase.h>
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"


#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h>
#include <RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h>
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "TMath.h"

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsConstraintFitProducer.h"


//
// constructors and destructor
//
HZZ4LeptonsConstraintFitProducer::HZZ4LeptonsConstraintFitProducer(const edm::ParameterSet& pset)
{

  vertexTag_     = pset.getParameter<edm::InputTag>("VertexLabel");
  RECOcollName   = pset.getParameter<edm::InputTag>("RECOcollName");
  nParticles     = pset.getParameter<uint>("nParticles");

  debug	         = pset.getUntrackedParameter<bool> ("debug", false);

  //now do what ever initialization is needed
  massSqr_ = 0.105658*0.1056583;

  std::string alias;

  produces<reco::VertexCollection>( alias = "KinematicFitVertex");
  produces<reco::VertexCollection>( alias = "StandardFitVertex");
  produces<edm::ValueMap<float> > ( alias = "RefittedMass");  
}


HZZ4LeptonsConstraintFitProducer::~HZZ4LeptonsConstraintFitProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HZZ4LeptonsConstraintFitProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   //ESHandle<SetupData> pSetup;
  //iSetup.get<SetupRecord>().get(pSetup);
  
  using namespace edm;
  using namespace std;
  using namespace reco;

  std::auto_ptr<reco::VertexCollection> KinFitVtx( new reco::VertexCollection );
  std::auto_ptr<reco::VertexCollection> StdFitVtx( new reco::VertexCollection );
  
  auto_ptr<edm::ValueMap<float> >       RefittedMassMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler          fillerMass(*RefittedMassMap);

  Handle<vector<Vertex> >  vertexs;
  iEvent.getByLabel(vertexTag_,vertexs);

  Vertex::Point pvPosition(vertexs->front().position());

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

  KalmanVertexFitter vtxFitter(true);

  //Creating a KinematicParticleFactory
  KinematicParticleFactoryFromTransientTrack factory;
 
  ParticleMass muon_mass = 0.1056583;
  ParticleMass electron_mass = 0.0005;
  //ParticleMass pion_mass = 0.139570;
  float muon_sigma = 0.0000000001;
  float electron_sigma = 0.0000000001;
  //float pion_sigma = 0.0000000001;

  //initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered 
  float chi;
  float ndof;
  
  vector<RefCountedKinematicParticle> Particles;
  vector<TransientTrack> t_tks;
  
  std::vector<reco::TransientTrack> refit_tks;

  // RECO 4lepton candidate
  edm::Handle<edm::View<Candidate> > Candidates;
  iEvent.getByLabel(RECOcollName, Candidates);

  int I_cand=0;
  size_t nCand = Candidates->size();
  std::vector<float> refittedmass(nCand);
  

  // Prepare Dummy Vertices for insertion as place holders for candidate array matching
  reco::Vertex::Error err;
  reco::Vertex::Point dummyPoint(-999,-999,-999);
  reco::Vertex dummyVertex_std(dummyPoint,err, -999*5,5,4);
      
  
  if(debug) cout << "Accessing " << nParticles << " lepton collection=" << RECOcollName.label() << " with size=" << Candidates->size() << endl;
  for( edm::View<Candidate>::const_iterator cand = Candidates->begin();cand != Candidates->end(); ++ cand ) { 

// first created integer to keep track of vertex push_backs
   int StdFit_push_back = 0;
   int KinFit_push_back = 0;
    
    cout << "Mass of candidate " << RECOcollName.label() << " is " << cand->mass() << endl;
    if (cand->mass() >12.){
      Particles.clear();
      t_tks.clear();
      
      //      bool isGoodMass_original=true;                                                          
      // bool isGoodMass_refitted=true;
      
      for (unsigned k = 0; k < cand->numberOfDaughters(); ++k ) {
	
	if (cand->daughter(k)->numberOfDaughters() >0){
	  for (unsigned l = 0; l < cand->daughter(k)->numberOfDaughters(); ++l ) {
	    chi = 0.; ndof = 0.;
	    if (cand->daughter(k)->daughter(l)->isGlobalMuon() || cand->daughter(k)->daughter(l)->isTrackerMuon()){
	      TransientTrack tt = theTTBuilder->build(cand->daughter(k)->daughter(l)->get<TrackRef>());
	      if(debug) cout << "Track momentum=" << tt.track().pt() << " chi2= " << tt.track().normalizedChi2() << endl;
	      Particles.push_back(factory.particle (tt,muon_mass,chi,ndof,muon_sigma));
	      t_tks.push_back(tt);
	    }
	    else if (cand->daughter(k)->daughter(l)->isElectron()){
	      TransientTrack tt = theTTBuilder->build(cand->daughter(k)->daughter(l)->get<GsfTrackRef>());
	      if(debug) cout << "Track momentum=" << tt.track().pt() << " chi2= " << tt.track().normalizedChi2() << endl;
	      Particles.push_back(factory.particle (tt,electron_mass,chi,ndof,electron_sigma));
	      t_tks.push_back(tt);
	    } 
	  }     	
	}
	else{
	  chi = 0.; ndof = 0.;
	  if (cand->daughter(k)->isGlobalMuon() || cand->daughter(k)->isTrackerMuon()){
	    TransientTrack tt = theTTBuilder->build(cand->daughter(k)->get<TrackRef>());
	    if(debug) cout << "Track momentum=" << tt.track().pt() << " chi2= " << tt.track().normalizedChi2() << endl;
	    Particles.push_back(factory.particle (tt,muon_mass,chi,ndof,muon_sigma));
	    t_tks.push_back(tt);
	  }
	  else if (cand->daughter(k)->isElectron()){
	    TransientTrack tt = theTTBuilder->build(cand->daughter(k)->get<GsfTrackRef>());
	    if(debug) cout << "Track momentum=" << tt.track().pt() << " chi2= " << tt.track().normalizedChi2() << endl;
	    Particles.push_back(factory.particle (tt,electron_mass,chi,ndof,electron_sigma));
	    t_tks.push_back(tt);
	  } 
	}
      }
      
      
      if(debug) cout << "Number of particle for constrain fitter= " << Particles.size()<< endl;
      if (Particles.size()>=nParticles){
	KinematicParticleVertexFitter fitter; 
	RefCountedKinematicTree myTree = fitter.fit(Particles); 
	
	if ( !myTree->isEmpty()) {
	  
	  //accessing the tree components
	  myTree->movePointerToTheTop();
	  
	  
	  RefCountedKinematicParticle allMuonsCand     = myTree->currentParticle();
	  RefCountedKinematicVertex allMuonsVertex     = myTree->currentDecayVertex();
	  
	  if(debug) cout << "m(" << cand->numberOfDaughters() << "muons): " << allMuonsCand->currentState().mass() << " +- "
			 << sqrt(allMuonsCand->currentState().kinematicParametersError().matrix()(6,6) ) << endl;
	  
	  
	  if (allMuonsVertex->vertexIsValid()) {	  
	    if(debug) cout << "kinematicFit vertex, ndof, chi2, prob: " 
			   << allMuonsVertex->position() << " , " 
			   << allMuonsVertex->degreesOfFreedom() << " , "
			   << allMuonsVertex->chiSquared()   << " , "
			   << TMath::Prob(allMuonsVertex->chiSquared(),allMuonsVertex->degreesOfFreedom()) << endl;
	    
// push back dummy vertex if not already stored:

             if(KinFit_push_back ==0){
	         KinFitVtx->push_back(*allMuonsVertex);
                 KinFit_push_back = 1;
             }
	  }

	  
	  //This is the standard vertex fit without constraint
	  TransientVertex myVertex = vtxFitter.vertex(t_tks);
	  reco::Vertex myRecoVertex = myVertex;
	  
	  if (myRecoVertex.isValid()) {
	    if(debug) cout << "standardFit vertex, ndof, chi2, prob: " 
			   << myRecoVertex.position() << " , " 
			   << myRecoVertex.ndof() << " , "
			   << myRecoVertex.chi2() << " , "
			   << TMath::Prob(myRecoVertex.chi2(),myRecoVertex.ndof()) <<  endl;
            // push back good vertex if place is free:
            if(StdFit_push_back ==0){
               StdFitVtx ->push_back(myRecoVertex);
               StdFit_push_back = 1;
            }
	    
	    
	    //////MM
	    if(debug) cout<<endl<<"Candidate "<<I_cand+1<<endl;
	    if(debug) cout<<endl<<"ORIGINAL TRACKS:"<<endl;
	    
	    cout<<"Original Track size = "<<t_tks.size()<<endl;
	    
	    double px_sum_original=0; 
	    double py_sum_original=0; 
	    double pz_sum_original=0; 
	    double e_sum_original=0;
	    
	    
	    for(unsigned int i=0; i< t_tks.size(); i++){
	      
	      
	      reco::Track myOriginalTrack= t_tks[i].track();
	      
	      px_sum_original+=myOriginalTrack.px();
	      py_sum_original+=myOriginalTrack.py();
	      pz_sum_original+=myOriginalTrack.pz();
	      
	      
	      double p=myOriginalTrack.p();
	      
	      double mass=Particles[i]->currentState().mass();
	      
	      if(debug) cout<<"\tmass["<<i+1<<"] ="<< mass;
	      if(debug) cout<<"\tp["<<i+1<<"] ="<< p<<endl;
	      
	      double energy= sqrt(mass*mass+p*p);
	      
	      e_sum_original+=energy;
	      
	    }
	    
	    
	    if(debug) cout<<endl<<"px="<<px_sum_original<<"\tpy="<<py_sum_original<<"\tpz="<<pz_sum_original<<"\t E="<<e_sum_original<<endl;
	    
	    double sq_M4l_vtx_original= 
	      e_sum_original*e_sum_original - 
	      (px_sum_original*px_sum_original + py_sum_original*py_sum_original + pz_sum_original*pz_sum_original);
	    
	    double M4l_vtx_original= -999.;
	    if(sq_M4l_vtx_original >= 0.) M4l_vtx_original= sqrt(sq_M4l_vtx_original);
	    
	    cout<<"M4l="<<M4l_vtx_original<<endl<<endl;
	    
	    // Refitted tracks
	    
	    if(debug) cout<<endl<<"REFITTED TRACKS:"<<endl;
	    
	    
	    double px_sum=0; 
	    double py_sum=0; 
	    double pz_sum=0; 
	    double e_sum=0;
	    
	    bool hasRefittedTracks = myVertex.hasRefittedTracks();
	    if(hasRefittedTracks){
	      
	      refit_tks= myVertex.refittedTracks(); 
	      
	      
	      if(debug) cout<<"RefittedTrack size = "<<refit_tks.size()<<endl;
	      
	      for(unsigned int i=0; i< refit_tks.size(); i++){
	      	
		reco::Track myRefittedTrack= refit_tks[i].track();
		
		cout << "Track momentum Refit=" << myRefittedTrack.pt() << endl;
		
		px_sum+=myRefittedTrack.px();
		py_sum+=myRefittedTrack.py();
		pz_sum+=myRefittedTrack.pz();
		
		double p=myRefittedTrack.p();
		
		double mass=Particles[i/*-fromI_tk*/]->currentState().mass();
		
		if(debug) cout<<"\tmass["<<i+1<<"] ="<< mass;
		if(debug) cout<<"\tp["<<i+1<<"] ="<< p<<endl;
		
		double energy= sqrt(mass*mass+p*p);
		
		e_sum+=energy;
		
	      }	    
	      
	    }
	    
	    if(debug) cout<<endl<<"px="<<px_sum<<"\tpy="<<py_sum<<"\tpz="<<pz_sum<<"\t E="<<e_sum<<endl;
	    
	    double sq_M4l_vtx= e_sum*e_sum - (px_sum*px_sum + py_sum*py_sum + pz_sum*pz_sum);
	    double M4l_vtx= -999.;
	    if(sq_M4l_vtx >= 0.) M4l_vtx = sqrt(sq_M4l_vtx);
	    
	    
	    if(debug) cout<<"M4l="<<M4l_vtx<<endl;
	    
	    cout<<"Original vs Refitted M4l: "<<setprecision(6)<<M4l_vtx_original<<"\t"<<setprecision(6)<<M4l_vtx<<endl<<endl;
	    refittedmass[I_cand]=float(M4l_vtx);	  	 	  	  	  	 	  
	    
	  }
	}
      }
    }

//  Increment Local Counter to keep track of candidate:
    I_cand++;

// FR:
// push back dummy vertex if not already stored: in StdFitVtx array ....
    if(StdFit_push_back ==0){
      StdFitVtx ->push_back(dummyVertex_std);
      StdFit_push_back = 1;
    }

// .... and in KinFit array 
    if(KinFit_push_back ==0){
      KinFitVtx ->push_back(dummyVertex_std);
      KinFit_push_back = 1;
    }


  } //end loop over candidates
  
  // filling map
  fillerMass.insert(Candidates, refittedmass.begin(), refittedmass.end());
  fillerMass.fill();
  
  // iEvent.put( KinFitVtx, "KinematicFitVertex" );
  // iEvent.put( StdFitVtx, "StandardFitVertex" );
  // iEvent.put( RefittedMassMap, "RefittedMass");

  iEvent.put(std::make_unique<reco::VertexCollection>(*KinFitVtx), "KinematicFitVertex" );
  iEvent.put(std::make_unique<reco::VertexCollection>(*StdFitVtx), "StandardFitVertex" );
  iEvent.put(std::make_unique<edm::ValueMap<float> >(*RefittedMassMap), "RefittedMass");
  
}

