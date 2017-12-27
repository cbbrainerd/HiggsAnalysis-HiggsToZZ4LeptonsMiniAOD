
/* Original author:  Nicola De Filippis - Politecnico and INFN - Bari
 *                   
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsCP.h"


// system include files
#include <memory>

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
// Distances
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
// Boosting particles handling
#include "CommonTools/CandUtils/interface/CenterOfMassBooster.h"
#include "CommonTools/CandUtils/interface/Booster.h"
// close decay tree
#include "CommonTools/CandUtils/interface/cloneDecayTree.h"



// namespaces
using namespace edm;
using namespace std;
using namespace reco;


double HZZ4LeptonsCP::Distance( const reco::Candidate & c1, const reco::Candidate & c2 ) {
	return  deltaR(c1,c2);
}

double HZZ4LeptonsCP::DistancePhi( const reco::Candidate & c1, const reco::Candidate & c2 ) {
	return  deltaPhi(c1.p4().phi(),c2.p4().phi());
}

// Constructor
HZZ4LeptonsCP::HZZ4LeptonsCP(const edm::ParameterSet& pset) {
	
	// Get the various input parameters
	typedef std::vector<edm::InputTag> vtag;
	RECOcollName  = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("RECOcollName"));
	decayChain_   = pset.getParameter<std::string>("decayChain");
	debug	      =	pset.getUntrackedParameter<bool> ("debug");

	string alias;
	produces<vector<double> >(alias = decayChain_+ "cosTheta1").setBranchAlias( alias );
	produces<vector<double> >(alias = decayChain_+ "cosTheta2").setBranchAlias( alias );
	produces<vector<double> >(alias = decayChain_+ "cosThetaStar").setBranchAlias( alias );
	produces<vector<double> >(alias = decayChain_+ "Phi").setBranchAlias( alias );
	produces<vector<double> >(alias = decayChain_+ "Phi1").setBranchAlias( alias );
	produces<vector<double> >(alias = decayChain_+ "Phi2").setBranchAlias( alias );
	produces<vector<double> >(alias = decayChain_+ "phi1RF").setBranchAlias( alias );
	produces<vector<double> >(alias = decayChain_+ "phi2RF").setBranchAlias( alias );
	produces<vector<double> >(alias = decayChain_+ "MELA").setBranchAlias( alias );
     	
}

// Destructor
HZZ4LeptonsCP::~HZZ4LeptonsCP() {
	
}


void HZZ4LeptonsCP::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	
	auto_ptr<vector<double> > cosTheta1 ( new vector<double> );
	auto_ptr<vector<double> > cosTheta2 ( new vector<double> );
	auto_ptr<vector<double> > cosThetaStar ( new vector<double> );
	auto_ptr<vector<double> > Phi ( new vector<double> );
	auto_ptr<vector<double> > Phi1 ( new vector<double> );
	auto_ptr<vector<double> > Phi2 ( new vector<double> );
	auto_ptr<vector<double> > phi1RF ( new vector<double> );
	auto_ptr<vector<double> > phi2RF ( new vector<double> );
	auto_ptr<vector<double> > MELA ( new vector<double> );
	
	cosTheta1->clear();
	cosTheta2->clear();
	cosThetaStar->clear();
	Phi->clear();
	Phi1->clear();
	Phi2->clear();
	phi1RF->clear();
	phi2RF->clear();
	MELA->clear();
	
	
	// ZZ candidates
	Handle<edm::View<Candidate> > CandidateH;
	iEvent.getByToken(RECOcollName,CandidateH);
	
	//if(debug) cout << "Candidate H " << CandidateH->size() << " of type=" << RECOcollName.label().c_str() << endl;
	
	// extract the 4 vectors
	// higgs
	TLorentzVector hP4,Z1P4,L11P4,L12P4,Z2P4,L21P4,L22P4;
	bool H2Z=false,Z12l=false,Z22l=false;

	double angle_costheta1=-999., angle_costheta2=-999., angle_phi=-999., angle_costhetastar=-999., 
	    angle_phistar1=-999., angle_phistar2=-999., angle_phistar12=-999., 
	    angle_phi1=-999., angle_phi2=-999.,D=-999.;

	for( edm::View<Candidate>::const_iterator candH = CandidateH->begin();candH != CandidateH->end(); ++ candH ){
	  if(debug) cout << "Candidate H with mass= " << candH->mass() << endl;
	  const reco::Candidate& theParticle1 = (*candH);
	  if (candH->numberOfDaughters()>=2){
	    H2Z=true;
	    hP4.SetPxPyPzE( theParticle1.px(), theParticle1.py(), theParticle1.pz(), theParticle1.energy() );
	    
	    // z1
	    Z1P4.SetPxPyPzE( candH->daughter(0)->px(), candH->daughter(0)->py(), candH->daughter(0)->pz(), candH->daughter(0)->energy() );
	    
	    if (candH->daughter(0)->numberOfDaughters()>=2){
	      Z12l=true;
	      if (candH->daughter(0)->daughter(0)->charge() < 0){
		L11P4.SetPxPyPzE( candH->daughter(0)->daughter(0)->px(), candH->daughter(0)->daughter(0)->py(), 
				  candH->daughter(0)->daughter(0)->pz(), candH->daughter(0)->daughter(0)->energy() );
		L12P4.SetPxPyPzE( candH->daughter(0)->daughter(1)->px(), candH->daughter(0)->daughter(1)->py(), 
				  candH->daughter(0)->daughter(1)->pz(), candH->daughter(0)->daughter(1)->energy() );
	      }
	      else{
		L12P4.SetPxPyPzE( candH->daughter(0)->daughter(0)->px(), candH->daughter(0)->daughter(0)->py(), 
				  candH->daughter(0)->daughter(0)->pz(), candH->daughter(0)->daughter(0)->energy() );
		L11P4.SetPxPyPzE( candH->daughter(0)->daughter(1)->px(), candH->daughter(0)->daughter(1)->py(), 
				  candH->daughter(0)->daughter(1)->pz(), candH->daughter(0)->daughter(1)->energy() );
	      }
	    }
	    else {
	      cout << "There are no 2 daughters of the Z1" << endl;
	      continue;
	    }
	    
	    // z2
	    Z2P4.SetPxPyPzE( candH->daughter(1)->px(), candH->daughter(1)->py(), candH->daughter(1)->pz(), candH->daughter(1)->energy() );
	    
	    if (candH->daughter(1)->numberOfDaughters()>=2){
	      Z22l=true;
	      if (candH->daughter(1)->daughter(0)->charge() < 0){
		L21P4.SetPxPyPzE( candH->daughter(1)->daughter(0)->px(), candH->daughter(1)->daughter(0)->py(), 
				  candH->daughter(1)->daughter(0)->pz(), candH->daughter(1)->daughter(0)->energy() );
		L22P4.SetPxPyPzE( candH->daughter(1)->daughter(1)->px(), candH->daughter(1)->daughter(1)->py(), 
				  candH->daughter(1)->daughter(1)->pz(), candH->daughter(1)->daughter(1)->energy() );
	      }
	      else{
		L22P4.SetPxPyPzE( candH->daughter(1)->daughter(0)->px(), candH->daughter(1)->daughter(0)->py(), 
				  candH->daughter(1)->daughter(0)->pz(), candH->daughter(1)->daughter(0)->energy() );		      		      
		L21P4.SetPxPyPzE( candH->daughter(1)->daughter(1)->px(), candH->daughter(1)->daughter(1)->py(), 
				  candH->daughter(1)->daughter(1)->pz(), candH->daughter(1)->daughter(1)->energy() );
	      }
	    }
	    else {
	      cout << "There are no 2 daughters of the Z2" << endl;
	      continue;
	    }
	  }
	  else {
	    cout << "There are no 2 daughters of the higgs" << endl;
	    continue;
	  }
			
	  // Initialize
	  angle_costheta1=-999., angle_costheta2=-999., angle_phi=-999., angle_costhetastar=-999., 
	  angle_phistar1=-999., angle_phistar2=-999., angle_phistar12=-999., 
	  angle_phi1=-999., angle_phi2=-999.,D=-999.;
	  
	  if (H2Z && Z12l && Z22l) {
	    calculateAngles( hP4, Z1P4, L11P4, L12P4, Z2P4, L21P4, L22P4, angle_costheta1, angle_costheta2, angle_phi, angle_costhetastar, angle_phistar1, angle_phistar2, angle_phistar12, angle_phi1, angle_phi2);
	    
	    	    
	    double z1mass=Z1P4.M(), z2mass=Z2P4.M();
	    //mzz=hP4.M();
	    
	    if (debug){
	      std::cout << "Phi: " << angle_phi << std::endl;
	      std::cout << "cos theta_1: " << angle_costheta1 << std::endl;
	      std::cout << "cos theta_2: " << angle_costheta2 << std::endl;
	      std::cout << "Phi1: " << angle_phistar1 << std::endl;
	      std::cout << "cos theta*: " << angle_costhetastar << std::endl;
	      std::cout << "mz1: " << z1mass << std::endl;
	      std::cout << "mz2: " << z2mass << std::endl;
	    }


	    // checkZorder<double>(z1mass,z2mass,angle_costhetastar,angle_costheta1,angle_costheta2,angle_phi,angle_phi1);
	    
	    
// 	    if(mzz>100. && mzz<800. && z2mass>4. ) {	    
// 	      //MELA LD
// 	      pair<double,double> P = likelihoodDiscriminant(mzz,z1mass,z2mass,angle_costhetastar,angle_costheta1,angle_costheta2,angle_phi,angle_phi1);
// 	      D=P.first/(P.first+P.second);	    
// 	      if (debug) cout << "MELA: " << D << endl;	    
// 	    }
	  }
	  		
	  cosTheta1->push_back(angle_costheta1);
	  cosTheta2->push_back(angle_costheta2);
	  cosThetaStar->push_back(angle_costhetastar);
	  Phi->push_back(angle_phi);
	  Phi1->push_back(angle_phistar1);
	  Phi2->push_back(angle_phistar2);
	  phi1RF->push_back(angle_phi1);
	  phi2RF->push_back(angle_phi2);
	  MELA->push_back(D);
	  
	}

	iEvent.put(cosTheta1, decayChain_ + "cosTheta1");
	iEvent.put(cosTheta2, decayChain_ + "cosTheta2");
	iEvent.put(cosThetaStar, decayChain_ + "cosThetaStar");
	iEvent.put(Phi, decayChain_ + "Phi");
	iEvent.put(Phi1, decayChain_ + "Phi1");
	iEvent.put(Phi2, decayChain_ + "Phi2");
	iEvent.put(phi1RF, decayChain_ + "phi1RF");
	iEvent.put(phi2RF, decayChain_ + "phi2RF");
	iEvent.put(MELA, decayChain_ + "MELA");
	
	
}

void HZZ4LeptonsCP::beginJob() {
  
  if(debug) cout << "in begin job of CP module!" << std::endl;
	
}

void HZZ4LeptonsCP::endJob() {
  
}


//=======================================================================

void HZZ4LeptonsCP::calculateAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4Lep11, TLorentzVector thep4Lep12, TLorentzVector thep4Z2, TLorentzVector thep4Lep21, TLorentzVector thep4Lep22, double& costheta1, double& costheta2, double& phi, double& costhetastar, double& phistar1, double& phistar2, double& phistar12, double& phi1, double& phi2){
	
	
	//std::cout << "In calculate angles..." << std::endl;
	
	double norm;
	
	TVector3 boostX = -(thep4H.BoostVector());
	TLorentzVector thep4Z1inXFrame( thep4Z1 );
	TLorentzVector thep4Z2inXFrame( thep4Z2 );	
	thep4Z1inXFrame.Boost( boostX );
	thep4Z2inXFrame.Boost( boostX );
	TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
	TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );
	
	// calculate phi1, phi2, costhetastar
	phi1 = theZ1X_p3.Phi();
	phi2 = theZ2X_p3.Phi();
	
	///////////////////////////////////////////////
	// check for z1/z2 convention, redefine all 4 vectors with convention
	///////////////////////////////////////////////	
	TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;

	/* old convention of choosing Z1 ------------------------------
	p4H = thep4H;
	if ((phi1 < 0)&&(phi1 >= -TMath::Pi())){
		p4Z1 = thep4Z2; p4M11 = thep4Lep21; p4M12 = thep4Lep22;
		p4Z2 = thep4Z1; p4M21 = thep4Lep11; p4M22 = thep4Lep12;		
		costhetastar = theZ2X_p3.CosTheta();
	}
	else{
		p4Z1 = thep4Z1; p4M11 = thep4Lep11; p4M12 = thep4Lep12;
		p4Z2 = thep4Z2; p4M21 = thep4Lep21; p4M22 = thep4Lep22;
		costhetastar = theZ1X_p3.CosTheta();
	} ---------------------------------------------- */

	p4Z1 = thep4Z1; p4M11 = thep4Lep11; p4M12 = thep4Lep12;
	p4Z2 = thep4Z2; p4M21 = thep4Lep21; p4M22 = thep4Lep22;
	costhetastar = theZ1X_p3.CosTheta();
	
	//std::cout << "phi1: " << phi1 << ", phi2: " << phi2 << std::endl;
	
	// now helicity angles................................
	// ...................................................
	TVector3 boostZ1 = -(p4Z1.BoostVector());
	TLorentzVector p4Z2Z1(p4Z2);
	p4Z2Z1.Boost(boostZ1);
	//find the decay axis
	/////TVector3 unitx_1 = -Hep3Vector(p4Z2Z1);
	TVector3 unitx_1( -p4Z2Z1.X(), -p4Z2Z1.Y(), -p4Z2Z1.Z() );
	norm = 1/(unitx_1.Mag());
	unitx_1*=norm;
	//boost daughters of z2
	TLorentzVector p4M21Z1(p4M21);
	TLorentzVector p4M22Z1(p4M22);
	p4M21Z1.Boost(boostZ1);
	p4M22Z1.Boost(boostZ1);
	//create z and y axes
	/////TVector3 unitz_1 = Hep3Vector(p4M21Z1).cross(Hep3Vector(p4M22Z1));
	TVector3 p4M21Z1_p3( p4M21Z1.X(), p4M21Z1.Y(), p4M21Z1.Z() );
	TVector3 p4M22Z1_p3( p4M22Z1.X(), p4M22Z1.Y(), p4M22Z1.Z() );
	TVector3 unitz_1 = p4M21Z1_p3.Cross( p4M22Z1_p3 );
	norm = 1/(unitz_1.Mag());
	unitz_1 *= norm;
	TVector3 unity_1 = unitz_1.Cross(unitx_1);
	
	//caculate theta1
	TLorentzVector p4M11Z1(p4M11);
	p4M11Z1.Boost(boostZ1);
	TVector3 p3M11( p4M11Z1.X(), p4M11Z1.Y(), p4M11Z1.Z() );
	TVector3 unitM11 = p3M11.Unit();
	double x_m11 = unitM11.Dot(unitx_1); double y_m11 = unitM11.Dot(unity_1); double z_m11 = unitM11.Dot(unitz_1);
	TVector3 M11_Z1frame(y_m11, z_m11, x_m11);
	costheta1 = M11_Z1frame.CosTheta();
	//std::cout << "theta1: " << M11_Z1frame.Theta() << std::endl;
	//////-----------------------old way of calculating phi---------------/////////
	phi = M11_Z1frame.Phi();
	
	//set axes for other system
	TVector3 boostZ2 = -(p4Z2.BoostVector());
	TLorentzVector p4Z1Z2(p4Z1);
	p4Z1Z2.Boost(boostZ2);
	TVector3 unitx_2( -p4Z1Z2.X(), -p4Z1Z2.Y(), -p4Z1Z2.Z() );
	norm = 1/(unitx_2.Mag());
	unitx_2*=norm;
	//boost daughters of z2
	TLorentzVector p4M11Z2(p4M11);
	TLorentzVector p4M12Z2(p4M12);
	p4M11Z2.Boost(boostZ2);
	p4M12Z2.Boost(boostZ2);
	TVector3 p4M11Z2_p3( p4M11Z2.X(), p4M11Z2.Y(), p4M11Z2.Z() );
	TVector3 p4M12Z2_p3( p4M12Z2.X(), p4M12Z2.Y(), p4M12Z2.Z() );
	TVector3 unitz_2 = p4M11Z2_p3.Cross( p4M12Z2_p3 );
	norm = 1/(unitz_2.Mag());
	unitz_2*=norm;
	TVector3 unity_2 = unitz_2.Cross(unitx_2);
	//calcuate theta2
	TLorentzVector p4M21Z2(p4M21);
	p4M21Z2.Boost(boostZ2);
	TVector3 p3M21( p4M21Z2.X(), p4M21Z2.Y(), p4M21Z2.Z() );
	TVector3 unitM21 = p3M21.Unit();
	double x_m21 = unitM21.Dot(unitx_2); double y_m21 = unitM21.Dot(unity_2); double z_m21 = unitM21.Dot(unitz_2);
	TVector3 M21_Z2frame(y_m21, z_m21, x_m21);
	costheta2 = M21_Z2frame.CosTheta();
	
	// calculate phi
	//calculating phi_n
	TLorentzVector n_p4Z1inXFrame( p4Z1 );
	TLorentzVector n_p4M11inXFrame( p4M11 );
	n_p4Z1inXFrame.Boost( boostX );
	n_p4M11inXFrame.Boost( boostX );        
	TVector3 n_p4Z1inXFrame_unit = n_p4Z1inXFrame.Vect().Unit();
	TVector3 n_p4M11inXFrame_unit = n_p4M11inXFrame.Vect().Unit();  
	TVector3 n_unitz_1( n_p4Z1inXFrame_unit );
	//// y-axis is defined by neg lepton cross z-axis
	//// the subtle part is here...
	//////////TVector3 n_unity_1 = n_p4M11inXFrame_unit.Cross( n_unitz_1 );
	TVector3 n_unity_1 = n_unitz_1.Cross( n_p4M11inXFrame_unit );
	TVector3 n_unitx_1 = n_unity_1.Cross( n_unitz_1 );
	
	TLorentzVector n_p4M21inXFrame( p4M21 );
	n_p4M21inXFrame.Boost( boostX );
	TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();
	//rotate into other plane
	TVector3 n_p4M21inXFrame_unitprime( n_p4M21inXFrame_unit.Dot(n_unitx_1), n_p4M21inXFrame_unit.Dot(n_unity_1), n_p4M21inXFrame_unit.Dot(n_unitz_1) );
	
	///////-----------------new way of calculating phi-----------------///////
	//double phi_n =  n_p4M21inXFrame_unitprime.Phi();
	/// and then calculate phistar1
	TVector3 n_p4PartoninXFrame_unit( 0.0, 0.0, 1.0 );
	TVector3 n_p4PartoninXFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_1), n_p4PartoninXFrame_unit.Dot(n_unity_1), n_p4PartoninXFrame_unit.Dot(n_unitz_1) );
	// negative sign is for arrow convention in paper
	phistar1 = (n_p4PartoninXFrame_unitprime.Phi());
	
	// and the calculate phistar2
	TLorentzVector n_p4Z2inXFrame( p4Z2 );
	n_p4Z2inXFrame.Boost( boostX );
	TVector3 n_p4Z2inXFrame_unit = n_p4Z2inXFrame.Vect().Unit();
	///////TLorentzVector n_p4M21inXFrame( p4M21 );
	//////n_p4M21inXFrame.Boost( boostX );        
	////TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();  
	TVector3 n_unitz_2( n_p4Z2inXFrame_unit );
	//// y-axis is defined by neg lepton cross z-axis
	//// the subtle part is here...
	//////TVector3 n_unity_2 = n_p4M21inXFrame_unit.Cross( n_unitz_2 );
	TVector3 n_unity_2 = n_unitz_2.Cross( n_p4M21inXFrame_unit );
	TVector3 n_unitx_2 = n_unity_2.Cross( n_unitz_2 );
	TVector3 n_p4PartoninZ2PlaneFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_2), n_p4PartoninXFrame_unit.Dot(n_unity_2), n_p4PartoninXFrame_unit.Dot(n_unitz_2) );
	phistar2 = (n_p4PartoninZ2PlaneFrame_unitprime.Phi());
	
	double phistar12_0 = phistar1 + phistar2;
	if (phistar12_0 > TMath::Pi()) phistar12 = phistar12_0 - 2*TMath::Pi();
	else if (phistar12_0 < (-1.)*TMath::Pi()) phistar12 = phistar12_0 + 2*TMath::Pi();
	else phistar12 = phistar12_0;
	
}


 
