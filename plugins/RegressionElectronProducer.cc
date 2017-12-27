/**\class RegressionElectronProducer.cc
 *
 * Original Author:  Nicola De Filippis
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/RegressionElectronProducer.h"

#include "FWCore/Framework/interface/ESHandle.h"

// Electrons:
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectronFwd.h>

#include "DataFormats/Common/interface/ValueMap.h"


using namespace edm;
using namespace std;
using namespace reco;


// constructor
RegressionElectronProducer::RegressionElectronProducer(const edm::ParameterSet& pset) {

  eleTag_                     = pset.getParameter<edm::InputTag>("eleCollection");
  eleRegressionEnergyErrorTag_= pset.getParameter<edm::InputTag>("eleRegressionEnergyErrorLabel");
  eleRegressionEnergyTag_     = pset.getParameter<edm::InputTag>("eleRegressionEnergyLabel");
        

  produces<reco::GsfElectronCollection>(); 

  counterelectron=0;	


}


RegressionElectronProducer::~RegressionElectronProducer() {

}

void RegressionElectronProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<reco::GsfElectronCollection> Gelec( new reco::GsfElectronCollection );

  // Get all pixel match GSF electron candidates
  edm::Handle<reco::GsfElectronCollection> electrons;
  iEvent.getByLabel(eleTag_, electrons);

   // electron regression
  edm::Handle<edm::ValueMap<double> > eleRegressionEnergymap;
  iEvent.getByLabel(eleRegressionEnergyTag_, eleRegressionEnergymap);

  edm::Handle<edm::ValueMap<double> > eleRegressionEnergyErrormap;
  iEvent.getByLabel(eleRegressionEnergyErrorTag_, eleRegressionEnergyErrormap);


  counterelectron=0;


  for (reco::GsfElectronCollection::const_iterator electron = electrons->begin(); electron != electrons->end(); ++electron) {
    
    edm::Ref<reco::GsfElectronCollection> eletrackref(electrons,counterelectron);

    // Electron Regression
    double newEnergy_=(*eleRegressionEnergymap)[eletrackref];
    double newEnergyError_=(*eleRegressionEnergyErrormap)[eletrackref];  
    math::XYZTLorentzVector newMomentum_ ;
    float errorTrackMomentum_ ;
    float finalMomentumError_ ;
    
    
    //float scEnergy = electron.ecalEnergy() ;
    float scEnergy = newEnergy_ ;
    int elClass = electron->classification() ;
    
    float trackMomentum  = electron->trackMomentumAtVtx().R() ;
    errorTrackMomentum_ = 999. ;
    
    // retreive momentum error 
    //MultiGaussianState1D qpState(MultiGaussianStateTransform::multiState1D(vtxTsos,0));
    //GaussianSumUtilities1D qpUtils(qpState);
    errorTrackMomentum_ = electron->trackMomentumError();

    
    float finalMomentum = electron->p4().t(); // initial
    float finalMomentumError = 999.;
    
    // first check for large errors
    
    if (errorTrackMomentum_/trackMomentum > 0.5 && newEnergyError_/scEnergy <= 0.5) {
      finalMomentum = scEnergy;    finalMomentumError = newEnergyError_;
    }
    else if (errorTrackMomentum_/trackMomentum <= 0.5 && newEnergyError_/scEnergy > 0.5){
      finalMomentum = trackMomentum;  finalMomentumError = errorTrackMomentum_;
    }
    else if (errorTrackMomentum_/trackMomentum > 0.5 && newEnergyError_/scEnergy > 0.5){
      if (errorTrackMomentum_/trackMomentum < newEnergyError_/scEnergy) {
	finalMomentum = trackMomentum; finalMomentumError = errorTrackMomentum_;
      }
      else{
	finalMomentum = scEnergy; finalMomentumError = newEnergyError_;
      }
    }
    
    //then apply the combination algorithm
    else {
      
      // calculate E/p and corresponding error
      float eOverP = scEnergy / trackMomentum;
      float errorEOverP = sqrt(
			       (newEnergyError_/trackMomentum)*(newEnergyError_/trackMomentum) +
			       (scEnergy*errorTrackMomentum_/trackMomentum/trackMomentum)*
			       (scEnergy*errorTrackMomentum_/trackMomentum/trackMomentum));
      
      //new comb
      
      bool eleIsNotInCombination = false ;
      if ( (eOverP  > 1 + 2.5*errorEOverP) || (eOverP  < 1 - 2.5*errorEOverP) || (eOverP < 0.8) || (eOverP > 1.3) )
	{ eleIsNotInCombination = true ; }
      
      if (eleIsNotInCombination)
	{
	  if (eOverP > 1)
	    { finalMomentum = scEnergy ; finalMomentumError = newEnergyError_ ; }
	  else
	    {
	      if (elClass == reco::GsfElectron::GOLDEN)
		{ finalMomentum = scEnergy; finalMomentumError = newEnergyError_; }
	      if (elClass == reco::GsfElectron::BIGBREM)
		{
		  if (scEnergy<36)
		    { finalMomentum = trackMomentum ; finalMomentumError = errorTrackMomentum_ ; }
		  else
		    { finalMomentum = scEnergy ; finalMomentumError = newEnergyError_ ; }
		}
	      if (elClass == reco::GsfElectron::BADTRACK)
		{ finalMomentum = scEnergy; finalMomentumError = newEnergyError_ ; }
	      if (elClass == reco::GsfElectron::SHOWERING)
		{
		  if (scEnergy<30)
		    { finalMomentum = trackMomentum ; finalMomentumError = errorTrackMomentum_; }
		  else
		    { finalMomentum = scEnergy; finalMomentumError = newEnergyError_;}
		}
	      if (elClass == reco::GsfElectron::GAP)
		{
		  if (scEnergy<60)
		    { finalMomentum = trackMomentum ; finalMomentumError = errorTrackMomentum_ ; }
		  else
		    { finalMomentum = scEnergy; finalMomentumError = newEnergyError_ ; }
		}
	    }
	}
      
      else
	{
	  // combination
	  finalMomentum = (scEnergy/newEnergyError_/newEnergyError_ + trackMomentum/errorTrackMomentum_/errorTrackMomentum_) /
	    (1/newEnergyError_/newEnergyError_ + 1/errorTrackMomentum_/errorTrackMomentum_);
	  float finalMomentumVariance = 1 / (1/newEnergyError_/newEnergyError_ + 1/errorTrackMomentum_/errorTrackMomentum_);
	  finalMomentumError = sqrt(finalMomentumVariance);
	}
    }         
    
    // }
    
    math::XYZTLorentzVector oldMomentum = electron->p4() ;
    
    newMomentum_ = math::XYZTLorentzVector
      ( oldMomentum.x()*finalMomentum/oldMomentum.t(),
	oldMomentum.y()*finalMomentum/oldMomentum.t(),
	oldMomentum.z()*finalMomentum/oldMomentum.t(),
	finalMomentum ) ;
    finalMomentumError_ =  finalMomentumError;  
    std::cout << "[ElectronEnergCorrector] old comb momentum " << oldMomentum.t() << " new comb momentum " << newMomentum_.t() << std::endl;
    std::cout << "[ElectronEnergCorrector] old comb transv. momentum " << oldMomentum.pt() << " new comb transv. momentum " << newMomentum_.pt() << std::endl;
    std::cout << "[ElectronEnergCorrector] new energy "  << float(newEnergy_) << " " << float(newEnergyError_)<< std::endl;
    std::cout << "[ElectronEnergCorrector] new error tracks and final momentum error "  << errorTrackMomentum_ << " " << finalMomentumError_<< std::endl;
    
    
    reco::GsfElectron clone =*electron;

    if (clone.ecalDrivenSeed()) {     
      clone.correctEcalEnergy(float(newEnergy_), float(newEnergyError_));
      clone.correctMomentum(newMomentum_,errorTrackMomentum_, finalMomentumError_);
    }
    
    Gelec->push_back( clone);
    counterelectron++;
  }
  
  
  const string iName = "";
  iEvent.put( Gelec, iName );

}


void RegressionElectronProducer::computeEpCombination
( const reco::GsfElectron & electron, double newEnergy_, double newEnergyError_ ,math::XYZTLorentzVector newMomentum_,float errorTrackMomentum_,float finalMomentumError_){
  
  //math::XYZTLorentzVector newMomentum_ ;
  //float errorTrackMomentum_ ;
  //float finalMomentumError_ ;

  //float scEnergy = electron.ecalEnergy() ;
  float scEnergy = newEnergy_ ;
  int elClass = electron.classification() ;

  float trackMomentum  = electron.trackMomentumAtVtx().R() ;
  errorTrackMomentum_ = 999. ;
  
  // retreive momentum error 
  //MultiGaussianState1D qpState(MultiGaussianStateTransform::multiState1D(vtxTsos,0));
  //GaussianSumUtilities1D qpUtils(qpState);
  errorTrackMomentum_ = electron.trackMomentumError();

  float finalMomentum = electron.p4().t(); // initial
  float finalMomentumError = 999.;
  
  // first check for large errors
 
  if (errorTrackMomentum_/trackMomentum > 0.5 && newEnergyError_/scEnergy <= 0.5) {
    finalMomentum = scEnergy;    finalMomentumError = newEnergyError_;
  }
  else if (errorTrackMomentum_/trackMomentum <= 0.5 && newEnergyError_/scEnergy > 0.5){
    finalMomentum = trackMomentum;  finalMomentumError = errorTrackMomentum_;
  }
  else if (errorTrackMomentum_/trackMomentum > 0.5 && newEnergyError_/scEnergy > 0.5){
    if (errorTrackMomentum_/trackMomentum < newEnergyError_/scEnergy) {
      finalMomentum = trackMomentum; finalMomentumError = errorTrackMomentum_;
    }
    else{
      finalMomentum = scEnergy; finalMomentumError = newEnergyError_;
    }
  }
  
  // then apply the combination algorithm
  else {

    // calculate E/p and corresponding error
    float eOverP = scEnergy / trackMomentum;
    float errorEOverP = sqrt(
			     (newEnergyError_/trackMomentum)*(newEnergyError_/trackMomentum) +
			     (scEnergy*errorTrackMomentum_/trackMomentum/trackMomentum)*
			     (scEnergy*errorTrackMomentum_/trackMomentum/trackMomentum));
    //old comb  
    //     if ( eOverP  > 1 + 2.5*errorEOverP )
    //       {
    // finalMomentum = scEnergy; finalMomentumError = newEnergyError_;
    // if ((elClass==reco::GsfElectron::GOLDEN) && electron.isEB() && (eOverP<1.15))
    //   {
    //     if (scEnergy<15) {finalMomentum = trackMomentum ; finalMomentumError = errorTrackMomentum_;}
    //   }
    //       }
    //     else if ( eOverP < 1 - 2.5*errorEOverP )
    //       {
    // finalMomentum = scEnergy; finalMomentumError = newEnergyError_;
    // if (elClass==reco::GsfElectron::SHOWERING)
    //   {
    //     if (electron.isEB())
    //       {
    // if(scEnergy<18) {finalMomentum = trackMomentum; finalMomentumError = errorTrackMomentum_;}
    //       }
    //     else if (electron.isEE())
    //       {
    // if(scEnergy<13) {finalMomentum = trackMomentum; finalMomentumError = errorTrackMomentum_;}
    //       }
    //     else
    //       { edm::LogWarning("ElectronMomentumCorrector::correct")<<"nor barrel neither endcap electron ?!" ; }
    //   }
    // else if (electron.isGap())
    //   {
    //     if(scEnergy<60) {finalMomentum = trackMomentum; finalMomentumError = errorTrackMomentum_;}
    //   }
    //       }
    //     else 
    //       {
    // // combination
    // finalMomentum = (scEnergy/newEnergyError_/newEnergyError_ + trackMomentum/errorTrackMomentum_/errorTrackMomentum_) /
    //   (1/newEnergyError_/newEnergyError_ + 1/errorTrackMomentum_/errorTrackMomentum_);
    // float finalMomentumVariance = 1 / (1/newEnergyError_/newEnergyError_ + 1/errorTrackMomentum_/errorTrackMomentum_);
    // finalMomentumError = sqrt(finalMomentumVariance);
    //       }
    //   }
    
    //new comb

    bool eleIsNotInCombination = false ;
    if ( (eOverP  > 1 + 2.5*errorEOverP) || (eOverP  < 1 - 2.5*errorEOverP) || (eOverP < 0.8) || (eOverP > 1.3) )
      { eleIsNotInCombination = true ; }
    if (eleIsNotInCombination)
      {
	if (eOverP > 1)
	  { finalMomentum = scEnergy ; finalMomentumError = newEnergyError_ ; }
	else
	  {
	    if (elClass == reco::GsfElectron::GOLDEN)
	      { finalMomentum = scEnergy; finalMomentumError = newEnergyError_; }
	    if (elClass == reco::GsfElectron::BIGBREM)
	      {
		if (scEnergy<36)
		  { finalMomentum = trackMomentum ; finalMomentumError = errorTrackMomentum_ ; }
		else
		  { finalMomentum = scEnergy ; finalMomentumError = newEnergyError_ ; }
	      }
	    if (elClass == reco::GsfElectron::BADTRACK)
	      { finalMomentum = scEnergy; finalMomentumError = newEnergyError_ ; }
	    if (elClass == reco::GsfElectron::SHOWERING)
	      {
		if (scEnergy<30)
		  { finalMomentum = trackMomentum ; finalMomentumError = errorTrackMomentum_; }
		else
		  { finalMomentum = scEnergy; finalMomentumError = newEnergyError_;}
	      }
	    if (elClass == reco::GsfElectron::GAP)
	      {
		if (scEnergy<60)
		  { finalMomentum = trackMomentum ; finalMomentumError = errorTrackMomentum_ ; }
		else
		  { finalMomentum = scEnergy; finalMomentumError = newEnergyError_ ; }
	      }
	  }
      }
 
    else
      {
	// combination
	finalMomentum = (scEnergy/newEnergyError_/newEnergyError_ + trackMomentum/errorTrackMomentum_/errorTrackMomentum_) /
	  (1/newEnergyError_/newEnergyError_ + 1/errorTrackMomentum_/errorTrackMomentum_);
	float finalMomentumVariance = 1 / (1/newEnergyError_/newEnergyError_ + 1/errorTrackMomentum_/errorTrackMomentum_);
	finalMomentumError = sqrt(finalMomentumVariance);
      }
  } 


  
  // }
  
  math::XYZTLorentzVector oldMomentum = electron.p4() ;
  newMomentum_ = math::XYZTLorentzVector
    ( oldMomentum.x()*finalMomentum/oldMomentum.t(),
      oldMomentum.y()*finalMomentum/oldMomentum.t(),
      oldMomentum.z()*finalMomentum/oldMomentum.t(),
      finalMomentum ) ;
  finalMomentumError_ =  finalMomentumError;  
  std::cout << "[ElectronEnergCorrector] old comb momentum " << oldMomentum.t() << " new comb momentum " << newMomentum_.t() << std::endl;
 
 
}

