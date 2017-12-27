#ifndef HZZ4LeptonsConstraintFitProducer_h
#define HZZ4LeptonsConstraintFitProducer_h

/**\class HZZ4LeptonsConstraintFitProducer 
 *
 * Original Author:  Boris Mangano  - UCSD
 * modified by N. De Filippis - Politecnico and INFN Bari
 * modified by M. Masciovecchio - Univ. Bologna
 * modified by J. Smith - UC Davis
 */


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class HZZ4LeptonsConstraintFitProducer : public edm::EDProducer {


 public:
  explicit HZZ4LeptonsConstraintFitProducer(const edm::ParameterSet&);
  ~HZZ4LeptonsConstraintFitProducer();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  edm::InputTag muonTag_, vertexTag_,RECOcollName;
  uint nParticles;

  // PG and FRC 06-07-11 try to reduce printout!
	bool debug;
  
  // ----------member data ---------------------------
  double massSqr_;
};

#endif
