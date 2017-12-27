#ifndef  RegressionElectronProducer_h
#define  RegressionElectronProducer_h

/**\class RegressionElectronProducer
 *
 *
 * Original Author:  Nicola De Filippis
 * Refine electron collection to begin with
 *
 */

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

class RegressionElectronProducer : public edm::EDProducer {
 public:
  explicit RegressionElectronProducer(const edm::ParameterSet& );
  ~RegressionElectronProducer();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  void computeEpCombination( const reco::GsfElectron &, double , double, math::XYZTLorentzVector, float, float);

  edm::InputTag eleTag_;
  edm::InputTag eleRegressionEnergyErrorTag_,eleRegressionEnergyTag_;
  int counterelectron;
 
};
  
#endif

