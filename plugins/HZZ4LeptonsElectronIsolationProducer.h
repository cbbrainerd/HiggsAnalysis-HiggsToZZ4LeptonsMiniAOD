#ifndef HZZ4LeptonsElectronIsolationProducer_h
#define HZZ4LeptonsElectronIsolationProducer_h

/**\class HZZ4LeptonsElectronIsolationProducer
 *
 *
 * Original Author: Nicola De Filippis
 *
 * Compute isolation for cones around electron candidates
 */
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"


class HZZ4LeptonsElectronIsolationProducer : public edm::EDProducer {

 public:
  explicit HZZ4LeptonsElectronIsolationProducer(const edm::ParameterSet&);
  ~HZZ4LeptonsElectronIsolationProducer();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::InputTag electronTag_,electronVetoTag_,muonsTag_;
  edm::InputTag tracksTag_;
  double isoCone;
  //  double ptMin;
  double isoVeto;
  // doudle maxDz;
  double isocut;

};

#endif
