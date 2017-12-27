#ifndef  HZZ4LeptonsPFfsrPhoton_h
#define  HZZ4LeptonsMuonnSelector_h

/**\class HZZ4LeptonsPFfsrPhoton
 *
 *
 * Original Author:  Nicola De Filippis
 *
 * Create muon reco collection from PF muons
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

class HZZ4LeptonsPFfsrPhoton : public edm::EDProducer {
 public:
  explicit HZZ4LeptonsPFfsrPhoton(const edm::ParameterSet& );
  ~HZZ4LeptonsPFfsrPhoton();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  bool isMuon;
  edm::InputTag pfLabel;

};

#endif
