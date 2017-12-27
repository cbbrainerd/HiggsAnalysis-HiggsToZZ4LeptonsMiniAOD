#ifndef  HZZ4LeptonsMuonCalibrator_h
#define  HZZ4LeptonsMuonCalibrator_h

/**\class HZZ4LeptonsMuonCalibrator
 *
 *
 * Original Author:  Nicola De Filippis
 *
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
#include "DataFormats/PatCandidates/interface/Muon.h"

class HZZ4LeptonsMuonCalibrator : public edm::EDProducer {
 public:
  explicit HZZ4LeptonsMuonCalibrator(const edm::ParameterSet& );
  ~HZZ4LeptonsMuonCalibrator();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  bool isData;	
  edm::EDGetTokenT<edm::View<pat::Muon> > muonLabel;

};

#endif
