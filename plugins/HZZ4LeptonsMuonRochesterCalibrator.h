#ifndef HZZ4LeptonsMuonRochesterCalibrator_h
#define HZZ4LeptonsMuonRochesterCalibrator_h

/**\class HZZ4LeptonsMuonRochesterCalibrator
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

class HZZ4LeptonsMuonRochesterCalibrator: public edm::EDProducer {
 public:
  explicit HZZ4LeptonsMuonRochesterCalibrator(const edm::ParameterSet& );
  ~HZZ4LeptonsMuonRochesterCalibrator();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  bool isData;	
  edm::EDGetTokenT<edm::View<pat::Muon> > muonLabel;  
  edm::EDGetTokenT<edm::Association<std::vector<reco::GenParticle> > > goodMuonMCMatch_;
  edm::EDGetTokenT<reco::CandidateCollection> myMuons_;
  bool MCTruth;
  
};

#endif
