#ifndef  HZZ4LeptonsMuonSelector_h
#define  HZZ4LeptonsMuonnSelector_h

/**\class HZZ4LeptonsMuonSelector
 *
 *
 * Original Author:  Dominique Fortin
 *
 * Refine muon collection to begin with
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
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"


class HZZ4LeptonsMuonSelector : public edm::EDProducer {
 public:
  explicit HZZ4LeptonsMuonSelector(const edm::ParameterSet& );
  ~HZZ4LeptonsMuonSelector();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  bool isGlobalMuon,isTrackerMuon;
  edm::EDGetTokenT<edm::View<pat::Muon> > muonLabel;
  float muonPtMin;
  float muonEtaMax;

};

#endif
