#ifndef  HZZ4LeptonsPFJetSelector_h
#define  HZZ4LeptonsPFJetSelector_h

/**\class HZZ4LeptonsPFJetIDSelector
 *
 *
 * Original Author:  Nicola De Filippis
 *
 * Refine PFJet collection to begin with
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
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"


class HZZ4LeptonsPFJetSelector : public edm::EDProducer {
 public:
  explicit HZZ4LeptonsPFJetSelector(const edm::ParameterSet& );
  ~HZZ4LeptonsPFJetSelector();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  bool isLoosePFJetID,isMediumPFJetID,isTightPFJetID;
  edm::EDGetTokenT<edm::View<pat::Jet> > pfjetsLabel;

  ////JEC
  //edm::EDGetTokenT<reco::JetCorrector> mJetCorrector;

};

#endif
