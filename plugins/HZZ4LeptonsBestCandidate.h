#ifndef HZZ4LeptonsBestCandidate_h
#define HZZ4LeptonsBestCandidate_h

/** \class HZZ4LeptonsBestCandidate
 *
 * Original Author:  Nicola De Filippis
 *
 */

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <vector>

class HZZ4LeptonsBestCandidate : public edm::EDProducer {
 public:
  explicit HZZ4LeptonsBestCandidate(const edm::ParameterSet&);
  ~HZZ4LeptonsBestCandidate();

 private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool find(const std::auto_ptr<reco::CandidateCollection>& c1Coll, const reco::Candidate& c2);

  // PG and FRC 06-07-11 try to reduce printout!
  bool debug;

  std::string decaychannel;
  std::string decayChain_;
  std::vector<std::string> valiasleptons,valiasbosons;
  std::vector<edm::InputTag> RECOcollName;
 
};

#endif
