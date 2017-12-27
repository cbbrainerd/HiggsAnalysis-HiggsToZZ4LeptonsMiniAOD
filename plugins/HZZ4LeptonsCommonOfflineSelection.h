#ifndef HZZ4LeptonsCommonOfflineSelection_h
#define HZZ4LeptonsCommonOfflineSelection_h

/* \class HZZ4LeptonsCommonOfflineSelection
 *
 *
 * Analysis selection:
 * - tight isolation
 * - vertexing
 *
 * author:  Nicola De Filippis - LLR-Ecole Polytechnique
 *
 */

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"


// Class declaration
class HZZ4LeptonsCommonOfflineSelection : public edm::EDProducer {
  
 public:
  // Constructor
  explicit HZZ4LeptonsCommonOfflineSelection(const edm::ParameterSet&);

  // Destructor
  ~HZZ4LeptonsCommonOfflineSelection();

 private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool match(double mass, double pt, int charge, const reco::CandidateCollection *c1Coll);

  std::string decaychannel;
  bool useBestCandidate;

  edm::InputTag electronTag_,muonTag_,BestCandidatesLeptonsTag_;
  edm::InputTag electronMapTag_,muonMapTag_;
  edm::InputTag isoVarTagElectrons;
  std::vector<double> isoVarCutElectrons;

  edm::InputTag isoVarTagMuons;
  std::vector<double> isoVarCutMuons;

  edm::InputTag electronTag_Vert,muonTag_Vert;
  edm::InputTag electronMapTag_Vert,muonMapTag_Vert;
  std::vector<double> vertVarCut;

};

#endif
