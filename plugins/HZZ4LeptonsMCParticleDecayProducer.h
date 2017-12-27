#ifndef HZZ4LeptonsMCParticleDecayProducer_h
#define HZZ4LeptonsMCParticleDecayProducer_h

/** \class HZZ4LeptonsMCParticleDecayProducer
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
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <vector>

class HZZ4LeptonsMCParticleDecayProducer : public edm::EDProducer {
 public:
  explicit HZZ4LeptonsMCParticleDecayProducer(const edm::ParameterSet&);
  ~HZZ4LeptonsMCParticleDecayProducer();

 private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::InputTag genCandidates_;
  int motherPdgId_;
  std::vector<int> firstdaughtersPdgId_,seconddaughtersPdgId_;
  std::string decayChain_;
  size_t daughtersize,firstdaughtersize,seconddaughtersize;
  std::vector<std::string> valiasleptons,valiasbosons;
  
};

#endif
