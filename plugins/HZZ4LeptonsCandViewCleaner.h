#ifndef HZZ4LeptonsCandViewCleaner_h
#define HZZ4LeptonsCandViewCleaner_h

////////////////////////////////////////////////////////////////////////////////
//
// HZZ4LeptonsCandViewCleaner
// --------------
// 
////////////////////////////////////////////////////////////////////////////////


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
  
#include <memory>
#include <vector>
#include <sstream>
  
  
////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////

class HZZ4LeptonsCandViewCleaner : public edm::EDProducer
{
public:
  // construction/destruction
  HZZ4LeptonsCandViewCleaner(const edm::ParameterSet& iConfig);
  virtual ~HZZ4LeptonsCandViewCleaner();

  // member functions
  void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
  void endJob();
  
private:
  // member data
  edm::InputTag              srcCands_;
  edm::InputTag              srcObjects_;

  std::string  moduleLabel_;
  unsigned int nCandidatesTot_;
  unsigned int nCandidatesClean_;
};



#endif
