
/* Original author:  Nicola De Filippis - LLR -Ecole Polytechnique 
 *
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/SimpleCounter.h"


// system include files
#include <memory>

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// namespaces
using namespace edm;
using namespace std;
using namespace reco;


// Constructor
SimpleCounter::SimpleCounter(const edm::ParameterSet& pset) {

  counter=0;
	
}

// Destructor
SimpleCounter::~SimpleCounter() {
  cout << " N_events_read  = " << counter << endl;
}


void SimpleCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  counter++;	

}

void SimpleCounter::beginJob() {
	
	
}

void SimpleCounter::endJob() {
	
}


