/* \class HZZ4LeptonsRunEventFilter
 *
 *
 * Run Event filter
 *
 * author:     Nicola De Filippis   - LLR-Ecole Polytechnique
 */


// system include files
#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsRunEventFilter.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <vector>

// namespaces
using namespace edm;
using namespace std;

// Constructor
HZZ4LeptonsRunEventFilter::HZZ4LeptonsRunEventFilter(const edm::ParameterSet& pset) {

   runlist = pset.getParameter<vector<int> >("RunList");
   lumislist = pset.getParameter<vector<int> >("LumiSectionList");
   eventlist = pset.getParameter<vector<int> >("EventList");                         
  
}


// Destructor
HZZ4LeptonsRunEventFilter::~HZZ4LeptonsRunEventFilter() {

}


// Filter Run Event
bool HZZ4LeptonsRunEventFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  // Dump Run and Event
  irun=iEvent.id().run();
  ievt=iEvent.id().event();
  ils=iEvent.luminosityBlock();

  for (unsigned int i=0; i<runlist.size(); i++){
    if ( irun==runlist.at(i) && ievt==eventlist.at(i) && ils==lumislist.at(i) ) {
      cout << "Found Run Event= " << irun << " " << ievt << " " << ils << endl;
          return true;          
        }
  }

  return false;

}

