#ifndef HZZ4LeptonsRunEventFilter_h
#define HZZ4LeptonsRunEventFilter_h

/* \class HZZ4LeptonsRunEventFilter
 *
 *
 * RunEventFilter
 *
 * author:  Nicola De Filippis - LLR-Ecole Polytechnique
 *
 */

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// User include files
#include <FWCore/ParameterSet/interface/ParameterSet.h>

// Class declaration
class HZZ4LeptonsRunEventFilter : public edm::EDFilter {
  
 public:
  // Constructor
  explicit HZZ4LeptonsRunEventFilter(const edm::ParameterSet&);

  // Destructor
  ~HZZ4LeptonsRunEventFilter();

  /// Get event properties to send to builder to fill seed collection
  virtual bool filter(edm::Event&, const edm::EventSetup& );

 private:

 int irun, ils, ievt;
 std::vector<int> runlist,lumislist,eventlist;

};

#endif
