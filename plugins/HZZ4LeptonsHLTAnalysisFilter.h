#ifndef HZZ4LeptonsHLTAnalysisFilter_h
#define HZZ4LeptonsHLTAnalysisFilter_h

/* \class HZZ4LeptonsHLTAnalysisFilter
 *
 *
 * HLTAnalysisFilter
 *
 * author:  Nicola De Filippis - Politecnico di Bari
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
class HZZ4LeptonsHLTAnalysisFilter : public edm::EDFilter {
  
 public:
  // Constructor
  explicit HZZ4LeptonsHLTAnalysisFilter(const edm::ParameterSet&);

  // Destructor
  ~HZZ4LeptonsHLTAnalysisFilter();

  /// Get event properties to send to builder to fill seed collection
  virtual bool filter(edm::Event&, const edm::EventSetup& );

  void respondToOpenInputFile(edm::FileBlock const& fb);

 private:
  std::string inputfileName;
  edm::EDGetTokenT<std::vector<std::string> > HLTInfoFired;


};

#endif
