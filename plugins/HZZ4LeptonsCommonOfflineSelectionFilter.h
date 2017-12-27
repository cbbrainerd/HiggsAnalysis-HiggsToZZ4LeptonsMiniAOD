#ifndef HZZ4LeptonsCommonOfflineSelectionFilter_h
#define HZZ4LeptonsCommonOfflineSelectionFilter_h

/* \class HZZ4LeptonsCommonOfflineSelectionFilter
 *
 *
 * TightIsolation for electrons and muons
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
#include <TH1.h>
#include <TFile.h>
#include <TDirectory.h>

// Class declaration
class HZZ4LeptonsCommonOfflineSelectionFilter : public edm::EDFilter {
  
 public:
  // Constructor
  explicit HZZ4LeptonsCommonOfflineSelectionFilter(const edm::ParameterSet&);

  // Destructor
  ~HZZ4LeptonsCommonOfflineSelectionFilter();

  /// Get event properties to send to builder to fill seed collection
  virtual bool filter(edm::Event&, const edm::EventSetup& );

  void endJob();

 private:

  std::string decaychannel,offselinst;
  std::vector<std::string> offseltags;
  std::string rootFileName_;
  TFile *theFile_ ;
  TH1I *Eff;
  char *locdir;
  TDirectory *adir;

  // Efficiency variables
  int nPresel;
  int nTightEle;
  int nTightMu;
  int nVert;
  std::vector<int> vcounter;

  std::string offSelectFileName;
  
 

};

#endif
