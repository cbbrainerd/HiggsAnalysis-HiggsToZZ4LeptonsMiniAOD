#ifndef HZZ4LeptonsPreselectionFilter_h
#define HZZ4LeptonsPreselectionFilter_h

/* \class HZZ4LeptonsPreselectionFilter
 *
 *
 * Analysis preselection:
 * m_ll > 12 GeV
 * m_H  > 100 GeV
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
class HZZ4LeptonsCommonPreselectionFilter : public edm::EDFilter {
  
 public:
  // Constructor
  explicit HZZ4LeptonsCommonPreselectionFilter(const edm::ParameterSet&);

  // Destructor
  ~HZZ4LeptonsCommonPreselectionFilter();

  /// Get event properties to send to builder to fill seed collection
  virtual bool filter(edm::Event&, const edm::EventSetup& );

  void endJob();

 private:
  
  std::string decaychannel,preselinst;
  std::vector<std::string> preseltags;
  std::string rootFileName_;
  TFile *theFile_ ;
  TH1I *Eff;
  char *locdir;
  TDirectory *adir;

  // Efficiency variables
  int nSkim;
  int nElec;
  int nEId;
  int nMu;
  int nMuId;
  int nZEE;
  int nZMM;
  int nH4leptons;
  int nLooseIsolEle;
  int nLooseIsolMu;
  std::vector<int> vcounter;

  std::string preSelectFileName;
  
  
 

};

#endif
