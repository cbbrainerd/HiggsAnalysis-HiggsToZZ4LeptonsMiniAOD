#ifndef HZZ4LeptonsCommonPreselection_h
#define HZZ4LeptonsCommonPreselection_h

/* \class HZZ4LeptonsCommonPreSelection
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
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

// Class declaration
class HZZ4LeptonsCommonPreselection : public edm::EDProducer {
  
 public:
  // Constructor
  explicit HZZ4LeptonsCommonPreselection(const edm::ParameterSet&);

  // Destructor
  ~HZZ4LeptonsCommonPreselection();

 private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  std::string decaychannel;
  edm::InputTag electronTag_,muonTag_;
  edm::InputTag zToEETag_,zToMuMuTag_,hTozzTo4leptonsTag_;
  edm::InputTag muonLooseIsolTag_, electronLooseIsolTag_; 

  int nfourlept,nElectron,nMuon,nLepton,nZEE,nZMuMu,nHiggs,nLooseIsolEle,nLooseIsolElepos,nLooseIsolEleneg,nLooseIsolMu,nLooseIsolMupos,nLooseIsolMuneg;

  int nEle_cut, nMu_cut;
  double eeMass_cut,mumuMass_cut;
  int  numberOfeeCombs_cut,numberOfmumuCombs_cut;
  double fourlMass_cut;
  int  numberOf4lCombs_cut;
  int nlooseEle_cut, nlooseMu_cut;

};

#endif
