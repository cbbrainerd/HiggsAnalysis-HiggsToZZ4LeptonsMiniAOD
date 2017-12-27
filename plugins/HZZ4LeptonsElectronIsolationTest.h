#ifndef  HZZ4LeptonsElectronIsolationTest_h
#define  HZZ4LeptonsElectronIsolationTest_h

/**\class HZZ4LeptonsElectronIsolationTest
 *
 *
 * Original Author: Nicola De Filippis
 * Refine electron collection to begin with
 *
 */

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectronFwd.h>
#include "DataFormats/Common/interface/ValueMap.h"

class HZZ4LeptonsElectronIsolationTest : public edm::EDProducer {
 public:
  explicit HZZ4LeptonsElectronIsolationTest(const edm::ParameterSet& );
  ~HZZ4LeptonsElectronIsolationTest();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::EDGetTokenT<edm::View<reco::GsfElectron> >  elecLabel;
  edm::EDGetTokenT<edm::ValueMap<double> >
    electronPFIsoValueChargedAllTag_,
    electronPFIsoValueChargedTag_,
    electronPFIsoValueNeutralTag_,
    electronPFIsoValueGammaTag_,
    electronPFIsoValuePUTag_;   
  
  float elecPtMin;
  float elecEtaMax;
  
  int counterelectron;
 
};
  
#endif

