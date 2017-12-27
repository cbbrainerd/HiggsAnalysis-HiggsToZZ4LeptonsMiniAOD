#ifndef  HZZ4LeptonsElectronOrdering_h
#define  HZZ4LeptonsElectronOrdering_h

/**\class HZZ4LeptonsElectronOrdering
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
#include "DataFormats/PatCandidates/interface/Electron.h"

class HZZ4LeptonsElectronOrdering : public edm::EDProducer {
 public:
  explicit HZZ4LeptonsElectronOrdering(const edm::ParameterSet& );
  ~HZZ4LeptonsElectronOrdering();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::EDGetTokenT<pat::ElectronCollection> elecLabel;

  int counterelectron;
 
};
  
#endif

