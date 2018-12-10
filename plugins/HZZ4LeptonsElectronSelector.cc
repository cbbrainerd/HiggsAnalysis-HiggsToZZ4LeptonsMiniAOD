/**\class HZZ4LeptonsElectronSelector.cc
 *
 * Original Author:  Nicola De Filippis
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsElectronSelector.h"

#include "FWCore/Framework/interface/ESHandle.h"

// Electrons:
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectronFwd.h>
#include "DataFormats/PatCandidates/interface/Electron.h"

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <Math/VectorUtil.h>
#include <memory>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;


// constructor
HZZ4LeptonsElectronSelector::HZZ4LeptonsElectronSelector(const edm::ParameterSet& pset) {

  elecLabel   = consumes<edm::View<pat::Electron> >(pset.getParameter<edm::InputTag>("electronCollection"));
  elecPtMin   = pset.getParameter<double>("electronPtMin");
  elecEtaMax  = pset.getParameter<double>("electronEtaMax");
   
  string iName = "selectedElectron";
  produces<pat::ElectronCollection>(); 

  counterelectron=0;	


}


HZZ4LeptonsElectronSelector::~HZZ4LeptonsElectronSelector() {

}

void HZZ4LeptonsElectronSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<pat::ElectronCollection> Gelec( new pat::ElectronCollection );

  // Get all pixel match slimmed electron candidates
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByToken(elecLabel, electrons);

    // Loop over slimmedElectrons
    for (unsigned int i = 0; i < electrons->size(); ++i) {
      Ref<edm::View<pat::Electron> > electronRef(electrons,i);
      cout << "Electron selector found with pT= " << electronRef->pt() << "and eta = "<< electronRef->eta()<<endl;
      if (electronRef->pt() >= elecPtMin && fabs(electronRef->eta()) < elecEtaMax){
	Gelec->push_back( *electronRef );
	++counterelectron;
      }
    } 

  const string iName = "";
  //iEvent.put( Gelec, iName );
  iEvent.put(std::make_unique<pat::ElectronCollection>(*Gelec), iName );

}

