/**\class HZZ4LeptonsElectronOrdering.cc
 *
 * Original Author:  Nicola De Filippis
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsElectronOrdering.h"

#include "FWCore/Framework/interface/ESHandle.h"

// Electrons:
//#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
//#include <DataFormats/EgammaCandidates/interface/GsfElectronFwd.h>
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


struct SortCandByDecreasingPt {
  bool operator()( const pat::Electron &c1, const pat::Electron &c2) const {
    return c1.pt() > c2.pt();
  }
};

// constructor
HZZ4LeptonsElectronOrdering::HZZ4LeptonsElectronOrdering(const edm::ParameterSet& pset) {

  elecLabel   = consumes<pat::ElectronCollection>(pset.getParameter<edm::InputTag>("electronCollection"));
 
  produces<pat::ElectronCollection>(); 

  counterelectron=0;	


}


HZZ4LeptonsElectronOrdering::~HZZ4LeptonsElectronOrdering() {

}

void HZZ4LeptonsElectronOrdering::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<pat::ElectronCollection> Gelec( new pat::ElectronCollection );

  // Get all pixel match electron candidates
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(elecLabel, electrons);

  if (electrons.isValid()){
    // Loop over Electrons
    for (unsigned int i = 0; i < electrons->size(); ++i) {
      Ref<pat::ElectronCollection> electronRef(electrons,i);
      Gelec->push_back( *electronRef );      
    } 
  }
  
  // Ordering leptons in pT;
  std::sort (Gelec->begin(),Gelec->end(),SortCandByDecreasingPt());

  for (pat::ElectronCollection::const_iterator cands=Gelec->begin(); cands!= Gelec->end(); ++cands) {
    cout << "Electron Ordered with pt= " << cands->pt() << endl;   
  }
  

  const string iName = "";
  //iEvent.put( Gelec, iName );
  iEvent.put(std::make_unique<pat::ElectronCollection>(*Gelec), iName );

}

