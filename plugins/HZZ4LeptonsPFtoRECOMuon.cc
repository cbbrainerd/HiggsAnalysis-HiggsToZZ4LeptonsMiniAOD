/**\class HZZ4LeptonsPFtoRECOMuon.cc
 *
 * Original Author:  Nicola De Filippis 
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsPFtoRECOMuon.h"

#include "FWCore/Framework/interface/ESHandle.h"

// PF candidates:
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"

// system include files
#include <memory>


using namespace edm;
using namespace std;
using namespace reco;

// constructor
HZZ4LeptonsPFtoRECOMuon::HZZ4LeptonsPFtoRECOMuon(const edm::ParameterSet& pset) {
  pfLabel        = consumes<edm::View<reco::PFCandidate> >(pset.getParameter<edm::InputTag>("pfCollection"));
  string alias;
  produces<reco::MuonCollection>(); 
}


HZZ4LeptonsPFtoRECOMuon::~HZZ4LeptonsPFtoRECOMuon() {
 
}


//
// member functions
//
void HZZ4LeptonsPFtoRECOMuon::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {



  // PF candidate to RECO Muon
  auto_ptr<reco::MuonCollection> muon( new reco::MuonCollection );

  edm::Handle<edm::View<PFCandidate> > particles;    
  iEvent.getByToken(pfLabel, particles);

  // Loop over PF candidates
  for (edm::View<reco::PFCandidate>::const_iterator mIter = particles->begin(); mIter != particles->end(); ++mIter ) {
    if ( abs(mIter->pdgId()) == 13  && mIter->muonRef().isNonnull()) {
      // cout << "PF muons found" << mIter->pt() << endl;
      muon->push_back( *mIter->muonRef() );   
    }
  }

  
  const string iName = "";
  iEvent.put( muon, iName );

}

