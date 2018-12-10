/**\class HZZ4LeptonsPFfsrPhoton.cc
 *
 * Original Author:  Nicola De Filippis 
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsPFfsrPhoton.h"

#include "FWCore/Framework/interface/ESHandle.h"

// PF candidates:
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"


// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/ValueMap.h"

// system include files
#include <memory>
#include "TMath.h"

using namespace edm;
using namespace std;
using namespace reco;

// constructor
HZZ4LeptonsPFfsrPhoton::HZZ4LeptonsPFfsrPhoton(const edm::ParameterSet& pset) {
  pfLabel        = pset.getParameter<edm::InputTag>("pfCollection");
  string alias;
  produces<reco::PhotonCollection>(); 
}


HZZ4LeptonsPFfsrPhoton::~HZZ4LeptonsPFfsrPhoton() {
 
}


//
// member functions
//
void HZZ4LeptonsPFfsrPhoton::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<reco::PhotonCollection> photon( new reco::PhotonCollection );

  // PF photon for FSR recovery
  edm::Handle<edm::View<PFCandidate> > particles;    
  iEvent.getByLabel(pfLabel.label(), particles);
  reco::Photon phfsr;

  // Loop over PF candidates for photons

  for (edm::View<reco::PFCandidate>::const_iterator pIter = particles->begin(); pIter != particles->end(); ++pIter ) {
    if ( pIter->charge()==0 && abs(pIter->pdgId()) == 22 && pIter->pt()>2 && fabs(pIter->eta())<2.4 ) {
      cout << "PF photons found with pt= " << pIter->pt() << " and eta=" << pIter->eta() << endl;

      phfsr.setP4(pIter->p4());
      phfsr.setCharge(0.);
      phfsr.setPdgId(22);

      photon->push_back(phfsr);   
    }    
  }

 


  /*
  // PF photon candidate from ecal energy of the muon
  reco::Photon mufsr; 


  for (edm::View<reco::PFCandidate>::const_iterator mIter = particles->begin(); mIter != particles->end(); ++mIter ) {
    if ( abs(mIter->pdgId()) == 13  && mIter->muonRef().isNonnull()) {
      // cout << "PF muons found" << mIter->pt() << endl;
      //muon->push_back( *mIter->photonRef() );   

      // Keep only those with a swallowed ECAL energy > 2 GeV
      if ( mIter->ecalEnergy() < 2.0 ) continue;
      // Request the pt to be in excess of 2 GeV/c
      double sintet = mIter->pt()/mIter->energy();
      double phpt = mIter->ecalEnergy() * sintet;
      if ( phpt < 2.0 ) continue;
      
      double ratio = mIter->ecalEnergy()/mIter->energy();
      mufsr.setP4(mIter->p4() * ratio);
      mufsr.setCharge(0.);
      mufsr.setMass(0.);
      mufsr.setPdgId(22);

      cout << "mu PF photons found with pt= " << mufsr.pt() << " and eta=" << mIter->p4().eta() << endl;   
      photon->push_back(mufsr);
    }
  }
  */
   // filling map
  
  const string iName = "";
  //iEvent.put( photon, iName );
  iEvent.put(std::make_unique<reco::PhotonCollection>(*photon), iName );

}

