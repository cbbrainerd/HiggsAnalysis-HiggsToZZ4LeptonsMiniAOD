/**\class HZZ4LeptonsElectronIsolationTest.cc
 *
 * Original Author:  Nicola De Filippis
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsElectronIsolationTest.h"

#include "FWCore/Framework/interface/ESHandle.h"

// Electrons:
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectronFwd.h>

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
HZZ4LeptonsElectronIsolationTest::HZZ4LeptonsElectronIsolationTest(const edm::ParameterSet& pset) {

  elecLabel   = consumes<edm::View<reco::GsfElectron> >(pset.getParameter<edm::InputTag>("electronCollection"));
  elecPtMin   = pset.getParameter<double>("electronPtMin");
  elecEtaMax  = pset.getParameter<double>("electronEtaMax");
   
  electronPFIsoValueChargedAllTag_= consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueChargedAll"));
  electronPFIsoValueChargedTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueCharged"));
  electronPFIsoValueNeutralTag_   = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueNeutral"));
  electronPFIsoValueGammaTag_     = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValueGamma"));
  electronPFIsoValuePUTag_        = consumes<edm::ValueMap<double> >(pset.getParameter<edm::InputTag>("ElectronPFIsoValuePU"));


  string iName = "selectedElectron";
  produces<reco::GsfElectronCollection>(); 

  counterelectron=0;	


}


HZZ4LeptonsElectronIsolationTest::~HZZ4LeptonsElectronIsolationTest() {

}

void HZZ4LeptonsElectronIsolationTest::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<reco::GsfElectronCollection> Gelec( new reco::GsfElectronCollection );

  // Get all pixel match GSF electron candidates
  edm::Handle<edm::View<GsfElectron> > electrons;
  iEvent.getByToken(elecLabel, electrons);

  // Particle Flow Isolation
  edm::Handle<edm::ValueMap<double> > isoPFChargedAllelemap;
  iEvent.getByToken(electronPFIsoValueChargedAllTag_, isoPFChargedAllelemap);
    
  edm::Handle<edm::ValueMap<double> > isoPFChargedelemap;
  iEvent.getByToken(electronPFIsoValueChargedTag_, isoPFChargedelemap);
  
  edm::Handle<edm::ValueMap<double> > isoPFNeutralelemap;
  iEvent.getByToken(electronPFIsoValueNeutralTag_, isoPFNeutralelemap);
  
  edm::Handle<edm::ValueMap<double> > isoPFGammaelemap;
  iEvent.getByToken(electronPFIsoValueGammaTag_, isoPFGammaelemap);
  
  edm::Handle<edm::ValueMap<double> > isoPFPUelemap;
  iEvent.getByToken(electronPFIsoValuePUTag_, isoPFPUelemap);

  if (electrons.isValid()){
    // Loop over GsfElectrons
    for (unsigned int i = 0; i < electrons->size(); ++i) {
      Ref<edm::View<reco::GsfElectron> > eletrackref(electrons,i);
      cout << "Electron selector found with pT= " << eletrackref->pt() << endl;
      
      // PF isolation
      cout << "chAllPart= "     << (*isoPFChargedAllelemap)[eletrackref] 
	   << " PFchHad= "       << (*isoPFChargedelemap)[eletrackref] 
	   << " PFneuHad= "      << (*isoPFNeutralelemap)[eletrackref] 
	   << " photon= "        << (*isoPFGammaelemap)[eletrackref] 
	   << " PFPUchAllPart= " << (*isoPFPUelemap)[eletrackref] << endl;

      if (eletrackref->pt() >= elecPtMin && fabs(eletrackref->eta()) < elecEtaMax){
	Gelec->push_back( *eletrackref );
	++counterelectron;
      }
    } 
  }

  const string iName = "";
  //iEvent.put( Gelec, iName );
  iEvent.put(std::make_unique<reco::GsfElectronCollection>(*Gelec), iName );

}

