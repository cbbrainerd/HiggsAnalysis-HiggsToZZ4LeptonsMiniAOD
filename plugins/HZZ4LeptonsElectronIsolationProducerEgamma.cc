/**\class HZZ4LeptonsElectronIsolationProducerEgamma
 *
 *
 * Original Author:  Roberto Salerno (Univ. Milano Bicocca)
 * Modified by: Adish P. Vartak
 * 
 * Compute isolation for cones around electron candidates
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsElectronIsolationProducerEgamma.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

// Electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

// ValueMap
#include "DataFormats/Common/interface/ValueMap.h"

#include <memory>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;

HZZ4LeptonsElectronIsolationProducerEgamma::HZZ4LeptonsElectronIsolationProducerEgamma(const edm::ParameterSet& pset): 
    electronTag(pset.getParameter<edm::InputTag>("ElectronsLabel")),
    eleTkIsoTag(pset.getParameter<edm::InputTag>("eleTkIso")),
    eleEcalIsoTag(pset.getParameter<edm::InputTag>("eleEcalIso")),
    eleHcalIsoTag(pset.getParameter<edm::InputTag>("eleHcalIso")),
    coeffTk(pset.getParameter<double>("coeffTk")),
    coeffEcal(pset.getParameter<double>("coeffEcal")),
    coeffHcal(pset.getParameter<double>("coeffHcal")),
    useRelativeIso(pset.getParameter<bool>("useRelativeIso")),
    threshold(pset.getParameter<double>("threshold"))
{
    produces<pat::ElectronCollection>();
    produces<pat::ElectronRefVector>();
}


HZZ4LeptonsElectronIsolationProducerEgamma::~HZZ4LeptonsElectronIsolationProducerEgamma() {
}


void HZZ4LeptonsElectronIsolationProducerEgamma::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    auto_ptr<pat::ElectronCollection> isolatedElectrons(new pat::ElectronCollection);
    auto_ptr<pat::ElectronRefVector> isolatedElectronRefs(new pat::ElectronRefVector);
    
    Handle<pat::ElectronCollection> electronsH;
    iEvent.getByLabel(electronTag, electronsH);
    pat::ElectronCollection electrons = *electronsH;
    
    edm::Handle<edm::ValueMap<double> > eleTkIsoH;
    iEvent.getByLabel(eleTkIsoTag, eleTkIsoH);
    const edm::ValueMap<double>& eleTkIso = *eleTkIsoH;

    edm::Handle<edm::ValueMap<double> > eleEcalIsoH;
    iEvent.getByLabel(eleEcalIsoTag, eleEcalIsoH);
    const edm::ValueMap<double>& eleEcalIso = *eleEcalIsoH;

    edm::Handle<edm::ValueMap<double> > eleHcalIsoH;
    iEvent.getByLabel(eleHcalIsoTag, eleHcalIsoH);
    const edm::ValueMap<double>& eleHcalIso = *eleHcalIsoH;

    for(std::size_t i = 0; i < electronsH->size(); i++) {
        
        pat::ElectronRef electronRef(electronsH,i);
        double tkIsoVal = eleTkIso[electronRef];
	double ecalIsoVal = eleEcalIso[electronRef];
	double hcalIsoVal = eleHcalIso[electronRef];
        double combinedIsoVal = coeffTk*tkIsoVal + coeffEcal*ecalIsoVal + coeffHcal*hcalIsoVal;
        if (useRelativeIso) combinedIsoVal /= electrons[i].pt();
        
         if (combinedIsoVal < threshold){
             isolatedElectrons->push_back(electrons[i]);
             isolatedElectronRefs->push_back(electronRef);      
         }
    }
    
    // iEvent.put(isolatedElectrons);
    // iEvent.put(isolatedElectronRefs);
    iEvent.put(std::make_unique<pat::ElectronCollection>(*isolatedElectrons));
    iEvent.put(std::make_unique<pat::ElectronRefVector>(*isolatedElectronRefs));
}

