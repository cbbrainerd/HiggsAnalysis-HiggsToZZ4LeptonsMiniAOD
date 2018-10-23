/**\class HZZ4LeptonsMuonIsolationProducerMu
 *
 *
 * Original Author:  Nicola De Filippis
 * 
 * Compute isolation for cones around electron candidates
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsMuonIsolationProducerMu.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

// Muons
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// ValueMap
#include "DataFormats/Common/interface/ValueMap.h"

#include <memory>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;

HZZ4LeptonsMuonIsolationProducerMu::HZZ4LeptonsMuonIsolationProducerMu(const edm::ParameterSet& pset): 
    muonTag(pset.getParameter<edm::InputTag>("MuonsLabel")),
    muTkIsoTag(pset.getParameter<edm::InputTag>("muTkIso")),
    muEcalIsoTag(pset.getParameter<edm::InputTag>("muEcalIso")),
    muHcalIsoTag(pset.getParameter<edm::InputTag>("muHcalIso")),
    coeffTk(pset.getParameter<double>("coeffTk")),
    coeffEcal(pset.getParameter<double>("coeffEcal")),
    coeffHcal(pset.getParameter<double>("coeffHcal")),
    useRelativeIso(pset.getParameter<bool>("useRelativeIso")),
    threshold(pset.getParameter<double>("threshold"))
{
    produces<reco::MuonCollection>();
    produces<reco::MuonRefVector>();
}


HZZ4LeptonsMuonIsolationProducerMu::~HZZ4LeptonsMuonIsolationProducerMu() {
}


void HZZ4LeptonsMuonIsolationProducerMu::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    auto_ptr<reco::MuonCollection> isolatedMuons(new reco::MuonCollection);
    auto_ptr<reco::MuonRefVector> isolatedMuonRefs(new reco::MuonRefVector);
    
    Handle<reco::MuonCollection> muonsH;
    iEvent.getByLabel(muonTag, muonsH);
    reco::MuonCollection muons = *muonsH;
    
    edm::Handle<edm::ValueMap<double> > muTkIsoH;
    iEvent.getByLabel(muTkIsoTag, muTkIsoH);
    const edm::ValueMap<double>& muTkIso = *muTkIsoH;

    edm::Handle<edm::ValueMap<double> > muEcalIsoH;
    iEvent.getByLabel(muEcalIsoTag, muEcalIsoH);
    const edm::ValueMap<double>& muEcalIso = *muEcalIsoH;

    edm::Handle<edm::ValueMap<double> > muHcalIsoH;
    iEvent.getByLabel(muHcalIsoTag, muHcalIsoH);
    const edm::ValueMap<double>& muHcalIso = *muHcalIsoH;

    for(std::size_t i = 0; i < muonsH->size(); i++) {
        
        reco::MuonRef muonRef(muonsH,i);
        double tkIsoVal = muTkIso[muonRef];
	double ecalIsoVal = muEcalIso[muonRef];
	double hcalIsoVal = muHcalIso[muonRef];
        double combinedIsoVal = coeffTk*tkIsoVal + coeffEcal*ecalIsoVal + coeffHcal*hcalIsoVal;
        if (useRelativeIso) combinedIsoVal /= muons[i].pt();
        
         if (combinedIsoVal < threshold){
             isolatedMuons->push_back(muons[i]);
             isolatedMuonRefs->push_back(muonRef);      
         }
    }
    
    // iEvent.put(isolatedMuons);
    // iEvent.put(isolatedMuonRefs);

    iEvent.put(std::make_unique<reco::MuonCollection>(*isolatedMuons));
    iEvent.put(std::make_unique<reco::MuonRefVector>(*isolatedMuonRefs));
}

