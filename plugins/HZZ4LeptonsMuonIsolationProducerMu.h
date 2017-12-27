#ifndef HZZ4LeptonsMuonIsolationProducerMu_h
#define HZZ4LeptonsMuonIsolationProducerMu_h

/**\class HZZ4LeptonsMuonIsolationProducerMu
 *
 *
 * Original Author: Niocla De Filippis
 *
 * Compute isolation for cones around muon candidates
 */
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <string>

class HZZ4LeptonsMuonIsolationProducerMu : public edm::EDProducer {
    public:
        explicit HZZ4LeptonsMuonIsolationProducerMu(const edm::ParameterSet&);
        ~HZZ4LeptonsMuonIsolationProducerMu();
    
    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);
        
        edm::InputTag muonTag;
        edm::InputTag muTkIsoTag;
        edm::InputTag muEcalIsoTag;
        edm::InputTag muHcalIsoTag;
        double coeffTk;
        double coeffEcal;
        double coeffHcal;
        bool useRelativeIso;
        double threshold;
};

#endif
