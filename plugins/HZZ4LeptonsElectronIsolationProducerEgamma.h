#ifndef HZZ4LeptonsElectronIsolationProducerEgamma_h
#define HZZ4LeptonsElectronIsolationProducerEgamma_h

/**\class HZZ4LeptonsElectronIsolationProducerEgamma
 *
 *
 * Original Author: Roberto Salerno 
 *
 * Compute isolation for cones around electron candidates
 */
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <string>

class HZZ4LeptonsElectronIsolationProducerEgamma : public edm::EDProducer {
    public:
        explicit HZZ4LeptonsElectronIsolationProducerEgamma(const edm::ParameterSet&);
        ~HZZ4LeptonsElectronIsolationProducerEgamma();
    
    private:
        virtual void produce(edm::Event&, const edm::EventSetup&);
        
        edm::InputTag electronTag;
        edm::InputTag eleTkIsoTag;
        edm::InputTag eleEcalIsoTag;
        edm::InputTag eleHcalIsoTag;
        double coeffTk;
        double coeffEcal;
        double coeffHcal;
        bool useRelativeIso;
        double threshold;
};

#endif
