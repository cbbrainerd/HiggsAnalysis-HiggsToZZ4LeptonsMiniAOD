#ifndef ConvValueMapProd_h
#define ConvValueMapProd_h

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

class ConvValueMapProd : public edm::EDProducer {
    public:
        explicit ConvValueMapProd(const edm::ParameterSet&);
        ~ConvValueMapProd();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::EDGetTokenT<std::vector<reco::GsfElectron> > gsfLabel_;
        edm::EDGetTokenT<std::vector<reco::Track> > tkLabel_;
};

#endif

