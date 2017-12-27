#ifndef HZZ4LeptonsIpToVtxProducer_h
#define HZZ4LeptonsIpToVtxProducer_h

/**\class HZZ4LeptonsIpToVtxProducer (from HZZ2e2muIpToVtxProducer)
 *
 * Original Author:  Alexis Pompili   - Bari
 *
 * Computes for each lepton the IP3D significance
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
// Beam Spot                                                                                                                                                                      
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
// Vertex utilities                                                                                                                                                               
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
// Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include "DataFormats/PatCandidates/interface/Muon.h"

// Electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

class HZZ4LeptonsIpToVtxProducer : public edm::EDProducer {

 public:

  explicit HZZ4LeptonsIpToVtxProducer(const edm::ParameterSet&);
  ~HZZ4LeptonsIpToVtxProducer();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  std::string decaychannel;
  edm::EDGetTokenT<edm::View<pat::Muon> >muonTag_;
  edm::EDGetTokenT<edm::View<pat::Electron> >electronTag_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > vertexTag_;
  edm::EDGetTokenT<reco::BeamSpot> offlineBeamSpot_;
  bool useBeamSpot_;
  bool debug;
  //
};

#endif
