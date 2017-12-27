#ifndef HZZ4LeptonsTipLipToVtxProducer_h
#define HZZ4LeptonsTipLipToVtxProducer_h

/**\class HZZ4LeptonsTipLipToVtxProducer (from HZZ2e2muIpToVtxProducer)
 *
 * Original Author:  Alexis Pompili   - Bari
 *
 * Computes for each lepton the Transverse & Longitudinal IP w.r.t to PrimaryVertex (Tip & Lip)
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
// Electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"


class HZZ4LeptonsTipLipToVtxProducer : public edm::EDProducer {

 public:

  explicit HZZ4LeptonsTipLipToVtxProducer(const edm::ParameterSet&);
  ~HZZ4LeptonsTipLipToVtxProducer();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
 
  std::string decaychannel;
  edm::EDGetTokenT<edm::View<reco::Muon> >muonTag_;
  edm::EDGetTokenT<edm::View<reco::GsfElectron> >electronTag_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > vertexTag_;
  edm::EDGetTokenT<reco::BeamSpot> offlineBeamSpot_;
  bool useBeamSpot_;
  bool debug;
  //
};

#endif
