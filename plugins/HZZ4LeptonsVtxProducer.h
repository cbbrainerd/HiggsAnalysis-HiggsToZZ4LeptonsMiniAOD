#ifndef HZZ4LeptonsVtxProducer_h
#define HZZ4LeptonsVtxProducer_h

/**\class HZZ4LeptonsVtxProducer 
 *
 * Original Author:  Alexis Pompili   - Bari
 *
 * Computes for each lepton the IP3D significance and other VTX stuff
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// the following seems to not work:
//#include "FWCore/ParameterSet/interface/InputTag.h" //so i use the following one
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/RefToBase.h"

//using namespace edm;
//using namespace reco;

class HZZ4LeptonsVtxProducer : public edm::EDProducer {

 public:

  explicit HZZ4LeptonsVtxProducer(const edm::ParameterSet&);
  ~HZZ4LeptonsVtxProducer();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  std::string decaychannel;
  edm::InputTag muonTag_, electronTag_;
  //
  float getMuDauIP3Dsignif(reco::TrackRef muDauTrack,
			  //const reco::Candidate* muDau,
                          const edm::ESHandle<TransientTrackBuilder> TrackBuilder,
                          reco::Vertex primaryVertex);
  //
  float getMuDauIP2Dsignif(reco::TrackRef muDauTrack,
			  //const reco::Candidate* muDau,
                          const edm::ESHandle<TransientTrackBuilder> TrackBuilder,
                          reco::Vertex primaryVertex);
  //
  /////
  //
  float getEleDauIP3Dsignif(reco::GsfTrackRef eleDauTrack,
                          const edm::ESHandle<TransientTrackBuilder> TrackBuilder,
                          reco::Vertex primaryVertex);
  //
  float getEleDauIP2Dsignif(reco::GsfTrackRef eleDauTrack,
                          const edm::ESHandle<TransientTrackBuilder> TrackBuilder,
                          reco::Vertex primaryVertex);

};

#endif
