#ifndef PhiPlusPlusMCGenProducer_h
#define PhiPlusPlusMCGenProducer_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <TH1.h>

// class to create TTree variables
#include <TFile.h>
#include <TTree.h>

using namespace std;
using namespace edm;
using namespace reco;

//the 21 final states
enum {data,eeee,eeem,eeet,eemm,eemt,eett,emem,emet,emmm,emmt,emtt,etet,etmm,etmt,ettt,mmmm,mmmt,mmtt,mtmt,mttt,tttt};

// class declaration

class PhiPlusPlusMCGenProducer : public edm::EDProducer {
 public:
  explicit PhiPlusPlusMCGenProducer(const edm::ParameterSet& iConfig);
  ~PhiPlusPlusMCGenProducer();

  static void FillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  virtual void beginRun(edm::Run&, edm::EventSetup const&);
  virtual void endRun(edm::Run&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

// Root definition

  TFile *theFile_;
  TTree *theTree_;

 edm::InputTag sourceLabel;

 int nevt;
 int d1,d2;
 double ee,em,et,mm,mt,tt;
 int fstate, modpoint;
 int dc1,dc2,dc3,dc4;
 bool debug;
};

#endif






