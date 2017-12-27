#ifndef GenAnalyzer_h
#define GenAnalyzer_h

#include <functional>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <utility>
#include <string>
#include <vector>
#include <memory>
#include <cmath>


#include <TH1.h>

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// Class to create TTree variables
#include <TFile.h> 
#include <TTree.h> 

using namespace std;
using namespace edm;
using namespace reco;

class GenAnalyzer : public edm::EDAnalyzer {
public:
  explicit GenAnalyzer(const edm::ParameterSet &params);
  virtual ~GenAnalyzer();
  
protected:
  virtual void analyze(const edm::Event &event,
		       const edm::EventSetup &es);
  std::string getParticleName(int id) const;
  
  void Initialize();

private:

  // ROOT definition
  TFile *theFile_ ;
  TTree *theTree_ ;


  edm::ESHandle<ParticleDataTable>        pdt_;
  edm::InputTag				sourceLabel;
   
  vector<float> leptonpt;
   
  float PT_h,PT_z[2],PT_l,PT_l_z[4],PT_l_stable[50],PT_l_1,PT_l_2,PT_l_3,PT_l_4,Mass_fourl[50],Mass_ZZ[50],PT_ZZ[50];
  double weightgen;
  int nevt;
  int irun,ievt,ils;

};

#endif
