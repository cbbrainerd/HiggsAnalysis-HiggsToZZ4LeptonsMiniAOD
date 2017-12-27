/* \class HZZ4leptonsMCGenFilter
 *
 *  Filter of h->zz->4leptons channel at the level of generation
 * author:  Simranjit Singh Chhibra & Nicola De Filippis & Gurpreet Singh
 *
 */

// system include files
#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsMCGenFilter.h"

// User include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// C++
#include <iostream>
#include <vector>

// Constructor
HZZ4LeptonsMCGenFilter::HZZ4LeptonsMCGenFilter(const edm::ParameterSet& pset) {
  
  // LeptonFlavour for H->ZZ final states
  // 0  = 4l including tau
  // 1  = 4 mu
  // 2  = 4 e
  // 3  = 2e 2mu
  // 4  = 2e 2tau
  // 5  = 2mu2tau
  // 6  = 4tau
  // 7  = 2e 2tau(LH or HL)
  // 8  = 2e 2tau(HH)
  // 9  = 2e 2tau(LL)
  // 10 = 2mu2tau(LH or HL)
  // 11 = 2mu2tau(HH)
  // 12 = 2mu2tau(LL) 
  // 13 = zlightjets  
  // 14 = zbbbar 
  // 15 = zccbar 
  
  // Local Debug flag
  gen                = pset.getParameter<edm::InputTag>("genParticles"); 
  debug              = pset.getParameter<bool>("DebugHZZ4LeptonsMCGenFilter");
  leptonFlavour      = pset.getParameter<int>("HZZ4LeptonsMCFilterLeptonFlavour");
  acceptance         = pset.getParameter<double>("acceptance");
  // zgen               = pset.getParameter<edm::InputTag>("Z"); 
  // zbbgen             = pset.getParameter<edm::InputTag>("llbBCands"); 
  // zccgen             = pset.getParameter<edm::InputTag>("llccCands"); 
  
  ikept   =0;
  evt     =0;
  lepdecay=0;
  haddecay=0;
  
}

// Destructor
HZZ4LeptonsMCGenFilter::~HZZ4LeptonsMCGenFilter() {
  
  std::cout << "number of events processed: " << evt << std::endl;
  std::cout << "number of events kept: " << ikept << std::endl;
  std::cout << "expected number of taus = " << taus << std::endl;
}

// Filter event
bool HZZ4LeptonsMCGenFilter::filter(edm::Event& event, const edm::EventSetup& setup ) {
  
  bool keepEvent   = false;
  evt++;
  
  bool FourL    = false;
  bool FourE    = false;
  bool FourM    = false;
  bool TwoETwoM = false;
  bool FourTau  = false;
  bool TauEM    = false;
  bool TwoETwoTau = false;
  bool TwoMTwoTau = false;
  bool lepton   = false;
  bool hadron   = false;
  
  //bool z        = false;
  //bool zbb=false, zbbbar   = false; 
  //bool zcc=false, zccbar   = false; 
  //bool zlightjets    = false;
  
  // get gen particle candidates 
  edm::Handle<reco::GenParticleCollection> genCandidates;	    
  event.getByLabel(gen, genCandidates);
  
  int nElec = 0, nMuon = 0, nTau = 0;
  
  for ( reco::GenParticleCollection::const_iterator mcIter=genCandidates->begin(); mcIter!=genCandidates->end(); ++mcIter ) {
    // Muons:
    std::cout<<mcIter->pdgId()<<" ";
    if ( mcIter->pdgId() == 13 || mcIter->pdgId() == -13) {
      // Mother is a Z and her mother is H
      
      if ( mcIter->mother()->pdgId() == 23 && mcIter->mother()->mother()->pdgId() == 25 && std::abs(mcIter->eta()) < acceptance ) {
       	nMuon++;std::cout<<"<==== muon"<<std::endl;
      }
    }
    // Electrons:
    if ( mcIter->pdgId() == 11 || mcIter->pdgId() == -11) {
        std::cout<<"ele"<<std::endl;
	// Mother is a Z and her mother is H
      if ( mcIter->mother()->pdgId() == 23 && mcIter->mother()->mother()->pdgId() == 25 && std::abs(mcIter->eta()) < acceptance ) {
	nElec++;
      }
    }
    // Taus:
    if ( mcIter->pdgId() == 15 || mcIter->pdgId() == -15) {
      // Mother is a Z and her mother is H
      if ( mcIter->mother()->pdgId() == 23 && mcIter->mother()->mother()->pdgId() == 25 && std::abs(mcIter->eta()) < acceptance ) {
        nTau++;
      }
    }
    
    if ( (abs(mcIter->pdgId())==11 || abs(mcIter->pdgId())==13) && abs(mcIter->mother(0)->pdgId())==15 ){
      lepdecay++;
      lepton=true;
      //cout <<"\n lepton found from tau decay with id="<<  mcIter->pdgId() << " " << getParticleName(int(mcIter->pdgId())) <<endl;    
    }
    
    // taus decay hadronically 
    else if ( 
	     ( abs(mcIter->pdgId())!=15 &&
	       abs(mcIter->pdgId())!=16 &&
	       abs(mcIter->pdgId())!=12 &&
	       abs(mcIter->pdgId())!=14 &&
	       abs(mcIter->pdgId())!=22 &&
	       abs(mcIter->pdgId())!=21 &&
	       mcIter->status()!=3)
	     &&
	     abs(mcIter->mother(0)->pdgId())==15
	     ){
      haddecay++;
      hadron=true;
      //cout <<"\n hadron found from tau decay with id="<< mcIter->pdgId() << " " << getParticleName(int(mcIter->pdgId()))  <<endl;
    }
  }
  
  // get Z Candidates
  /*edm::Handle<edm::View<reco::Candidate> >  zCandidates;
    event.getByLabel(zgen, zCandidates);
    if (zCandidates.isValid()) {  
    for (edm::View<reco::Candidate>::const_iterator mcIter=zCandidates->begin(); mcIter!=zCandidates->end(); ++mcIter ) {
    z = true;
    } 
    }*/
  
  //  get Zbbbar candidates                                                    
  /*edm::Handle<edm::View<reco::Candidate> >  ZbbbarCandidates;                                                             
  event.getByLabel(zbbgen, ZbbbarCandidates); 
  if (ZbbbarCandidates.isValid()) {  
    for (edm::View<reco::Candidate>::const_iterator mcIter=ZbbbarCandidates->begin(); mcIter!=ZbbbarCandidates->end(); ++mcIter ) { 
      zbb = true;
    } 
    }*/
  
  //  get Zccbar candidates
  /*edm::Handle<edm::View<reco::Candidate> >  ZccbarCandidates;
  event.getByLabel(zccgen, ZccbarCandidates);
  if (ZccbarCandidates.isValid()) {
    for (edm::View<reco::Candidate>::const_iterator mcIter=ZccbarCandidates->begin(); mcIter!=ZccbarCandidates->end(); ++mcIter ) {
      zcc = true;
    }
    }*/
  
  if (nElec == 4) FourE = true;
  if (nMuon == 4) FourM = true;
  if (nMuon == 2 && nElec == 2 ) TwoETwoM = true;
  if (nMuon == 2 && nTau  == 2  ) TwoMTwoTau = true;
  if (nElec == 2 && nTau  == 2  ) TwoETwoTau = true;
  if (nTau  == 4) FourTau = true;
  if ((nTau  == 2 && nElec == 2 ) || (nTau  == 2 && nMuon == 2 )) TauEM = true;
  if ( FourE || FourM || TwoETwoM || FourTau || TauEM ) FourL = true;
  
  //if (z == true && zbb ==true) zbbbar = true;
  //if (z == true && zbb ==false && zcc ==true) zccbar = true;
  //if (zbbbar == false && zccbar == false) zlightjets = true;
  
  
  if ( leptonFlavour ==  0 && FourL      ) keepEvent = true;    
  if ( leptonFlavour ==  1 && FourM      ) keepEvent = true;    
  if ( leptonFlavour ==  2 && FourE      ) keepEvent = true;    
  if ( leptonFlavour ==  3 && TwoETwoM   ) keepEvent = true;
  if ( leptonFlavour ==  4 && TwoETwoTau ) keepEvent = true;
  if ( leptonFlavour ==  5 && TwoMTwoTau ) keepEvent = true;
  if ( leptonFlavour ==  6 && FourTau    ) keepEvent = true;
  if ( leptonFlavour ==  7 && TwoETwoTau && hadron && lepton  ) keepEvent = true;
  if ( leptonFlavour ==  8 && TwoETwoTau && hadron && !lepton ) keepEvent = true;
  if ( leptonFlavour ==  9 && TwoETwoTau && !hadron && lepton ) keepEvent = true;
  if ( leptonFlavour == 10 && TwoMTwoTau && hadron && lepton  ) keepEvent = true;
  if ( leptonFlavour == 11 && TwoMTwoTau && hadron && !lepton ) keepEvent = true;
  if ( leptonFlavour == 12 && TwoMTwoTau && !hadron && lepton ) keepEvent = true;
  //if ( leptonFlavour == 13 && zlightjets )  keepEvent  = true;
  //if ( leptonFlavour == 14 && zbbbar ) keepEvent = true;
  //if ( leptonFlavour == 15 && zccbar ) keepEvent = true;
  
  
  if (keepEvent ) ikept++;
  taus= 2*ikept;
  
  return keepEvent;
  
}


void HZZ4LeptonsMCGenFilter::endJob(){
  
  std::cout << "Number of evetns filtered in the acceptance= " << ikept << std::endl;
  
}
