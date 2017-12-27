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
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/ZccGenAnalyzer.h"

using namespace std;
using namespace edm;
using namespace reco;

ZccGenAnalyzer::ZccGenAnalyzer(const edm::ParameterSet &params) :
  sourceLabel(params.getParameter<edm::InputTag>("src"))
{
  edm::Service<TFileService> fs;
  histocpt         = fs->make<TH1F>("histocpt","histocpt",200, 0, 200.);
  histocptcm       = fs->make<TH1F>("histocptcm","histocptcm",200, 0, 200.);
  histoPt_c        = fs->make<TH1F>("histoPt_c","histoPt_c",200, 0, 200.);
  histoPt_dc       = fs->make<TH1F>("histoPt_dc","histoPt_dc",200, 0, 200.);
  histoPt_dcbar    = fs->make<TH1F>("histoPt_dcbar","histoPt_dcbar",200, 0, 200.);
  
  histoMass_c      = fs->make<TH1F>("histoMass_c","histoMass_c",200, 0, 3.);
  histophi_c       = fs->make<TH1F>("histophi_c","histophi_c",100, -4, 4.);
  histoeta_c       = fs->make<TH1F>("histoeta_c","histoeta_c",100, -6, 6.);
  histoMass_cbar   = fs->make<TH1F>("histoMass_cbar","histoMass_cbar",200, 0, 5.);
  histoPt_cbar     = fs->make<TH1F>("histoPt_cbar","histoPt_cbar",200, 0, 200.);
  histophi_cbar    = fs->make<TH1F>("histophi_cbar","histophi_cbar",100, -4, 4.);
  histoeta_cbar    = fs->make<TH1F>("histoeta_cbar","histoeta_cbar",100, -6, 6.);
  histoMass_l_c    = fs->make<TH1F>("histoMass_l_c","histoMass_l_c",200, 0, 5.);
  histoPt_l_c      = fs->make<TH1F>("histoPt_l_c","histoPt_l_c",100, 0, 100.);
  histophi_l_c     = fs->make<TH1F>("histophi_l_c","histophi_l_c",100, -4, 4.);
  histoeta_l_c     = fs->make<TH1F>("histoeta_l_c","histoeta_l_c",100, -7, 7.);
  histoPt_z        = fs->make<TH1F>("histoPt_z","histoPt_z",300, 0, 300.);
  histoeta_z       = fs->make<TH1F>("histoeta_z","histoeta_z",100, -7, 7.);
  histophi_z       = fs->make<TH1F>("histophi_z","histophi_z",100, -4, 4.);
  histoMass_z      = fs->make<TH1F>("histoMass_z","histoMass_z",200, 0, 200.);
  histoZmass       = fs->make<TH1F>("histoZmass","histoZmass",200, 0, 200.);
  histoMass_zdileptons   = fs->make<TH1F>("histoMass_zdileptons","histoMass_zdileptons",300, 0, 200.);
  
  histoPt_l_z       = fs->make<TH1F>("histoPt_l_z","histoPt_l_z",200, 0, 200.);
  histoeta_l_z      = fs->make<TH1F>("histoeta_l_z","histoeta_l_z",100, -7, 7.);
  histophi_l_z      = fs->make<TH1F>("histophi_l_z","histophi_l_z",100, -4, 4.);
  histoMass_l_z     = fs->make<TH1F>("histoMass_l_z","histoMass_l_z",200, 0, 30.);
  histoPt_l_stable  = fs->make<TH1F>("histoPt_l_stable","histoPt_l_stable",100, 0, 100.);	
  histoeta_l_stable = fs->make<TH1F>("histoeta_l_stable","histoeta_l_stable",100, -8, 8.);	
  histophi_l_stable = fs->make<TH1F>("histophi_l_stable","histophi_l_stable",100, -4, 4.);
  histoPt_1         = fs->make<TH1F>("histoPt_1","histoPt_1",200, 0, 200.);
  histoPt_2         = fs->make<TH1F>("histoPt_2","histoPt_2",100, 0, 100.);
  histoPt_3         = fs->make<TH1F>("histoPt_3","histoPt_3",100, 0, 100.);
  histoPt_4         = fs->make<TH1F>("histoPt_4","histoPt_4",100, 0, 100.);
  
  histoMass_ccbarnc = fs->make<TH1F>("histoMass_ccbarnc","histoMass_ccbarnc",500, 0, 500.);
  histoMass_ccbar   = fs->make<TH1F>("histoMass_ccbar","histoMass_ccbar",500, 0, 500.);
  histoPt_ccbar     = fs->make<TH1F>("histoPt_ccbar","histoPt_ccbar",300, 0, 300.);
  histoeta_ccbar    = fs->make<TH1F>("histoeta_ccbar","histoeta_ccbar",100, -7, 7.);
  histophi_ccbar    = fs->make<TH1F>("histophi_ccbar","histophi_ccbar",100, -4, 4.);
  histoMass_Zccbar  = fs->make<TH1F>("histoMass_Zccbar","histoMass_Zccbar",600, 0, 600.);
  histoPt_Zccbar    = fs->make<TH1F>("histoPt_Zccbar","histoPt_Zccbar",200, 0, 200.);
  histoeta_Zccbar   = fs->make<TH1F>("histoeta_Zccbar","histoeta_Zccbar",100, -8, 8.);	
  histophi_Zccbar   = fs->make<TH1F>("histophi_Zccbar","histophi_Zccbar",100, -4, 4.);
  histoMass_fourl   = fs->make<TH1F>("histoMass_fourl","histoMass_fourl",500, 0, 500.);
  histoPt_fourl     = fs->make<TH1F>("histoPt_fourl","histoPt_fourl",300, 0, 300.);
  histoPt_D1        = fs->make<TH1F>("histoPt_D1","histoPt_D1",200, 0, 200.);	
  histoPt_D2        = fs->make<TH1F>("histoPt_D2","histoPt_D2",100, 0, 100.);	
  histoPt_D3        = fs->make<TH1F>("histoPt_D3","histoPt_D3",100, 0, 100.);
  histoPt_D4        = fs->make<TH1F>("histoPt_D4","histoPt_D4",100, 0, 100.);	
  histoeta_fourl    = fs->make<TH1F>("histoeta_fourl","histoeta_fourl",100, -7, 7.);	
  histoeta_D1       = fs->make<TH1F>("histoeta_D1","histoeta_D1",100, -7, 7.);
  histoeta_D2       = fs->make<TH1F>("histoeta_D2","histoeta_D2",100, -7, 7.);
  histoeta_D3       = fs->make<TH1F>("histoeta_D3","histoeta_D3",100, -7, 7.);
  histoeta_D4       = fs->make<TH1F>("histoeta_D4","histoeta_D4",100, -7, 7.);
  histophi_fourl    = fs->make<TH1F>("histophi_fourl","histophi_fourl",100, -4, 4.);
  histophi_D1       = fs->make<TH1F>("histophi_D1","histophi_D1",100, -4, 4.);
  histophi_D2       = fs->make<TH1F>("histophi_D2","histophi_D2",100, -4, 4.);
  histophi_D3       = fs->make<TH1F>("histophi_D3","histophi_D3",100, -4, 4.);
  histophi_D4       = fs->make<TH1F>("histophi_D4","histophi_D4",100, -4, 4.);
  histoMass_trilepton = fs->make<TH1F>("histoMass_trilepton","histoMass_trilepton",500, 0, 500.);	
  histophi_trilepton  = fs->make<TH1F>("histophi_trilepton","histophi_trilepton",100, -4, 4.);	
  histoeta_trilepton  = fs->make<TH1F>("histoeta_trilepton","histoeta_trilepton",100, -7, 7.);	
  histoPt_trilepton   = fs->make<TH1F>("histoPt_trilepton","histoPt_trilepton",200, 0, 200.);
  histonlept          = fs->make<TH1I>("histonlept","histonlept",20, 0, 20); 
}

ZccGenAnalyzer::~ZccGenAnalyzer()
{
}

void ZccGenAnalyzer::analyze(const edm::Event &event, const edm::EventSetup &es)
{
  // get gen particle candidates 
  
  edm::Handle<reco::GenParticleCollection> genCandidates;	    
  event.getByLabel(sourceLabel, genCandidates);
  es.getData( pdt_ );
  leptonpt.clear();
  
  for ( GenParticleCollection::const_iterator mcIter=genCandidates->begin(); mcIter!=genCandidates->end(); ++mcIter ) {
    
    // Find Z and Z->lepton
    
    if (abs(mcIter->pdgId()==23) && (mcIter->status()==3)){
      std::cout << " Found Z with mass= " << mcIter->mass() << " and charge= " << mcIter->charge() << " and daughters= " << mcIter->daughter(0)->pdgId() << " " << mcIter->daughter(1)->pdgId() << std::endl;
      histoPt_z->Fill(mcIter->pt());     
      histoeta_z->Fill(mcIter->eta());
      histophi_z->Fill(mcIter->phi());
      histoMass_z->Fill(mcIter->mass());
    }	
    if ((( abs(mcIter->pdgId())==11) || (abs(mcIter->pdgId())==13)) && 
	(abs(mcIter->mother(0)->pdgId())==23 && (mcIter->status()==3))) {	    
      std::cout << " Found l<-Z with mass= " << mcIter->mass() << " pdgId= " << mcIter->pdgId() << "status= " << mcIter->status() <<  std::endl;
      histoPt_l_z->Fill(mcIter->pt()); 
      histoeta_l_z->Fill(mcIter->eta()); 
      histophi_l_z->Fill(mcIter->phi());
      histoMass_l_z->Fill(mcIter->mass());
    }
    
    // leptons from direct c decay, generally we don't get lepton  directly from c and we consider leptons from decay chain of c but not here. Read MC info of events for more detail 
    
    if (mcIter->pdgId()==4 && (mcIter->status()==3)) {
      if(abs(mcIter->daughter(0)->pdgId())==4 && (mcIter->daughter(0)->status()==2)){
	histoPt_dc->Fill(mcIter->daughter(0)->pt());
      }
      histoPt_c->Fill(mcIter->pt());
      histoeta_c->Fill(mcIter->eta());
      histophi_c->Fill(mcIter->phi());
      histoMass_c->Fill(mcIter->mass());	    
      std::cout << " Found c with pdgId= " << mcIter->pdgId() << " mass= " << mcIter->mass() << " and charge= " << mcIter->charge()  << " Mother= " << mcIter->mother(0)->pdgId()<< " and daughters= " << mcIter->daughter(0)->pdgId() << std::endl;
    }
    
    if (mcIter->pdgId()==-4 && (mcIter->status()==3)) {
      if(abs(mcIter->daughter(0)->pdgId())==4 && (mcIter->daughter(0)->status()==2)){
	histoPt_dcbar->Fill(mcIter->daughter(0)->pt());
      }
      histoPt_cbar->Fill(mcIter->pt()); 
      histoeta_cbar->Fill(mcIter->eta());
      histophi_cbar->Fill(mcIter->phi());	    
      histoMass_cbar->Fill(mcIter->mass());
      std::cout << " Found c_bar with pdgId= " << mcIter->pdgId() << " mass= " << mcIter->mass() << " and charge= " << mcIter->charge() << " Mother= " << mcIter->mother(0)->pdgId() << " and daughters= " << mcIter->daughter(0)->pdgId() << std::endl;
    }
    
    // Lets find  lepton from c decay, if any 
    if (((abs(mcIter->pdgId())==11) || (abs(mcIter->pdgId())==13)) &&
	((abs(mcIter->mother(0)->pdgId())==4) &&  ((mcIter->status()==2) || (mcIter->status()==3)))){
      cout << " This c has n daughters= " << mcIter->mother(0)->numberOfDaughters() << endl;
      for (unsigned int d=0;d<mcIter->mother(0)->numberOfDaughters(); d++){
	cout << mcIter->mother(0)->daughter(d)->pdgId() << endl;
      }
      std::cout << " Found l<-c with mass= " << mcIter->mass() << " pdgId= " << mcIter->pdgId() << " status= " << mcIter->status() << std::endl;
      histoPt_l_c->Fill(mcIter->pt());
      histoeta_l_c->Fill(mcIter->eta());
      histophi_l_c->Fill(mcIter->phi());
      histoMass_l_c->Fill(mcIter->mass());
    } 
    
    // stable leptons
    if ((( abs(mcIter->pdgId())==11) || (abs(mcIter->pdgId())==13)) && (mcIter->status()==1)){
      histoPt_l_stable->Fill(mcIter->pt()); 
      histoeta_l_stable->Fill(mcIter->eta()); 
      histophi_l_stable->Fill(mcIter->phi());
      leptonpt.push_back(mcIter->pt());
      
      std::cout << " \n Found lepton with Id= " << mcIter->pdgId() << "  in the final state with mother= " ;
      // if (mcIter->numberOfMothers()>0 &&  mcIter->mother(0)->status()!=3){
      if (mcIter->numberOfMothers()>0){
	std::cout <<  (mcIter->mother(0)->pdgId()) << " << "; 
	if (mcIter->mother(0)->status()==3) continue;
	if (mcIter->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->status()!=3 ){
	  std::cout <<  (mcIter->mother(0)->mother(0)->pdgId()) << " << ";
	  if (mcIter->mother(0)->mother(0)->status()==3) continue;
	  if (mcIter->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->status()!=3 ){
	    std::cout<<  (mcIter->mother(0)->mother(0)->mother(0)->pdgId()) << " << ";
	    if (mcIter->mother(0)->mother(0)->mother(0)->status()==3) continue;
	    if (mcIter->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
	      std::cout<<  (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->pdgId()) << " << ";
	      if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
	      if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		std::cout<<  (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->pdgId()) << " << ";
		if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		  std::cout <<  (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->pdgId()) << " << " << std::endl;
		  if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  // c candidates
  edm::Handle<edm::View<Candidate> > cCandidates;
  event.getByLabel("c",cCandidates );
  for (edm::View<Candidate>::const_iterator mcIter=cCandidates->begin(); mcIter!=cCandidates->end(); ++mcIter ){
    cout << "\n c found from cCandidates having pdgId " << mcIter->pdgId() <<  " with status= " << mcIter->status() << endl;  
    histocpt->Fill(mcIter->pt());
    if ((abs(mcIter->mother(0)->pdgId())==4) && (mcIter->mother(0)->status()==3)){
      histocptcm->Fill(mcIter->pt());
    }
  }
  
  // ccbar without cleaning
  edm::Handle<edm::View<Candidate> > ccbarCand;
  event.getByLabel("ccCands",ccbarCand );
  if (ccbarCand.isValid()) {
    for (edm::View<Candidate>::const_iterator mcIter=ccbarCand->begin(); mcIter!=ccbarCand->end(); ++mcIter ) {
      std:: cout << " No CLEANING Charge c cbar= " << mcIter->daughter(0)->pdgId() << " " << mcIter->daughter(1)->pdgId() << " Mother c cbar= " << mcIter->daughter(0)->mother(0)->pdgId() << " " << mcIter->daughter(1)->mother(0)->pdgId() << " No CLEANING ccbar  Mass= " << mcIter->mass() << std::endl;
      histoMass_ccbarnc->Fill(mcIter->mass());
    }
  }
  
  // ccbar candidates cleaned 
  edm::Handle<edm::View<Candidate> > ccbarCandidates;
  event.getByLabel("ccCandscleaned",ccbarCandidates );
  if (ccbarCandidates.isValid()) { 
    for (edm::View<Candidate>::const_iterator mcIter=ccbarCandidates->begin(); mcIter!=ccbarCandidates->end(); ++mcIter ) {
      std:: cout << " CLEANED  Charge c cbar= " << mcIter->daughter(0)->pdgId() << " &  " << mcIter->daughter(1)->pdgId() << " Mother c cbar= " << mcIter->daughter(0)->mother(0)->pdgId() << " " << mcIter->daughter(1)->mother(0)->pdgId() << " CLEANED ccbar Mass= " << mcIter->mass() << std::endl;
      
      histoMass_ccbar->Fill(mcIter->mass()); 
      histoPt_ccbar->Fill(mcIter->pt()); 
      histoeta_ccbar->Fill(mcIter->eta()); 
      histophi_ccbar->Fill(mcIter->phi()); 
    } 
  }
  
  // just Z
  edm::Handle<edm::View<Candidate> > Z;
  event.getByLabel("Z",Z);
  for (edm::View<Candidate>::const_iterator mcIter=Z->begin(); mcIter!=Z->end(); ++mcIter){
    histoZmass->Fill(mcIter->mass());
  }
  
  // get Zccbarcandidates
  edm::Handle<edm::View<Candidate> >  ZccbarCandidates;
  event.getByLabel("llccCands", ZccbarCandidates);
  if (ZccbarCandidates.isValid()) {
    for (edm::View<Candidate>::const_iterator mcIter=ZccbarCandidates->begin(); mcIter!=ZccbarCandidates->end(); ++mcIter ) {
      cout << " Zccbar Mass= " << mcIter->mass()
	   << " pdgIds= " 
	   << mcIter->daughter(0)->pdgId() << " "
	   << mcIter->daughter(1)->daughter(0)->pdgId() << " " 
	   << mcIter->daughter(1)->daughter(1)->pdgId() << " " 
	   << endl;
      histoMass_Zccbar->Fill(mcIter->mass()); 
      histoPt_Zccbar->Fill(mcIter->pt()); 
      histoeta_Zccbar->Fill(mcIter->eta()); 
      histophi_Zccbar->Fill(mcIter->phi());
    }
  }
  
  // dileptons from z candidates (stable leptons mu and e) 
  edm::Handle<edm::View<Candidate> > dileptonzCandidate;
  event.getByLabel("dileptons",dileptonzCandidate );
  if (dileptonzCandidate.isValid()) {
    for (edm::View<Candidate>::const_iterator mcIter=dileptonzCandidate->begin(); mcIter!=dileptonzCandidate->end(); ++mcIter ){
      if ((abs(mcIter->daughter(0)->mother(0)->mother(0)->pdgId())==23) && (abs(mcIter->daughter(1)->mother(0)->mother(0)->pdgId())==23) && (abs(mcIter->daughter(0)->pdgId()) ==  (abs(mcIter->daughter(1)->pdgId())))) {
	histoMass_zdileptons->Fill(mcIter->mass());
	cout << "Mass of dileptons=" << mcIter->mass() << " and charge= " << mcIter->charge() << " and daughters= " << mcIter->daughter(0)->pdgId() << " " << mcIter->daughter(1)->pdgId() << endl;
      }
    }
  }
  
  // sort lepton pt and apply skimming if needed !
  sort(leptonpt.rbegin(),leptonpt.rend());
  std::cout  << "\n Number of stable leptons= " << leptonpt.size() << std::endl;
  histonlept->Fill(leptonpt.size());
  if (leptonpt.size()>=4){
    std::cout << "\n Sorted vector in decreasing order= " << leptonpt.at(0) << " " << leptonpt.at(1) << " " << leptonpt.at(2) << " " << leptonpt.at(3) << std::endl;
    histoPt_1->Fill(leptonpt.at(0));
    histoPt_2->Fill(leptonpt.at(1));
    histoPt_3->Fill(leptonpt.at(2));
    histoPt_4->Fill(leptonpt.at(3));
  }
  
  //  get 3l by assigning the possible mothers to daughters of 3l,these possible mothers are end product of decay chains started from c || cbar 
  edm::Handle<edm::View<Candidate> >  trileptonCandidates;
  event.getByLabel("trileptonscleaned", trileptonCandidates);
  if (trileptonCandidates.isValid()) { 
    std::cout <<"\n size of trileptonscleaned ===========================   " <<trileptonCandidates->size()<< std::endl; 
    for (edm::View<Candidate>::const_iterator mcIter=trileptonCandidates->begin(); mcIter!=trileptonCandidates->end(); ++mcIter ) {
      if (
	  (mcIter->daughter(0)->mother(0)->mother(0)->pdgId()==423  ||
	   mcIter->daughter(0)->mother(0)->mother(0)->pdgId()==23  ||
	   mcIter->daughter(0)->mother(0)->mother(0)->pdgId()==413 ||
	   mcIter->daughter(0)->mother(0)->mother(0)->pdgId()==10421 ||
	   mcIter->daughter(0)->mother(0)->mother(0)->pdgId()==10411 ||
	   (abs(mcIter->daughter(0)->mother(0)->mother(0)->pdgId())==411) || 
	   (abs(mcIter->daughter(0)->mother(0)->mother(0)->pdgId())==421) || 
	   mcIter->daughter(0)->mother(0)->mother(0)->pdgId()==521 ||
	   mcIter->daughter(0)->mother(0)->mother(0)->pdgId()==10511 ||
	   mcIter->daughter(0)->mother(0)->mother(0)->pdgId()==10521 ||
	   mcIter->daughter(0)->mother(0)->mother(0)->pdgId()==513 ||
	   mcIter->daughter(0)->mother(0)->mother(0)->pdgId()==523 ||
	   mcIter->daughter(0)->mother(0)->mother(0)->pdgId()==511 
	   ) &&  
	  (mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==423 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==23  ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==413 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==10421 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==10411 ||
	   (abs(mcIter->daughter(1)->mother(0)->mother(0)->pdgId())==411) || 
	   (abs(mcIter->daughter(1)->mother(0)->mother(0)->pdgId())==421) || 
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==521 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==10511 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==10521 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==513 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==523 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==511
	   ) &&
	  (mcIter->daughter(2)->mother(0)->mother(0)->pdgId()==423 ||
	   mcIter->daughter(2)->mother(0)->mother(0)->pdgId()==23  ||
	   mcIter->daughter(2)->mother(0)->mother(0)->pdgId()==413 ||
	   mcIter->daughter(2)->mother(0)->mother(0)->pdgId()==10421 ||
	   mcIter->daughter(2)->mother(0)->mother(0)->pdgId()==10411 ||
	   (abs(mcIter->daughter(2)->mother(0)->mother(0)->pdgId())==411) || 
	   (abs(mcIter->daughter(2)->mother(0)->mother(0)->pdgId())==421) || 
	   mcIter->daughter(2)->mother(0)->mother(0)->pdgId()==521 ||
	   mcIter->daughter(2)->mother(0)->mother(0)->pdgId()==10511 ||
	   mcIter->daughter(2)->mother(0)->mother(0)->pdgId()==10521 ||
	   mcIter->daughter(2)->mother(0)->mother(0)->pdgId()==513 ||
	   mcIter->daughter(2)->mother(0)->mother(0)->pdgId()==523 ||
	   mcIter->daughter(2)->mother(0)->mother(0)->pdgId()==511))
	{
	  std::cout << "3l Mass= " << mcIter->mass() << " "
		    << "pdgId= " 
		    << mcIter->daughter(0)->pdgId() << " " 
		    << mcIter->daughter(1)->pdgId() << " "
		    << mcIter->daughter(2)->pdgId() << " " 
		    <<std::endl; 
	  histoMass_trilepton->Fill(mcIter->mass()); 
	  histoPt_trilepton->Fill(mcIter->pt()); 
	  histoeta_trilepton->Fill(mcIter->eta()); 
	  histophi_trilepton->Fill(mcIter->phi());
	}
    }
  }
  
  // get clean 4l candidates by assigning the possible mothers to daughters of 4l 
  edm::Handle<edm::View<Candidate> >  fourlCandidates;
  event.getByLabel("fourleptonscleaned", fourlCandidates);
  if (fourlCandidates.isValid()) { 
    std::cout <<"\n size of fourleptonscleaned =========================== " <<fourlCandidates->size()<< std::endl; 
    for (edm::View<Candidate>::const_iterator mcIter=fourlCandidates->begin(); mcIter!=fourlCandidates->end(); ++mcIter) {
      if ( 
	  (mcIter->daughter(0)->daughter(0)->mother(0)->mother(0)->pdgId()==423  ||
	   mcIter->daughter(0)->daughter(0)->mother(0)->mother(0)->pdgId()==23  ||
	   mcIter->daughter(0)->daughter(0)->mother(0)->mother(0)->pdgId()==413  ||
	   mcIter->daughter(0)->daughter(0)->mother(0)->mother(0)->pdgId()==10421  ||
	   mcIter->daughter(0)->daughter(0)->mother(0)->mother(0)->pdgId()==10411  ||
	   (abs(mcIter->daughter(0)->daughter(0)->mother(0)->mother(0)->pdgId())==411)  ||
	   (abs(mcIter->daughter(0)->daughter(0)->mother(0)->mother(0)->pdgId())==421)  ||
	   mcIter->daughter(0)->daughter(0)->mother(0)->mother(0)->pdgId()==521  ||
	   mcIter->daughter(0)->daughter(0)->mother(0)->mother(0)->pdgId()==10511  ||
	   mcIter->daughter(0)->daughter(0)->mother(0)->mother(0)->pdgId()==10521  ||
	   mcIter->daughter(0)->daughter(0)->mother(0)->mother(0)->pdgId()==513  ||
	   mcIter->daughter(0)->daughter(0)->mother(0)->mother(0)->pdgId()==523  ||
	   mcIter->daughter(0)->daughter(0)->mother(0)->mother(0)->pdgId()==511
	   )
	  &&
	  (mcIter->daughter(0)->daughter(1)->mother(0)->mother(0)->pdgId()==423  ||
	   mcIter->daughter(0)->daughter(1)->mother(0)->mother(0)->pdgId()==23  ||
	   mcIter->daughter(0)->daughter(1)->mother(0)->mother(0)->pdgId()==413  ||
	   mcIter->daughter(0)->daughter(1)->mother(0)->mother(0)->pdgId()==10421  ||
	   mcIter->daughter(0)->daughter(1)->mother(0)->mother(0)->pdgId()==10411  ||
	   (abs(mcIter->daughter(0)->daughter(1)->mother(0)->mother(0)->pdgId())==411)  ||
	   (abs(mcIter->daughter(0)->daughter(1)->mother(0)->mother(0)->pdgId())==421)  ||
	   mcIter->daughter(0)->daughter(1)->mother(0)->mother(0)->pdgId()==521  ||
	   mcIter->daughter(0)->daughter(1)->mother(0)->mother(0)->pdgId()==10511  ||
	   mcIter->daughter(0)->daughter(1)->mother(0)->mother(0)->pdgId()==10521  ||
	   mcIter->daughter(0)->daughter(1)->mother(0)->mother(0)->pdgId()==513  ||
	   mcIter->daughter(0)->daughter(1)->mother(0)->mother(0)->pdgId()==523  ||
	   mcIter->daughter(0)->daughter(1)->mother(0)->mother(0)->pdgId()==511
	   ) &&
	  (mcIter->daughter(0)->daughter(2)->mother(0)->mother(0)->pdgId()==423  ||
	   mcIter->daughter(0)->daughter(2)->mother(0)->mother(0)->pdgId()==23  ||
	   mcIter->daughter(0)->daughter(2)->mother(0)->mother(0)->pdgId()==413  ||
	   mcIter->daughter(0)->daughter(2)->mother(0)->mother(0)->pdgId()==10421  ||
	   mcIter->daughter(0)->daughter(2)->mother(0)->mother(0)->pdgId()==10411  ||
	   (abs(mcIter->daughter(0)->daughter(2)->mother(0)->mother(0)->pdgId())==411)  ||
	   (abs(mcIter->daughter(0)->daughter(2)->mother(0)->mother(0)->pdgId())==421)  ||
	   mcIter->daughter(0)->daughter(2)->mother(0)->mother(0)->pdgId()==521  ||
	   mcIter->daughter(0)->daughter(2)->mother(0)->mother(0)->pdgId()==10511  ||
	   mcIter->daughter(0)->daughter(2)->mother(0)->mother(0)->pdgId()==10521  ||
	   mcIter->daughter(0)->daughter(2)->mother(0)->mother(0)->pdgId()==513  ||
	   mcIter->daughter(0)->daughter(2)->mother(0)->mother(0)->pdgId()==523  ||
	   mcIter->daughter(0)->daughter(2)->mother(0)->mother(0)->pdgId()==511
	   ) &&
	  (mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==423 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==23  ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==413 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==10421 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==10411 ||
	   (abs(mcIter->daughter(1)->mother(0)->mother(0)->pdgId())==411) || 
	   (abs(mcIter->daughter(1)->mother(0)->mother(0)->pdgId())==421) || 
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==521 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==10511 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==10521 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==513 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==523 ||
	   mcIter->daughter(1)->mother(0)->mother(0)->pdgId()==511
	   )
	  )
	{
	  std::cout << "4l Mass= " << mcIter->mass() << " "
		    << "pdgId= " 
		    << mcIter->daughter(0)->daughter(0)->pdgId() << " " 
		    << mcIter->daughter(0)->daughter(1)->pdgId() << " "
		    << mcIter->daughter(0)->daughter(2)->pdgId() << " " 
		    << mcIter->daughter(1)->pdgId() << " " 
		    <<"Pt= "
		    << mcIter->daughter(0)->daughter(0)->pt() << " " 
		    << mcIter->daughter(0)->daughter(1)->pt() << " "
		    << mcIter->daughter(0)->daughter(2)->pt() << " " 
		    << mcIter->daughter(1)->pt() << " "
		    <<"Phi= "
		    << mcIter->daughter(0)->daughter(0)->phi() << " " 
		    << mcIter->daughter(0)->daughter(1)->phi() << " "
		    << mcIter->daughter(0)->daughter(2)->phi() << " " 
		    << mcIter->daughter(1)->phi() << " " 
		    <<"eta= "
		    << mcIter->daughter(0)->daughter(0)->eta() << " " 
		    << mcIter->daughter(0)->daughter(1)->eta() << " "
		    << mcIter->daughter(0)->daughter(2)->eta() << " " 
		    << mcIter->daughter(1)->eta() << " "
		    <<std::endl; 
	  
	  histoMass_fourl->Fill(mcIter->mass()); 
	  histoPt_fourl->Fill(mcIter->pt());
	  histoPt_D1->Fill(mcIter->daughter(0)->daughter(0)->pt());
	  histoPt_D2->Fill(mcIter->daughter(0)->daughter(1)->pt());
	  histoPt_D3->Fill(mcIter->daughter(0)->daughter(2)->pt());
	  histoPt_D4->Fill(mcIter->daughter(1)->pt());
	  histoeta_D1->Fill(mcIter->daughter(0)->daughter(0)->eta());
	  histoeta_D2->Fill(mcIter->daughter(0)->daughter(1)->eta());
	  histoeta_D3->Fill(mcIter->daughter(0)->daughter(2)->eta());
	  histoeta_D4->Fill(mcIter->daughter(1)->eta());
	  histoeta_fourl->Fill(mcIter->eta()); 
	  histophi_D1->Fill(mcIter->daughter(0)->daughter(0)->phi());
	  histophi_D2->Fill(mcIter->daughter(0)->daughter(1)->phi());
	  histophi_D3->Fill(mcIter->daughter(0)->daughter(2)->phi());
	  histophi_D4->Fill(mcIter->daughter(1)->phi());
	  histophi_fourl->Fill(mcIter->phi());
	}
    }
  }
}

std::string ZccGenAnalyzer::getParticleName(int id) const
{
  const ParticleData * pd = pdt_->particle( id );
  if (!pd) {
    std::ostringstream ss;
    ss << "P" << id;
    return ss.str();
  } else
    return pd->name();
}


   
