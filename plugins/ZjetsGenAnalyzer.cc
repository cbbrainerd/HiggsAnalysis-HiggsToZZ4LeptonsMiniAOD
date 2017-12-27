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
#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/ZjetsGenAnalyzer.h"

using namespace std;
using namespace edm;
using namespace reco;

ZjetsGenAnalyzer::ZjetsGenAnalyzer(const edm::ParameterSet &params) :
  sourceLabel(params.getParameter<edm::InputTag>("src"))
{
  edm::Service<TFileService> fs;
  histojpt         = fs->make<TH1F>("histojpt","histojpt",200, 0, 200.);
  histojptjm       = fs->make<TH1F>("histojptjm","histojptjm",200, 0, 200.);
  histoPt_j       = fs->make<TH1F>("histoPt_j","histoPt_j",200, 0, 200.);
  histoPt_dj       = fs->make<TH1F>("histoPt_dj","histoPt_dj",200, 0, 200.);
  histoPt_djbar    = fs->make<TH1F>("histoPt_djbar","histoPt_djbar",200, 0, 200.);
  
  histoMass_j     = fs->make<TH1F>("histoMass_j","histoMass_j",200, 0, 0.01);
  histophi_j       = fs->make<TH1F>("histophi_j","histophi_j",100, -4, 4.);
  histoeta_j       = fs->make<TH1F>("histoeta_j","histoeta_j",100, -8, 8.);
  histoMass_jbar   = fs->make<TH1F>("histoMass_jbar","histoMass_jbar",200, 0, 5.);
  histoPt_jbar     = fs->make<TH1F>("histoPt_jbar","histoPt_jbar",200, 0, 200.);
  histophi_jbar    = fs->make<TH1F>("histophi_jbar","histophi_jbar",100, -4, 4.);
  histoeta_jbar    = fs->make<TH1F>("histoeta_jbar","histoeta_jbar",100, -8, 8.);
  histoMass_l_j    = fs->make<TH1F>("histoMass_l_j","histoMass_l_j",200, 0, 5.);
  histoPt_l_j      = fs->make<TH1F>("histoPt_l_j","histoPt_l_j",100, 0, 100.);
  histophi_l_j     = fs->make<TH1F>("histophi_l_j","histophi_l_j",100, -4, 4.);
  histoeta_l_j     = fs->make<TH1F>("histoeta_l_j","histoeta_l_j",100, -8, 8.);
  histoPt_z        = fs->make<TH1F>("histoPt_z","histoPt_z",300, 0, 300.);
  histoeta_z       = fs->make<TH1F>("histoeta_z","histoeta_z",100, -8, 8.);
  histophi_z       = fs->make<TH1F>("histophi_z","histophi_z",100, -4, 4.);
  histoMass_z      = fs->make<TH1F>("histoMass_z","histoMass_z",200, 0, 200.);
  histoZmass       = fs->make<TH1F>("histoZmass","histoZmass",200, 0, 200.);
  histoMass_dileptons   = fs->make<TH1F>("histoMass_dileptons","histoMass_dileptons",300, 0, 200.);
  
  histoPt_l_z       = fs->make<TH1F>("histoPt_l_z","histoPt_l_z",200, 0, 200.);
  histoeta_l_z      = fs->make<TH1F>("histoeta_l_z","histoeta_l_z",100, -8, 8.);
  histophi_l_z      = fs->make<TH1F>("histophi_l_z","histophi_l_z",100, -4, 4.);
  histoMass_l_z     = fs->make<TH1F>("histoMass_l_z","histoMass_l_z",200, 0, 30.);
  histoPt_l_stable  = fs->make<TH1F>("histoPt_l_stable","histoPt_l_stable",100, 0, 100.);	
  histoeta_l_stable = fs->make<TH1F>("histoeta_l_stable","histoeta_l_stable",100, -8, 8.);	
  histophi_l_stable = fs->make<TH1F>("histophi_l_stable","histophi_l_stable",100, -4, 4.);
  histoPt_1         = fs->make<TH1F>("histoPt_1","histoPt_1",200, 0, 200.);
  histoPt_2         = fs->make<TH1F>("histoPt_2","histoPt_2",100, 0, 100.);
  histoPt_3         = fs->make<TH1F>("histoPt_3","histoPt_3",100, 0, 100.);
  histoPt_4         = fs->make<TH1F>("histoPt_4","histoPt_4",100, 0, 100.);
  
  histoMass_jjbarnc = fs->make<TH1F>("histoMass_jjbarnc","histoMass_jjbarnc",500, 0, 500.);
  histoMass_jjbar   = fs->make<TH1F>("histoMass_jjbar","histoMass_jjbar",500, 0, 500.);
  histoPt_jjbar     = fs->make<TH1F>("histoPt_jjbar","histoPt_jjbar",300, 0, 300.);
  histoeta_jjbar    = fs->make<TH1F>("histoeta_jjbar","histoeta_jjbar",100, -8, 8.);
  histophi_jjbar    = fs->make<TH1F>("histophi_jjbar","histophi_jjbar",100, -4, 4.);
  histoMass_Zjjbar  = fs->make<TH1F>("histoMass_Zjjbar","histoMass_Zjjbar",600, 0, 600.);
  histoPt_Zjjbar    = fs->make<TH1F>("histoPt_Zjjbar","histoPt_Zjjbar",200, 0, 200.);
  histoeta_Zjjbar   = fs->make<TH1F>("histoeta_Zjjbar","histoeta_Zjjbar",100, -8, 8.);	
  histophi_Zjjbar   = fs->make<TH1F>("histophi_Zjjbar","histophi_Zjjbar",100, -4, 4.);
  histoMass_fourl   = fs->make<TH1F>("histoMass_fourl","histoMass_fourl",500, 0, 500.);
  histoPt_fourl     = fs->make<TH1F>("histoPt_fourl","histoPt_fourl",300, 0, 300.);
  histoPt_D1        = fs->make<TH1F>("histoPt_D1","histoPt_D1",200, 0, 200.);	
  histoPt_D2        = fs->make<TH1F>("histoPt_D2","histoPt_D2",100, 0, 100.);	
  histoPt_D3        = fs->make<TH1F>("histoPt_D3","histoPt_D3",100, 0, 100.);
  histoPt_D4        = fs->make<TH1F>("histoPt_D4","histoPt_D4",100, 0, 100.);	
  histoeta_fourl    = fs->make<TH1F>("histoeta_fourl","histoeta_fourl",100, -8, 8.);	
  histoeta_D1       = fs->make<TH1F>("histoeta_D1","histoeta_D1",100, -8, 8.);
  histoeta_D2       = fs->make<TH1F>("histoeta_D2","histoeta_D2",100, -8, 8.);
  histoeta_D3       = fs->make<TH1F>("histoeta_D3","histoeta_D3",100, -8, 8.);
  histoeta_D4       = fs->make<TH1F>("histoeta_D4","histoeta_D4",100, -8, 8.);
  histophi_fourl    = fs->make<TH1F>("histophi_fourl","histophi_fourl",100, -4, 4.);
  histophi_D1       = fs->make<TH1F>("histophi_D1","histophi_D1",100, -4, 4.);
  histophi_D2       = fs->make<TH1F>("histophi_D2","histophi_D2",100, -4, 4.);
  histophi_D3       = fs->make<TH1F>("histophi_D3","histophi_D3",100, -4, 4.);
  histophi_D4       = fs->make<TH1F>("histophi_D4","histophi_D4",100, -4, 4.);
  histoMass_trilepton = fs->make<TH1F>("histoMass_trilepton","histoMass_trilepton",500, 0, 500.);	
  histophi_trilepton  = fs->make<TH1F>("histophi_trilepton","histophi_trilepton",100, -4, 4.);	
  histoeta_trilepton  = fs->make<TH1F>("histoeta_trilepton","histoeta_trilepton",100, -8, 8.);	
  histoPt_trilepton   = fs->make<TH1F>("histoPt_trilepton","histoPt_trilepton",200, 0, 200.);
  histonlept          = fs->make<TH1I>("histonlept","histonlept",20, 0, 20); 
}

ZjetsGenAnalyzer::~ZjetsGenAnalyzer()
{
}

void ZjetsGenAnalyzer::analyze(const edm::Event &event, const edm::EventSetup &es)
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
      std::cout << " Found l<-Z with mass= " << mcIter->mass() << " pdgId= " << mcIter->pdgId() << std::endl;
      histoPt_l_z->Fill(mcIter->pt()); 
      histoeta_l_z->Fill(mcIter->eta()); 
      histophi_l_z->Fill(mcIter->phi());
      histoMass_l_z->Fill(mcIter->mass());
    }
    
    // Light jets decay 
    
    if (((mcIter->pdgId()==1)|| (mcIter->pdgId()==2) || (mcIter->pdgId()==3)) && (mcIter->status()==3)) {
      if(((abs(mcIter->daughter(0)->pdgId())==1) || (abs(mcIter->daughter(0)->pdgId())==2) || (abs(mcIter->daughter(0)->pdgId())==3)) && (mcIter->daughter(0)->status()==2)){
	histoPt_dj->Fill(mcIter->daughter(0)->pt());
      }
      histoPt_j->Fill(mcIter->pt());
      histoeta_j->Fill(mcIter->eta());
      histophi_j->Fill(mcIter->phi());
      histoMass_j->Fill(mcIter->mass());	    
      std::cout << " Found light jet with pdgId= " << mcIter->pdgId() << " mass= " << mcIter->mass() << " and charge= " << mcIter->charge()  << " Mother= " << mcIter->mother(0)->pdgId()<< " and daughters= " << mcIter->daughter(0)->pdgId() << std::endl;
    }
    
    // anti particls of light jets
    if (((mcIter->pdgId()==-1)|| (mcIter->pdgId()==-2) || (mcIter->pdgId()==-3)) && (mcIter->status()==3)) { 
      if(((abs(mcIter->daughter(0)->pdgId())==1) || (abs(mcIter->daughter(0)->pdgId())==2) || (abs(mcIter->daughter(0)->pdgId())==3)) && (mcIter->daughter(0)->status()==2)){
	histoPt_djbar->Fill(mcIter->daughter(0)->pt());
      }
      histoPt_jbar->Fill(mcIter->pt()); 
      histoeta_jbar->Fill(mcIter->eta());
      histophi_jbar->Fill(mcIter->phi());	    
      histoMass_jbar->Fill(mcIter->mass());
      std::cout << " Found anti light jet with pdgId= " << mcIter->pdgId() << " mass= " << mcIter->mass() << " and charge= " << mcIter->charge() << " Mother= " << mcIter->mother(0)->pdgId() << " and daughters= " << mcIter->daughter(0)->pdgId() << std::endl;
    }
    
    // find  lepton from light jets decay, if any 
    if (((abs(mcIter->pdgId())==11) || (abs(mcIter->pdgId())==13)) && 
	(((abs(mcIter->mother(0)->pdgId())==1) || (abs(mcIter->mother(0)->pdgId())==2) || (abs(mcIter->mother(0)->pdgId())==3)) && 
	 ((mcIter->status()==2) || (mcIter->status()==3)))) {
      cout << " This light jet has n daughters= " << mcIter->mother(0)->numberOfDaughters() << endl;
      for (unsigned int d=0;d<mcIter->mother(0)->numberOfDaughters(); d++){
	cout << mcIter->mother(0)->daughter(d)->pdgId() << endl;
      }
      std::cout << " Found l<- j with mass= " << mcIter->mass() << " pdgId= " << mcIter->pdgId() << " status= " << mcIter->status() << std::endl;
      histoPt_l_j->Fill(mcIter->pt());
      histoeta_l_j->Fill(mcIter->eta());
      histophi_l_j->Fill(mcIter->phi());
      histoMass_l_j->Fill(mcIter->mass());
    } 
    
    //  stable leptons in the event 
    if (((abs(mcIter->pdgId())==11) || (abs(mcIter->pdgId())==13)) && (mcIter->status()==1)){
      histoPt_l_stable->Fill(mcIter->pt()); 
      histoeta_l_stable->Fill(mcIter->eta()); 
      histophi_l_stable->Fill(mcIter->phi());
      leptonpt.push_back(mcIter->pt());
      
      std::cout << " \n Found lepton with Id= " << mcIter->pdgId() << "  in the final state with mother= " ;
      // if (mcIter->numberOfMothers()>0 &&  mcIter->mother(0)->status()!=3){
      if (mcIter->numberOfMothers()>0){
	std::cout << getParticleName(int(mcIter->mother(0)->pdgId())) << " << "; 
	if (mcIter->mother(0)->status()==3) continue;
	if (mcIter->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->status()!=3 ){
	  std::cout << getParticleName(int(mcIter->mother(0)->mother(0)->pdgId())) << " << ";
	  if (mcIter->mother(0)->mother(0)->status()==3) continue;
	  if (mcIter->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->status()!=3 ){
	    std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
	    if (mcIter->mother(0)->mother(0)->mother(0)->status()==3) continue;
	    if (mcIter->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
	      std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
	      if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
	      if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
		if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		  std::cout << getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << " << std::endl;
		  if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		}
	      }
	    }
	  }
	}
      }      
    }     
  }
  
  // light jet (lj) candidates
  edm::Handle<edm::View<Candidate> > ljCandidates;
  event.getByLabel("lj",ljCandidates );
  for (edm::View<Candidate>::const_iterator mcIter=ljCandidates->begin(); mcIter!=ljCandidates->end(); ++mcIter ){
    cout << "\n light jet found from ljCandidates having pdgId " << mcIter->pdgId() <<  " with status= " << mcIter->status() << endl;  
    histojpt->Fill(mcIter->pt());
    if (((abs(mcIter->mother(0)->pdgId())==1) || (abs(mcIter->mother(0)->pdgId())==2) || (abs(mcIter->mother(0)->pdgId())==3)) && (mcIter->mother(0)->status()==3)){
      histojptjm->Fill(mcIter->pt());
    }
  }
  
  // jjbar (light jet - anti light jet) without cleaning
  edm::Handle<edm::View<Candidate> > jjbarCand;
  event.getByLabel("jjCands",jjbarCand );
  if (jjbarCand.isValid()) {
    for (edm::View<Candidate>::const_iterator mcIter=jjbarCand->begin(); mcIter!=jjbarCand->end(); ++mcIter ) {
      std:: cout << " No CLEANING Charge j jbar= " << mcIter->daughter(0)->pdgId() << " " << mcIter->daughter(1)->pdgId() << " Mother j jbar= " << mcIter->daughter(0)->mother(0)->pdgId() << " " << mcIter->daughter(1)->mother(0)->pdgId() << " No CLEANING jjbar  Mass= " << mcIter->mass() << std::endl;
      histoMass_jjbarnc->Fill(mcIter->mass());
    }
  }
  
  // jjbar candidates cleaned 
  edm::Handle<edm::View<Candidate> > jjbarCandidates;
  event.getByLabel("jjCandscleaned",jjbarCandidates );
  if (jjbarCandidates.isValid()) { 
    for (edm::View<Candidate>::const_iterator mcIter=jjbarCandidates->begin(); mcIter!=jjbarCandidates->end(); ++mcIter ) {
      std:: cout << " CLEANED  Charge j jbar= " << mcIter->daughter(0)->pdgId() << " &  " << mcIter->daughter(1)->pdgId() << " Mother j jbar= " << mcIter->daughter(0)->mother(0)->pdgId() << " " << mcIter->daughter(1)->mother(0)->pdgId() << " CLEANED jjbar Mass= " << mcIter->mass() << std::endl;
      
      histoMass_jjbar->Fill(mcIter->mass()); 
      histoPt_jjbar->Fill(mcIter->pt()); 
      histoeta_jjbar->Fill(mcIter->eta()); 
      histophi_jjbar->Fill(mcIter->phi()); 
    } 
  }
  
  // just Z
  edm::Handle<edm::View<Candidate> > Z;
  event.getByLabel("Z",Z);
  for (edm::View<Candidate>::const_iterator mcIter=Z->begin(); mcIter!=Z->end(); ++mcIter){
    histoZmass->Fill(mcIter->mass());
  }
  
  // get Zjjbarcandidates (Z+ light jets)
  edm::Handle<edm::View<Candidate> >  ZjjbarCandidates;
  event.getByLabel("lljjCands", ZjjbarCandidates);
  if (ZjjbarCandidates.isValid()) {
    for (edm::View<Candidate>::const_iterator mcIter=ZjjbarCandidates->begin(); mcIter!=ZjjbarCandidates->end(); ++mcIter ) {
      cout << " Zjjbar Mass= " << mcIter->mass()
	   << " pdgIds= " 
	   << mcIter->daughter(0)->pdgId() << " "
	   << mcIter->daughter(1)->daughter(0)->pdgId() << " " 
	   << mcIter->daughter(1)->daughter(1)->pdgId() << " " 
	   << endl;
      histoMass_Zjjbar->Fill(mcIter->mass()); 
      histoPt_Zjjbar->Fill(mcIter->pt()); 
      histoeta_Zjjbar->Fill(mcIter->eta()); 
      histophi_Zjjbar->Fill(mcIter->phi());
    }
  }
  
  // dileptons from Z candidates (stable leptons mu and e) 
  edm::Handle<edm::View<Candidate> > dileptonzCandidate;
  event.getByLabel("dileptons",dileptonzCandidate );
  if (dileptonzCandidate.isValid()) {
    for (edm::View<Candidate>::const_iterator mcIter=dileptonzCandidate->begin(); mcIter!=dileptonzCandidate->end(); ++mcIter ){
      if (((abs(mcIter->daughter(0)->mother(0)->mother(0)->pdgId())==23) && (abs(mcIter->daughter(1)->mother(0)->mother(0)->pdgId())==23)) && (abs(mcIter->daughter(0)->pdgId()) ==  abs(mcIter->daughter(1)->pdgId()))) {
	histoMass_dileptons->Fill(mcIter->mass());
	cout << "Mass of dileptons=" << mcIter->mass() << " and charge= " << mcIter->charge() << " and daughters= " << mcIter->daughter(0)->pdgId() << " " << mcIter->daughter(1)->pdgId() << endl;
      }
    }
  }
  
  // sort lepton pt and apply skimm if needed !
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
  
  //  trileptons cleaned candidates, no assignment of possible mothers because heavy flavors has been already filtered 
  edm::Handle<edm::View<Candidate> >  trileptonCandidates;
  event.getByLabel("trileptonscleaned", trileptonCandidates);
  if (trileptonCandidates.isValid()) { 
    std::cout <<"\n size of trileptonscleaned ===========================   " <<trileptonCandidates->size()<< std::endl; 
    for (edm::View<Candidate>::const_iterator mcIter=trileptonCandidates->begin(); mcIter!=trileptonCandidates->end(); ++mcIter ) {
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
  
  // get clean 4l candidates
  edm::Handle<edm::View<Candidate> >  fourlCandidates;
  event.getByLabel("fourleptonscleaned", fourlCandidates);
  if (fourlCandidates.isValid()) { 
    std::cout <<"\n size of fourleptonscleaned ===========================   " <<fourlCandidates->size()<< std::endl; 
    for (edm::View<Candidate>::const_iterator mcIter=fourlCandidates->begin(); mcIter!=fourlCandidates->end(); ++mcIter) {
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

std::string ZjetsGenAnalyzer::getParticleName(int id) const
{
  const ParticleData * pd = pdt_->particle( id );
  if (!pd) {
    std::ostringstream ss;
    ss << "P" << id;
    return ss.str();
  } else
    return pd->name();
}

   
