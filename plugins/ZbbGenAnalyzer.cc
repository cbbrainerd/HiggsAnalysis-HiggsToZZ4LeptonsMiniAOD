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
#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/ZbbGenAnalyzer.h"

using namespace std;
using namespace edm;
using namespace reco;

ZbbGenAnalyzer::ZbbGenAnalyzer(const edm::ParameterSet &params) :
  sourceLabel(params.getParameter<edm::InputTag>("src"))
{
  edm::Service<TFileService> fs;
  histobpt         = fs->make<TH1F>("histobpt","histobpt",200, 0, 200.);
  histobptbm       = fs->make<TH1F>("histobptbm","histobptbm",200, 0, 200.);
  histoPt_b        = fs->make<TH1F>("histoPt_b","histoPt_b",200, 0, 200.);
  histoPt_db       = fs->make<TH1F>("histoPt_db","histoPt_db",200, 0, 200.);
  histoPt_dbbar    = fs->make<TH1F>("histoPt_dbbar","histoPt_dbbar",200, 0, 200.);
  
  histoMass_b      = fs->make<TH1F>("histoMass_b","histoMass_b",200, 0, 5.);
  histophi_b       = fs->make<TH1F>("histophi_b","histophi_b",100, -4, 4.);
  histoeta_b       = fs->make<TH1F>("histoeta_b","histoeta_b",100, -6, 6.);
  histoMass_bbar   = fs->make<TH1F>("histoMass_bbar","histoMass_bbar",200, 0, 5.);
  histoPt_bbar     = fs->make<TH1F>("histoPt_bbar","histoPt_bbar",200, 0, 200.);
  histophi_bbar    = fs->make<TH1F>("histophi_bbar","histophi_bbar",100, -4, 4.);
  histoeta_bbar    = fs->make<TH1F>("histoeta_bbar","histoeta_bbar",100, -6, 6.);
  histoMass_l_b    = fs->make<TH1F>("histoMass_l_b","histoMass_l_b",200, 0, 5.);
  histoPt_l_b      = fs->make<TH1F>("histoPt_l_b","histoPt_l_b",100, 0, 100.);
  histophi_l_b     = fs->make<TH1F>("histophi_l_b","histophi_l_b",100, -4, 4.);
  histoeta_l_b     = fs->make<TH1F>("histoeta_l_b","histoeta_l_b",100, -7, 7.);
  histoPt_z        = fs->make<TH1F>("histoPt_z","histoPt_z",300, 0, 300.);
  histoeta_z       = fs->make<TH1F>("histoeta_z","histoeta_z",100, -7, 7.);
  histophi_z       = fs->make<TH1F>("histophi_z","histophi_z",100, -4, 4.);
  histoMass_z      = fs->make<TH1F>("histoMass_z","histoMass_z",200, 0, 200.);
  histoZmass       = fs->make<TH1F>("histoZmass","histoZmass",200, 0, 200.);
  histoMass_dileptons   = fs->make<TH1F>("histoMass_dileptons","histoMass_dileptons",300, 0, 200.);
  
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
  
  histoMass_bbbarnc = fs->make<TH1F>("histoMass_bbbarnc","histoMass_bbbarnc",500, 0, 500.);
  histoMass_bbbar   = fs->make<TH1F>("histoMass_bbbar","histoMass_bbbar",500, 0, 500.);
  histoPt_bbbar     = fs->make<TH1F>("histoPt_bbbar","histoPt_bbbar",300, 0, 300.);
  histoeta_bbbar    = fs->make<TH1F>("histoeta_bbbar","histoeta_bbbar",100, -7, 7.);
  histophi_bbbar    = fs->make<TH1F>("histophi_bbbar","histophi_bbbar",100, -4, 4.);
  histoMass_Zbbbar  = fs->make<TH1F>("histoMass_Zbbbar","histoMass_Zbbbar",600, 0, 600.);
  histoPt_Zbbbar    = fs->make<TH1F>("histoPt_Zbbbar","histoPt_Zbbbar",200, 0, 200.);
  histoeta_Zbbbar   = fs->make<TH1F>("histoeta_Zbbbar","histoeta_Zbbbar",100, -8, 8.);	
  histophi_Zbbbar   = fs->make<TH1F>("histophi_Zbbbar","histophi_Zbbbar",100, -4, 4.);
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

ZbbGenAnalyzer::~ZbbGenAnalyzer()
{
}

void ZbbGenAnalyzer::analyze(const edm::Event &event, const edm::EventSetup &es)
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
    if (((abs(mcIter->pdgId())==11) || (abs(mcIter->pdgId())==13)) && 
	(abs(mcIter->mother(0)->pdgId())==23 && (mcIter->status()==3))) {	    
      std::cout << " Found l<-Z with mass= " << mcIter->mass() << " pdgId= " << mcIter->pdgId() << std::endl;
      histoPt_l_z->Fill(mcIter->pt()); 
      histoeta_l_z->Fill(mcIter->eta()); 
      histophi_l_z->Fill(mcIter->phi());
      histoMass_l_z->Fill(mcIter->mass());
    }
    
    // Try to get leptons from direct b decay, generally we don't get lepton  directly from b and we consider leptons from decay chain of b but not here. Read MC info of events for more detail 
    
    if (mcIter->pdgId()==5 && (mcIter->status()==3)) {
      if(abs(mcIter->daughter(0)->pdgId())==5 && (mcIter->daughter(0)->status()==2)){
	histoPt_db->Fill(mcIter->daughter(0)->pt());
      }
      histoPt_b->Fill(mcIter->pt());
      histoeta_b->Fill(mcIter->eta());
      histophi_b->Fill(mcIter->phi());
      histoMass_b->Fill(mcIter->mass());	    
      std::cout << " Found b with pdgId= " << mcIter->pdgId() << " mass= " << mcIter->mass() << " and charge= " << mcIter->charge()  << " Mother= " << mcIter->mother(0)->pdgId()<< " and daughters= " << mcIter->daughter(0)->pdgId() << std::endl;
    }
    
    if (mcIter->pdgId()==-5 && (mcIter->status()==3)) {
      if(abs(mcIter->daughter(0)->pdgId())==5 && (mcIter->daughter(0)->status()==2)){
	histoPt_dbbar->Fill(mcIter->daughter(0)->pt());
      }
      histoPt_bbar->Fill(mcIter->pt()); 
      histoeta_bbar->Fill(mcIter->eta());
      histophi_bbar->Fill(mcIter->phi());	    
      histoMass_bbar->Fill(mcIter->mass());
      std::cout << " Found b_bar with pdgId= " << mcIter->pdgId() << " mass= " << mcIter->mass() << " and charge= " << mcIter->charge() << " Mother= " << mcIter->mother(0)->pdgId() << " and daughters= " << mcIter->daughter(0)->pdgId() << std::endl;
    }
    
    if (((abs(mcIter->pdgId())==11) || (abs(mcIter->pdgId())==13)) && 
	(abs(mcIter->mother(0)->pdgId())==5 && ((mcIter->status()==2) || (mcIter->status()==3))))
      {
	cout << " This b has n daughters= " << mcIter->mother(0)->numberOfDaughters() << endl;
	for (unsigned int d=0;d<mcIter->mother(0)->numberOfDaughters(); d++){
	  cout << mcIter->mother(0)->daughter(d)->pdgId() << endl;
	}
	std::cout << " Found l<-b with mass= " << mcIter->mass() << " pdgId= " << mcIter->pdgId() << " status= " << mcIter->status() << std::endl;
	histoPt_l_b->Fill(mcIter->pt());
	histoeta_l_b->Fill(mcIter->eta());
	histophi_l_b->Fill(mcIter->phi());
	histoMass_l_b->Fill(mcIter->mass());
      } 
    
    //  stable leptons
    if ((( abs(mcIter->pdgId())==11) || (abs(mcIter->pdgId())==13)) && (mcIter->status()==1)){
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
  
  // b candidates: 
  edm::Handle<edm::View<Candidate> > bCandidates;
  event.getByLabel("b",bCandidates );
  for (edm::View<Candidate>::const_iterator mcIter=bCandidates->begin(); mcIter!=bCandidates->end(); ++mcIter ){
    cout << "\n b found from bcandidates having pdgId " << mcIter->pdgId() <<  " with status= " << mcIter->status() << endl;  
    histobpt->Fill(mcIter->pt());
    if ((abs(mcIter->mother(0)->pdgId())==5) && (mcIter->mother(0)->status()==3)){
      histobptbm->Fill(mcIter->pt());
    }
  }
  
  // bbbar without cleaning
  edm::Handle<edm::View<Candidate> > bbbarCand;
  event.getByLabel("bBCands",bbbarCand );
  if (bbbarCand.isValid()) {
    for (edm::View<Candidate>::const_iterator mcIter=bbbarCand->begin(); mcIter!=bbbarCand->end(); ++mcIter ) {
      std:: cout << " No CLEANING Charge b bbar= " << mcIter->daughter(0)->pdgId() << " " << mcIter->daughter(1)->pdgId() << " Mother b bbar= " << mcIter->daughter(0)->mother(0)->pdgId() << " " << mcIter->daughter(1)->mother(0)->pdgId() << " NCbbbar  Mass= " << mcIter->mass() << std::endl;
      histoMass_bbbarnc->Fill(mcIter->mass());
    }
  }
  
  // bbbar candidates cleaned 
  edm::Handle<edm::View<Candidate> > bbbarCandidates;
  event.getByLabel("bBCandscleaned",bbbarCandidates );
  if (bbbarCandidates.isValid()) { 
    for (edm::View<Candidate>::const_iterator mcIter=bbbarCandidates->begin(); mcIter!=bbbarCandidates->end(); ++mcIter ) {
      std:: cout << " CLEANED  Charge b bbar= " << mcIter->daughter(0)->pdgId() << " &  " << mcIter->daughter(1)->pdgId() << " Mother b bbar= " << mcIter->daughter(0)->mother(0)->pdgId() << " " << mcIter->daughter(1)->mother(0)->pdgId() << " CLEANED bbbar Mass= " << mcIter->mass() << std::endl;
      histoMass_bbbar->Fill(mcIter->mass()); 
      histoPt_bbbar->Fill(mcIter->pt()); 
      histoeta_bbbar->Fill(mcIter->eta()); 
      histophi_bbbar->Fill(mcIter->phi()); 
    } 
  }
  
  // just Z
  edm::Handle<edm::View<Candidate> > Z;
  event.getByLabel("Z",Z);
  for (edm::View<Candidate>::const_iterator mcIter=Z->begin(); mcIter!=Z->end(); ++mcIter){
    histoZmass->Fill(mcIter->mass());
  }
  
  // get Zbbbarcandidates
  edm::Handle<edm::View<Candidate> >  ZbbbarCandidates;
  event.getByLabel("llbBCands", ZbbbarCandidates);
  if (ZbbbarCandidates.isValid()) {
    for (edm::View<Candidate>::const_iterator mcIter=ZbbbarCandidates->begin(); mcIter!=ZbbbarCandidates->end(); ++mcIter ) {
      cout << " Zbbbar Mass= " << mcIter->mass()
	   << " pdgIds= " 
	   << mcIter->daughter(0)->pdgId() << " "
	   << mcIter->daughter(1)->daughter(0)->pdgId() << " " 
	   << mcIter->daughter(1)->daughter(1)->pdgId() << " " 
	   << endl;
      histoMass_Zbbbar->Fill(mcIter->mass()); 
      histoPt_Zbbbar->Fill(mcIter->pt()); 
      histoeta_Zbbbar->Fill(mcIter->eta()); 
      histophi_Zbbbar->Fill(mcIter->phi());
    }
  }
  
  // same flavor dilepton  candidates from Z (stable leptons mu and e)
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
  
  // sort lepton pt and apply skimming if needed !
  sort(leptonpt.rbegin(),leptonpt.rend());
  std::cout  << " \n Number of stable leptons= " << leptonpt.size() << std::endl;
  histonlept->Fill(leptonpt.size());
  if (leptonpt.size()>=4){
    std::cout << " \n Sorted vector in decreasing order= " << leptonpt.at(0) << " " << leptonpt.at(1) << " " << leptonpt.at(2) << " " << leptonpt.at(3) << std::endl;
    histoPt_1->Fill(leptonpt.at(0));
    histoPt_2->Fill(leptonpt.at(1));
    histoPt_3->Fill(leptonpt.at(2));
    histoPt_4->Fill(leptonpt.at(3));
  }
  
  //  trileptons cleaned candidates
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
	   mcIter->daughter(2)->mother(0)->mother(0)->pdgId()==411 || 
	   mcIter->daughter(2)->mother(0)->mother(0)->pdgId()==421 || 
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
  
  // get clean 4l candidates
  edm::Handle<edm::View<Candidate> >  fourlCandidates;
  event.getByLabel("fourleptonscleaned", fourlCandidates);
  if (fourlCandidates.isValid()) { 
    std::cout <<"\n size of fourleptonscleaned ===========================   " <<fourlCandidates->size()<< std::endl; 
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

std::string ZbbGenAnalyzer::getParticleName(int id) const
{
  const ParticleData * pd = pdt_->particle( id );
  if (!pd) {
    std::ostringstream ss;
    ss << "P" << id;
    return ss.str();
  } else
    return pd->name();
}
