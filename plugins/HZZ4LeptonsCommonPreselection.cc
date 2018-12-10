/* \class HZZ4LeptonsCommonPreselection
 *
 *
 * H->ZZ->4l analysis preselection:
 * skim input
 * electron and muon selection
 * m_ll, m4l constraints
 * loose isolation on electrons and muons
 *
 * author:     Nicola De Filippis   - LLR-Ecole Plytechnique
 */


// system include files
#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsCommonPreselection.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectronFwd.h>
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <memory>

// namespaces
using namespace edm;
using namespace std;
using namespace reco;

struct SortCandByDecreasingPt {
  bool operator()( const Candidate &c1, const Candidate &c2) const {
    return c1.pt() > c2.pt();
  }
};

// Constructor
HZZ4LeptonsCommonPreselection::HZZ4LeptonsCommonPreselection(const edm::ParameterSet& pset) {

  // Decay Channel
  decaychannel                                                            = pset.getParameter<std::string>("decaychannel");
  if (decaychannel=="2e2mu" || decaychannel=="4mu") muonTag_              = pset.getParameter<edm::InputTag>("MuonsLabel");
  if (decaychannel=="2e2mu" || decaychannel=="4e")  electronTag_          = pset.getParameter<edm::InputTag>("ElectronsLabel");
  if (decaychannel=="2e2mu" || decaychannel=="4e")  zToEETag_             = pset.getParameter<edm::InputTag>("ZEELabel");
  if (decaychannel=="2e2mu" || decaychannel=="4mu") zToMuMuTag_           = pset.getParameter<edm::InputTag>("ZMMLabel");
  hTozzTo4leptonsTag_                                                     = pset.getParameter<edm::InputTag>("HLabel");
  if (decaychannel=="2e2mu" || decaychannel=="4mu") muonLooseIsolTag_     = pset.getParameter<edm::InputTag>("MuonsLooseIsolLabel");
  if (decaychannel=="2e2mu" || decaychannel=="4e")  electronLooseIsolTag_ = pset.getParameter<edm::InputTag>("ElectronsLooseIsolLabel");

  edm::ParameterSet cutsconf                                              = pset.getParameter<edm::ParameterSet>("cuts");
  if (decaychannel=="2e2mu" || decaychannel=="4e")  nEle_cut              = cutsconf.getParameter<int>("nEle");
  if (decaychannel=="2e2mu" || decaychannel=="4mu") nMu_cut               = cutsconf.getParameter<int>("nMu");
  if (decaychannel=="2e2mu" || decaychannel=="4e")  eeMass_cut            = cutsconf.getParameter<double>("eeMass");
  if (decaychannel=="2e2mu" || decaychannel=="4mu") mumuMass_cut          = cutsconf.getParameter<double>("mumuMass");
  fourlMass_cut                                                           = cutsconf.getParameter<double>("fourleptMass");
  if (decaychannel=="2e2mu" || decaychannel=="4e")  numberOfeeCombs_cut   = cutsconf.getParameter<int>("numberOfeeCombs");
  if (decaychannel=="2e2mu" || decaychannel=="4mu") numberOfmumuCombs_cut = cutsconf.getParameter<int>("numberOfmumuCombs");
  numberOf4lCombs_cut                                                     = cutsconf.getParameter<int>("numberOf4lCombs");
  if (decaychannel=="2e2mu" || decaychannel=="4e")  nlooseEle_cut         = cutsconf.getParameter<int>("nlooseEle");
  if (decaychannel=="2e2mu" || decaychannel=="4mu") nlooseMu_cut          = cutsconf.getParameter<int>("nlooseMu");

  cout << "Starting preselection for channel " << decaychannel << endl;

  string alias;
  if (decaychannel=="2e2mu"){
    produces<bool> (alias = decaychannel + "PreselAtleast2Ele"  ).setBranchAlias( alias );
    produces<bool> (alias = decaychannel + "PreselAtleast2Mu"   ).setBranchAlias( alias );
    produces<bool> (alias = decaychannel + "PreselAtleast1ZEE"  ).setBranchAlias( alias );
    produces<bool> (alias = decaychannel + "PreselAtleast1ZMuMu").setBranchAlias( alias );
    produces<bool> (alias = decaychannel + "PreselLoose2IsolEle").setBranchAlias( alias );
    produces<bool> (alias = decaychannel + "PreselLoose2IsolMu" ).setBranchAlias( alias );
  }
  else if (decaychannel=="4e"){
    produces<bool> (alias = decaychannel + "PreselAtleast4Ele"  ).setBranchAlias( alias );
    produces<bool> (alias = decaychannel + "PreselAtleast2ZEE"  ).setBranchAlias( alias );
    produces<bool> (alias = decaychannel + "PreselLoose4IsolEle").setBranchAlias( alias );
  }
  else if (decaychannel=="4mu"){
    produces<bool> (alias = decaychannel + "PreselAtleast4Mu"   ).setBranchAlias( alias );
    produces<bool> (alias = decaychannel + "PreselAtleast2ZMuMu").setBranchAlias( alias );
    produces<bool> (alias = decaychannel + "PreselLoose4IsolMu" ).setBranchAlias( alias );
  }

  produces<bool> (alias = decaychannel + "PreselAtleast1H"      ).setBranchAlias( alias );
  produces<bool> (alias = decaychannel + "Presel"               ).setBranchAlias( alias );

}


// Destructor
HZZ4LeptonsCommonPreselection::~HZZ4LeptonsCommonPreselection() {

}


// Filter event (event preselection)
void HZZ4LeptonsCommonPreselection::produce(edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  //2e2mu
  auto_ptr<bool> atleast2Ele ( new bool );
  auto_ptr<bool> atleast2Mu ( new bool);
  auto_ptr<bool> atleast1ZEE ( new bool );
  auto_ptr<bool> atleast1ZMuMu ( new bool );
  auto_ptr<bool> Loose2IsoEle ( new bool );
  auto_ptr<bool> Loose2IsoMu ( new bool );
  
  *atleast2Ele   = false;
  *atleast2Mu    = false;
  *atleast1ZEE   = false;
  *atleast1ZMuMu = false;
  *Loose2IsoEle  = false;
  *Loose2IsoMu   = false;
  
  //4e
  auto_ptr<bool> atleast4Ele ( new bool );
  auto_ptr<bool> atleast2ZEE ( new bool ); 
  auto_ptr<bool> Loose4IsoEle ( new bool );
  *atleast4Ele   = false;
  *atleast2ZEE   = false;
  *Loose4IsoEle  = false;
  
  //4mu
  auto_ptr<bool> atleast4Mu ( new bool );
  auto_ptr<bool> atleast2ZMuMu ( new bool );
  auto_ptr<bool> Loose4IsoMu ( new bool );
  *atleast4Mu    = false;
  *atleast2ZMuMu = false;
  *Loose4IsoMu   = false;
    
  auto_ptr<bool> atleast1H ( new bool );
  *atleast1H     = false;

  // Preselection flag
  auto_ptr<bool> boolPresel ( new bool );
  *boolPresel    = false;

   // Selected Electrons
  int posEle=0,negEle=0;
  if (decaychannel=="2e2mu" || decaychannel=="4e"){
    edm::Handle<edm::View<GsfElectron> > electronsHandle;
    edm::View<GsfElectron>::const_iterator eIter;
    iEvent.getByLabel(electronTag_.label(), electronsHandle);
    for (eIter = electronsHandle->begin(); eIter != electronsHandle->end(); ++eIter ) {
      cout << "Electron with pt= " <<  eIter->pt() << "  eta= " << eIter->eta() << "  p= " <<  eIter->p() << endl;
      if ( eIter->charge()==1){
        posEle++;
      }
      else if ( eIter->charge()==-1){
        negEle++;
    }
    }
  }

  // Selected Muons
  int posMu=0,negMu=0;
  if (decaychannel=="2e2mu" || decaychannel=="4mu"){
    Handle<edm::View<Muon> > muonsHandle;
    edm::View<Muon>::const_iterator muIter;
    iEvent.getByLabel(muonTag_.label(), muonsHandle);
    for (muIter = muonsHandle->begin(); muIter != muonsHandle->end(); ++muIter ) {
      cout << "Muon with pt= " <<  muIter->pt() << "  eta= " << muIter->eta() << "  p= " <<  muIter->p() << endl;
      if ( muIter->charge()==1 ){
        posMu++;
      }
      else if ( muIter->charge()==-1 ){
        negMu++;
      }
    }
  }
        
  /// channel conditions
  nElectron=0;
  nMuon=0;
  if (decaychannel=="2e2mu") {
    if ( posEle>=1 && negEle>=1 ) {
      nElectron=posEle+negEle;
    }
    if ( posMu>=1 && negMu>=1 ) {
      nMuon=posMu+negMu;
    }
  }
  else if (decaychannel=="4e") {
    if ( posEle>=2 && negEle>=2 ) {
      nElectron=posEle+negEle;
    }
  }
  else if (decaychannel=="4mu") {
    if ( posMu>=2 && negMu>=2 ) {
      nMuon=posMu+negMu;
    }  
  }
       

  // Pairs of EE MuMu 
  int nZEE=0,nZMuMu=0;
  if (decaychannel=="2e2mu"){
    Handle<edm::View<Candidate> >  zEECandidates;
    iEvent.getByLabel(zToEETag_.label(), zEECandidates);    
    for ( edm::View<Candidate>::const_iterator zIter=zEECandidates->begin(); zIter!= zEECandidates->end(); ++zIter ) {
      cout << "Zee mass= " << double(zIter->p4().mass()) << endl;
      if ( double(zIter->p4().mass())> eeMass_cut ){ 
	nZEE++;
      }
    }
    Handle<edm::View<Candidate> >  zMuMuCandidates;
    iEvent.getByLabel(zToMuMuTag_.label(), zMuMuCandidates);
    for ( edm::View<Candidate>::const_iterator zIter=zMuMuCandidates->begin(); zIter!= zMuMuCandidates->end(); ++zIter ) {
      cout << "Zmumu mass= " << double(zIter->p4().mass()) << endl;
      if ( zIter->p4().mass()> mumuMass_cut ){
	nZMuMu++;
      }
    }    
  }
  
  // exclusive couples ZEE and ZMUMU
  if (decaychannel=="4e" ){
    Handle<edm::View<Candidate> >  higgsCandidates;
    iEvent.getByLabel(hTozzTo4leptonsTag_.label(), higgsCandidates);
    for ( edm::View<Candidate>::const_iterator hIter=higgsCandidates->begin(); hIter!= higgsCandidates->end(); ++hIter ) {
      if (nZEE<2) nZEE=0;
      for (size_t ndau=0; ndau<hIter->numberOfDaughters();ndau++){
	if ( hIter->daughter(ndau)->p4().mass()> eeMass_cut){  
	  nZEE++;
	}
      }
    }
  }
  //
  if (decaychannel=="4mu" ){
    Handle<edm::View<Candidate> >  higgsCandidates;
    iEvent.getByLabel(hTozzTo4leptonsTag_.label(), higgsCandidates);
    for ( edm::View<Candidate>::const_iterator hIter=higgsCandidates->begin(); hIter!= higgsCandidates->end(); ++hIter ) {
      if (nZMuMu<2) nZMuMu=0;
      for (size_t ndau=0; ndau<hIter->numberOfDaughters();ndau++){
	if ( hIter->daughter(ndau)->p4().mass()> mumuMass_cut ){  
	  nZMuMu++;
	}
      }
    }
  }


  // 4 lepton combinations
  Handle<edm::View<Candidate> >  higgsCandidates;
  iEvent.getByLabel(hTozzTo4leptonsTag_.label(), higgsCandidates);
   int nHiggs=0;
  for ( edm::View<Candidate>::const_iterator hIter=higgsCandidates->begin(); hIter!= higgsCandidates->end(); ++hIter ) {
    if ( hIter->p4().mass()> fourlMass_cut ){  
      nHiggs++;
    }
  }    

  // Loose isolation for electrons
  if (decaychannel=="2e2mu" || decaychannel=="4e"){
    Handle<edm::View<GsfElectron> > electrons;
    edm::View<GsfElectron>::const_iterator eIter;
    iEvent.getByLabel(electronLooseIsolTag_.label(), electrons);
    nLooseIsolEle=0;
    nLooseIsolElepos=0;
    nLooseIsolEleneg=0;
    for (eIter = electrons->begin(); eIter != electrons->end(); ++eIter ) {
      cout << "Loose isolated electrons with pt and charge= " << eIter->pt() << " " << eIter->charge() <<endl;
      if (eIter->charge()>0) nLooseIsolElepos++;
      if (eIter->charge()<0) nLooseIsolEleneg++;
    }
    
    if (decaychannel=="2e2mu" && nLooseIsolElepos >=1 && nLooseIsolEleneg>=1 )  nLooseIsolEle=nLooseIsolElepos+nLooseIsolEleneg;
    if (decaychannel=="4e" && nLooseIsolElepos >=2 && nLooseIsolEleneg>=2 )  nLooseIsolEle=nLooseIsolElepos+nLooseIsolEleneg;
    cout <<"Number of loose isolated electrons (matching charge)= " << nLooseIsolEle << endl;
  }


  // Loose isolation for muons
  if (decaychannel=="2e2mu" || decaychannel=="4mu"){
    Handle<edm::View<Muon> > muons;
    edm::View<Muon>::const_iterator muIter;
    iEvent.getByLabel(muonLooseIsolTag_.label(), muons);
    nLooseIsolMu=0;
    nLooseIsolMupos=0;
    nLooseIsolMuneg=0;
    for (muIter = muons->begin(); muIter != muons->end(); ++muIter ) {
      cout << "Loose isolated muons with pt and charge= " << muIter->pt() << " " << muIter->charge() <<endl;
      if (muIter->charge()>0) nLooseIsolMupos++;
      if (muIter->charge()<0) nLooseIsolMuneg++;
    }
    
    if (decaychannel=="2e2mu" && nLooseIsolMupos >=1 && nLooseIsolMuneg>=1 )  nLooseIsolMu=nLooseIsolMupos+nLooseIsolMuneg;
    if (decaychannel=="4mu" && nLooseIsolMupos >=2 && nLooseIsolMuneg>=2 )  nLooseIsolMu=nLooseIsolMupos+nLooseIsolMuneg;
    cout <<"Number of loose isolated muons (matching charge)= " << nLooseIsolMu << endl;
  } 

  // 4Leptons channel
  if (decaychannel=="2e2mu"){
    if (nElectron >= nEle_cut ) *atleast2Ele     = true;
    if (nMuon     >= nMu_cut ) *atleast2Mu      = true;
    if (nZEE      >= numberOfeeCombs_cut )   *atleast1ZEE     = true;
    if (nZMuMu    >= numberOfmumuCombs_cut ) *atleast1ZMuMu   = true;
    if (nLooseIsolEle >= nlooseEle_cut )        *Loose2IsoEle = true;
    if (nLooseIsolMu  >= nlooseMu_cut)          *Loose2IsoMu  = true;
    
    // iEvent.put(atleast2Ele,   decaychannel + "PreselAtleast2Ele");
    // iEvent.put(atleast2Mu,    decaychannel + "PreselAtleast2Mu"); 
    // iEvent.put(atleast1ZEE,   decaychannel + "PreselAtleast1ZEE");
    // iEvent.put(atleast1ZMuMu, decaychannel + "PreselAtleast1ZMuMu");
    // iEvent.put(Loose2IsoEle,  decaychannel + "PreselLoose2IsolEle");
    // iEvent.put(Loose2IsoMu,   decaychannel + "PreselLoose2IsolMu"); 
    
    
    iEvent.put(std::make_unique<bool>(*atleast2Ele),   decaychannel + "PreselAtleast2Ele");
    iEvent.put(std::make_unique<bool>(*atleast2Mu),    decaychannel + "PreselAtleast2Mu"); 
    iEvent.put(std::make_unique<bool>(*atleast1ZEE),   decaychannel + "PreselAtleast1ZEE");
    iEvent.put(std::make_unique<bool>(*atleast1ZMuMu), decaychannel + "PreselAtleast1ZMuMu");
    iEvent.put(std::make_unique<bool>(*Loose2IsoEle),  decaychannel + "PreselLoose2IsolEle");
    iEvent.put(std::make_unique<bool>(*Loose2IsoMu),   decaychannel + "PreselLoose2IsolMu"); 
    
    
    if ( (nElectron >= nEle_cut)  && (nMuon >= nMu_cut) &&
         (nZEE >= numberOfeeCombs_cut) && (nZMuMu >= numberOfmumuCombs_cut) && (nHiggs >= numberOf4lCombs_cut ) &&
         (nLooseIsolEle >= nlooseEle_cut) && (nLooseIsolMu >= nlooseMu_cut) ) *boolPresel=true;
  }
  else if (decaychannel=="4e"){
    if (nElectron >= nEle_cut )  *atleast4Ele    = true;
    if (nZEE >= numberOfeeCombs_cut)        *atleast2ZEE    = true;
    if (nLooseIsolEle >= nlooseEle_cut)     *Loose4IsoEle   = true;
    
    // iEvent.put(atleast4Ele,   decaychannel + "PreselAtleast4Ele");
    // iEvent.put(atleast2ZEE,   decaychannel + "PreselAtleast2ZEE");
    // iEvent.put(Loose4IsoEle,  decaychannel + "PreselLoose4IsolEle");
    
    iEvent.put(std::make_unique<bool>(*atleast4Ele),   decaychannel + "PreselAtleast4Ele");
    iEvent.put(std::make_unique<bool>(*atleast2ZEE),   decaychannel + "PreselAtleast2ZEE");
    iEvent.put(std::make_unique<bool>(*Loose4IsoEle),  decaychannel + "PreselLoose4IsolEle");
    
    if ( (nElectron >= nEle_cut) && (nZEE >= numberOfeeCombs_cut) && (nHiggs >= numberOf4lCombs_cut) &&
         (nLooseIsolEle >= nlooseEle_cut) ) *boolPresel=true;
  }
  else if (decaychannel=="4mu"){
    if (nMuon >= nMu_cut)      *atleast4Mu     = true;
    if (nZMuMu >= numberOfmumuCombs_cut )   *atleast2ZMuMu = true;
    if (nLooseIsolMu >= nlooseMu_cut )      *Loose4IsoMu   = true;
    
    // iEvent.put(atleast4Mu,    decaychannel + "PreselAtleast4Mu");   
    // iEvent.put(atleast2ZMuMu, decaychannel + "PreselAtleast2ZMuMu");
    // iEvent.put(Loose4IsoMu,   decaychannel + "PreselLoose4IsolMu"); 
    
    iEvent.put(std::make_unique<bool>(*atleast4Mu),    decaychannel + "PreselAtleast4Mu");   
    iEvent.put(std::make_unique<bool>(*atleast2ZMuMu), decaychannel + "PreselAtleast2ZMuMu");
     iEvent.put(std::make_unique<bool>(*Loose4IsoMu),   decaychannel + "PreselLoose4IsolMu");
     
     if ( (nMuon >= nMu_cut) && (nZMuMu >= numberOfmumuCombs_cut) && (nHiggs >= numberOf4lCombs_cut) &&
	  (nLooseIsolMu >= nlooseMu_cut)) *boolPresel=true;
  }
         
  if (nHiggs >= numberOf4lCombs_cut)      *atleast1H         = true;
  
  // iEvent.put(atleast1H,  decaychannel + "PreselAtleast1H"); 
  //iEvent.put(boolPresel, decaychannel + "Presel");
  
  iEvent.put(std::make_unique<bool>(*atleast1H),  decaychannel + "PreselAtleast1H");
  iEvent.put(std::make_unique<bool>(*boolPresel), decaychannel + "Presel");
  
  
}

void HZZ4LeptonsCommonPreselection::beginJob() {
  cout << "Starting preselection" << endl;
}

void HZZ4LeptonsCommonPreselection::endJob() {
  cout << "Create preselection variables" << endl;
}


