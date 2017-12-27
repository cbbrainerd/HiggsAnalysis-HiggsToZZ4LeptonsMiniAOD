/* \class HZZ4LeptonsCommonOfflineSelection
 *
 *
 * Tight isolation: electron and muons
 *
 * author:     Nicola De Filippis   - LLR-Ecole Plytechnique
 */


// system include files
#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsCommonOfflineSelection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Common/interface/ValueMap.h"

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

bool cmp( float a, float b ) {
  return a > b;
}
 
// Constructor
HZZ4LeptonsCommonOfflineSelection::HZZ4LeptonsCommonOfflineSelection(const edm::ParameterSet& pset) {

  // Decay Channel
  decaychannel                                                            = pset.getParameter<std::string>("decaychannel");
  useBestCandidate                                                        = pset.getParameter<bool> ("useBestCandidate" );
  BestCandidatesLeptonsTag_                                               = pset.getParameter<edm::InputTag>("BestCandidatesLeptons");

  // Tight isolation
  if (decaychannel=="2e2mu" || decaychannel=="4mu") muonTag_              = pset.getParameter<edm::InputTag>("MuonsLabel");
  if (decaychannel=="2e2mu" || decaychannel=="4mu") muonMapTag_           = pset.getParameter<edm::InputTag>("MuonsMapLabel");
  if (decaychannel=="2e2mu" || decaychannel=="4e")  electronTag_          = pset.getParameter<edm::InputTag>("ElectronsLabel");
  if (decaychannel=="2e2mu" || decaychannel=="4e")  electronMapTag_       = pset.getParameter<edm::InputTag>("ElectronsMapLabel");
  // tight isolation cuts
  if (decaychannel=="2e2mu" || decaychannel=="4e")  isoVarTagElectrons    = pset.getParameter<edm::InputTag>  ("isoVarTagElectrons");
  if (decaychannel=="2e2mu" || decaychannel=="4e")  isoVarCutElectrons    = pset.getParameter<vector<double> >("isoVarCutElectrons");
  if (decaychannel=="2e2mu" || decaychannel=="4mu") isoVarTagMuons        = pset.getParameter<edm::InputTag>  ("isoVarTagMuons");
  if (decaychannel=="2e2mu" || decaychannel=="4mu") isoVarCutMuons        = pset.getParameter<vector<double> >("isoVarCutMuons");

  // vertexing
  if (decaychannel=="2e2mu" || decaychannel=="4mu") muonTag_Vert          = pset.getParameter<edm::InputTag>("MuonsLabelVert");
  if (decaychannel=="2e2mu" || decaychannel=="4mu") muonMapTag_Vert       = pset.getParameter<edm::InputTag>("MuonsMapLabelVert");
  if (decaychannel=="2e2mu" || decaychannel=="4e")  electronTag_Vert      = pset.getParameter<edm::InputTag>("ElectronsLabelVert");
  if (decaychannel=="2e2mu" || decaychannel=="4e")  electronMapTag_Vert   = pset.getParameter<edm::InputTag>("ElectronsMapLabelVert");
  // vertexing cuts
  vertVarCut                                                              = pset.getParameter<vector<double> >("vertVarCut");

  cout << "Starting Offline selection for channel " << decaychannel << endl;

  string alias;
  if (decaychannel=="2e2mu" || decaychannel=="4e"){
    produces<bool> (alias = decaychannel + "OffselTightCombIsolEle").setBranchAlias( alias );
  }
  if (decaychannel=="2e2mu" || decaychannel=="4mu"){
    produces<bool> (alias = decaychannel + "OffselTightCombIsolMu").setBranchAlias( alias );
  }

  produces<bool> (alias = decaychannel + "OffselVertComb").setBranchAlias( alias );
  produces<bool> (alias = decaychannel + "Offsel"  ).setBranchAlias( alias );

}


// Destructor
HZZ4LeptonsCommonOfflineSelection::~HZZ4LeptonsCommonOfflineSelection() {

}


// Filter event (event preselection)
void HZZ4LeptonsCommonOfflineSelection::produce(edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  bool matched=false;
  
  // flags
  auto_ptr<bool> TightCombIsoEle ( new bool );
  *TightCombIsoEle  = false;
  auto_ptr<bool> TightCombIsoMu ( new bool );
  *TightCombIsoMu   = false;
  
  auto_ptr<bool> VertComb ( new bool );
  *VertComb   = false;
  
  // Offline selection flag
  auto_ptr<bool> boolOffsel ( new bool );
  *boolOffsel    = false;  
  
  //Tight isolation for electrons
  int nTightCombIsolEle=0;
  if (decaychannel=="2e2mu" || decaychannel=="4e"){
    
    Handle<edm::View<GsfElectron> > EleCandidates;
    iEvent.getByLabel(electronTag_, EleCandidates);
    
    Handle<edm::ValueMap<float> > isoelemap;
    iEvent.getByLabel(electronMapTag_, isoelemap);
    
    
    vector<float> isoelevector;
    int index=0;
    for (edm::View<reco::GsfElectron>::const_iterator cand = EleCandidates->begin(); 
	 cand != EleCandidates->end(); ++cand) {
      edm::Ref<edm::View<reco::GsfElectron> > eletrackref(EleCandidates,index);
      cout << "Isolation value from electron map " << (*isoelemap)[eletrackref] << endl;
      
      if (useBestCandidate==true) {
	Handle<CandidateCollection> bestleptonsCands;
	iEvent.getByLabel(BestCandidatesLeptonsTag_, bestleptonsCands);
	const reco::CandidateCollection* bestleptons = bestleptonsCands.product () ;
	matched=match(cand->mass(),cand->pt(),cand->charge(),bestleptons);
      }
      else {
	matched=true;
      }
      
      if (matched) isoelevector.push_back((*isoelemap)[eletrackref]);
      index ++;
    }
    // sorting in decreasing order
    sort(isoelevector.begin(),isoelevector.end(),cmp);
    // cut on the two least isolated electrons; sum of iso; if the least are tight isolated --> all are tight isolated
    
    if (isoelevector.size()>=2){
      if ( (isoelevector.at(0)+isoelevector.at(1))< isoVarCutElectrons.at(0) ) nTightCombIsolEle++;
    }
    else {
      nTightCombIsolEle=-999;
      cout << "Warning: there are no 2 values for isolation variables" << endl;
    }     
  }
  
  // Tight isolation for muons
  int nTightCombIsolMu=0;
  if (decaychannel=="2e2mu" || decaychannel=="4mu"){    
    
    Handle<edm::View<Muon> > MuCandidates;
    iEvent.getByLabel(muonTag_, MuCandidates);
    
    Handle<edm::ValueMap<float> > isomumap;
    iEvent.getByLabel(muonMapTag_, isomumap);
    
    vector<float> isomuvector;
    int indexbis=0;
    for (edm::View<reco::Muon>::const_iterator cand = MuCandidates->begin(); 
	 cand != MuCandidates->end(); ++cand) {
      edm::Ref<edm::View<reco::Muon> > mutrackref(MuCandidates,indexbis);
      cout << "Isolation value from muon map " << (*isomumap)[mutrackref] << endl;
      
      if (useBestCandidate==true) {
	Handle<CandidateCollection> bestleptonsCands;
	iEvent.getByLabel(BestCandidatesLeptonsTag_, bestleptonsCands);
	const reco::CandidateCollection* bestleptons = bestleptonsCands.product () ;
	matched=match(cand->mass(),cand->pt(),cand->charge(),bestleptons);
      }
      else {
	matched=true;
      }
      
      if (matched) isomuvector.push_back((*isomumap)[mutrackref]);
      indexbis++;
    }
    // sorting in decreasing order
    sort(isomuvector.begin(),isomuvector.end(),cmp);
    // cut on the two least isolated muons; sum of iso; if the least are tight isolated --> all are tight isolated
    if (isomuvector.size()>=2){
      if ( (isomuvector.at(0)+isomuvector.at(1))< isoVarCutMuons.at(0) ) nTightCombIsolMu++;
    }
    else {
      nTightCombIsolMu=-999;
      cout << "Warning: there are no 2 values for isolation variables" << endl;
    }     
  }

  // Vertexing      
  vector<float> vertexleptvector;
  int nVertComb=0;
  
  if (decaychannel=="2e2mu" || decaychannel=="4mu"){
    Handle<edm::View<Muon> > VertMuCandidates;
    iEvent.getByLabel(muonTag_Vert, VertMuCandidates);
    
    Handle<edm::ValueMap<float> > vertexmumap;
    iEvent.getByLabel(muonMapTag_Vert, vertexmumap);
    
    int indexvertbis=0;
    for (edm::View<reco::Muon> ::const_iterator cand = VertMuCandidates->begin(); 
	 cand != VertMuCandidates->end(); ++cand) {
      edm::Ref<edm::View<reco::Muon> > mutrackref(VertMuCandidates,indexvertbis);
      cout << "Vertexing value from muon map " << (*vertexmumap)[mutrackref] << endl;

      if (useBestCandidate==true) {
	Handle<CandidateCollection> bestleptonsCands;
	iEvent.getByLabel(BestCandidatesLeptonsTag_, bestleptonsCands);
	const reco::CandidateCollection* bestleptons = bestleptonsCands.product();
	matched=match(cand->mass(),cand->pt(),cand->charge(),bestleptons);
      }
      else {
	matched=true;
      }

      if (matched) vertexleptvector.push_back(fabs((*vertexmumap)[mutrackref]));
      indexvertbis++;
    }    
  }
  
  if (decaychannel=="2e2mu" || decaychannel=="4e"){
    Handle<edm::View<GsfElectron> > VertEleCandidates;
    iEvent.getByLabel(electronTag_Vert, VertEleCandidates);
    
    Handle<edm::ValueMap<float> > vertexelemap;
    iEvent.getByLabel(electronMapTag_Vert, vertexelemap);

    int indexvert=0;
    for (edm::View<reco::GsfElectron>::const_iterator cand = VertEleCandidates->begin(); 
         cand != VertEleCandidates->end(); ++cand) {
      edm::Ref<edm::View<reco::GsfElectron> > eletrackref(VertEleCandidates,indexvert);
      cout << "Vertexing value from electron map " << (*vertexelemap)[eletrackref] << endl;
    
      if (useBestCandidate==true) {
	Handle<CandidateCollection> bestleptonsCands;
	iEvent.getByLabel(BestCandidatesLeptonsTag_, bestleptonsCands);
	const reco::CandidateCollection* bestleptons = bestleptonsCands.product();
	matched=match(cand->mass(),cand->pt(),cand->charge(),bestleptons);
      }
      else {
	matched=true;
      }

      if (matched) vertexleptvector.push_back(fabs((*vertexelemap)[eletrackref]));
      indexvert++;
    }
  }
  
  // sorting in decreasing order
  sort(vertexleptvector.begin(),vertexleptvector.end(),cmp);
  // cut on the two least IP leptons;  
  if (vertexleptvector.size()>=2){
    if ( vertexleptvector.at(0)<vertVarCut.at(0) && vertexleptvector.at(1)<vertVarCut.at(1) ) nVertComb++;
  }
  else {
    nVertComb=-999;
    cout << "Warning: there are no 2 values for vertex significance" << endl;
  }
  
  
  //  2e2mu, 4e and 4mu channels
  if (decaychannel=="2e2mu" || decaychannel=="4e"){
    if (nTightCombIsolEle >=1) *TightCombIsoEle = true;
    iEvent.put(TightCombIsoEle,    decaychannel + "OffselTightCombIsolEle");
  }
  if (decaychannel=="2e2mu" || decaychannel=="4mu"){
    if (nTightCombIsolMu  >=1) *TightCombIsoMu  = true;
    iEvent.put(TightCombIsoMu,  decaychannel + "OffselTightCombIsolMu" );
  }
  
  if (nVertComb>=1) *VertComb=true; 
  iEvent.put(VertComb,    decaychannel + "OffselVertComb");

  if      (decaychannel=="2e2mu"){
    if ( nTightCombIsolEle>=1 && nTightCombIsolMu>=1 && nVertComb>=1 ) *boolOffsel=true;
  }
  else if (decaychannel=="4e"   ){
    if ( nTightCombIsolEle>=1 && nVertComb >=1 ) *boolOffsel=true;
  }
  else if (decaychannel=="4mu"  ){
    if ( nTightCombIsolMu>=1 && nVertComb >=1 )  *boolOffsel=true;
  }
  
  iEvent.put(boolOffsel, decaychannel + "Offsel");
  
}

void HZZ4LeptonsCommonOfflineSelection::beginJob() {
}

void HZZ4LeptonsCommonOfflineSelection::endJob() {
  cout << "Created offline selection variables" << endl;
}


bool HZZ4LeptonsCommonOfflineSelection::match(double mass, double pt, int charge, const reco::CandidateCollection *c1Coll){

  bool found=false;

  for( CandidateCollection::const_iterator pp = c1Coll->begin();pp != c1Coll->end(); ++ pp ) {

    if ((abs(pp->p4().mass()-mass)  <0.001 ) &&
        (abs(pp->p4().pt()  -pt)    <0.001 ) &&
        (abs(pp->charge()   -charge)<0.001 )  ){
      found=true;
      cout << "Found lepton in the best leptons collection" << endl;
    }
  }
  return found;
}

