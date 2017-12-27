/**\class HZZ4LeptonsBestCandidate
 *
 * Original Author:  Nicola De Filippis
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsBestCandidate.h"

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"



// system include files
#include <memory>
#include <sstream>

using namespace edm;
using namespace std;
using namespace reco;


struct SortCandByDecreasingPt {
  bool operator()( const Candidate &c1, const Candidate &c2) const {
    return c1.p4().pt() > c2.p4().pt();
  }
};

// constructors
HZZ4LeptonsBestCandidate::HZZ4LeptonsBestCandidate(const edm::ParameterSet& iConfig) {
  decaychannel = iConfig.getParameter<std::string>("decaychannel");
  typedef std::vector<edm::InputTag> vtag;
  RECOcollName = iConfig.getParameter<vtag>("RECOcollName");
  decayChain_  = iConfig.getParameter<std::string>("decayChain");

	// PG and FRC 06-07-11
	debug	=	iConfig.getUntrackedParameter<bool> ("debug", false);

  string alias;
  produces<CandidateCollection >(alias = decayChain_+ "Mother").setBranchAlias( alias );

  for (unsigned int j = 0; j < 2; ++ j ){
    ostringstream index,collection;
    index << j;
    collection << decayChain_ << "Boson" << index.str();
    valiasbosons.push_back(collection.str());
    produces<CandidateCollection >(valiasbosons.at(j)).setBranchAlias( valiasbosons.at(j) );
  }

  produces<CandidateCollection >(alias = decayChain_+ "Leptons").setBranchAlias( alias );

}


// destructor
HZZ4LeptonsBestCandidate::~HZZ4LeptonsBestCandidate() {

}


void HZZ4LeptonsBestCandidate::beginJob() {
  
}


void HZZ4LeptonsBestCandidate::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<CandidateCollection> mothercands_(new CandidateCollection);
  auto_ptr<CandidateCollection> daughterscands_Z(new CandidateCollection);
  auto_ptr<CandidateCollection> daughterscands_Zstar(new CandidateCollection);
  auto_ptr<CandidateCollection> leptonscands_(new CandidateCollection);
  auto_ptr<CandidateCollection> bestleptonscands_(new CandidateCollection);
  
  mothercands_         ->clear();
  daughterscands_Z     ->clear();
  daughterscands_Zstar ->clear();
  leptonscands_        ->clear();
  bestleptonscands_    ->clear();

  // Higgs candidate in input
  Handle<edm::View<Candidate> > Candidates;
  iEvent.getByLabel(RECOcollName.at(0), Candidates);

  float ZNomMass  = 91.1876;
  //  float zcandMass = 0.;
  float deltaZ    = 9999999;

  const Candidate *bestZshell=NULL;

  if (Candidates->size()>0){

    for ( edm::View<Candidate>::const_iterator hIter=Candidates->begin(); hIter!= Candidates->end(); ++hIter ){
      if(debug) cout << "Input higgs candidate with mass= " << hIter->p4().mass() << endl;
      for (size_t j=0; j<hIter->numberOfDaughters();j++){ 
	// in case of muons want to build Z on-shell with GM
	if (hIter->daughter(j)->daughter(0)->isMuon() && hIter->daughter(j)->daughter(1)->isMuon()){
	  if(debug) cout << "Z mass= " << hIter->daughter(j)->p4().mass()
	       << " made of Mu with pt= " << hIter->daughter(j)->daughter(0)->p4().pt() << " " << hIter->daughter(j)->daughter(1)->p4().pt() 
	       << " and isGM= " << hIter->daughter(j)->daughter(0)->isGlobalMuon() << " " << hIter->daughter(j)->daughter(1)->isGlobalMuon() << endl;
	  
	  if ( fabs(hIter->daughter(j)->p4().mass()-ZNomMass ) < deltaZ && hIter->daughter(j)->daughter(0)->isGlobalMuon() && hIter->daughter(j)->daughter(1)->isGlobalMuon() ) {
	    deltaZ    = fabs(hIter->daughter(j)->p4().mass()-ZNomMass);
	    if(debug) cout << "Delta Z= " << deltaZ << endl;
	    //zcandMass = hIter->daughter(j)->p4().mass();  
	    bestZshell=hIter->daughter(j)->clone();
	  }
	}
	else {
	  if(debug) cout << "Z mass= " << hIter->daughter(j)->p4().mass()
	       << " made of Mu with pt= " << hIter->daughter(j)->daughter(0)->p4().pt() << " " << hIter->daughter(j)->daughter(1)->p4().pt() << endl;
	  
	  if ( fabs(hIter->daughter(j)->p4().mass()-ZNomMass ) < deltaZ  ) {
	    deltaZ    = fabs(hIter->daughter(j)->p4().mass()-ZNomMass);
	    if(debug) cout << "Delta Z= " << deltaZ << endl;
	    //zcandMass = hIter->daughter(j)->p4().mass();  
	    bestZshell=hIter->daughter(j)->clone();
	  }
	}

	if (!find(leptonscands_,*hIter->daughter(j)->daughter(0)->clone()) ){	  
	  leptonscands_->push_back(hIter->daughter(j)->daughter(0)->clone());
	  if(debug) cout << "Saving lepton" << endl;
	}
	if (!find(leptonscands_,*hIter->daughter(j)->daughter(1)->clone())) {	  
	  leptonscands_->push_back(hIter->daughter(j)->daughter(1)->clone());
	  if(debug) cout << "Saving lepton" << endl;
	}
      }
    }      
    
    if (bestZshell && bestZshell->numberOfDaughters()==2){
      if(debug) cout << "On shell Z is= "    << bestZshell->p4().mass() 
	   << "  with leptons pt/charge=" << bestZshell->daughter(0)->p4().pt() << "/" << bestZshell->daughter(0)->charge()
	   << " and "              << bestZshell->daughter(1)->p4().pt() << "/" << bestZshell->daughter(1)->charge() << endl;
      daughterscands_Z->push_back(bestZshell->clone());
      bestleptonscands_->push_back(bestZshell->daughter(0)->clone());
      bestleptonscands_->push_back(bestZshell->daughter(1)->clone());
      
      // Build the best Z* from the high pt leptons
      CompositeCandidate bestZstar;
      int tmpcharge_l=0,pairfound=0;
      //      int flavour=0;
      int tmpflavour=0;
      bool isGM=true;
      
      // Ordering leptons in pT; the highest pt lepton is the first to build Z*
      leptonscands_->sort(SortCandByDecreasingPt());
      if(debug) cout << "Size of lepton collection from Zs is= " << leptonscands_->size() << endl;
      
      if (leptonscands_->size() >= 4) {
	
	for ( CandidateCollection::const_iterator leptonIter=leptonscands_->begin(); leptonIter!= leptonscands_->end(); ++leptonIter ) {
	  // cout << "lepton from z mass= " << leptonIter->mother(0)->p4().mass()<< endl;
	  if ( leptonIter->p4()==bestZshell->daughter(0)->p4() ||  leptonIter->p4()==bestZshell->daughter(1)->p4() ) {
	    if(debug) cout << "Leptons from on shell Z already used...skipped pdgId, pT= " << leptonIter->pdgId() << " " << leptonIter->p4().pt() << " and charge= " << leptonIter->charge() << endl;
	    //flavour=leptonIter->pdgId();
	  }
	  else {
	    // cout << "id, charge of lepton: " << leptonIter->pdgId() << " " << leptonIter->charge() << "  and pt= " << leptonIter->p4().pt() << endl;
	    if (decaychannel=="2e2mu") tmpflavour=bestZshell->daughter(0)->pdgId();
	    if (leptonIter->isMuon() ) isGM=leptonIter->isGlobalMuon();
	    if ( abs(leptonIter->pdgId())!=abs(tmpflavour) && leptonIter->charge()!=tmpcharge_l && pairfound<2 && isGM) {
	      bestZstar.addDaughter(*leptonIter);
	      tmpcharge_l=leptonIter->charge();
	      if(debug) cout << "Lepton used for Z* building: pdgId, pt, eta, phi= " 
		   << leptonIter->pdgId() << " " 
		   << leptonIter->pt() << " " 
		   << leptonIter->eta() << " " 
		   << leptonIter->phi() << endl;
	      pairfound++;
	    }
	  }
	}
	
	
	// Z* shuold have two daughters
	if (bestZstar.numberOfDaughters()==2){
	  AddFourMomenta addP4;
	  addP4.set( bestZstar );
	  if(debug) cout << "bestZstar mass=" << bestZstar.p4().mass() << endl;
	  // cout << "bestZstar leptons pt=" << bestZstar.daughter(0)->p4().pt() << " " << bestZstar.daughter(1)->p4().pt() << endl;
	  daughterscands_Zstar->push_back(bestZstar.clone());  
	  bestleptonscands_->push_back(bestZstar.daughter(0)->clone());
	  bestleptonscands_->push_back(bestZstar.daughter(1)->clone());
	  
	  // Build the best H from Z and Z*
	  CompositeCandidate bestH;
	  bestH.addDaughter(*bestZshell);
	  bestH.addDaughter(bestZstar);
	  
	  AddFourMomenta addP4_H;
	  addP4_H.set( bestH );
	  if(debug) cout << "bestH " << bestH.p4().mass() << endl;
	  mothercands_->push_back(bestH.clone());
	}
	else {
	  if(debug) cout << "No 2 daughters for Z* -> no Z* satisfying conditions" << endl;
	  if(debug) cout << "No Best Higgs candidate built" << endl;
	}
      }
    }
    else {
      if(debug) cout << "No 2 daughters for Z on shell satisfying conditions" << endl;
      if(debug) cout << "No Best Higgs candidate built" << endl;
    }
  }
  
  iEvent.put(mothercands_, decayChain_ + "Mother");
  iEvent.put(daughterscands_Z, valiasbosons.at(0));
  iEvent.put(daughterscands_Zstar, valiasbosons.at(1));
  iEvent.put(bestleptonscands_, decayChain_ + "Leptons");

}


void HZZ4LeptonsBestCandidate::endJob() {

}

bool HZZ4LeptonsBestCandidate::find(const std::auto_ptr<reco::CandidateCollection>& c1Coll, const reco::Candidate& c2){
  
  bool found=false;
  
  for( CandidateCollection::const_iterator pp = c1Coll->begin();pp != c1Coll->end(); ++ pp ) {
    if ((abs(pp->p4().mass() - c2.p4().mass())< 0.01 ) &&
	(abs(pp->p4().pt()   - c2.p4().pt())  < 0.01 ) &&
        (abs(pp->p4().phi()  - c2.p4().phi()) < 0.01 ) &&
	(abs(pp->charge()    - c2.charge())   < 0.01 )  ){
      found=true;
      if(debug) cout << "Found already registered candidate" << endl;
    }
  }
  return found;
}


