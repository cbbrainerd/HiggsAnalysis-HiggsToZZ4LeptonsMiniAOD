/**\class HZZ4LeptonsMCParticleDecayProducer
 *
 * Original Author:  Nicola De Filippis
 * 
 * Class to store all the particles of a double decay like H->ZZ->2e2mu
 * Final leptons could be ordered n pT
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsMCParticleDecayProducer.h"

// Candidate handling
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// system include files
#include <memory>
#include <sstream>

using namespace edm;
using namespace std;
using namespace reco;


struct SortCandByDecreasingPt {
  bool operator()( const Candidate &c1, const Candidate &c2) const {
    return c1.pt() > c2.pt();
  }
};

// constructors
HZZ4LeptonsMCParticleDecayProducer::HZZ4LeptonsMCParticleDecayProducer(const edm::ParameterSet& iConfig) {

  genCandidates_=iConfig.getParameter<InputTag>("src");
  decayChain_=iConfig.getParameter<std::string>("decayChain");
  motherPdgId_=iConfig.getParameter<int >("motherPdgId");
  firstdaughtersPdgId_=iConfig.getParameter<vector<int> >("firstdaughtersPdgId"); 
  seconddaughtersPdgId_=iConfig.getParameter<vector<int> >("seconddaughtersPdgId"); 

  string alias;
  produces<CandidateCollection >(alias = decayChain_+ "Mother").setBranchAlias( alias );
  produces<CompositeCandidateCollection >(alias = decayChain_+ "CompositeMother").setBranchAlias( alias );

  firstdaughtersize = firstdaughtersPdgId_.size();
  daughtersize = seconddaughtersPdgId_.size();

  for (unsigned int j = 0; j < firstdaughtersize; ++ j ){
    ostringstream index,collection;
    index << j;
    collection << decayChain_ << "Boson" << index.str();
    valiasbosons.push_back(collection.str());
    produces<CandidateCollection >(valiasbosons.at(j)).setBranchAlias( valiasbosons.at(j) );
  }


  for (unsigned int j = 0; j < daughtersize; ++ j ){
    ostringstream index,collection;
    index << j;
    collection << decayChain_ << "Lepton" << index.str();
    valiasleptons.push_back(collection.str());
    produces<CandidateCollection >(valiasleptons.at(j)).setBranchAlias( valiasleptons.at(j) );
  }
}


// destructor
HZZ4LeptonsMCParticleDecayProducer::~HZZ4LeptonsMCParticleDecayProducer() {

}

void HZZ4LeptonsMCParticleDecayProducer::beginJob() {
}


void HZZ4LeptonsMCParticleDecayProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {


  ESHandle<HepPDT::ParticleDataTable> pdt_;
  iSetup.getData( pdt_ );

  const HepPDT::ParticleData *pd;

  pd = pdt_->particle( motherPdgId_ );
  cout << "Producing EDM particles:" <<endl;
  cout << "\t Mother name= " <<  pd->name() << " and code= " << motherPdgId_  <<endl;

  for (unsigned int j = 0; j < firstdaughtersPdgId_.size(); ++ j) {
    pd = pdt_->particle( firstdaughtersPdgId_[j]);
    cout << "\t Daughter name= " <<  pd->name() << " and code= " << firstdaughtersPdgId_[j]  <<endl;
  }
  for (unsigned int j = 0; j < seconddaughtersPdgId_.size(); ++ j) {
    pd = pdt_->particle( seconddaughtersPdgId_[j]);
    cout << "\t Final decay products name= " << pd->name() << " and code= " << seconddaughtersPdgId_[j]  <<endl;
  }


  // get gen particle candidates
  edm::Handle<GenParticleCollection> genCandidatesCollection;
  iEvent.getByLabel(genCandidates_, genCandidatesCollection);

  auto_ptr<CandidateCollection> mothercands_(new CandidateCollection);
  auto_ptr<CompositeCandidateCollection> motherCompositecands_(new CompositeCandidateCollection);
  auto_ptr<CandidateCollection> firstdaughterscands_(new CandidateCollection);
  auto_ptr<CandidateCollection> daughterscands_(new CandidateCollection);
  

  mothercands_->clear();
  motherCompositecands_->clear();
  firstdaughterscands_->clear();
  daughterscands_->clear();
  
  for( GenParticleCollection::const_iterator p = genCandidatesCollection->begin();p != genCandidatesCollection->end(); ++ p ) {
    if (p->pdgId() == motherPdgId_  && p->status() == 2 ){
      CompositeCandidate H;
      mothercands_->push_back(p->clone());
      size_t nfirstdau = p->numberOfDaughters();
      
      for( size_t i = 0; i < nfirstdau; ++ i ){
	signed int tmpid=0;
	for (size_t j = 0; j < firstdaughtersize && firstdaughtersPdgId_[j]!=tmpid ; ++ j ){
	  if (p->daughter(i)->pdgId()==firstdaughtersPdgId_[j] && p->daughter(i)->status()==2 ){
	    size_t tmpindex=j;
	    size_t nseconddau= p->daughter(i)->numberOfDaughters();
	    for (size_t k = 0; k < nseconddau; ++ k ){	
	      signed int tmpid2=0;
	      for (size_t l = 0; l < daughtersize && seconddaughtersPdgId_[l]!=tmpid2 ; ++ l ){
		if (p->daughter(i)->daughter(k)->pdgId()==seconddaughtersPdgId_[l] && p->daughter(i)->daughter(k)->status()==1 ){
		  if (tmpindex==j){
		    firstdaughterscands_->push_back(p->daughter(i)->clone());
		    H.addDaughter(*(p->daughter(i)->clone()));
		    tmpid=firstdaughtersPdgId_[j];
		    tmpindex=999;
		  }
		  daughterscands_->push_back(p->daughter(i)->daughter(k)->clone());
		  tmpid2=seconddaughtersPdgId_[l];
		}
	      }	      
            }
	  }
	}
      } 
      AddFourMomenta addP4_H;
      addP4_H.set( H );
      motherCompositecands_->push_back(*(H.clone()));
    }    
  }

  // iEvent.put(mothercands_, decayChain_ + "Mother");
  // iEvent.put(motherCompositecands_, decayChain_ + "CompositeMother");
  iEvent.put(std::make_unique<reco::CandidateCollection>(*mothercands_), decayChain_ + "Mother");
  iEvent.put(std::make_unique<reco::CompositeCandidateCollection>(*motherCompositecands_), decayChain_ + "CompositeMother");
  
  for (unsigned int row = 0; row < firstdaughtersize; ++ row ){
    auto_ptr<CandidateCollection> bosonscands_(new CandidateCollection);
    if (! firstdaughterscands_->empty()) bosonscands_->push_back((firstdaughterscands_->begin()+row)->clone());
    //iEvent.put(bosonscands_, valiasbosons.at(row));
    iEvent.put(std::make_unique<reco::CandidateCollection>(*bosonscands_), valiasbosons.at(row));
  }

  // daughterscands_->sort(SortCandByDecreasingPt());
  for (unsigned int row = 0; row < daughtersize; ++ row ){
    auto_ptr<CandidateCollection> leptonscands_(new CandidateCollection);
    if (! daughterscands_->empty()) leptonscands_->push_back((daughterscands_->begin()+row)->clone());
    //iEvent.put(leptonscands_, valiasleptons.at(row));
    iEvent.put(std::make_unique<reco::CandidateCollection>(*leptonscands_), valiasleptons.at(row));
  }

}


void HZZ4LeptonsMCParticleDecayProducer::endJob() {

}


