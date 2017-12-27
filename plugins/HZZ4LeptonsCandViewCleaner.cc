////////////////////////////////////////////////////////////////////////////////
//
// HZZ4LeptonsCandViewCleaner
// --------------
//
////////////////////////////////////////////////////////////////////////////////


#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsCandViewCleaner.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
HZZ4LeptonsCandViewCleaner::HZZ4LeptonsCandViewCleaner(const edm::ParameterSet& iConfig)
  : srcCands_    (iConfig.getParameter<edm::InputTag>        ("srcObject"))
  , srcObjects_ (iConfig.getParameter<edm::InputTag>         ("srcObjectsToRemove"))
  , moduleLabel_(iConfig.getParameter<string>                ("@module_label"))
  , nCandidatesTot_(0)
  , nCandidatesClean_(0)
{
  produces<edm::RefToBaseVector<reco::Candidate> >();
}


//______________________________________________________________________________
HZZ4LeptonsCandViewCleaner::~HZZ4LeptonsCandViewCleaner()
{
  
}



////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void HZZ4LeptonsCandViewCleaner::produce(edm::Event& iEvent,const edm::EventSetup& iSetup)
{
  auto_ptr<edm::RefToBaseVector<reco::Candidate> >
    cleanCandidates(new edm::RefToBaseVector<reco::Candidate>());

  edm::Handle<reco::CandidateView> candidates;
  iEvent.getByLabel(srcCands_,candidates);
  
  bool* isClean = new bool[candidates->size()];
  for (unsigned int iCandidate=0;iCandidate<candidates->size();iCandidate++) isClean[iCandidate] = true;
  
  edm::Handle<reco::CandidateView> objects;
  iEvent.getByLabel(srcObjects_,objects);

  //std::cout << "Size=" << candidates->size() << " " << objects->size()<< std::endl;
  for (unsigned int iCandidate=0;iCandidate<candidates->size();iCandidate++) {
    int nfound=0;
    const reco::Candidate& candidate = candidates->at(iCandidate);
    //cout << "counter cand" << iCandidate << " with mass= " << candidate.mass() << endl;
    for (unsigned int iObj=0;iObj<objects->size();iObj++) {
      const reco::Candidate& obj = objects->at(iObj);
      //cout << "obj cand with mass= " << obj.mass() << endl;
      if (candidate.mass()==obj.mass() && candidate.pt()==obj.pt() )  nfound++;
      //std::cout << "nfound=" << nfound <<std::endl;
      if (nfound>=2) {
	//std::cout  << "found a double particle" << std::endl;
	isClean[iObj] = false;
	break;
      }
    }
  }
  
  for (unsigned int iCandidate=0;iCandidate<candidates->size();iCandidate++)
    if (isClean[iCandidate]) cleanCandidates->push_back(candidates->refAt(iCandidate));
  
  nCandidatesTot_  +=candidates->size();
  nCandidatesClean_+=cleanCandidates->size();

  delete [] isClean;  
  iEvent.put(cleanCandidates);
}


//______________________________________________________________________________
void HZZ4LeptonsCandViewCleaner::endJob()
{
  stringstream ss;
  ss<<"nCandidatesTot="<<nCandidatesTot_<<" nCandidatesClean="<<nCandidatesClean_
    <<" fCandidatesClean="<<100.*(nCandidatesClean_/(double)nCandidatesTot_)<<"%\n";
  cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++"
      <<"\n"<<moduleLabel_<<"(HZZ4LeptonsCandViewCleaner) SUMMARY:\n"<<ss.str()
      <<"++++++++++++++++++++++++++++++++++++++++++++++++++"
      <<endl;
}


