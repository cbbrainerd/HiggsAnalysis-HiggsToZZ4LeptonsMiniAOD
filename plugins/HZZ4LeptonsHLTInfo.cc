/** \class HZZ4LeptonsHLTInfo
 *
 * See header file for documentation
 *
 *
 *  \author Nicola De Filippis
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsHLTInfo.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "HLTrigger/HLTfilters/interface/HLTHighLevel.h"

// Muons:
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>

// Electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <cassert>

using namespace std;
using namespace edm;
using namespace reco;

//
// constructors and destructor
//
HZZ4LeptonsHLTInfo::HZZ4LeptonsHLTInfo(const edm::ParameterSet& iConfig)
{
  // get names from module parameters, then derive slot numbers

  inputTag_           = consumes<edm::TriggerResults >(iConfig.getParameter<edm::InputTag> ("TriggerResultsTag"));
  n_                  = 0;
  firstevent_         = true;  
  produces<vector<std::string> >();

  // PG and FRC 06-07-11
  debug	=	iConfig.getUntrackedParameter<bool> ("debug", false);

}

HZZ4LeptonsHLTInfo::~HZZ4LeptonsHLTInfo()
{
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void HZZ4LeptonsHLTInfo::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  auto_ptr<vector<std::string> > HLTPathsFired( new vector<std::string> );

   const string invalid("@@invalid@@");

   // get hold of TriggerResults Object
   Handle<TriggerResults> trh;
   iEvent.getByToken(inputTag_,trh);
   
   if (trh.isValid()) {
     if(debug) cout << "TriggerResults found, number of HLT paths: " << trh->size() << endl;
     
     // get hold of trigger names - based on TriggerResults object!
     //triggerNames_.init(*trh);
     triggerNames_=iEvent.triggerNames(*trh);
     if (firstevent_) {
       for (unsigned int i=0; i<triggerNames_.size(); i++) {
	 if(debug) cout << "Found the trigger path= " << triggerNames_.triggerName(i) << endl;
	 firstevent_=false;
       }
     }
     
     // for empty input vectors (n==0), default to all HLT trigger paths!
     n_=trh->size();
     HLTPathsByName_.resize(n_);
     HLTPathsByIndex_.resize(n_);
     for (unsigned int i=0; i!=n_; i++) {
       HLTPathsByName_[i]=triggerNames_.triggerName(i);
       HLTPathsByIndex_[i]=i;
     }
          
     // count number of requested HLT paths which have fired
     for (unsigned int i=0; i!=n_; i++) {
       if (HLTPathsByIndex_[i]<trh->size()) {
	 if (trh->accept(HLTPathsByIndex_[i])) {
	   if(debug) cout << "Fired HLT path= " << HLTPathsByName_[i] << endl;
	   HLTPathsFired->push_back(HLTPathsByName_[i]);
	 }
       }
     }
     
   }

   iEvent.put(HLTPathsFired);     
   
   
}

void HZZ4LeptonsHLTInfo::endJob() {
}



