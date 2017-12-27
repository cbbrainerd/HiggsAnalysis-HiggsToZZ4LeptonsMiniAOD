#ifndef HZZ4LeptonsHLTInfo_h
#define HZZ4LeptonsHLTInfo_h

/** \class HZZ4LeptonsHLTInfo
 *
 *  
 *  This class is an HLTProducer creating a flag for each HLT pattern satisfied
 *
 *
 *  \author Nicola De Filippis
 *
 */
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include<vector>
#include<string>

//
// class declaration
//

class HZZ4LeptonsHLTInfo : public edm::EDProducer {

  public:

    explicit HZZ4LeptonsHLTInfo(const edm::ParameterSet&);
    ~HZZ4LeptonsHLTInfo();

 private:
    virtual void produce(edm::Event& event, const edm::EventSetup& eventSetup);
    virtual void endJob();


    /// HLT TriggerResults EDProduct
    edm::EDGetTokenT<edm::TriggerResults> inputTag_;

    /// HLT trigger names
    edm::TriggerNames triggerNames_;

    /// number of HLT trigger paths requested in configuration
    unsigned int n_;
    bool firstevent_;

  // PG and FRC 06-07-11 try to reduce printout!
	bool debug;
	
    /// list of required HLT triggers by HLT name
    std::vector<std::string > HLTPathsByName_;
    /// list of required HLT triggers by HLT index
    std::vector<unsigned int> HLTPathsByIndex_;

};

#endif //HZZ4LeptonsHLTInfo_h
