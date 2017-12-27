/* \class HZZ4LeptonsCommonOfflineSelectionFilter
 *
 *
 * Tight isolation for electron and muons
 *
 * author:     Nicola De Filippis    - LLR-Ecole Polytechnique
 * modific. R.Casagrande, J.Zablocki - Purdue University
 */


// system include files
#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsCommonOfflineSelectionFilter.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <vector>

// ROOT histo
#include <TH1.h>

// namespaces
using namespace edm;
using namespace std;
using namespace reco;

// Constructor
HZZ4LeptonsCommonOfflineSelectionFilter::HZZ4LeptonsCommonOfflineSelectionFilter(const edm::ParameterSet& pset) {


  // Deacay Channel
  decaychannel     = pset.getParameter<std::string>("decaychannel");
  offselinst       = pset.getParameter<std::string>("offselinst");
  offseltags       = pset.getParameter<std::vector<std::string> >("offseltags");

  // Off selection file
  offSelectFileName= pset.getParameter<std::string>("offSelectFileName");

  //
  //  rootFileName_ = pset.getParameter<string>("rootFileName");
  // theFile_ = new TFile(rootFileName_.c_str(),"RECREATE");
//   locdir="OfflineselectionHistos";
//   TDirectory *savdir = gDirectory;
//   TDirectory *adir = savdir->mkdir(locdir);
//   adir->cd();

  nPresel   = 0;
  nTightEle = 0;
  nTightMu  = 0;
  nVert     = 0;
  vcounter.clear();
  vcounter.resize(offseltags.size());	

}


// Destructor
HZZ4LeptonsCommonOfflineSelectionFilter::~HZZ4LeptonsCommonOfflineSelectionFilter() {

  ofstream offsel_file;
  
  offsel_file.open( offSelectFileName.c_str() );
  
  offsel_file << "*********************** "             << std::endl;
  offsel_file << "Offline isolation efficiency "        << std::endl;
  offsel_file << "*********************** "             << std::endl;
  offsel_file <<   "nPresel           : " << nPresel    << std::endl;
  
  int i=-1; 
  if (vcounter.size()) {
    if ( decaychannel=="2e2mu" || decaychannel=="4e" ){
      i++;
      offsel_file << "nTightEle         : " << vcounter.at(i) << std::endl;
    }
    if ( decaychannel=="2e2mu" || decaychannel=="4mu" ) {
      i++;
      offsel_file << "nTightMu          : " << vcounter.at(i) << std::endl;
    }
    i++;
    offsel_file   << "nVert             : " << vcounter.at(i) << std::endl;
  }
  offsel_file.close();
  
}


// Filter event (event offlineselection)
bool HZZ4LeptonsCommonOfflineSelectionFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  
    
  // input number of events
  nPresel++;
  
  
  for (size_t i=0;i<offseltags.size();i++){
    string str1 = decaychannel.c_str();
    string str2 = offseltags.at(i).c_str();
    str1.append(str2);
    Handle<bool> OffselHandle;
    iEvent.getByLabel(offselinst,str1,OffselHandle);
    cout << "String=" << str1 << "  and value=" << *OffselHandle.product() <<endl;
    if ( *OffselHandle.product()==1 ) {
      vcounter.at(i)++;
    }
    else {
      return false;
    }    
  }
  
  return true;
  
}

void HZZ4LeptonsCommonOfflineSelectionFilter::endJob(){

//   theFile_->cd(locdir);
  
//   TH1I *Eff=new TH1I("","",1,0,1);
//   Eff->SetName("Presel");
//   Eff->SetBinContent(1,nPresel);
//   //Eff->Write();
//   delete Eff;

//   for (size_t i=0;i<offseltags.size();i++){
//     TH1I *Eff=new TH1I("","",1,0,1);
//     Eff->SetName(offseltags.at(i).c_str());
//     Eff->SetBinContent(1,vcounter.at(i));
//     //Eff->Write();
//     //cout << "Counter= " << vcounter.at(i) << endl;
//     delete Eff;
//   }

//   theFile_->Write() ;
//   theFile_->Close();


}

