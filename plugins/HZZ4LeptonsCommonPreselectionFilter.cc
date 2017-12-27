/* \class HZZ4LeptonsPreselectionFilter
 *
 *
 * H->ZZ->4l analysis preselection:
 * skim input 
 * electron and muon selection
 * m_ll, m4l constraints
 * loose isolation on electrons and muons
 *
 * author:     Nicola De Filippis   - LLR-Ecole Polytechnique
 */


// system include files
#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsCommonPreselectionFilter.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <vector>

// namespaces
using namespace edm;
using namespace std;
using namespace reco;

// ROOT histo
#include <TH1.h>

// Constructor
HZZ4LeptonsCommonPreselectionFilter::HZZ4LeptonsCommonPreselectionFilter(const edm::ParameterSet& pset) {


  // Decay Channel
  decaychannel     = pset.getParameter<std::string>("decaychannel");
  preselinst       = pset.getParameter<std::string>("preselinst");
  preseltags       = pset.getParameter<std::vector<std::string> >("preseltags");

  // Preselection file
  preSelectFileName= pset.getParameter<std::string>("preSelectFileName");

  //
  //  rootFileName_ = pset.getParameter<string>("rootFileName"); 
  //theFile_ = new TFile(rootFileName_.c_str(), "RECREATE");
//   locdir="PreselectionHistos";
//   TDirectory *savdir = gDirectory;
//   TDirectory *adir = savdir->mkdir(locdir);
//   adir->cd();

  nSkim          = 0;
  nElec          = 0;
  nMu            = 0;
  nZEE           = 0;
  nZMM           = 0;
  nH4leptons     = 0;
  nLooseIsolEle = 0;
  nLooseIsolMu  = 0;
  vcounter.clear();
  vcounter.resize(preseltags.size());
}


// Destructor
HZZ4LeptonsCommonPreselectionFilter::~HZZ4LeptonsCommonPreselectionFilter() {

  ofstream presel_file;

  presel_file.open( preSelectFileName.c_str() );
                                                     
  presel_file   << "*********************** "             << std::endl;
  presel_file   << "Preselection efficiency "             << std::endl;
  presel_file   << "*********************** "             << std::endl;
  presel_file   << "nSkim           : " << nSkim          << std::endl;
					      
  int i=-1;

  if (vcounter.size()) {
    if ( decaychannel=="2e2mu" || decaychannel=="4e" ){
      i++;  
      presel_file << "nElec           : " << vcounter.at(i) << std::endl;
    }
    if ( decaychannel=="2e2mu" || decaychannel=="4mu" ) {
      i++;
      presel_file << "nMuon           : " << vcounter.at(i) << std::endl;
    }
    if ( decaychannel=="2e2mu" || decaychannel=="4e" ) {
      i++;
      presel_file << "Z->EE           : " << vcounter.at(i) << std::endl;
    }
    if ( decaychannel=="2e2mu" || decaychannel=="4mu" ) {
      i++;
      presel_file << "Z->MuMu         : " << vcounter.at(i) << std::endl;
    }
    i++;
    presel_file   << "H->ZZ           : " << vcounter.at(i) << std::endl;
    
    if ( decaychannel=="2e2mu" || decaychannel=="4e" ) {
      i++;
      presel_file << "loose IsolEle   : " << vcounter.at(i) << std::endl;
    }
    if ( decaychannel=="2e2mu" || decaychannel=="4mu" ) {
      i++;
      presel_file << "loose IsolMu    : " << vcounter.at(i) << std::endl;
    }    
  }
  presel_file.close();

}


// Filter event (event preselection)
bool HZZ4LeptonsCommonPreselectionFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  // input number of events
  nSkim++;

  for (size_t i=0;i<preseltags.size();i++){
    string str1 = decaychannel.c_str();
    string str2 = preseltags.at(i).c_str();
    str1.append(str2);
    Handle<bool> PreselHandle;
    iEvent.getByLabel(preselinst.c_str(),str1,PreselHandle);
    cout << "String=" << str1 << "  and value=" << *PreselHandle.product() <<endl;
    if ( *PreselHandle.product()==1 ) {
      vcounter.at(i)++;
    }
    else {
      return false;
    }    
  }

  //for (size_t i=0;i<preseltags.size();i++){
  //  cout << "Counter= " << vcounter.at(i) << endl;
  //}
  
  return true;
  

}

void HZZ4LeptonsCommonPreselectionFilter::endJob(){

  //theFile_->cd(locdir);

 //  Eff=new TH1I("","",1,0,1);
//   Eff->SetName("Skim");
//   Eff->SetBinContent(1,nSkim);
//   Eff->Write();
//   delete Eff;

//   for (size_t i=0;i<preseltags.size();i++){
//     Eff=new TH1I("","",1,0,1);
//     Eff->SetName(preseltags.at(i).c_str());
//     Eff->SetBinContent(1,vcounter.at(i));
//     Eff->Write();
//     //cout << "Counter= " << vcounter.at(i) << endl;
//     delete Eff;
//   }
  
  //theFile_->Write() ;
  //theFile_->Close();
 
}
