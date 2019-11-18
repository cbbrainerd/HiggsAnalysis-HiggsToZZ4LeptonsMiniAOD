// system include files
#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsHLTAnalysisFilter.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "FWCore/Framework/interface/FileBlock.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// namespaces
using namespace edm;
using namespace std;

// Constructor
HZZ4LeptonsHLTAnalysisFilter::HZZ4LeptonsHLTAnalysisFilter(const edm::ParameterSet& pset) {

  HLTInfoFired = consumes<std::vector<std::string> >(pset.getParameter<edm::InputTag>("HLTInfoFired"));
  for(auto &x : pset.getParameter<std::vector<std::string> >("pass_triggers")) {
    pass_triggers.insert(x);       
  }
  for(auto &x : pset.getParameter<std::vector<std::string> >("veto_triggers")) {
    veto_triggers.insert(x);
  }
  std::cout << "Trigger setup finished. Passing triggers:";
  for(auto &x : pass_triggers) std::cout << x << ",";
  std::cout << "\nVeto triggers:";
  for(auto &x : veto_triggers) std::cout << x << ",";
  std::cout << "\n";
}

#include <algorithm>
//Algorithm that searches two sorted containers for elements such that element1 is a prefix of element2 (i.e. HLT_Ele35_v in the first container matches HLT_Ele35_v1)
template <class ForwardIt>
std::pair<ForwardIt,ForwardIt> search_for_prefix(ForwardIt first_prefix,ForwardIt last_prefix,ForwardIt first_search,ForwardIt last_search) {
    for(;first_prefix!=last_prefix;++first_prefix) {
        first_search=std::lower_bound(first_search,last_search,*first_prefix);
        if(first_search==last_search) break;
        if((first_prefix->size() <= first_search->size()) && std::mismatch(first_prefix->begin(),first_prefix->end(),first_search->begin()).first==first_prefix->end()) return std::make_pair(first_prefix,first_search);
    }
    return std::make_pair(last_prefix,last_search);
}

bool check_trigger(const std::set<std::string>& to_check,const std::set<std::string>& fired) {
    return search_for_prefix(to_check.begin(),to_check.end(),fired.begin(),fired.end()).first!=to_check.end();
}

// Destructor
HZZ4LeptonsHLTAnalysisFilter::~HZZ4LeptonsHLTAnalysisFilter() {

}


// Filter Run Event
bool HZZ4LeptonsHLTAnalysisFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  edm::Handle<vector<std::string> > HLTfired_;
  iEvent.getByToken(HLTInfoFired,HLTfired_);
  std::set<std::string> HLTFired;
  for(auto const cand : *HLTfired_) {
    HLTFired.insert(cand);
  }
  return check_trigger(pass_triggers,HLTFired) && !check_trigger(veto_triggers,HLTFired);
/*  
  vector<string> HLTimported;
  string tmpstring="";
  
  for (vector<string>::const_iterator cand=HLTfired_->begin(); cand!=HLTfired_->end(); ++cand){
    unsigned int i=cand-HLTfired_->begin();
    HLTimported.push_back(cand->c_str());
    string newstr=HLTimported.at(i) + ":" + tmpstring;
    tmpstring=newstr;
  }
  
  cout << "HLTFiredString= " << tmpstring.c_str() << endl;

  char HLTPathsFired[20000];
  sprintf(HLTPathsFired,tmpstring.c_str());
  
  stringstream ss (stringstream::in | stringstream::out);
  ss << HLTPathsFired;
  TString hlt(ss.str());

  TString out = inputfileName;

  bool debug=true;
  cout << "Filename is= " << out.Data() << endl;
*/
/*
  if(
     (out.Contains("2017") && out.Contains("data") && !out.Contains("Moriond17")) && 
     (out.Contains("2017") && out.Contains("data") && !out.Contains("Spring16")) && 
     (out.Contains("2017") && out.Contains("data") && !out.Contains("Fall17")) ){
      
    if( out.Contains("DoubleEG")){
      
      if( debug ){ cout << "\n ** Step 2 (Trigger): 2017 DoubleElectron"<< endl ;
	
	cout << "This is HLT in data" << endl;
	cout<<" HLTPathsFired... "<<hlt<<endl;
      }
      
      if(
	 !hlt.Contains("HLT_Ele35_WPTight_Gsf_v") && // single-ele
	 !hlt.Contains("HLT_Ele38_WPTight_Gsf_v") &&
	 !hlt.Contains("HLT_Ele40_WPTight_Gsf_v") &&     // single-ele
	 
	 !hlt.Contains("HLT_IsoMu27") &&   // single-muon
	 	 
	 !hlt.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_") && // Di-Ele
	 !hlt.Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL") && //Di-Ele
	 !hlt.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL") && // Tri-Ele
	 
	 !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8") &&  //Di-Muon
	 !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8") && // Di-Muon
	 !hlt.Contains("HLT_TripleMu_12_10_5") && //Tri-Muon
	 !hlt.Contains("HLT_TripleMu_10_5_5_D2") && //Tri-Muon 
	 
	 !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") && //Mu-Ele
	 !hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
	 !hlt.Contains("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
	 !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
	 
	 !hlt.Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ") && //Di-Muon Ele
	 !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL") && //Muon-DiEle
	 !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ") //Muon-DiEle
	 ) {
	if( debug )cout << "Event not passing the HLT trigger paths" << endl;
	return false;
      }
      
      if(	 
	 hlt.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_") || // Di-Ele 
	 hlt.Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL") || // Double-Ele
	 hlt.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL") //Tri-Ele
	 ) {
	if( debug )cout << "Event passing the HLT trigger vetos for DoubleEG PD" << endl;
	return true;
      }
    }
    else if( out.Contains("DoubleMuon") ){
      
      if( debug ){ cout << "\n ** Step 2 (Trigger): 2017 DiMuon"<< endl ;
	
	cout << "This is HLT in data" << endl;
	cout<<" HLTPathsFired... "<<hlt<<endl;
      }
      
      if(
	 
       !hlt.Contains("HLT_Ele35_WPTight_Gsf_v") && // single-ele
       !hlt.Contains("HLT_Ele38_WPTight_Gsf_v") &&
       !hlt.Contains("HLT_Ele40_WPTight_Gsf_v") &&     // single-ele
       
       !hlt.Contains("HLT_IsoMu27_v") &&   // single-muon
       
       !hlt.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL") && // Di-Ele //+ HLT same name with _DZ
       !hlt.Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL") && //Di-Ele
       !hlt.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL") && // Tri-Ele
       
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8") &&  //Di-Muon
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8") && // Di-Muon
       !hlt.Contains("HLT_TripleMu_12_10_5") && //Tri-Muon
       !hlt.Contains("HLT_TripleMu_10_5_5_D") && //Tri-Muon 
       
       !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") && //Mu-Ele
       !hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       !hlt.Contains("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       
       !hlt.Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ") && //Di-Muon Ele
       !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL") && //Muon-DiEle
       !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ") //Muon-DiEle
	 
		) {
	if( debug )cout << "Event not passing the HLT trigger paths" << endl;
	return false;
      }
      
      if(
	 (
	  hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8") ||  //Di-Muon
	  hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8") || // Di-Muon
	  hlt.Contains("HLT_TripleMu_12_10_5") || //Tri-Muon
	  hlt.Contains("HLT_TripleMu_10_5_5_D") //Tri-Muon
	  ) &&
	 (
	  !hlt.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_") && // Di-Ele
	  !hlt.Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL") && //Di-Ele
	  !hlt.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL")  // Tri-Ele
	  )
	 ) {
	if( debug )cout << "Event passing the HLT trigger vetos for DoubleMuon PD" << endl;
	return true;
      }
    }  
    else if( out.Contains("MuonEG") ){
      if( debug ){ cout << "\n ** Step 2 (Trigger): 2017 MuEle"<< endl ;
	
	cout << "This is HLT in data" << endl;
	cout<<" HLTPathsFired... "<<hlt<<endl;
      }
      
      if(
	 
       !hlt.Contains("HLT_Ele35_WPTight_Gsf_v") && // single-ele
       !hlt.Contains("HLT_Ele38_WPTight_Gsf_v") &&
       !hlt.Contains("HLT_Ele40_WPTight_Gsf_v") &&     // single-ele
       
       !hlt.Contains("HLT_IsoMu27_v") &&   // single-muon
       
       !hlt.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL") && // Di-Ele //+ HLT same name with _DZ
       !hlt.Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL") && //Di-Ele
       !hlt.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL") && // Tri-Ele
       
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8") &&  //Di-Muon
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8") && // Di-Muon
       !hlt.Contains("HLT_TripleMu_12_10_5") && //Tri-Muon
       !hlt.Contains("HLT_TripleMu_10_5_5_D") && //Tri-Muon 
       
       !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") && //Mu-Ele
       !hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       !hlt.Contains("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       
       !hlt.Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ") && //Di-Muon Ele
       !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL") && //Muon-DiEle
       !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ") //Muon-DiEle
	 
	 ) {
	if( debug )cout << "Event not passing the HLT trigger paths" << endl;
	return false;	      
      }
      
      if(
	 (	  
	  hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") || //Mu-Ele //new 2017
	  hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") || //Mu-Ele
	  hlt.Contains("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") || //Mu-Ele //new 2017
	  hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") || //Mu-Ele //new 2017

	  hlt.Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ") || //Di-Muon Ele
	  hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL") || //Muon-DiEle 2017
	  hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ") //Muon-DiEle 2017	  
	  ) &&
	 (

	  !hlt.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_") && // Di-Ele
	  !hlt.Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL") && //Di-Ele
	  !hlt.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL") && // Tri-Ele
	  
	  !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8") &&  //Di-Muon
	  !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8") && // Di-Muon
	  !hlt.Contains("HLT_TripleMu_12_10_5") && //Tri-Muon
	  !hlt.Contains("HLT_TripleMu_10_5_5_D") //Tri-Muon
	  )
	 ) {
	   if( debug )cout << "Event passing the HLT trigger vetos for MuEG PD" << endl;
	return true;
      }
    }

    else if( out.Contains("SingleElectron") ){
      if( debug ){ cout << "\n ** Step 2 (Trigger): 2017 SingleElectron"<< endl ;
	
	cout << "This is HLT in data" << endl;
	cout<<" HLTPathsFired... "<<hlt<<endl;
      }
      
      if(	 
       !hlt.Contains("HLT_Ele35_WPTight_Gsf_v") && // single-ele
       !hlt.Contains("HLT_Ele38_WPTight_Gsf_v") &&
       !hlt.Contains("HLT_Ele40_WPTight_Gsf_v") &&     // single-ele
       
       !hlt.Contains("HLT_IsoMu27_v") &&   // single-muon
       
       !hlt.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL") && // Di-Ele //+ HLT same name with _DZ
       !hlt.Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL") && //Di-Ele
       !hlt.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL") && // Tri-Ele
       
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8") &&  //Di-Muon
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8") && // Di-Muon
       !hlt.Contains("HLT_TripleMu_12_10_5") && //Tri-Muon
       !hlt.Contains("HLT_TripleMu_10_5_5_D") && //Tri-Muon 
       
       !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") && //Mu-Ele
       !hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       !hlt.Contains("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       
       !hlt.Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ") && //Di-Muon Ele
       !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL") && //Muon-DiEle
       !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ") //Muon-DiEle	 
	 ) {
	if( debug )cout << "Event not passing the HLT trigger paths" << endl;
	return false;	      
      }
      
      if(
	 (
	  hlt.Contains("HLT_Ele35_WPTight_Gsf_v") ||    // single-ele
	  hlt.Contains("HLT_Ele38_WPTight_Gsf_v") ||
	  hlt.Contains("HLT_Ele40_WPTight_Gsf_v")  // single-ele	  
	  ) &&
	 (
	  
	 !hlt.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_") && // Di-Ele
	 !hlt.Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL") && //Di-Ele
	 !hlt.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL") && // Tri-Ele
	  
	 !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8") &&  //Di-Muon
	 !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8") && // Di-Muon
	 !hlt.Contains("HLT_TripleMu_12_10_5") && //Tri-Muon
	 !hlt.Contains("HLT_TripleMu_10_5_5_D") && //Tri-Muon
	 
	 !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") && //Mu-Ele
	 !hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
	 !hlt.Contains("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
	 !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
	  
	 !hlt.Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ") && //Di-Muon Ele
	 !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL") && //Muon-DiEle
	 !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ") //Muon-DiEle
	 
	  )
	 ){
	if( debug )cout << "Event passing the HLT trigger vetos for SingleElectron PD" << endl;
	return true;
      }
    }
    /////////   
    else if( out.Contains("SingleMuon") ){
      if( debug ){ cout << "\n ** Step 2 (Trigger): 2017 SingleElectron"<< endl ;
	
	cout << "This is HLT in data" << endl;
	cout<<" HLTPathsFired... "<<hlt<<endl;
      }
      
      if(	 
       !hlt.Contains("HLT_Ele35_WPTight_Gsf_v") && // single-ele
       !hlt.Contains("HLT_Ele38_WPTight_Gsf_v") &&
       !hlt.Contains("HLT_Ele40_WPTight_Gsf_v") &&     // single-ele
       
       !hlt.Contains("HLT_IsoMu27_v") &&   // single-muon
       
       !hlt.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL") && // Di-Ele //+ HLT same name with _DZ
       !hlt.Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL") && //Di-Ele
       !hlt.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL") && // Tri-Ele
       
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8") &&  //Di-Muon
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8") && // Di-Muon
       !hlt.Contains("HLT_TripleMu_12_10_5") && //Tri-Muon
       !hlt.Contains("HLT_TripleMu_10_5_5_D") && //Tri-Muon 
       
       !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") && //Mu-Ele
       !hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       !hlt.Contains("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       
       !hlt.Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ") && //Di-Muon Ele
       !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL") && //Muon-DiEle
       !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ") //Muon-DiEle	 
	 ) {
	if( debug )cout << "Event not passing the HLT trigger paths" << endl;
	return false;	      
      }
      
      if(
	 (
	  hlt.Contains("HLT_IsoMu27")   // single-muon
 
	  ) &&
	 (
	  
	 !hlt.Contains("HLT_Ele35_WPTight_Gsf_v") && // single-ele
	 !hlt.Contains("HLT_Ele38_WPTight_Gsf_v") &&
	 !hlt.Contains("HLT_Ele40_WPTight_Gsf_v") &&     // single-ele
	  
	 !hlt.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_") && // Di-Ele
	 !hlt.Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL") && //Di-Ele
	 !hlt.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL") && // Tri-Ele
	  
	 !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8") &&  //Di-Muon
	 !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8") && // Di-Muon
	 !hlt.Contains("HLT_TripleMu_12_10_5") && //Tri-Muon
	 !hlt.Contains("HLT_TripleMu_10_5_5_D") && //Tri-Muon 
	 
	 !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") && //Mu-Ele
	 !hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
	 !hlt.Contains("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
	 !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
	  
	 !hlt.Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ") && //Di-Muon Ele
	 !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL") && //Muon-DiEle
	 !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ") //Muon-DiEle
	 
	  )
	 ){
	if( debug )cout << "Event passing the HLT trigger vetos for SingleElectron PD" << endl;
	return true;
      }
    }    
  }
  else if( out.Contains("Spring15") ){
    if( debug ){ cout << "\n ** Step 2 (Trigger): "<< endl ;

      cout << "This is HLT in MC" << endl;
      cout<<" HLTPathsFired... "<<hlt<<endl;
    }
    
    if(
       !hlt.Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") &&    // di-electron trigger
       !hlt.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v") && // Triele
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") &&     // di-muon trigger
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") &&   // di-muon trigger
       !hlt.Contains("HLT_TripleMu_12_10_5_v") &&  // Trimuon
       !hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v")   && // MuEle
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") &&    // MuEle
       !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v") &&  // Mu-DiEle
       !hlt.Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v") && //DiMu-Ele
       !hlt.Contains("HLT_Ele27_WP85_Gsf_v") //Single-Ele
       ) {
      if( debug )cout << "Event not passing the HLT trigger paths" << endl;
      return false;
    }
    else return true; 
    
  }
  else if( out.Contains("Fall15")){
    if( debug ){ cout << "\n ** Step 2 (Trigger): "<< endl ;

      cout << "This is HLT in MC" << endl;
      cout<<" HLTPathsFired... "<<hlt<<endl;
    }

    if(
       !hlt.Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") &&    // di-electron trigger                       
       !hlt.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v") && // Triele                                          
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") &&     // di-muon trigger                               
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") &&   // di-muon trigger                               
       !hlt.Contains("HLT_TripleMu_12_10_5_v") &&  // Trimuon                                                                                                                         
       !hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v")   && // MuEle                                
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") &&    // MuEle                              
       !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v") &&  // Mu-DiEle                                            
       !hlt.Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v") && //DiMu-Ele 
       !hlt.Contains("HLT_Ele23_WPLoose_Gsf_v") //Single-Ele  
       ) {
      if( debug )cout << "Event not passing the HLT trigger paths" << endl;
      return false;
    }
    else return true;

  }
  else if(out.Contains("Fall17")){
    if( debug ){ cout << "\n ** Step 2 (Trigger): "<< endl ;
      
      cout << "This is HLT in MC" << endl;
      cout<<" HLTPathsFired... "<<hlt<<endl;
    }

    if(
       !hlt.Contains("HLT_Ele35_WPTight_Gsf_v") && // single-ele
       !hlt.Contains("HLT_Ele38_WPTight_Gsf_v") &&
       !hlt.Contains("HLT_Ele40_WPTight_Gsf_v") &&     // single-ele
       
       !hlt.Contains("HLT_IsoMu27_v") &&   // single-muon
       
       !hlt.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL") && // Di-Ele //+ HLT same name with _DZ
       !hlt.Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL") && //Di-Ele
       !hlt.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL") && // Tri-Ele
       
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8") &&  //Di-Muon
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8") && // Di-Muon
       !hlt.Contains("HLT_TripleMu_12_10_5") && //Tri-Muon
       !hlt.Contains("HLT_TripleMu_10_5_5_D") && //Tri-Muon 
       
       !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") && //Mu-Ele
       !hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       !hlt.Contains("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") && //Mu-Ele
       
       !hlt.Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ") && //Di-Muon Ele
       !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL") && //Muon-DiEle
       !hlt.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ") //Muon-DiEle	 
       ) {
      // if( debug )
      cout << "Event not passing the HLT trigger paths" << endl;
      return false;
    }
    else {
      cout << "Event passing the HLT trigger paths for MC" << endl;
      return true; }
  }//
  
  cout << " Event not recognized"<< endl;
  return false;
*/
}

void HZZ4LeptonsHLTAnalysisFilter::respondToOpenInputFile(edm::FileBlock const& fb) {
  inputfileName = fb.fileName();
  cout << "Input Filename is=" << inputfileName.c_str() << endl; 
}
