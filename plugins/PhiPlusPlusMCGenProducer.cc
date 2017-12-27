// -*- C++ -*-
//
// Package:    PhiPlusPlusMCGenProducer
// Class:      PhiPlusPlusMCGenProducer
// 
/**\class PhiPlusPlusMCGenProducer PhiPlusPlusMCGenProducer.cc liliana_code/PhiPlusPlusMCGenProducer/src/PhiPlusPlusMCGenProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Liliana Losurdo - Università di Bari 
//         Created:  Mon Dec  5 13:49:12 CET 2011
// $Id: PhiPlusPlusMCGenProducer.cc,v 1.1 2011/12/29 16:51:05 ndefilip Exp $
//
//

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/PhiPlusPlusMCGenProducer.h"

#include <memory>
#include <iostream>



//
// constructors and destructor
//
PhiPlusPlusMCGenProducer::PhiPlusPlusMCGenProducer(const edm::ParameterSet& iConfig) : 
  sourceLabel(iConfig.getParameter<edm::InputTag>("src"))
{
  debug=true;
  produces<vector<std::string> >("ExoticFired").setBranchAlias("ExoticFired");
  
}


PhiPlusPlusMCGenProducer::~PhiPlusPlusMCGenProducer()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PhiPlusPlusMCGenProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  auto_ptr<vector<std::string> > FiredExoticHiggs( new vector<std::string> );


  std::vector<int> dc; // will contain the Higgs decay products

  //if (iEvent.isRealData()) return data;

  //Get generator particles
  Handle<vector<GenParticle> >genCandidates;
  iEvent.getByLabel(sourceLabel,genCandidates);

  //Clear the previous state vector
  dc.clear();
  nevt++;
  d1=d2=0;
  dc1=dc2=dc3=dc4=999;

  //Iterate over all generator particles and write down the H++ and H-- decay products
  for (vector<GenParticle>::const_iterator mcIter=genCandidates->begin(); mcIter!=genCandidates->end(); ++mcIter){
    if ((abs(mcIter->pdgId())==9900041 || abs(mcIter->pdgId())==37) && mcIter->status()==3) {
      d1 = mcIter->daughterRef(0)->pdgId();
      d2 = mcIter->daughterRef(1)->pdgId();    

      dc.push_back(d1);
      dc.push_back(d2);
    }
  }
    
  // Sort them (positive particles to the front, negative to the back) cioè stamparle in ordine crescente secondo il loro pdgId
  sort(dc.begin(), dc.end());

  if (dc.size()==4){
    //if (debug) cout << "sorted dc content: " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
    
    dc1=dc[0]; dc2=dc[1]; dc3=dc[2]; dc4=dc[3];
    
    if (dc1 == -11 && dc2 == -11 && dc3 == 11 && dc4 == 11) {
      FiredExoticHiggs->push_back("eeee");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: eeee"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -13 && dc2 == -11 && dc3 == 11 && dc4 == 11 ) || ( dc1 == -11 && dc2 == -11 && dc3 == 11 && dc4 == 13 ) ) {
      FiredExoticHiggs->push_back("eeem");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: eeem"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -15 && dc2 == -11 && dc3 == 11 && dc4 == 11 ) || ( dc1 == -11 && dc2 == -11 && dc3 == 11 && dc4 == 15 ) ) {
      FiredExoticHiggs->push_back("eeet");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: eeet"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -13 && dc2 == -13 && dc3 == 11 && dc4 == 11 ) || ( dc1 == -11 && dc2 == -11 && dc3 == 13 && dc4 == 13 ) ) {
      FiredExoticHiggs->push_back("eemm");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: eemm"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -15 && dc2 == -13 && dc3 == 11 && dc4 == 11 ) || ( dc1 == -11 && dc2 == -11 && dc3 == 13 && dc4 == 15 ) ) {
      FiredExoticHiggs->push_back("eemt");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: eemm"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -15 && dc2 == -13 && dc3 == 11 && dc4 == 11 ) || ( dc1 == -11 && dc2 == -11 && dc3 == 13 && dc4 == 15 ) ) {
      FiredExoticHiggs->push_back("eemt");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: eemt"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -15 && dc2 == -15 && dc3 == 11 && dc4 == 11 ) || ( dc1 == -11 && dc2 == -11 && dc3 == 15 && dc4 == 15 ) ) {
      FiredExoticHiggs->push_back("eett");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: eett"<<endl;
      cout<<" "<<endl;
    }
    if ( dc1 == -13 && dc2 == -11 && dc3 == 11 && dc4 == 13 ) {
      FiredExoticHiggs->push_back("emem");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: emem"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -13 && dc2 == -11 && dc3 == 11 && dc4 == 15 ) || ( dc1 == -15 && dc2 == -11 && dc3 == 11 && dc4 == 13 ) ) {
      FiredExoticHiggs->push_back("emet");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: emet"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -13 && dc2 == -11 && dc3 == 13 && dc4 == 13 ) || ( dc1 == -13 && dc2 == -13 && dc3 == 11 && dc4 == 13 ) ) {
      FiredExoticHiggs->push_back("emmm");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: emmm"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -13 && dc2 == -11 && dc3 == 13 && dc4 == 15 ) || ( dc1 == -15 && dc2 == -13 && dc3 == 11 && dc4 == 13 ) ) {
      FiredExoticHiggs->push_back("emmt");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: emmt"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -13 && dc2 == -11 && dc3 == 15 && dc4 == 15 ) || ( dc1 == -15 && dc2 == -15 && dc3 == 11 && dc4 == 13 ) ) {
      FiredExoticHiggs->push_back("emtt");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: emtt"<<endl;
      cout<<" "<<endl;
    }
    if (dc1 == -15 && dc2 == -11 && dc3 == 11 && dc4 == 15) {
      FiredExoticHiggs->push_back("etet");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: etet"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -15 && dc2 == -11 && dc3 == 13 && dc4 == 13 ) || ( dc1 == -13 && dc2 == -13 && dc3 == 11 && dc4 == 15 ) ) {
      FiredExoticHiggs->push_back("etmm");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: etmm"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -15 && dc2 == -11 && dc3 == 13 && dc4 == 15 ) || ( dc1 == -15 && dc2 == -13 && dc3 == 11 && dc4 == 15 ) ) {
      FiredExoticHiggs->push_back("etmt");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: etmt"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -15 && dc2 == -11 && dc3 == 15 && dc4 == 15 ) || ( dc1 == -15 && dc2 == -15 && dc3 == 11 && dc4 == 15 ) ) {
      FiredExoticHiggs->push_back("ettt");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: ettt"<<endl;
      cout<<" "<<endl;
    }
    if (dc1 == -13 && dc2 == -13 && dc3 == 13 && dc4 == 13) {
      FiredExoticHiggs->push_back("mmmm");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: mmmm"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -13 && dc2 == -13 && dc3 == 13 && dc4 == 15 ) || ( dc1 == -15 && dc2 == -13 && dc3 == 13 && dc4 == 13 ) ) {
      FiredExoticHiggs->push_back("mmmt");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: mmmt"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -13 && dc2 == -13 && dc3 == 15 && dc4 == 15 ) || ( dc1 == -15 && dc2 == -15 && dc3 == 13 && dc4 == 13 ) ) {
      FiredExoticHiggs->push_back("mmtt");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: mmtt"<<endl;
      cout<<" "<<endl;
    }
    if (dc1 == -15 && dc2 == -13 && dc3 == 13 && dc4 == 15) {
      FiredExoticHiggs->push_back("mtmt");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: mtmt"<<endl;
      cout<<" "<<endl;
    }
    if ( ( dc1 == -15 && dc2 == -13 && dc3 == 15 && dc4 == 15 ) || ( dc1 == -15 && dc2 == -15 && dc3 == 13 && dc4 == 15 ) ) {
      FiredExoticHiggs->push_back("mttt");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: mttt"<<endl;
      cout<<" "<<endl;
    }
    if (dc1 == -15 && dc2 == -15 && dc3 == 15 && dc4 == 15) {
      FiredExoticHiggs->push_back("tttt");
      cout << " " << endl;
      cout << "H++H-- case -> 4-leptons " << endl;
      cout << "Four-leptons = " << dc[0] << "," << dc[1] << "," << dc[2] << "," << dc[3] << endl;
      cout << " " << endl;
      cout<<" final state: tttt"<<endl;
      cout<<" "<<endl;
    }

    bool trilepton=false;
    vector<int> dctri;
    dctri.clear();

    for (unsigned int i=0; i<dc.size(); i++){
      if ( abs(dc.at(i))==12 || abs(dc.at(i))==14 || abs(dc.at(i))==16) {
	trilepton=true;       
	continue;
      }
      dctri.push_back(dc.at(i));
    }
    
    // cout << "Tri-leptons = " << dctri.at(0) << ", " << dctri.at(1) << ", " << dctri.at(2) << endl;
    if (trilepton==true) {
      //cout << "Tri-leptons = " << dctri.at(0) << ", " << dctri.at(1) << ", " << dctri.at(2) << endl;
      cout << " " << endl;
      cout << "H++H- case -> 3 leptons " << endl;
      cout << "Tri-leptons = " << dctri.at(0) << ", " << dctri.at(1) << ", " << dctri.at(2) << endl;
      cout << " " << endl;
      dc1=dctri[0]; dc2=dctri[1]; dc3=dctri[2];

      if ( (dc1 == -11 && dc2 == -11 && dc3 == 11) || (dc1 == 11 && dc2 == 11 && dc3 ==-11)) {
	FiredExoticHiggs->push_back("eee");
	cout << " final state: eee" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -11 && dc2 == -11 && dc3 == 13) || (dc1 == 11 && dc2 == 11 && dc3 ==-13)) {
	FiredExoticHiggs->push_back("eem");
	cout << " final state: eem" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -11 && dc2 == -11 && dc3 == 15) || (dc1 == 11 && dc2 == 11 && dc3 ==-15)) {
	FiredExoticHiggs->push_back("eet");
	cout << " final state: eet" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -13 && dc2 == -13 && dc3 == 11) || (dc1 == 13 && dc2 == 13 && dc3 ==-11)) {
	FiredExoticHiggs->push_back("mme");
	cout << " final state: mme" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -13 && dc2 == -13 && dc3 == 13) || (dc1 == 13 && dc2 == 13 && dc3 ==-13)) {
	FiredExoticHiggs->push_back("mmm");
	cout << " final state: mmm" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -13 && dc2 == -13 && dc3 == 15) || (dc1 == 13 && dc2 == 13 && dc3 ==-15)) {
	FiredExoticHiggs->push_back("mmt");
	cout << " final state: mmt" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -15 && dc2 == -15 && dc3 == 11) || (dc1 == 15 && dc2 == 15 && dc3 ==-11)) {
	FiredExoticHiggs->push_back("tte");
	cout << " final state: tte" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -15 && dc2 == -15 && dc3 == 13) || (dc1 == 15 && dc2 == 15 && dc3 ==-13)) {
	FiredExoticHiggs->push_back("ttm");
	cout << " final state: ttm" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -15 && dc2 == -15 && dc3 == 15) || (dc1 == 15 && dc2 == 15 && dc3 ==-15)) {
	FiredExoticHiggs->push_back("ttt");
	cout << " final state: ttt" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -11 && dc2 == -13 && dc3 == 11) || (dc1 == 11 && dc2 == 13 && dc3 ==-11)) {
	FiredExoticHiggs->push_back("eme");
	cout << " final state: eme" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -11 && dc2 == -13 && dc3 == 13) || (dc1 == 11 && dc2 == 13 && dc3 ==-13)) {
	FiredExoticHiggs->push_back("emm");
	cout << " final state: emm" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -11 && dc2 == -13 && dc3 == 15) || (dc1 == 11 && dc2 == 13 && dc3 ==-15)) {
	FiredExoticHiggs->push_back("emt");
	cout << " final state: emt" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -11 && dc2 == -15 && dc3 == 11) || (dc1 == 11 && dc2 == 15 && dc3 ==-11)) {
	FiredExoticHiggs->push_back("ete");
	cout << " final state: ete" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -11 && dc2 == -15 && dc3 == 13) || (dc1 == 11 && dc2 == 15 && dc3 ==-13)) {
	FiredExoticHiggs->push_back("etm");
	cout << " final state: etm" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -11 && dc2 == -15 && dc3 == 15) || (dc1 == 11 && dc2 == 15 && dc3 ==-15)) {
	FiredExoticHiggs->push_back("ett");
	cout << " final state: ett" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -13 && dc2 == -15 && dc3 == 11) || (dc1 == 13 && dc2 == 15 && dc3 ==-11)) {
	FiredExoticHiggs->push_back("mte");
	cout << "final state: mte" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -13 && dc2 == -15 && dc3 == 13) || (dc1 == 13 && dc2 == 15 && dc3 ==-13)) {
	FiredExoticHiggs->push_back("mtm");
	cout << "final state: mtm" << endl;
	cout << " " << endl;
      }
      if ( (dc1 == -13 && dc2 == -15 && dc3 == 15) || (dc1 == 13 && dc2 == 15 && dc3 ==-15)) {
	FiredExoticHiggs->push_back("mtt");
	cout << " final state: mtt" << endl;
	cout << " " << endl;
      }      
    }
   
    
    
    //if (modpoint > 4 || modpoint < 1) return -1;
    //   if (modpoint == 1) {
    //     // BR=1/6 model point
    //     ee=1./6;
    //     em=1./6;
    //     et=1./6;
    //     mm=1./6;
    //     mt=1./6;
    //     tt=1./6;
    //   }
    //   if (modpoint == 2) {
    //     // Normal hierarchy
    //     ee=0.0032;
    //     em=0.0064;
    //     et=0.0064;
    //     mm=0.3017;
    //     mt=0.3807;
    //     tt=0.3017;
    //   }
    //   if (modpoint == 3) {
    //     // Inverse hierarchy
    //     ee=0.4974;
    //     em=0;
    //     et=0;
    //     mm=0.1256;
    //     mt=0.2513;
    //     tt=0.1256;
    //   }
    //   if (modpoint == 4) {
    //     // Degenerate state
    //     ee=0.34;
    //     em=0;
    //     et=0;
    //     mm=0.3299;
    //     mt=0.0002;
    //     tt=0.3299;
    //   }
    //   switch (fstate) {
    //   case 1:
    //     return ee*ee;
    //   case 2:
    //     return 2*ee*em;
    //   case 3:
    //     return 2*ee*et;
    //   case 4:
    //     return 2*ee*mm;
    //   case 5:
    //     return 2*ee*mt;
    //   case 6:
    //     return 2*ee*tt;
    //   case 7:
    //     return em*em;
    //   case 8:
    //     return 2*em*et;
    //   case 9:
    //     return 2*em*mm;
    //   case 10:
    //     return 2*em*mt;
    //   case 11:
    //     return 2*em*tt;
    //   case 12:
    //     return et*et;
    //   case 13:
    //     return 2*et*mm;
    //   case 14:
    //     return 2*et*mt;
    //   case 15: 
    //     return 2*et*tt;
    //   case 16:
    //     return mm*mm;
    //   case 17:
    //     return 2*mm*mt;
    //   case 18:
    //     return 2*mm*tt;
    //   case 19:
    //     return mt*mt;
    //   case 20:
    //     return 2*mt*tt;
    //   case 21:
    //     return tt*tt;
    //   default:
    //     return 1;
    //   }
    //   return -1;
 
  }
  
  iEvent.put(FiredExoticHiggs,"ExoticFired");     
  

}


// ------------ method called once each job just before starting event loop  ------------
void 
PhiPlusPlusMCGenProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhiPlusPlusMCGenProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
PhiPlusPlusMCGenProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PhiPlusPlusMCGenProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PhiPlusPlusMCGenProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PhiPlusPlusMCGenProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}


