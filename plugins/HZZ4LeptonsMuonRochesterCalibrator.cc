/**\class HZZ4LeptonsMuonRochesterCalibrator.cc
 *
 * Original Author:  Nicola De Filippis 
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsMuonRochesterCalibrator.h"

#include "FWCore/Framework/interface/ESHandle.h"

// Muons:
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/ValueMap.h"


#include "RoccoR.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

// system include files
#include <Math/VectorUtil.h>
#include <memory>
#include <vector>


using namespace edm;
using namespace std;
using namespace reco;
using namespace math;

// constructor
HZZ4LeptonsMuonRochesterCalibrator::HZZ4LeptonsMuonRochesterCalibrator(const edm::ParameterSet& pset) {
  isData           = pset.getParameter<bool>("isData");
  muonLabel        = consumes<edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("muonCollection"));
  MCTruth               = pset.getParameter<bool>("MCTruth");
 
 // Gen Match information
 
  goodMuonMCMatch_      = consumes<edm::Association<std::vector<reco::GenParticle> > >(pset.getParameter<edm::InputTag>("goodMuonMCMatch"));
  myMuons_              = consumes<reco::CandidateCollection>(pset.getParameter<edm::InputTag>("myMuons"));
 
  string alias;
  produces<pat::MuonCollection>(); 
  produces<edm::ValueMap<float> >("CorrPtError");
}


HZZ4LeptonsMuonRochesterCalibrator::~HZZ4LeptonsMuonRochesterCalibrator() {
  
}

//
// member functions
//
void HZZ4LeptonsMuonRochesterCalibrator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  // muons
  auto_ptr<pat::MuonCollection> Gmuon( new pat::MuonCollection );
 
  edm::Handle<edm::View<pat::Muon> > muons;
  edm::View<pat::Muon>::const_iterator mIter;    
  iEvent.getByToken(muonLabel, muons);
  
  std::vector<float> pterror;
  pterror.clear();
  
  auto_ptr<edm::ValueMap<float> > CorrPtErrorMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerCorrPtError(*CorrPtErrorMap);	
  
  
  // To add matching informations Reham
  // Matching

    edm::Handle<edm::Association<vector<reco::GenParticle> > > GenParticlesMatchMu;
    iEvent.getByToken(goodMuonMCMatch_, GenParticlesMatchMu);
 
    edm::Handle<reco::CandidateCollection > CollMu;
    iEvent.getByToken(myMuons_, CollMu);
  
  //  if (GenParticlesMatchMu.isValid()){
  // // cout << endl<< "Muons:"<<endl<<"The reco collection to be matched has size= " <<  CollMu->size() << endl;
  // // cout << "The matched map collection has size= " <<  GenParticlesMatchMu->size() << endl;
  //  }
  
    int jj=0;

  // Loop over muons
  for (mIter = muons->begin(); mIter != muons->end(); ++mIter ) {

    edm::Ref<edm::View<pat::Muon>>muref(muons,jj);
    
    pat::Muon* calibmu = mIter->clone();

    //
    vector<double> vcorrPt, vcorrPtError;
    double oldPt=0., oldPtError=0.;
    double  smearedPt=0., smearedPtError=0.;
    double correction=1, correction_error=0;
    TLorentzVector p4;
    int nl;
    TRandom3* rgen_;
    rgen_ = new TRandom3(0);
    double u = rgen_->Rndm();
    RoccoR* calibrator; //muon calibrator
    
    bool Match = false;
    double gen_Mu_pt= 0.;

    // cout<<"######### Muon Rochester correction ############"<<endl;  
    // cout<<"PAT Muon PT = "<<mIter->pt()<<endl;
    
    oldPt=mIter->pt();
    oldPtError=mIter->muonBestTrack()->ptError();
  
    cout<<"Muon PT = "<<oldPt<<endl;
    
    // To add matching informations Reham
      // Matching
      if (MCTruth==true){
      // 	int i=0;
      // 	for ( reco::CandidateCollection::const_iterator hIter=CollMu->begin(); hIter!= CollMu->end(); ++hIter ){
      // 	   cout << "Reco Muon with pT= " << hIter->pt() << " and mass="<< hIter->mass()<< endl;

      // 	  if (fabs(hIter->pt()- mIter->pt())<0.01){
      // 	    cout<<"Match found "<<endl;
      // 	    i=hIter-(CollMu->begin());
      // 	    cout<<"i = "<<i<<endl;
      // 	    CandidateRef Ref( CollMu, i );
	             edm::Ref<std::vector<reco::GenParticle> > genrefMu = (*GenParticlesMatchMu)[muref];
		     if (!genrefMu.isNull()){
		       cout<<"Matched found "<<endl;
		       cout << "GenMuon with pT= " << genrefMu->p4().pt() << " and mass="<< genrefMu->p4().mass()<< endl;
		       Match = true;
		       gen_Mu_pt = genrefMu->p4().pt();
		       cout<<"Calibrator Match Gen Muon PT = "<<gen_Mu_pt<<endl;
		     } 
		     else {
       	      cout << "Calibrator There is no reference to a genMuon" << endl;
       	      Match = false;
		     }
      // 	  } //if  
      // 	}//for   
      }  //end of matching


      // cout<<"######## Muon correction ########"<<endl;
 
    edm::FileInPath corrPath("roccor_Run2_v2/data/RoccoR2017.txt");
    calibrator = new RoccoR(corrPath.fullPath());

    // cout<<"open the txt file for muon corrections"<<endl;
     
     cout<<" muon pt = "<<mIter->pt()<<" and best track type = "<<mIter->muonBestTrackType()<<endl;
     
      if (mIter->muonBestTrackType() == 1 && mIter->pt()<=200.){

	nl =  mIter->track()->hitPattern().trackerLayersWithMeasurement();

    	if (isData && nl>5 ){ //on data correction oly

	   cout<<"Data Muon correction "<<endl;

    	  if (mIter->pt()>2.0 && fabs(mIter->eta())<2.4){

	    //RoccoR Calibrator;
	    correction = calibrator->kScaleDT( mIter->charge(), mIter->pt(), mIter->eta(), mIter->phi(), 0 ,0); 
	    correction_error = calibrator->kScaleDTerror( mIter->charge(), mIter->pt(), mIter->eta(), mIter->phi()); }//end eta > 2.4
	    else {// with pt<2 or not in the acceptance
	       correction =1;
	       correction_error =0;}//end else
	}//end data
   
	 	else if(!isData && nl>5 ){ // isMC - calibration from data + smearing

		   cout<<"Calibrator This is MC Muon correction"<<endl;

		    if (  Match == true){// when matched gen muon available

		       cout<<"Calibrator MC and matched gen muon "<<endl;
		      cout<<"Calibrator gen muon pt = "<<gen_Mu_pt<<endl;

	     correction = calibrator->kSpreadMC( mIter->charge(), mIter->pt(), mIter->eta(), mIter->phi(), gen_Mu_pt , 0 ,0 ); 
	    correction_error = calibrator->kSpreadMCerror( mIter->charge(), mIter->pt(), mIter->eta(), mIter->phi(), gen_Mu_pt);
	    cout<<"Calibrator correction = "<<correction<<endl;
	     }

		    else { //when matched gen muon not available

		       cout<<"Calibrator MC and No matched gen muon "<<endl;		 

		    correction = calibrator->kSmearMC( mIter->charge(), mIter->pt(), mIter->eta(), mIter->phi(), nl , u, 0 ,0 ); 
		    correction_error = calibrator->kSmearMCerror( mIter->charge(), mIter->pt(), mIter->eta(), mIter->phi(), nl ,u );

	    cout<<"Calibrator correction = "<<correction<<endl;}

		}//ends if !isdata
      }//end if tracker , pt<200
   
    smearedPt = oldPt*correction;
    smearedPtError = oldPtError+correction;

    cout<<"In calibrator Old Muon Pt = "<<oldPt<<" , correction = "<<correction<<" smeared muon pt = "<<smearedPt<<endl;
     cout<<"Old pt error (from muon best track pt error) = "<<oldPtError<<", correction_error = "<<correction_error<<"smeared pt Error = "<< smearedPtError<<"+/- "<<correction_error<<endl;
     
    pterror.push_back(smearedPtError);    
    p4.SetPtEtaPhiM(smearedPt, mIter->eta(), mIter->phi(), mIter->mass());
    calibmu->setP4(reco::Particle::PolarLorentzVector(p4.Pt(), p4.Eta(), p4.Phi(), mIter->mass())); 
    Gmuon->push_back( *calibmu );
    jj++;
  }
  
  fillerCorrPtError.insert(muons,pterror.begin(),pterror.end());
  fillerCorrPtError.fill();
  
  const string iName = ""; 
  iEvent.put(std::make_unique<pat::MuonCollection>(*Gmuon), iName );
  iEvent.put(std::make_unique<edm::ValueMap<float>>(*CorrPtErrorMap), "CorrPtError");
  
}

