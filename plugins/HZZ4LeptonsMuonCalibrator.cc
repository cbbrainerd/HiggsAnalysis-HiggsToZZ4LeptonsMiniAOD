/**\class HZZ4LeptonsMuonCalibrator.cc
 *
 * Original Author:  Nicola De Filippis 
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsMuonCalibrator.h"

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

#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"

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
HZZ4LeptonsMuonCalibrator::HZZ4LeptonsMuonCalibrator(const edm::ParameterSet& pset) {
  isData           = pset.getParameter<bool>("isData");
  muonLabel        = consumes<edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("muonCollection"));

  string alias;
  produces<pat::MuonCollection>(); 
  produces<edm::ValueMap<float> >("CorrPtError");
}


HZZ4LeptonsMuonCalibrator::~HZZ4LeptonsMuonCalibrator() {
 
}


//
// member functions
//
void HZZ4LeptonsMuonCalibrator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {



  // muons
  auto_ptr<pat::MuonCollection> Gmuon( new pat::MuonCollection );
  edm::Handle<edm::View<pat::Muon> > muons;
  edm::View<pat::Muon>::const_iterator mIter;    
  iEvent.getByToken(muonLabel, muons);

		 
  std::vector<float> pterror;
  pterror.clear();

  auto_ptr<edm::ValueMap<float> > CorrPtErrorMap(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler fillerCorrPtError(*CorrPtErrorMap);	

  double corrPt=0.,corrPtError=0.;
  double smearedPt=0., smearedPtError=0.;
    
  TLorentzVector p4;

  // Loop over muons
  for (mIter = muons->begin(); mIter != muons->end(); ++mIter ) {
    
    pat::Muon* calibmu = mIter->clone();
    
    smearedPt=mIter->pt();
    smearedPtError=mIter->muonBestTrack()->ptError();

    /*   if (mIter->muonBestTrackType() == 1 && mIter->pt()<=200.){

      if (isData){
	if (mIter->pt()>2.0 && fabs(mIter->eta())<2.4){
	  KalmanMuonCalibrator calibrator("DATA_80X_13TeV");
	  corrPt = calibrator.getCorrectedPt(mIter->pt(), mIter->eta(), mIter->phi(), mIter->charge());
	  corrPtError = corrPt * calibrator.getCorrectedError(corrPt, mIter->eta(), mIter->bestTrack()->ptError()/corrPt );
	  smearedPt=corrPt; // no smearing on data
	  smearedPtError=corrPtError; // no smearing on data
	}
      }
      else { // isMC - calibration from data + smearing
	KalmanMuonCalibrator calibrator("MC_80X_13TeV");
	corrPt = calibrator.getCorrectedPt(mIter->pt(), mIter->eta(), mIter->phi(), mIter->charge());
	corrPtError = corrPt * calibrator.getCorrectedError(corrPt, mIter->eta(), mIter->bestTrack()->ptError()/corrPt );
	// smearedPt = calibrator.smearForSync(corrPt, mIter->eta()); // for synchronization
	smearedPt = calibrator.smear(corrPt, mIter->eta());
	smearedPtError = corrPtError;
	//smearedPt * calibrator.getCorrectedErrorAfterSmearing(smearedPt, mIter->eta(), corrPtError /smearedPt );
      }
      
      }*/

    pterror.push_back(smearedPtError);
    // cout << "Muon pT= " << calibmu->pt() << " Corrected Muon pT= " << smearedPt << " and pT error= " << smearedPtError << " GlobalMuon = " << calibmu->isGlobalMuon() <<" number of matches= " << calibmu->numberOfMatches(reco::Muon::SegmentArbitration) << "Tracknonnull= " << calibmu->reco::Muon::innerTrack().isNonnull() << endl;

    p4.SetPtEtaPhiM(smearedPt, mIter->eta(), mIter->phi(), mIter->mass());
    calibmu->setP4(reco::Particle::PolarLorentzVector(p4.Pt(), p4.Eta(), p4.Phi(), mIter->mass()));
  
    //reco::Candidate::PolarLorentzVector p4Polar_;
    //p4Polar_ = reco::Candidate::PolarLorentzVector(smearedPt, mIter->eta(), mIter->phi(), mIter->mass());
    //calibmu->setP4(p4Polar_);
    
    Gmuon->push_back( *calibmu );
  }

  fillerCorrPtError.insert(muons,pterror.begin(),pterror.end());
  fillerCorrPtError.fill();
  
  const string iName = "";
  // iEvent.put( Gmuon, iName );
  // iEvent.put( CorrPtErrorMap, "CorrPtError");
  
  iEvent.put(std::make_unique<pat::MuonCollection>(*Gmuon), iName );
  iEvent.put(std::make_unique<edm::ValueMap<float>>(*CorrPtErrorMap), "CorrPtError");
  
}

