#ifndef HZZ4LeptonsZZMassErrors_h
#define HZZ4LeptonsZZMassErrors_h

/** \class HZZ4LeptonsZZMassErrors
 *
 *  No description available.
 *
 *  $Date: 2012/10/05 10:25:22 $
 *  $Revision: 1.1 $
 *  \author N. Amapane - CERN
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/PatCandidates/interface/Electron.h"


class HZZ4LeptonsZZMassErrors {
public:
  /// Constructor
  HZZ4LeptonsZZMassErrors(const edm::EventSetup& eventSetup);

  /// Destructor
  virtual ~HZZ4LeptonsZZMassErrors(){};
  
  // Operations
  float massError(const reco::Candidate* Z1Lp, 
		  const reco::Candidate* Z1Lm, 
		  const reco::Candidate* Z2Lp, 
		  const reco::Candidate* Z2Lm);

  double pError(const pat::Electron* ele);

 private:
  edm::ESHandle<MagneticField> magfield;
  
};
#endif

