#ifndef SimpleCounter_h
#define SimpleCounter_h

/** \class SimpleCounter
 *
 * Original Author:  Nicola De Filippis
 *
 */

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/ESHandle.h"


class SimpleCounter : public edm::EDAnalyzer {
public:

        SimpleCounter(const edm::ParameterSet&);
	~SimpleCounter();

	void analyze(const edm::Event& e, const edm::EventSetup& c);
	void beginJob();
	void endJob();
	
private:
	
	int counter;
		
	};
	
#endif
