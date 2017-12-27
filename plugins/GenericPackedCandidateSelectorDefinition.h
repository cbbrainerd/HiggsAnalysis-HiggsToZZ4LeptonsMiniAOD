#ifndef CommonTools_ParticleFlow_GenericPackedCandidateSelectorDefinition
#define CommonTools_ParticleFlow_GenericPackedCandidateSelectorDefinition

/**
   \class    pf2pat::GenericPFCandidateSelectorDefinition GenericPFCandidateSelectorDefinition.h "CommonTools/ParticleFlow/interface/GenericPFCandidateSelectorDefinition.h"
   \brief    Selects PFCandidates basing on cuts provided with string cut parser

   \author   Giovanni Petrucciani
   \version  $Id: GenericPFCandidateSelectorDefinition.h,v 1.1 2011/01/28 20:56:44 srappocc Exp $
*/

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "PackedCandidateSelectorDefinition.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

namespace pf2pat {

  struct GenericPackedCandidateSelectorDefinition : public PackedCandidateSelectorDefinition {

    GenericPackedCandidateSelectorDefinition ( const edm::ParameterSet & cfg, edm::ConsumesCollector && iC ) :
      selector_( cfg.getParameter< std::string >( "cut" ) ) { }

    void select( const HandleToCollection & hc,
		 const edm::Event & e,
		 const edm::EventSetup& s) {
      selected_.clear();

      unsigned key=0;
      for( collection::const_iterator pfc = hc->begin();
	   pfc != hc->end(); ++pfc, ++key) {

	if( selector_(*pfc) ) {
	  selected_.push_back( pat::PackedCandidate(*pfc) );
//	  pat::PackedCandidatePtr ptrToMother( hc, key );
//	  selected_.back().setSourceCandidatePtr( ptrToMother );

	}
      }
    }

    private:
    StringCutObjectSelector<pat::PackedCandidate> selector_;
  };
}

#endif
