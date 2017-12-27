#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "GenericPackedCandidateSelectorDefinition.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

typedef ObjectSelector<pf2pat::GenericPackedCandidateSelectorDefinition> GenericPackedCandidateSelector;

DEFINE_FWK_MODULE(GenericPackedCandidateSelector);
