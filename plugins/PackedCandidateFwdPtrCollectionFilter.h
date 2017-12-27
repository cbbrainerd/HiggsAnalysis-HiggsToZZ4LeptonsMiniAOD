#include "CommonTools/UtilAlgos/interface/FwdPtrCollectionFilter.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/PdgIdSelector.h"
#include "CommonTools/ParticleFlow/interface/PackedCandidateWithSrcPtrFactory.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

typedef edm::FwdPtrCollectionFilter< pat::PackedCandidate, 
                                     reco::StringCutObjectSelectorHandler<pat::PackedCandidate,false>, 
                                     reco::PackedCandidateWithSrcPtrFactory >  PackedCandidateFwdPtrCollectionStringFilter;
typedef edm::FwdPtrCollectionFilter< pat::PackedCandidate, reco::PdgIdSelectorHandler, 
                                     reco::PackedCandidateWithSrcPtrFactory >  PackedCandidateFwdPtrCollectionPdgIdFilter;

