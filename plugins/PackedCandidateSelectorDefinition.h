#ifndef CommonTools_ParticleFlow_PackedCandidateSelectorDefinition
#define CommonTools_ParticleFlow_PackedCandidateSelectorDefinition

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "boost/iterator/transform_iterator.hpp"
#include <functional>

namespace pf2pat {

  class PackedCandidateSelectorDefinition {

  public:
    typedef pat::PackedCandidateCollection collection;
    typedef edm::Handle< collection > HandleToCollection;
    typedef std::vector<pat::PackedCandidate>  container;
    
    struct Pointer : public std::unary_function<pat::PackedCandidate,const pat::PackedCandidate *> { 
      const pat::PackedCandidate * operator()(const pat::PackedCandidate &c) const { return &c; } 
    };
    
    typedef boost::transform_iterator<Pointer,container::const_iterator> const_iterator;
    
    PackedCandidateSelectorDefinition () {}
    
    const_iterator begin() const { return const_iterator(selected_.begin()); }
    
    const_iterator end() const { return const_iterator(selected_.end()); }
    
    size_t size() const { return selected_.size(); }
    
    const container& selected() const {return selected_;}
    
  protected:
    container selected_;
  };
}

#endif
