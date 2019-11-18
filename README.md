#!/bin/bash
set -e
. /cvmfs/cms.cern.ch/cmsset_default.sh
CMSSW_VERSION="10_2_15"
PREFIX="HZZ"
scram p -n "${PREFIX}_${CMSSW_VERSION}" CMSSW "CMSSW_${CMSSW_VERSION}"
cd "${PREFIX}_${CMSSW_VERSION}"/src
eval `scramv1 runtime -sh`
git cms-init
git clone https://github.com/cbbrainerd/HiggsAnalysis-HiggsToZZ4LeptonsMiniAOD HiggsAnalysis/HiggsToZZ4Leptons -b 2018
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
git cms-merge-topic cms-egamma:EgammaPostRecoTools
#git cms-merge-topic cms-egamma:EgammaID_949
pushd ZZMatrixElement
bash setup.sh
popd
git clone https://github.com/mkovac/MuonMVAReader.git
git clone https://github.com/usarica/MelaAnalytics.git
