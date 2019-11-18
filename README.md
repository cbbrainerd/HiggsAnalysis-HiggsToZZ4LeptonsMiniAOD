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
git cms-merge-topic mkovac:Electron_XGBoost_MVA_2016_and_2018_CMSSW_10_2_15
#git cms-merge-topic cms-egamma:EgammaID_949
pushd ZZMatrixElement
bash setup.sh
popd
git cms-addpkg EgammaAnalysis/ElectronTools
rm EgammaAnalysis/ElectronTools/data -rf
git clone git@github.com:cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
git clone https://github.com/mkovac/MuonMVAReader.git
git clone https://github.com/usarica/MelaAnalytics.git
