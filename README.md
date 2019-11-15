# HiggsAnalysis-HiggsToZZ4LeptonsMiniAOD 
Package for  H->ZZ->4l anaysis for Run2 miniAOD  

```
#!/bin/bash
CMSSW_VERSION="CMSSW_10_2_15"
cmsrel "$CMSSW_VERSION"
cd "$CMSSW_VERSION"/src
eval `scramv1 runtime -sh`
git cms-init
git clone https://github.com/cbbrainerd/HiggsAnalysis-HiggsToZZ4LeptonsMiniAOD HiggsAnalysis/HiggsToZZ4Leptons -b 2018
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
git cms-merge-topic cms-egamma:EgammaPostRecoTools
git cms-merge-topic cms-egamma:EgammaID_949
pushd ZZMatrixElement
bash setup.sh
popd
git clone https://github.com/mkovac/MuonMVAReader.git
git clone https://github.com/usarica/MelaAnalytics.git
#Download Rochester corrections from https://twiki.cern.ch/twiki/bin/view/CMS/RochcorMuon
```
