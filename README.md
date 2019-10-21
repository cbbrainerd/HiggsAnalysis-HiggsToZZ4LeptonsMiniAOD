# HiggsAnalysis-HiggsToZZ4LeptonsMiniAOD 
Package for  H->ZZ->4l anaysis for Run2 miniAOD  

```
cmsrel CMSSW_9_4_10
cd CMSSW_9_4_10/src
cmsenv
git cms-init
git clone https://github.com/cbbrainerd/HiggsAnalysis-HiggsToZZ4LeptonsMiniAOD HiggsAnalysis/HiggsToZZ4Leptons -b synchronization
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
git cms-merge-topic cms-egamma:EgammaPostRecoTools
git cms-merge-topic cms-egamma:EgammaID_949
pushd ZZMatrixElement
bash setup.sh
popd
#Download Rochester corrections from https://twiki.cern.ch/twiki/bin/view/CMS/RochcorMuon
```
