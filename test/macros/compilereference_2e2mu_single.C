#if !defined(__CINT__) || defined(__MAKECINT__)

#include "HZZ4LeptonsAnalysis_2e2mu.h"
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
#include <string>
#include <iostream>
#include <TSystem.h>
#include <TH2.h>
#include "TChain.h"
#include <stdlib.h>
#include <TDCacheFile.h>

#endif

using namespace std;

int main(int argc, char ** argv){
  
  string site=argv[1];
  cout << site << " configuration" <<endl;

  string dataconf=argv[2];
  cout << dataconf << " dataconfiguration" <<endl;

  string mcconf=argv[3];
  cout << mcconf << " mc configuration" <<endl;

  
  Char_t nome[300];

  if (site.find("CERN")<5){
    sprintf(nome,"/castor/cern.ch/user/n/ndefilip/DAS2013/bkg/roottree_leptons_ZZTo2e2mu_8TeV-powheg-pythia6.root"); 
  }
  else if (site.find("FNAL")<5){
    sprintf(nome,"dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/cmsdas/2013/HZZ4lExercise/sig/roottree_leptons_GluGluToHToZZTo4L_M-1000_8TeV-powheg-pythia6.root");
  }
  else {
    //sprintf(nome,"/localdata/Syncr13TeV/roottree_leptons_sync_Fall15_HiggsToZZ_76x.root");
    sprintf(nome,"roottree_leptons_2015.root");
  }


  TFile *file3;
  file3 = TFile::Open(nome);
  
  //if (site.find("FNAL")<5){
  //  file3 = new  TDCacheFile (nome,"READ","ROOT file with trees",0);
  //}
  //else {
  //  file3 = TFile::Open(nome);
  //}

  cout << "Read file with name: " << nome << endl;
  TTree *tree3 = (TTree*)file3->Get("HZZ4LeptonsAnalysis");
  HZZ4LeptonsAnalysis make3(tree3,1.,dataconf,mcconf);
  //HZZ4LeptonsAnalysis make3(tree3);

  sprintf(nome,"output_2015_Higgs.root");
  make3.Loop(nome);

  cout << "Create file with name: " << nome << endl;
  delete tree3;
  file3 -> Close();

  return 0; 

}

