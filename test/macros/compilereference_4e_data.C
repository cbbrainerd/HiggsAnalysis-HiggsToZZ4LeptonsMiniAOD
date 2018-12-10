#if !defined(__CINT__) || defined(__MAKECINT__)

#include "HZZ4LeptonsAnalysis_4e.h"
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
#include <iostream>
#include <TSystem.h>
#include <TH2.h>
#include "TChain.h"
#include <stdlib.h>
#include <TDCacheFile.h>

#endif

using namespace std;


int main (int argc, char ** argv){
  
  cout << "Ciao" << endl;
  string site=argv[7];
  cout << site << " configuration" <<endl;
  
  string dataconf=argv[8];
  

  if (dataconf.find("2012_postICHEP")<15) dataconf="2012";
  cout << dataconf << " data configuration" <<endl;
  
  string mcconf=argv[9];
  cout << mcconf << " mc configuration" <<endl;
  
  // data
  ifstream fdata;
  fdata.open(argv[5]);
  
  cout << "Data samples" << endl;
  
  const int ndata = atoi(argv[6]);
  string datasamples[ndata];
  
  for(int i=0;i<ndata;i++){
    fdata >> datasamples[i];
    cout << "Sample=" << datasamples[i] << endl;
  }
  
  
  // 
  
  float lumifb=4710./1000.;
  
  Double_t mH=150.;
  cout << "mH= " << mH << endl;
  
  
  // data
  bool runondata=true;
  
  if (runondata==true){
    
    bool itera=false;
    
    //if (mH<=160.){
    //  if (i==0) itera=true;  
    //}
    //else {
    itera=true;
    //}
    
    for(int i=0;i<ndata && itera==true;i++){
      
      string name= "roottree_leptons_crab_";
      name.append(datasamples[i]);
      name.append(".root");
      TString dirInput;
      
      if ( site.find("CERN")<5 ){
	dirInput="/castor/cern.ch/user/n/ndefilip/Paper/data";    // to run at CERN
      }
      else if (site.find("DESY")<5){
	dirInput="/nfs/dust/test/cmsdas/school16/HZZ4lExercise/data"; //to run at DESY
      }
      else if ( site.find("FNAL")<5 && name.find("2011")<38){
        dirInput="dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/cmsdas/2014/HZZ4lExercise/data/2011";
      }
      else if ( site.find("FNAL")<5 && name.find("2012")<38){
        dirInput="dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/cmsdas/2014/HZZ4lExercise/data/2012";
      }
      else if (site.find("FNAL")<5 && dataconf.find("2016")<5 ){
	dirInput="root://cmseos.fnal.gov///store/user/cmsdas/2017/long_exercises/MonoHiggsHZZ/Data2016_MonoHiggs_13TeV_merged";
      }
      else {
        if ( name.find("2015")<100){
          dirInput="/lustre/cms/store/user/gminiell/MonoHiggs/Data2015_MonoHiggs_13TeV_76X_merged";
        }
        else if ( name.find("2016")<38){
          dirInput="/lustre/cms/store/user/defilip/MonoHiggs/Data2016_MonoHiggs_13TeV_merged_Syst/";
        }
      }

      cout << "Input directory is:" << dirInput<< endl;

      TString dataFile=name;
      if (dataFile.Contains("miniaod")) dirInput="/lustre/cms/store/user/defilip/MonoHiggs/Data2016_MonoHiggs_13TeV_merged_Syst/";  

      Char_t nome[300];
      sprintf(nome,"%s/%s",dirInput.Data(),name.c_str());
     

      TFile *file3;
      TTree *tree3;
      
      //if (dataFile.Contains("SingleElectron_Run2015D-16Dec2015-v1")){
      //  TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
      //  chain->Add("/lustre/cms/store/user/gminiell/MonoHiggs/Data2015_MonoHiggs_13TeV_76X_merged/roottree_leptons_crab_SingleElectron_Run2015D-16Dec2015-v1_part_1.root");
      //  chain->Add("/lustre/cms/store/user/gminiell/MonoHiggs/Data2015_MonoHiggs_13TeV_76X_merged/roottree_leptons_crab_SingleElectron_Run2015D-16Dec2015-v1_part_1_1.root");
      //  chain->Add("/lustre/cms/store/user/gminiell/MonoHiggs/Data2015_MonoHiggs_13TeV_76X_merged/roottree_leptons_crab_SingleElectron_Run2015D-16Dec2015-v1_part_2.root");
      //  chain->Add("/lustre/cms/store/user/gminiell/MonoHiggs/Data2015_MonoHiggs_13TeV_76X_merged/roottree_leptons_crab_SingleElectron_Run2015D-16Dec2015-v1_part_2_1.root");
      //  chain->Add("/lustre/cms/store/user/gminiell/MonoHiggs/Data2015_MonoHiggs_13TeV_76X_merged/roottree_leptons_crab_SingleElectron_Run2015D-16Dec2015-v1_part_3.root");
      //  chain->Add("/lustre/cms/store/user/gminiell/MonoHiggs/Data2015_MonoHiggs_13TeV_76X_merged/roottree_leptons_crab_SingleElectron_Run2015D-16Dec2015-v1_part_3_1.root");
      //  chain->Add("/lustre/cms/store/user/gminiell/MonoHiggs/Data2015_MonoHiggs_13TeV_76X_merged/roottree_leptons_crab_SingleElectron_Run2015D-16Dec2015-v1_part_3_2.root");
      //  chain->Add("/lustre/cms/store/user/gminiell/MonoHiggs/Data2015_MonoHiggs_13TeV_76X_merged/roottree_leptons_crab_SingleElectron_Run2015D-16Dec2015-v1_part_4.root");
      //  chain->Add("/lustre/cms/store/user/gminiell/MonoHiggs/Data2015_MonoHiggs_13TeV_76X_merged/roottree_leptons_crab_SingleElectron_Run2015D-16Dec2015-v1_part_4_1.root");
      //  chain->Add("/lustre/cms/store/user/gminiell/MonoHiggs/Data2015_MonoHiggs_13TeV_76X_merged/roottree_leptons_crab_SingleElectron_Run2015D-16Dec2015-v1_part_4_2.root");
      //  tree3 = chain;
      //}
      //{
	file3 = TFile::Open(nome);
	cout << "Read file with name: " << nome << endl;
	tree3 = (TTree*)file3->Get("HZZ4LeptonsAnalysis");
	//}

      HZZ4LeptonsAnalysis make3(tree3,1.,dataconf,mcconf);      

      char *path=NULL;
      size_t size=300;
      path=getcwd(path,size);

      sprintf(nome,"%s/output_%s.root",path,datasamples[i].c_str());
      make3.Loop(nome);
      
      cout << "Create file with name: " << nome << endl;
      delete tree3;
      file3 -> Close();
      
    }
    
  }
  
  return 0; 

}

