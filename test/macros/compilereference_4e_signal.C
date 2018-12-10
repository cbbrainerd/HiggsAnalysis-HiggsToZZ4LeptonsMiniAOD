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

#endif

using namespace std;


int main (int argc, char ** argv){

  string site=argv[7];
  cout << site << " configuration" <<endl;
  
  string dataconf=argv[8];
  cout << dataconf << " data configuration" <<endl;

  string mcconf=argv[9];
  cout << mcconf << " mc configuration" <<endl;

  ifstream fsig;
  fsig.open(argv[1]);

  const int signlines = atoi(argv[2]);

  string sigsamples[signlines];
  float signinput[signlines];
  float signhlt[signlines];
  float signskim[signlines];
  float sigxsection[signlines];
  
  // Signal
  cout << "Signal samples" << endl;
  for(int i=0;i<signlines;i++){
    fsig >> sigsamples[i] >> signinput[i] >> signhlt[i] >> signskim[i] >> sigxsection[i];
    cout << "Sample=" << sigsamples[i] << " Ninput=" << signinput[i] << " NHLT=" << signhlt[i] << " NSkim=" << signskim[i] << " Xsection(pb)=" << sigxsection[i] << endl;
  }


  float lumifb=0.;

  if (mcconf.find("Spring16")<5) lumifb=35.867;
  //if (mcconf.find("Spring16")<5) lumifb=10.693;

  // running on each signal and all bkg
  for(int i=0;i<signlines;i++){

    string name= "roottree_leptons_crab_"+sigsamples[i]+".root";
    TString dirInput;
    if (site.find("CERN")<5){
      if (mcconf.find("Fall11")<5) dirInput="/castor/cern.ch/user/n/ndefilip/Paper/MCFall11";    // to run at CERN
      else dirInput="/castor/cern.ch/user/n/ndefilip/Paper/MCSummer12";    // to run at CERN
    }
    else if (site.find("DESY")<5){
      dirInput="/nfs/dust/test/cmsdas/school16/HZZ4lExercise/sig"; //to run at DESY
    }
    else if (site.find("FNAL")<5 && mcconf.find("Fall11")<5){
      dirInput="dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/cmsdas/2014/HZZ4lExercise/sig/Fall11";
    }
    else if (site.find("FNAL")<5 && mcconf.find("Summer12")<5 ){
      dirInput="dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/cmsdas/2014/HZZ4lExercise/sig/Summer12";
    }
    else if (site.find("FNAL")<5 && mcconf.find("Spring16")<5 ){
      dirInput="root://cmseos.fnal.gov///store/user/cmsdas/2017/long_exercises/MonoHiggsHZZ/Spring16_merged";
    }
    else if (mcconf.find("Fall11")<5){
      dirInput="/lustre/cms/store/user/defilip/Fall11_445_paper_step_analysis_merged"; 
    }
    else if (mcconf.find("Summer12")<5){
      dirInput="/lustre/cms/store/user/defilip/Summer12_53X_paper_step_analysis_merged";  
    }
    else if (mcconf.find("Phys14")<5){
       dirInput="/lustre/cms/store/user/defilip/MonoHiggs/Phys14_720_merged";  
    }
    else if (mcconf.find("Fall15")<5){
       dirInput="/lustre/cms/store/user/gminiell/MonoHiggs/Signal_MonoHiggs_76X_merged";  
    }
    else if (mcconf.find("Spring16")<5){
      dirInput="/lustre/cms/store/user/dburns/MonoHiggs/Spring16_merged";
    }


    TString sigFile=name;
    sigFile=sigFile.ReplaceAll("-4e","");
 
    //if (sigFile.Contains("2HDM") ) dirInput="/lustre/cms/store/user/gminiell/MonoHiggs/Spring16_merged";
    if (sigFile.Contains("2HDM") || sigFile.Contains("Baryon") ) dirInput="/lustre/cms/store/user/defilip/MonoHiggs/Moriond17_merged_Syst";    
    if (sigFile.Contains("HHTo2B4L_madgraph_pythia8") ) dirInput="/lustre/cms/store/user/defilip/MonoHiggs/Moriond17_merged_Syst";

    Char_t nome[300];
    sprintf(nome,"%s/%s",dirInput.Data(),sigFile.Data());

    TString mhstring=sigFile.ReplaceAll("GluGluToHToZZTo4L_M-","");
    if (mcconf.find("Fall11")<5) { 
      mhstring=mhstring.ReplaceAll("_7TeV-powheg-pythia6","");
    }
    else if (mcconf.find("Summer12")<5) {
      mhstring=mhstring.ReplaceAll("_8TeV-powheg-pythia6","");
    }
    else if (mcconf.find("Phys14")<5) {
      mhstring=mhstring.ReplaceAll("Phys14_Higgs_Zprime_nohdecay_","");
      mhstring=mhstring.ReplaceAll("GeV_13TeV","");
    }
    
    mhstring=mhstring.ReplaceAll("roottree_leptons_crab","");
    mhstring=mhstring.ReplaceAll(".root","");
  
    Double_t mH=0.;
    mH=atof(mhstring.Data());
    
    
    float weight= lumifb*(sigxsection[i]*1000.*signskim[i]/signinput[i])/signskim[i];


    cout << "mH= " << mH << " weight= " << weight << endl;

    bool runonsig=true;

    if (runonsig==true){
      TFile *file3;
      file3 = TFile::Open(nome);
      cout << "Read file with name: " << nome << endl;
      TTree *tree3 = (TTree*)file3->Get("HZZ4LeptonsAnalysis");

      HZZ4LeptonsAnalysis make3(tree3,weight,dataconf,mcconf);

      char *path=NULL;
      size_t size=300;
      path=getcwd(path,size);

      sprintf(nome,"%s/output_%s.root",path,sigsamples[i].c_str());
      make3.Loop(nome);
      
      cout << "Create file with name: " << nome << endl;
      delete tree3;
      file3 -> Close();
    }

    
    
  }
  
  return 0; 

}

