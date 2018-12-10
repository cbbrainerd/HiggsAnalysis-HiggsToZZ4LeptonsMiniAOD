#if !defined(__CINT__) || defined(__MAKECINT__)

#include "HZZ4LeptonsAnalysis_4e.h"
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
//#include <string>
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


  // bkg
  ifstream fbkg;
  fbkg.open(argv[3]);

  const int nlines = atoi(argv[4]);

  string bkgsamples[nlines];
  float bkgninput[nlines];
  float bkgnhlt[nlines];
  float bkgnskim[nlines];
  float bkgxsection[nlines];
  
  cout << "Bkg samples" << endl;

  for(int i=0;i<nlines;i++){
    fbkg >> bkgsamples[i] >> bkgninput[i] >> bkgnhlt[i] >> bkgnskim[i] >> bkgxsection[i];
    cout << "Sample=" << bkgsamples[i] << " Ninput=" << bkgninput[i] << " NHLT=" << bkgnhlt[i] << " NSkim=" << bkgnskim[i] << " Xsection(pb)=" << bkgxsection[i] << endl;
  }

  
  // 
  float lumifb=0.;

 if (mcconf.find("Fall17")<5) lumifb=41.37;  
 // if (mcconf.find("Spring16")<5) lumifb=35.867;
  //if (mcconf.find("Spring16")<5) lumifb=26.958;

  // Run on bkg
  for(int i=0;i<nlines;i++){
    
    string name= "roottree_leptons_crab_"+bkgsamples[i]+".root";
    
    TString dirInput;
    if (site.find("CERN")<5){
      if (mcconf.find("Fall11")<5) dirInput="/castor/cern.ch/user/n/ndefilip/Paper/MCFall11";    // to run at CERN
      else dirInput="/castor/cern.ch/user/n/ndefilip/Paper/MCSummer12";    // to run at CERN
    }
    else if (site.find("DESY")<5){
      dirInput="/nfs/dust/test/cmsdas/school16/HZZ4lExercise/bkg"; //to run at DESY
    }
    else if (site.find("FNAL")<5){
      dirInput="root://cmseos.fnal.gov///store/user/cmsdas/2017/long_exercises/MonoHiggsHZZ/Spring16_merged";
    }
    else if (mcconf.find("Spring16")<5){
      dirInput="/lustre/cms/store/user/mmelodea/RunIISummer16";
    }
    else if (mcconf.find("Fall17")<5){
      dirInput="/lustre/cms/store/user/reham/RunIIFall17";
    }
    
    TString bkgFile=name;
    cout << bkgFile.Data() << endl;
    
    if (bkgFile.Contains("HHTo2B4L")  ){
      TString datasetName=bkgFile.ReplaceAll("_","_");
      bkgFile=datasetName;
      dirInput="/lustre/cms/store/user/defilip/MonoHiggs/Moriond17_merged_Syst";
    }
    
    //if (bkgFile.Contains("QCD")  ){
    //  TString datasetName=bkgFile.ReplaceAll("_","_");
    //  bkgFile=datasetName;
    //  dirInput="/lustre/cms/store/user/gminiell/MonoHiggs/Spring16_merged";
    //}
    
    
    Char_t nome[300];
    sprintf(nome,"%s/%s",dirInput.Data(),bkgFile.Data());
    
    float weight= lumifb*(bkgxsection[i]*1000.*bkgnskim[i]/bkgninput[i])/bkgnskim[i];
    cout << "weight is " << weight << endl;	
    
    TFile *file3;
    TTree *tree3;
    
    if (bkgFile.Contains("DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")){
      TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
      chain->Add("/lustre/cms/store/user/mmelodea/RunIISummer16/roottree_leptons_crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0.root");
      chain->Add("/lustre/cms/store/user/mmelodea/RunIISummer16/roottree_leptons_crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1.root");
      chain->Add("/lustre/cms/store/user/mmelodea/RunIISummer16/roottree_leptons_crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2.root");
      chain->Add("/lustre/cms/store/user/mmelodea/RunIISummer16/roottree_leptons_crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_3.root");
      chain->Add("/lustre/cms/store/user/mmelodea/RunIISummer16/roottree_leptons_crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_4.root");
      chain->Add("/lustre/cms/store/user/mmelodea/RunIISummer16/roottree_leptons_crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_5.root");
      chain->Add("/lustre/cms/store/user/mmelodea/RunIISummer16/roottree_leptons_crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_6.root");
      chain->Add("/lustre/cms/store/user/mmelodea/RunIISummer16/roottree_leptons_crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_7.root");
      chain->Add("/lustre/cms/store/user/mmelodea/RunIISummer16/roottree_leptons_crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_8.root");
      chain->Add("/lustre/cms/store/user/mmelodea/RunIISummer16/roottree_leptons_crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_9.root");
      tree3 = chain;
    }
    else if (bkgFile.Contains("WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")){
      TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
      chain->Add("/lustre/cms/store/user/mmelodea/RunIISummer16/roottree_leptons_crab_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
      chain->Add("/lustre/cms/store/user/mmelodea/RunIISummer16/roottree_leptons_crab_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1.root");
      tree3 = chain;
    }  
    /*
	if (bkgFile.Contains("QCD_Pt_10to15_TuneCUETP8M1_13TeV_pythia8")){
	  TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_10to15_TuneCUETP8M1_13TeV_pythia8.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_10to15_TuneCUETP8M1_13TeV_pythia8_1.root");
	  tree3 = chain;
	}
	else if (bkgFile.Contains("QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8")){
	  TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8_1.root");
	  tree3 = chain;
	}
	else if (bkgFile.Contains("QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8")){
	  TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8_1.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8_2.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8_3.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8_4.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8_5.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8_6.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8_7.root");
	  tree3 = chain;
	}
	else if (bkgFile.Contains("QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8")){
	  TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8_1.root");	  
	  tree3 = chain;
	}	  
	else if (bkgFile.Contains("QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8")){
	  TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8_1.root");	  
	  tree3 = chain;
	}	  
	else if (bkgFile.Contains("QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8")){
	  TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8_1.root");	  
	  tree3 = chain;
	}
	else if (bkgFile.Contains("QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8")){
	  TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8_1.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8_2.root");
	  tree3 = chain;
	}
	else if (bkgFile.Contains("QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8")){
	  TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8_1.root");
	  tree3 = chain;
	}
	else if (bkgFile.Contains("QCD_Pt_5to10_TuneCUETP8M1_13TeV_pythia8")){
	  TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_5to10_TuneCUETP8M1_13TeV_pythia8.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_5to10_TuneCUETP8M1_13TeV_pythia8_1.root");
	  tree3 = chain;
	}
	else if (bkgFile.Contains("QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8")){
	  TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8.root");
	  chain->Add("/lustre/cms/store/user/dburns/MonoHiggs/Fall15_25ns_merged/roottree_leptons_crab_QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8_1.root");
	  tree3 = chain;
	}
        */
    else {	          
      file3 = TFile::Open(nome);
      cout << "Read file with name: " << nome << endl;
      tree3 = (TTree*)file3->Get("HZZ4LeptonsAnalysis");	  
    }
    
    
    HZZ4LeptonsAnalysis make3(tree3,weight,dataconf,mcconf);
    
    char *path=NULL;
    size_t size=300;
    path=getcwd(path,size);
    
    // cout << "Path is " << path << endl;
    sprintf(nome,"%s/output_%s.root",path,bkgsamples[i].c_str());
    
    // cout << "Output file name is " << nome << endl;
    make3.Loop(nome);
    
    cout << "Create file with name: " << nome << endl;
    delete tree3;
    file3 -> Close();
    
  }
  
  return 0; 

}

