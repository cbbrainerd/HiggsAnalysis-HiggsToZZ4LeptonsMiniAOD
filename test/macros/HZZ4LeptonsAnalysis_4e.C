#define HZZ4LeptonsAnalysis_cxx
#include "HZZ4LeptonsAnalysis_4e.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TSpline.h>

#include <TMath.h>
#include <fstream>
#include <iostream>
#include <sstream>

#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <libgen.h>

//#include "ZZMatrixElement/MELA/src/computeAngles.h"
//#include "ZZMatrixElement/MELA/src/computeAngles.cc"
#include "ZZMatrixElement/MELA/interface/TUtil.hh"
#include "ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h"


using namespace std;
// using namespace RooFit;
// using namespace meMCFM;
using namespace MEMNames;

const double Zmass = 91.188; // nominal Z boson mass
const double mPI = 3.141592654; 

void HZZ4LeptonsAnalysis::Loop(Char_t *output)
{

   if (fChain == 0) return;

   // Declare MEM class
   MEMs combinedMEM(13,125,"CTEQ6L");     
         
   // BNN
   Char_t datasetChar[500],bnnOUT[500],eventsOUT[500];
  
   cout << "The output file is " << output << endl;
   TString out = output;
   TString datasetName=out.ReplaceAll(".root","");
   sprintf(datasetChar,"%s",datasetName.Data());
   sprintf(bnnOUT,"%s_bnn.txt",datasetName.Data());
   sprintf(eventsOUT,"%s_bnn.root",datasetName.Data());
   cout << "bnnOUT= " << bnnOUT << endl;
   bnn_file.open(bnnOUT);

   bnn_file << "weight int_weight lept1_pt lept1_eta lept1_phi lept1_charge lept1_pfx lept1_sip lept1_mvaid lept2_pt lept2_eta lept2_phi lept2_charge lept2_pfx lept2_sip lept2_mvaid lept3_pt lept3_eta lept3_phi lept3_charge lept3_pfx lept3_sip lept3_mvaid lept4_pt lept4_eta lept4_phi lept4_charge lept4_pfx lept4_sip lept4_mvaid iso_max sip_max Z1mass Z2mass angle_costhetastar angle_costheta1 angle_costheta2 angle_phi angle_phistar1 KD psKD gravKD pt4l mass4l pfmet jet1_pt jet1_eta jet1_phi jet1_et jet2_pt jet2_eta jet2_phi jet2_et deltaetajj massjj VD njets pfmet" << endl;

   // Book txt file for candidate events
   Char_t txtOUT[500];
   sprintf(txtOUT,"%s_txt.txt",datasetName.Data());
   cout << "Opening a txt file with candidate events " << txtOUT << endl;
   ofstream output_txt;
   output_txt.open(txtOUT);

   Char_t txtOUT_VBF[500];
   sprintf(txtOUT_VBF,"%s_txt_vbf.txt",datasetName.Data());
   cout << "Opening a txt file with candidate events " << txtOUT_VBF << endl;
   ofstream output_txt_vbf;
   output_txt_vbf.open(txtOUT_VBF);


   // isSignal DM
   //TString signal="2HDM";
   Bool_t isSignal=false;
   //isSignal=(datasetName.Contains(signal));
   isSignal=(datasetName.Contains("2HDM") || datasetName.Contains("Baryonic"));
   cout << "isSignal= " << isSignal << endl;

   char *basename(char *path);
   cout << "Basename " << basename(datasetChar) << endl;
   TString datasetBase=basename(datasetChar);

   // pT reweigthing
   if ( datasetBase.Contains("2HDM") && datasetBase.Contains("Target") ){
     datasetBase=datasetBase.ReplaceAll("_reweighted","");
     datasetBase=datasetBase.ReplaceAll("output_ZprimeToA0hToA0chichihZZTo4l_2HDM_13TeV-madgraph_","");
     datasetBase=datasetBase.ReplaceAll("MZp-","");
     datasetBase=datasetBase.ReplaceAll("MA0-","");
     cout << "DatasetBase is= " << datasetBase.Data() << endl;
   }
   else if ( datasetBase.Contains("Baryonic") && datasetBase.Contains("Target") ){
     datasetBase=datasetBase.ReplaceAll("_reweighted","");
     datasetBase=datasetBase.ReplaceAll("output_MonoHZZ4l_ZpBaryonic_13TeV-madgraph_","");
     datasetBase=datasetBase.ReplaceAll("MZp-","");
     datasetBase=datasetBase.ReplaceAll("MChi-","");
     cout << "DatasetBase is= " << datasetBase.Data() << endl;
   }

   TString ggH="GluGluHToZZTo4L";
   Bool_t isggH=false;
   TString ggHnew="SMHiggsToZZTo4L";
   isggH=(datasetName.Contains(ggH) || datasetName.Contains(ggHnew));
   cout << "isggH " << isggH << endl;
   

   TString powheg15="powheg15";
   Bool_t ispowheg15=false;
   ispowheg15=datasetName.Contains(powheg15);
   cout << "ispowheg15= " << int(ispowheg15) << endl;

   TString vbf="VBF";
   Bool_t isvbf=false;
   isvbf=datasetName.Contains(vbf);
   cout << "isvbf= " << int(isvbf) << endl;

   TString mhstring=datasetBase.ReplaceAll("GluGluToHToZZTo4L_M-","");
   mhstring=mhstring.ReplaceAll("VBF_HToZZTo4L_M-","");
   mhstring=mhstring.ReplaceAll("SMHiggsToZZTo4L_M-","");
   mhstring=mhstring.ReplaceAll("_7TeV-powheg-pythia6","");
   mhstring=mhstring.ReplaceAll("_8TeV-powheg-pythia6","");
   mhstring=mhstring.ReplaceAll("_7TeV-powheg15-pythia6","");
   mhstring=mhstring.ReplaceAll("_8TeV-powheg15-pythia6","");
   mhstring=mhstring.ReplaceAll("_8TeV-powheg15-JHUgenV3-pythia6","");
   mhstring=mhstring.ReplaceAll("_7TeV-powheg15-JHUgenV3-pythia6","");
   mhstring=mhstring.ReplaceAll("output_","");
   cout << "mhstring=" << mhstring << endl;
  
   float mHgen=0.;
   mHgen=atof(mhstring.Data());
   cout << "mHgen= " << mHgen << endl;
   
   //pT gen reweighting
   char histotitle[500];
   sprintf(histotitle,"%s",datasetBase.Data());
   TH1F *pTweighthisto;

   TFile *f2HDM=TFile::Open("HISTOShapes2HDM_READ_ext.root");
   TFile *fBaryonic=TFile::Open("HISTOShapesZpB_READ.root");

   if ( datasetBase.Contains("2HDM") && datasetBase.Contains("Target")){
     pTweighthisto=(TH1F*)f2HDM->Get(histotitle);
   }
   else if ( datasetBase.Contains("Baryonic") && datasetBase.Contains("Target")){
     pTweighthisto=(TH1F*)fBaryonic->Get(histotitle);
   }

   // Pileup reweighting 2016 data vs Spring16 MC in 80x                                                                                                                                                  

   TFile *_filePU;
   //_filePU= TFile::Open("pileup_MC_80x_271036-276811_69200.root");
   //TH1D *puweight = (TH1D*)_filePU->Get("puweight");
   _filePU= TFile::Open("PU_Reweight_2017.root");     
   TH1F *puweight = (TH1F*)_filePU->Get("PU_Ratio"); 

   /////////////Lepton Efficiency Scale Factrons/////////////
   // Load histograms
   //
   // TFile *ele_scale_mvafactors2016 = new TFile("MVA_vs_RECO.root");
   // TH2F *ele_scale_mva2016 = (TH2F*)gDirectory->Get("hScaleFactors_ID"); 
   // TFile *ele_scale_mvafactors2016_Crack = new TFile("MVA_vs_RECO_Cracks.root");
   // TH2F *ele_scale_mva2016_Cracks = (TH2F*)gDirectory->Get("hScaleFactors_ID_Cracks"); 
   
   //2016  
   //TFile *ele_scale_factors_v3 = new TFile("ele_scale_factors_v3.root");
   //TH2F *ele_scale_factors2016     = (TH2F*)gDirectory->Get("ele_scale_factors"); 
   //TH2F *ele_scale_factors_gap2016 = (TH2F*)gDirectory->Get("ele_scale_factors_gap"); 

   //2017
   //electron Id,Iso efficiency
    TFile *ele_scale_factors2017 = new TFile("egammaEffi_txt_EGM2D_Moriond2018v1.root");
    TH2F *ele_scale_2017 = (TH2F*)gDirectory->Get("EGamma_SF2D"); 
    TFile *ele_scale_factors2017_gap = new TFile("egammaEffi_txt_EGM2D_Moriond2018v1_gap.root");
    TH2F *ele_scale_2017_gap = (TH2F*)gDirectory->Get("EGamma_SF2D");

    //electron Reconstruction efficiency
    
    TFile *ele_RecoEff_2017_lowEt = new TFile("egammaEffi_txt_EGM2D_runBCDEF_passingRECO_lowEt.root"); //Low Et<20
    TH2F *ele_Reco_eff_2017_lowEt = (TH2F*)gDirectory->Get("EGamma_SF2D");

    TFile *ele_RecoEff_2017_highEt = new TFile("egammaEffi_txt_EGM2D_runBCDEF_passingRECO.root"); //high Et>20
    TH2F *ele_Reco_eff_2017_highEt = (TH2F*)gDirectory->Get("EGamma_SF2D");
    
   // TFile *ele_scale_Sipfactors2016 = new TFile("SIP.root");
   // TH2F *ele_scale_Sip2016 = (TH2F*)gDirectory->Get("hScaleFactors_SIP"); 
   // TFile *ele_scale_Sipfactors2016_Crack = new TFile("SIP_Cracks.root");
   // TH2F *ele_scale_Sip2016_Cracks = (TH2F*)gDirectory->Get("hScaleFactors_SIP_Cracks"); 
   
   // correction to the error
   /*
   TFile *ebe_corr2012= new TFile("ebeOverallCorrections.Legacy2013.v0.root");
   TH2F *ebe_mu_reco53x= (TH2F*)gDirectory->Get("mu_reco53x");
   TH2F *ebe_mu_mc53x= (TH2F*)gDirectory->Get("mu_mc53x");

   TFile *ebe_corr2011= new TFile("ebeOverallCorrections.LegacyPaper.42x.root");
   TH2F *ebe_mu_reco42x= (TH2F*)gDirectory->Get("mu_reco42x");
   TH2F *ebe_mu_mc42x= (TH2F*)gDirectory->Get("mu_mc42x");
   */

   // kfactor_ggZZ(float GENmassZZ, int finalState)     
   TString strSystTitle[9] ={
     "Nominal",
     "PDFScaleDn",
     "PDFScaleUp",
     "QCDScaleDn",
     "QCDScaleUp",
     "AsDn",
     "AsUp",
     "PDFReplicaDn",
     "PDFReplicaUp"
   };
   
   TFile* fin = TFile::Open("Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
   // Open the files
   TSpline3* ggZZ_kf[9]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
   for(int f=0;f<9;f++){
     ggZZ_kf[f] = (TSpline3*)fin->Get(Form("sp_kfactor_%s", strSystTitle[f].Data()));
   }   
   fin->Close();
   

   // Book root file (for output):
   TFile * theFile = new TFile(output,"RECREATE");

   //TString Vars("Weight:Run:Event:LumiSection:massZ1:massZ2:mass4l:Iso_max:Sip_max:MELA:FSR");
   //TNtuple * thePlots=new TNtuple("Candidates","Candidates",Vars);

    // Clone tree for Z1
   //TTree *z1tree = fChain->CloneTree(0);
   

   double DELTAPHI( double , double ) ; //call the function  
   double invmass (float M1, float PT1, float ETA1, float PHI1, float M2, float PT2, float ETA2, float PHI2 );
   
   // Book relevant variables -- counters:

   int N_0 = 0;	  // MC truth & acceptance
   int N_01 = 0;
   int N_02 = 0;

   int N_1 = 0;	  // Skim
   int N_2 = 0;

   int N_3a = 0;
   int N_3_FSR = 0;
   int N_3b = 0;

   int N_4a = 0;
   int N_4b = 0;
   int N_4c = 0;
   int N_4d = 0;

   int N_5 = 0;
   int N_6 = 0;
   int N_7 = 0;
   int N_7_PFMET = 0;
   int N_8 = 0;
   int N_8_a = 0;
   int N_8_b = 0;
   int N_9 = 0;

   int N_9_1FSR = 0;
   int N_9_2FSR = 0;

   int N_9PS = 0;
   int N_9GRAV = 0;
   
   int N_9a_VBF = 0;
   int N_9b_VBF = 0;
   int N_9_PFMET = 0;
   int N_VBF = 0;

   int N_10 = 0;
   int N_11 = 0;

   // counter weighted
   double N_0_w = 0;	  // MC truth & acceptance
   double N_01_w = 0;
   double N_02_w = 0;

   double N_1_w = 0;	  // Skim
   double N_2_w = 0;

   double N_3a_w = 0;
   double N_3_FSR_w = 0;
   double N_3b_w = 0;

   double N_4a_w = 0;
   double N_4b_w = 0;
   double N_4c_w = 0;
   double N_4d_w = 0;

   double N_5_w = 0;
   double N_6_w = 0;
   double N_7_w = 0;
   double N_7_PFMET_w = 0;

   double N_8_w = 0;
   double N_8_a_w = 0;
   double N_8_b_w = 0;
   double N_9_w = 0;

   double N_9_1FSR_w = 0;
   double N_9_2FSR_w = 0;

   double N_9PS_w = 0;
   double N_9GRAV_w = 0;
   
   double N_9a_VBF_w = 0;
   double N_9b_VBF_w = 0;
   double N_9_PFMET_w = 0;

   double N_VBF_w = 0;

   double N_10_w = 0;
   double N_11_w = 0;

   // Book Histos ***
   TH1D *nEvent_4l_w = new TH1D("nEvent_4l_w", "nEventComplete Weightd", 23, 0., 23.);
   TH1D *nEvent_4l = new TH1D("nEvent_4l", "nEventComplete", 23, 0., 23.);

   TH1D *nEvent_red = new TH1D("nEvent_red", "nEventCompleteReduced", 7, 0., 7.);
   TH1D *nEvent_ZZ = new TH1D("nEvent_ZZ", "nEventCompleteFromZZ", 5, 0., 5.);

   TH1F *hweight = new TH1F("hweight", "hweight",10000,0.,10.);

   TH1F *Gen_H_MASS              = new TH1F("Gen_H_MASS", "Gen_H_MASS",8000,0.,2000.);  
   TH1F *Gen_H_MASS_ReWeighted   = new TH1F("Gen_H_MASS_ReWeighted", "Gen_H_MASS_ReWeighted",8000,0.,2000.);  
   TH1F *Gen_H_MASS_ReWeightedP  = new TH1F("Gen_H_MASS_ReWeightedP", "Gen_H_MASS_ReWeightedP",8000,0.,2000.);  
   TH1F *Gen_H_MASS_ReWeightedM  = new TH1F("Gen_H_MASS_ReWeightedM", "Gen_H_MASS_ReWeightedM",8000,0.,2000.);  

   // Pileup reweighting
   TH1F *hPUvertices             = new TH1F("hPUvertices", "hPUvertices",70,0.,70.);  
   TH1F *hPUvertices_ReWeighted  = new TH1F("hPUvertices_ReWeighted", "hPUvertices_ReWeighted",70,0.,70.);  

   // pfmet with logaritmic bin width
   const int NMBINS = 10;
   const double MMIN = 1.0, MMAX = 1000.0;
   double logMbins[NMBINS+1];
   float binNormNr=0.;
   for (int ibin = 0; ibin <= NMBINS; ibin++) {
     logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);
     cout << logMbins[ibin] << endl;
   }

   const int NMOBINS = 5;
   const double MOMIN = 0.0, MOMAX = 1000.0;
   double loglinMbins[6]={0.,25.,50.,200., 500.,1000.};

   //step 3
   TH1F * hMZ_3 = new TH1F("hMZ_3", "Mass of Z after selection step 3", 200 , -0.5 , 199.5 );
   hMZ_3->SetXTitle("mass_Z  (GeV)");
   TH1F * hMZBB_3 = new TH1F("hMZBB_3", "Mass of Z after selection step 3 barrel barrel", 200 , -0.5 , 199.5 );
   hMZBB_3->SetXTitle("mass_Z barrel barrel  (GeV)");
   TH1F * hMZEE_3 = new TH1F("hMZEE_3", "Mass of Z after selection step 3 endcap endcap", 200 , -0.5 , 199.5 );
   hMZEE_3->SetXTitle("mass_Z endcap endcap  (GeV)");
   TH1F * hPtZ_3 = new TH1F("hPtZ_3", "Pt of Z after selection step 3", 200 , -0.5 , 199.5 );
   hPtZ_3->SetXTitle("pt_Z  (GeV)");
   TH1F * hYZ_3 = new TH1F("hYZ_3", "Y of Z after selection step 3", 500 , -5. , 5.);
   hYZ_3->SetXTitle("Y_Z");
   
   TH1F * hPtLep_3 = new TH1F("hPtLep_3", "Pt of Lep after selection step 3", 500 , -0.5 , 499.5 );
   hPtLep_3->SetXTitle("pt_Lep  (GeV)");
   TH1F * hYLep_3 = new TH1F("hEtaLep_3", "Y of Lep after selection step 3", 500 , -5. , 5. );
   hYLep_3->SetXTitle("Y of Lep");
   TH1F * hIsoLep_3 = new TH1F("hIsoLep_3", "Isolation of Lep after selection step 3", 2000 , -10., 10.);
   hIsoLep_3->SetXTitle("Iso");
   TH1F * hSipLep_3 = new TH1F("hSipLep_3", "Sip of Lep after selection step 3",  1000 , -20. , 40. );
   hSipLep_3->SetXTitle("Sip");
   TH1F * hIpLep_3 = new TH1F("hIpLep_3", "Ip of Lep after selection step 3",  1000 , -20. , 40. );
   hIpLep_3->SetXTitle("Ip");
   TH1F * hIpErLep_3 = new TH1F("hIpErLep_3", "IpEr of Lep after selection step 3",  1000 , 0. , 10. );
   hIpErLep_3->SetXTitle("IpEr");

   TH1F * hPtLep2_3 = new TH1F("hPtLep2_3", "Pt of Lep2 after selection step 3", 200 , -0.5 , 199.5 );
   hPtLep2_3->SetXTitle("pt_Lep2  (GeV)");
   TH1F * hYLep2_3 = new TH1F("hEtaLep2_3", "Y of Lep2 after selection step 3", 500 , -5. , 5. );
   hYLep2_3->SetXTitle("Y of Lep2");
   TH1F * hIsoLep2_3 = new TH1F("hIsoLep2_3", "Isolation of Lep2 after selection step 3", 2000 , -10., 10. );
   hIsoLep2_3->SetXTitle("Iso");
   TH1F * hSipLep2_3 = new TH1F("hSipLep2_3", "Sip of Lep2 after selection step 3",  1000 , -20. , 40. );
   hSipLep2_3->SetXTitle("Sip");
   TH1F * hIpLep2_3 = new TH1F("hIpLep2_3", "Ip of Lep2 after selection step 3",  1000 , -20. , 40. );
   hIpLep2_3->SetXTitle("Ip");
   TH1F * hIpErLep2_3 = new TH1F("hIpErLep2_3", "IpEr of Lep2 after selection step 3",  1000 , 0. , 10. );
   hIpErLep2_3->SetXTitle("IpEr");

   TH1F * hIso_3 = new TH1F("hIso_3", "Isolation maxima after selection step 3", 2000 , -10. , 10. );
   hIso_3->SetXTitle("Iso");
   TH1F * hSip_3 = new TH1F("hSip_3", "Sip maxima after selection step 3",  1000 , -20. , 40. );
   hSip_3->SetXTitle("Sip");
   TH1F * hIp_3 = new TH1F("hIp_3", "Ip maxima after selection step 3",  1000 , -20. , 40. );
   hIp_3->SetXTitle("Ip");

   TH1F * hDjj_3 = new TH1F("hDjj_3", "Delta jets vbf selection step 1", 200, -19.5, 19.5 );
   hDjj_3->SetXTitle("Delta jets");
   TH1F * hMjj_3 = new TH1F("hMjj_3", "Mass jets vbf selection step 1", 200, -0.5, 499.5 );
   hMjj_3->SetXTitle("Mass jets");
   TH1F * hVD_3 = new TH1F("hVD_3", "Discriminant vbf selection step 1", 200, -0.5, 9.5 );
   hMjj_3->SetXTitle("Discriminant");
   
   TH1F * hPFMET_3 = new TH1F("hPFMET_3", "PF MET after selection step 3", 4000 , 0., 4000.);
   hPFMET_3->SetXTitle("PF MET");
   
   // Build the histo with constant log bin width   
   TH1F *hLogXPFMET_3             = new TH1F("hLogXPFMET_3","hLogXPFMET_3",NMBINS, logMbins);   
   hLogXPFMET_3->Sumw2();
   TH1F *hLogLinXPFMET_3             = new TH1F("hLogLinXPFMET_3","hLogLinXPFMET_3",NMOBINS, loglinMbins);
   hLogLinXPFMET_3->Sumw2();

   //step 5
   TH1F * hM4l_5 = new TH1F("hM4l_5", "Mass of four leptons after selection step 5", 1200, 4.5,1204.5 );
   hM4l_5->SetXTitle("4 lepton mass  (GeV)");


   TH1F * hMZ1_5 = new TH1F("hMZ1_5", "Mass of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ1_5->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_5 = new TH1F("hPtZ1_5", "Pt of Z1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ1_5->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_5 = new TH1F("hYZ1_5", "Y of Z1 after selection step 5", 500 , -5. , 5.);
   hYZ1_5->SetXTitle("Y_Z1");

   TH1F * hMZ2_5 = new TH1F("hMZ2_5", "Mass of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hMZ2_5->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_5 = new TH1F("hPtZ2_5", "Pt of Z2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtZ2_5->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_5 = new TH1F("hYZ2_5", "Y of Z2 after selection step 5", 500 , -5. , 5. );
   hYZ2_5->SetXTitle("Y_Z2");

   TH1F * hPtLep1_5 = new TH1F("hPtLep1_5", "Pt of Lep1 after selection step 5", 200 , -0.5 , 199.5 );
   hPtLep1_5->SetXTitle("pt_Lep1  (GeV)");
   TH1F * hYLep1_5 = new TH1F("hEtaLep1_5", "Y of Lep1 after selection step 5", 500 , -5. , 5. );
   hYLep1_5->SetXTitle("Y of Lep2");
   TH1F * hIsoLep1_5 = new TH1F("hIsoLep1_5", "Isolation of Lep1 after selection step 5", 2000 , -10. , 10. );
   hIsoLep1_5->SetXTitle("Iso");
   TH1F * hSipLep1_5 = new TH1F("hSipLep1_5", "Sip of Lep1 after selection step 5",  1000 , -20. , 40. );
   hSipLep1_5->SetXTitle("Sip");
   TH1F * hIpLep1_5 = new TH1F("hIpLep1_5", "Ip of Lep1 after selection step 5",  1000 , -20. , 40. );
   hIpLep1_5->SetXTitle("Ip");
   TH1F * hIpErLep1_5 = new TH1F("hIpErLep1_5", "IpEr of Lep1 after selection step 5",  1000 , 0. , 10. );
   hIpErLep1_5->SetXTitle("IpEr");

   TH1F * hPtLep2_5 = new TH1F("hPtLep2_5", "Pt of Lep2 after selection step 5", 200 , -0.5 , 199.5 );
   hPtLep2_5->SetXTitle("pt_Lep2  (GeV)");
   TH1F * hYLep2_5 = new TH1F("hEtaLep2_5", "Y of Lep2 after selection step 5", 500 , -5. , 5. );
   hYLep2_5->SetXTitle("Y of Lep2");
   TH1F * hIsoLep2_5 = new TH1F("hIsoLep2_5", "Isolation of Lep2 after selection step 5", 2000 , -10. , 10. );
   hIsoLep2_5->SetXTitle("Iso");
   TH1F * hSipLep2_5 = new TH1F("hSipLep2_5", "Sip of Lep2 after selection step 5",  1000 , -20. , 40. );
   hSipLep2_5->SetXTitle("Sip");
   TH1F * hIpLep2_5 = new TH1F("hIpLep2_5", "Ip of Lep2 after selection step 5",  1000 , -20. , 40. );
   hIpLep2_5->SetXTitle("Ip");
   TH1F * hIpErLep2_5 = new TH1F("hIpErLep2_5", "IpEr of Lep2 after selection step 5",  1000 , 0. , 10. );
   hIpErLep2_5->SetXTitle("IpEr");
   
   TH1F * hPtLep3_5 = new TH1F("hPtLep3_5", "Pt of Lep3 after selection step 5", 200 , -0.5 , 199.5 );
   hPtLep3_5->SetXTitle("pt_Lep3  (GeV)");
   TH1F * hYLep3_5 = new TH1F("hEtaLep3_5", "Y of Lep3 after selection step 5", 500 , -5. , 5. );
   hYLep3_5->SetXTitle("Y of Lep2");
   TH1F * hIsoLep3_5 = new TH1F("hIsoLep3_5", "Isolation of Lep3 after selection step 5", 2000 , -10. , 10. );
   hIsoLep3_5->SetXTitle("Iso");
   TH1F * hSipLep3_5 = new TH1F("hSipLep3_5", "Sip of Lep3 after selection step 5",  1000 , -20. , 40. );
   hSipLep3_5->SetXTitle("Sip");
   TH1F * hIpLep3_5 = new TH1F("hIpLep3_5", "Ip of Lep3 after selection step 5",  1000 , -20. , 40. );
   hIpLep3_5->SetXTitle("Ip");
   TH1F * hIpErLep3_5 = new TH1F("hIpErLep3_5", "IpEr of Lep3 after selection step 5",  1000 , 0. , 10. );
   hIpErLep3_5->SetXTitle("IpEr");

   TH1F * hPtLep4_5 = new TH1F("hPtLep4_5", "Pt of Lep4 after selection step 5", 200 , -0.5 , 199.5 );
   hPtLep4_5->SetXTitle("pt_Lep4  (GeV)");
   TH1F * hYLep4_5 = new TH1F("hEtaLep4_5", "Y of Lep4 after selection step 5", 50 , -5. , 5. );
   hYLep4_5->SetXTitle("Y of Lep2");
   TH1F * hIsoLep4_5 = new TH1F("hIsoLep4_5", "Isolation of Lep4 after selection step 5", 2000 , -10. , 10. );
   hIsoLep4_5->SetXTitle("Iso");
   TH1F * hSipLep4_5 = new TH1F("hSipLep4_5", "Sip of Lep4 after selection step 5",  1000 , -20. , 40. );
   hSipLep4_5->SetXTitle("Sip");
   TH1F * hIpLep4_5 = new TH1F("hIpLep4_5", "Ip of Lep4 after selection step 5",  1000 , -20. , 40. );
   hIpLep4_5->SetXTitle("Ip");
   TH1F * hIpErLep4_5 = new TH1F("hIpErLep4_5", "IpEr of Lep4 after selection step 5",  1000 , 0. , 10. );
   hIpErLep4_5->SetXTitle("IpEr");


   TH1F * hIso_5 = new TH1F("hIso_5", "Isolation maxima after selection step 5", 2000 , -10. , 10. );
   hIso_5->SetXTitle("Iso");
   TH1F * hSip_5 = new TH1F("hSip_5", "Sip maxima after selection step 5",  1000 , -20. , 40. );
   hSip_5->SetXTitle("Sip");
   TH1F * hIp_5 = new TH1F("hIp_5", "Ip maxima after selection step 5",  1000 , -20. , 40. );
   hIp_5->SetXTitle("Ip");
   

   //step 6
   
   TH1F * hminMll_6 = new TH1F("hminMll_6", "minMll at selection step 6", 400 , 0. , 200.);
   hminMll_6->SetXTitle("minMll  (GeV)");
   

   //step 7
   TH1F * hM4l_7 = new TH1F("hM4l_7", "Mass of four leptons after selection step 7", 1200, 4.5,1204.5 );
   hM4l_7->SetXTitle("4 lepton mass  (GeV)");

   TH1F * hMZ1_7 = new TH1F("hMZ1_7", "Mass of Z1 after selection step 7", 200 , -0.5 , 199.5);
   hMZ1_7->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_7 = new TH1F("hPtZ1_7", "Pt of Z1 after selection step 7", 200 , -0.5 , 199.5);
   hPtZ1_7->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_7 = new TH1F("hYZ1_7", "Y of Z1 after selection step 7", 500 , -5. , 5.);
   hYZ1_7->SetXTitle("Y_Z1");

   TH1F * hMZ2_7 = new TH1F("hMZ2_7", "Mass of Z2 after selection step 7", 200 , -0.5 , 199.5);
   hMZ2_7->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_7 = new TH1F("hPtZ2_7", "Pt of Z2 after selection step 7", 200 , -0.5 , 199.5);
   hPtZ2_7->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_7 = new TH1F("hYZ2_7", "Y of Z2 after selection step 7", 500 , -5. , 5.);
   hYZ2_7->SetXTitle("Y_Z2");

   TH1F * hPtLep1_7 = new TH1F("hPtLep1_7", "Pt of Lep1 after selection step 7", 200 , -0.5 , 199.5 );
   hPtLep1_7->SetXTitle("pt_Lep1  (GeV)");
   TH1F * hYLep1_7 = new TH1F("hEtaLep1_7", "Y of Lep1 after selection step 7", 500 , -5. , 5. );
   hYLep1_7->SetXTitle("Y of Lep2");
   TH1F * hIsoLep1_7 = new TH1F("hIsoLep1_7", "Isolation of Lep1 after selection step 7", 2000 , -10. , 10. );
   hIsoLep1_7->SetXTitle("Iso");
   TH1F * hSipLep1_7 = new TH1F("hSipLep1_7", "Sip of Lep1 after selection step 7",  1000 , -20. , 40. );
   hSipLep1_7->SetXTitle("Sip");
   TH1F * hIpLep1_7 = new TH1F("hIpLep1_7", "Ip of Lep1 after selection step 7",  1000 , -20. , 40. );
   hIpLep1_7->SetXTitle("Ip");
   TH1F * hIpErLep1_7 = new TH1F("hIpErLep1_7", "IpEr of Lep1 after selection step 7",  1000 , 0. , 10. );
   hIpErLep1_7->SetXTitle("IpEr");


   TH1F * hPtLep2_7 = new TH1F("hPtLep2_7", "Pt of Lep2 after selection step 7", 200 , -0.5 , 199.5 );
   hPtLep2_7->SetXTitle("pt_Lep2  (GeV)");
   TH1F * hYLep2_7 = new TH1F("hEtaLep2_7", "Y of Lep2 after selection step 7", 500 , -5. , 5. );
   hYLep2_7->SetXTitle("Y of Lep2");
   TH1F * hIsoLep2_7 = new TH1F("hIsoLep2_7", "Isolation of Lep2 after selection step 7", 2000 , -10. , 10. );
   hIsoLep2_7->SetXTitle("Iso");
   TH1F * hSipLep2_7 = new TH1F("hSipLep2_7", "Sip of Lep2 after selection step 7",  1000 , -20. , 40. );
   hSipLep2_7->SetXTitle("Sip");
   TH1F * hIpLep2_7 = new TH1F("hIpLep2_7", "Ip of Lep2 after selection step 7",  1000 , -20. , 40. );
   hIpLep2_7->SetXTitle("Ip");
   TH1F * hIpErLep2_7 = new TH1F("hIpErLep2_7", "IpEr of Lep2 after selection step 7",  1000 , 0. , 10. );
   hIpErLep2_7->SetXTitle("IpEr");

   TH1F * hPtLep3_7 = new TH1F("hPtLep3_7", "Pt of Lep3 after selection step 7", 200 , -0.5 , 199.5 );
   hPtLep3_7->SetXTitle("pt_Lep3  (GeV)");
   TH1F * hYLep3_7 = new TH1F("hEtaLep3_7", "Y of Lep3 after selection step 7", 500 , -5. , 5. );
   hYLep3_7->SetXTitle("Y of Lep2");
   TH1F * hIsoLep3_7 = new TH1F("hIsoLep3_7", "Isolation of Lep3 after selection step 7", 2000 , -10. , 10. );
   hIsoLep3_7->SetXTitle("Iso");
   TH1F * hSipLep3_7 = new TH1F("hSipLep3_7", "Sip of Lep3 after selection step 7",  1000 , -20. , 40. );
   hSipLep3_7->SetXTitle("Sip");
   TH1F * hIpLep3_7 = new TH1F("hIpLep3_7", "Ip of Lep3 after selection step 7",  1000 , -20. , 40. );
   hIpLep3_7->SetXTitle("Ip");
   TH1F * hIpErLep3_7 = new TH1F("hIpErLep3_7", "IpEr of Lep3 after selection step 7",  1000 , 0. , 10. );
   hIpErLep3_7->SetXTitle("IpEr");

   TH1F * hPtLep4_7 = new TH1F("hPtLep4_7", "Pt of Lep4 after selection step 7", 200 , -0.5 , 199.5 );
   hPtLep4_7->SetXTitle("pt_Lep4  (GeV)");
   TH1F * hYLep4_7 = new TH1F("hEtaLep4_7", "Y of Lep4 after selection step 7", 500 , -5. , 5. );
   hYLep4_7->SetXTitle("Y of Lep2");
   TH1F * hIsoLep4_7 = new TH1F("hIsoLep4_7", "Isolation of Lep4 after selection step 7", 2000 , -10. , 10. );
   hIsoLep4_7->SetXTitle("Iso");
   TH1F * hSipLep4_7 = new TH1F("hSipLep4_7", "Sip of Lep4 after selection step 7",  1000 , -20. , 40. );
   hSipLep4_7->SetXTitle("Sip");
   TH1F * hIpLep4_7 = new TH1F("hIpLep4_7", "Ip of Lep4 after selection step 7",  1000 , -20. , 40. );
   hIpLep4_7->SetXTitle("Ip");
   TH1F * hIpErLep4_7 = new TH1F("hIpErLep4_7", "IpEr of Lep4 after selection step 7",  1000 , 0. , 10. );
   hIpErLep4_7->SetXTitle("IpEr");

   TH1F * hIso_7 = new TH1F("hIso_7", "Isolation maxima after selection step 7", 2000 , -10. , 10. );
   hIso_7->SetXTitle("Iso");
   TH1F * hSip_7 = new TH1F("hSip_7", "Sip maxima after selection step 7",  1000 , -20. , 40. );
   hSip_7->SetXTitle("Sip");
   TH1F * hIp_7 = new TH1F("hIp_7", "Ip maxima after selection step 7",  1000 , -20. , 40. );
   hIp_7->SetXTitle("Ip");
   TH1F * hPFMET_7 = new TH1F("hPFMET_7", "PF MET after selection step 7", 1000 , 0., 1000.);
   hPFMET_7->SetXTitle("PF MET (GeV)"); 

   TH1F * hDPHI_7 = new TH1F("DPHI_7", "polar angle between 4l and E_{T,miss}", 1000, 0., 5. );
   hDPHI_7->SetXTitle("#DELTA#phi(4l,E_{T,miss})");


   //step 8
   TH1F * hM4l_8 = new TH1F("hM4l_8", "Mass of four leptons after selection step 8", 1200, 4.5, 1204.5 );
   hM4l_8->SetXTitle("4 lepton mass  (GeV)");
   TH1F * hM4l_8_100_800 = new TH1F("hM4l_8_100_800", "Mass of four leptons after selection step 8", 1200, 4.5,1204.5 );
   hM4l_8_100_800->SetXTitle("4 lepton mass  (GeV)");
   
   TH1F * hMZ1_8 = new TH1F("hMZ1_8", "Mass of Z1 after selection step 8", 200 , -0.5 , 199.5 );
   hMZ1_8->SetXTitle("mass_Z1  (GeV)");
   TH1F * hPtZ1_8 = new TH1F("hPtZ1_8", "Pt of Z1 after selection step 8", 200 , -0.5 , 199.5 );
   hPtZ1_8->SetXTitle("pt_Z1  (GeV)");
   TH1F * hYZ1_8 = new TH1F("hYZ1_8", "Y of Z1 after selection step 8", 500 , -5. , 5.);
   hYZ1_8->SetXTitle("Y_Z1");

   TH1F * hMZ2_8 = new TH1F("hMZ2_8", "Mass of Z2 after selection step 8", 200 , -0.5 , 199.5 );
   hMZ2_8->SetXTitle("mass_Z2  (GeV)");
   TH1F * hPtZ2_8 = new TH1F("hPtZ2_8", "Pt of Z2 after selection step 8", 200 , -0.5 , 199.5 );
   hPtZ2_8->SetXTitle("pt_Z2  (GeV)");
   TH1F * hYZ2_8 = new TH1F("hYZ2_8", "Y of Z2 after selection step 8", 500 , -5. , 5. );
   hYZ2_8->SetXTitle("Y_Z2");

   TH1F * hPtLep1_8 = new TH1F("hPtLep1_8", "Pt of Lep1 after selection step 8", 200 , -0.5 , 199.5 );
   hPtLep1_8->SetXTitle("pt_Lep1  (GeV)");
   TH1F * hYLep1_8 = new TH1F("hEtaLep1_8", "Y of Lep1 after selection step 8", 500 , -5. , 5. );
   hYLep1_8->SetXTitle("Y of Lep2");
   TH1F * hIsoLep1_8 = new TH1F("hIsoLep1_8", "Isolation of Lep1 after selection step 8", 2000 , -10. , 10. );
   hIsoLep1_8->SetXTitle("Iso");
   TH1F * hSipLep1_8 = new TH1F("hSipLep1_8", "Sip of Lep1 after selection step 8",  1000 , -20. , 40. );
   hSipLep1_8->SetXTitle("Sip");
   TH1F * hIpLep1_8 = new TH1F("hIpLep1_8", "Ip of Lep1 after selection step 8",  1000 , -20. , 40. );
   hIpLep1_8->SetXTitle("Ip");
   TH1F * hIpErLep1_8 = new TH1F("hIpErLep1_8", "IpEr of Lep1 after selection step 8",  1000 , 0. , 10. );
   hIpErLep1_8->SetXTitle("IpEr");

   
   TH1F * hPtLep2_8 = new TH1F("hPtLep2_8", "Pt of Lep2 after selection step 8", 200 , -0.5 , 199.5 );
   hPtLep2_8->SetXTitle("pt_Lep2  (GeV)");
   TH1F * hYLep2_8 = new TH1F("hEtaLep2_8", "Y of Lep2 after selection step 8", 500 , -5. , 5. );
   hYLep2_8->SetXTitle("Y of Lep2");
   TH1F * hIsoLep2_8 = new TH1F("hIsoLep2_8", "Isolation of Lep2 after selection step 8", 2000 , -10. , 10. );
   hIsoLep2_8->SetXTitle("Iso");
   TH1F * hSipLep2_8 = new TH1F("hSipLep2_8", "Sip of Lep2 after selection step 8",  1000 , -20. , 40. );
   hSipLep2_8->SetXTitle("Sip");
   TH1F * hIpLep2_8 = new TH1F("hIpLep2_8", "Ip of Lep2 after selection step 8",  1000 , -20. , 40. );
   hIpLep2_8->SetXTitle("Ip");
   TH1F * hIpErLep2_8 = new TH1F("hIpErLep2_8", "IpEr of Lep2 after selection step 8",  1000 , 0. , 10. );
   hIpErLep2_8->SetXTitle("IpEr");

   TH1F * hPtLep3_8 = new TH1F("hPtLep3_8", "Pt of Lep3 after selection step 8", 200 , -0.5 , 199.5 );
   hPtLep3_8->SetXTitle("pt_Lep3  (GeV)");
   TH1F * hYLep3_8 = new TH1F("hEtaLep3_8", "Y of Lep3 after selection step 8", 500 , -5. , 5. );
   hYLep3_8->SetXTitle("Y of Lep2");
   TH1F * hIsoLep3_8 = new TH1F("hIsoLep3_8", "Isolation of Lep3 after selection step 8", 2000 , -10. , 10. );
   hIsoLep3_8->SetXTitle("Iso");
   TH1F * hSipLep3_8 = new TH1F("hSipLep3_8", "Sip of Lep3 after selection step 8",  1000 , -20. , 40. );
   hSipLep3_8->SetXTitle("Sip");
   TH1F * hIpLep3_8 = new TH1F("hIpLep3_8", "Ip of Lep3 after selection step 8",  1000 , -20. , 40. );
   hIpLep3_8->SetXTitle("Ip");
   TH1F * hIpErLep3_8 = new TH1F("hIpErLep3_8", "IpEr of Lep3 after selection step 8",  1000 , 0. , 10. );
   hIpErLep3_8->SetXTitle("IpEr");

   TH1F * hPtLep4_8 = new TH1F("hPtLep4_8", "Pt of Lep4 after selection step 8", 200 , -0.5 , 199.5 );
   hPtLep4_8->SetXTitle("pt_Lep4  (GeV)");
   TH1F * hYLep4_8 = new TH1F("hEtaLep4_8", "Y of Lep4 after selection step 8", 500 , -5. , 5. );
   hYLep4_8->SetXTitle("Y of Lep2");
   TH1F * hIsoLep4_8 = new TH1F("hIsoLep4_8", "Isolation of Lep4 after selection step 8", 2000 , -10. , 10. );
   hIsoLep4_8->SetXTitle("Iso");
   TH1F * hSipLep4_8 = new TH1F("hSipLep4_8", "Sip of Lep4 after selection step 8",  1000 , -20. , 40. );
   hSipLep4_8->SetXTitle("Sip");
   TH1F * hIpLep4_8 = new TH1F("hIpLep4_8", "Ip of Lep4 after selection step 8",  1000 , -20. , 40. );
   hIpLep4_8->SetXTitle("Ip");
   TH1F * hIpErLep4_8 = new TH1F("hIpErLep4_8", "IpEr of Lep4 after selection step 8",  1000 , 0. , 10. );
   hIpErLep4_8->SetXTitle("IpEr");
   

   TH1F * hIso_8 = new TH1F("hIso_8", "Isolation maxima after selection step 8", 2000 , -10. , 10. );
   hIso_8->SetXTitle("Iso");
   TH1F * hSip_8 = new TH1F("hSip_8", "Sip maxima after selection step 8",  1000 , -20. , 40. );
   hSip_8->SetXTitle("Sip");
   TH1F * hIp_8 = new TH1F("hIp_8", "Ip maxima after selection step 8",  1000 , -20. , 40. );
   hIp_8->SetXTitle("Ip");
   
   TH1F * hDjj_8 = new TH1F("hDjj_8", "Delta jets vbf selection step 8", 200, -19.5, 19.5 );
   hDjj_8->SetXTitle("Delta jets");
   TH1F * hMjj_8 = new TH1F("hMjj_8", "Mass jets vbf selection step 8", 200, -0.5, 499.5 );
   hMjj_8->SetXTitle("Mass jets");
   TH1F * hVD_8 = new TH1F("hVD_8", "Discriminant vbf selection step 8", 200, -0.5, 9.5 );
   hMjj_8->SetXTitle("Discriminant");
   TH2F * hMjj_Djj_8 = new TH2F("hMjj_Djj_8", "Mjj vs Djj selection step 8", 200, -0.5,499.5,200,-19.5, 19.5 );
   hMjj_Djj_8->SetXTitle("Mjj_vs_Djj");

   //no FSR   
   TH1F * hMZ1_noFSR_8 = new TH1F("hMZ1_noFSR_8", "Mass of Z1 after selection step 8 _noFSR", 200 , -0.5 , 199.5 );
   hMZ1_8->SetXTitle("mass_Z1  (GeV)");
   TH1F * hMZ2_noFSR_8 = new TH1F("hMZ2_noFSR_8", "Mass of Z2 after selection step 8 _noFSR", 200 , -0.5 , 199.5 );
   hMZ2_8->SetXTitle("mass_Z2  (GeV)");
   
   TH1F * hPFMET_8 = new TH1F("hPFMET_8", "PF MET after selection step 8", 1000 , 0., 1000.);
   hPFMET_8->SetXTitle("PF MET");
   
   TH1F * hM4l_T_8 = new TH1F("hM4l_T_8", "Transverse Mass of four leptons after full selection + MET", 1200, 4.5, 1204.5 );
   hM4l_T_8->SetXTitle("m_{T} + PF MET (GeV)");
   TH1F *hLogLinXM4l_T_8          = new TH1F("hLogLinXM4l_T_8","hLogLinXM4l_T_8",NMOBINS, loglinMbins);
   hLogLinXM4l_T_8->Sumw2();

   TH1F * hDPHI_8 = new TH1F("DPHI_8", "polar angle between 4l and E_{T,miss}", 1000, 0., 5. );
   hDPHI_8->SetXTitle("#DELTA#phi(4l,E_{T,miss})");
   
   TH1D * hNgood_8 = new TH1D("hNgood", "Number of good leptons", 10, -0.5, 9.5);
   hNgood_8->SetXTitle("# good leptons");
   TH1D * hNbjets_8 = new TH1D("hNbjets", "Number of b jets", 10, -0.5, 9.5);
   hNbjets_8->SetXTitle("# b-jets");    

   // Build the histo with constant log bin width   
   TH1F *hLogXPFMET_8             = new TH1F("hLogXPFMET_8","hLogXPFMET_8",NMBINS, logMbins);   
   hLogXPFMET_8->Sumw2();
   TH1F *hLogLinXPFMET_8          = new TH1F("hLogLinXPFMET_8","hLogLinXPFMET_8",NMOBINS, loglinMbins);
   hLogLinXPFMET_8->Sumw2();

    // Step 8a with m4l window cut
   TH1F * hPFMET_8a = new TH1F("hPFMET_8a", "PF MET after selection step 8a", 1000 , 0., 1000.);
   hPFMET_8a->SetXTitle("PF MET (GeV)");   
   TH1F *hLogLinXPFMET_8a          = new TH1F("hLogLinXPFMET_8a","hLogLinXPFMET_8a",NMOBINS, loglinMbins);
   hLogLinXPFMET_8a->Sumw2();

   TH1F * hM4l_T_8a = new TH1F("hM4l_T_8a", "Transverse Mass of four leptons after full selection + MET", 1200, 4.5, 1204.5 );
   hM4l_T_8a->SetXTitle("m_{T} + PF MET (GeV)");
   
   TH1F * hDPHI_8a = new TH1F("DPHI_8a", "polar angle between 4l and E_{T,miss}", 1000, 0., 5. );
   hDPHI_8a->SetXTitle("#DELTA#phi(4l,E_{T,miss})");

   TH1D * hNgood_8a = new TH1D("hNgood_8a", "Number of good leptons", 10, -0.5, 9.5);
   hNgood_8a->SetXTitle("# good leptons");
   TH1D * hNbjets_8a = new TH1D("hNbjets_8a", "Number of b jets", 10, -0.5, 9.5);
   hNbjets_8a->SetXTitle("# b-jets");

   // Step 8b with mZ window cut
   TH1F * hPFMET_8b = new TH1F("hPFMET_8b", "PF MET after selection step 8b", 1000 , 0., 1000.);
   hPFMET_8b->SetXTitle("PF MET (GeV)");   
   TH1F *hLogLinXPFMET_8b = new TH1F("hLogLinXPFMET_8b","hLogLinXPFMET_8b",NMOBINS, loglinMbins);
   hLogLinXPFMET_8b->Sumw2();   
   
   // Step 9 with PFMET cut   
   TH1F * hM4l_9 = new TH1F("hM4l_9", "Mass of four leptons after selection step 9", 1200, 4.5, 1204.5 );
   hM4l_9->SetXTitle("4 lepton mass  (GeV)");
   
   TH1F * hPFMET_9 = new TH1F("hPFMET_9", "PF MET after selection step 9", 1000 , 0., 1000.);
   hPFMET_9->SetXTitle("PF MET (GeV)");   
   
   TH1F * hM4l_T_9 = new TH1F("hM4l_T_9", "Transverse Mass of four leptons after full selection + MET", 1200, 4.5, 1204.5 );
   hM4l_T_9->SetXTitle("m_{T} + PF MET (GeV)");
   
   TH1F * hDPHI_9 = new TH1F("DPHI_9", "polar angle between 4l and E_{T,miss}", 1000, 0., 5. );
   hDPHI_9->SetXTitle("#DELTA#phi(4l,E_{T,miss})");

   TH1D * hNgood_9 = new TH1D("hNgood_9", "Number of good leptons", 10, -0.5, 9.5);
   hNgood_9->SetXTitle("# good leptons");
   TH1D * hNbjets_9 = new TH1D("hNbjets_9", "Number of b jets", 10, -0.5, 9.5);
   hNbjets_9->SetXTitle("# b-jets");


   // Step 10 with signal region cuts
   TH1F * hPFMET_10 = new TH1F("hPFMET_10", "PF MET after selection step 10", 1000 , 0., 1000.);
   hPFMET_10->SetXTitle("PF MET (GeV)");  
   hPFMET_10->Sumw2();

   TH1F * hPFMET_JetEnUp_10= new TH1F("hPFMET_JetEnUp_10", "PF MET after selection step 10", 1000 , 0., 1000.);
   hPFMET_JetEnUp_10->Sumw2();
   TH1F * hPFMET_JetEnDn_10= new TH1F("hPFMET_JetEnDn_10", "PF MET after selection step 10", 1000 , 0., 1000.);
   hPFMET_JetEnDn_10->Sumw2();
   TH1F * hPFMET_ElectronEnUp_10= new TH1F("hPFMET_ElectronEnUp_10", "PF MET after selection step 10", 1000 , 0., 1000.);
   hPFMET_ElectronEnUp_10->Sumw2();
   TH1F * hPFMET_ElectronEnDn_10= new TH1F("hPFMET_ElectronEnDn_10", "PF MET after selection step 10", 1000 , 0., 1000.);
   hPFMET_ElectronEnDn_10->Sumw2();
   TH1F * hPFMET_MuonEnUp_10= new TH1F("hPFMET_MuonEnUp_10", "PF MET after selection step 10", 1000 , 0., 1000.);
   hPFMET_MuonEnUp_10->Sumw2();
   TH1F * hPFMET_MuonEnDn_10= new TH1F("hPFMET_MuonEnDn_10", "PF MET after selection step 10", 1000 , 0., 1000.);
   hPFMET_MuonEnDn_10->Sumw2();
   TH1F * hPFMET_JetResUp_10= new TH1F("hPFMET_JetResUp_10", "PF MET after selection step 10", 1000 , 0., 1000.);
   hPFMET_JetResUp_10->Sumw2();
   TH1F * hPFMET_JetResDn_10= new TH1F("hPFMET_JetResDn_10", "PF MET after selection step 10", 1000 , 0., 1000.);
   hPFMET_JetResDn_10->Sumw2();
   TH1F * hPFMET_UnclusteredEnUp_10= new TH1F("hPFMET_UnclusteredEnUp_10", "PF MET after selection step 10", 1000 , 0., 1000.);
   hPFMET_UnclusteredEnUp_10->Sumw2();
   TH1F * hPFMET_UnclusteredEnDn_10= new TH1F("hPFMET_UnclusteredEnDn_10", "PF MET after selection step 10", 1000 , 0., 1000.);
   hPFMET_UnclusteredEnDn_10->Sumw2();
   TH1F * hPFMET_PhotonEnUp_10= new TH1F("hPFMET_PhotonEnUp_10", "PF MET after selection step 10", 1000 , 0., 1000.);
   hPFMET_PhotonEnUp_10->Sumw2();
   TH1F * hPFMET_PhotonEnDn_10= new TH1F("hPFMET_PhotonEnDn_10", "PF MET after selection step 10", 1000 , 0., 1000.);
   hPFMET_PhotonEnDn_10->Sumw2();

   TH1F *hLogLinXPFMET_10             = new TH1F("hLogLinXPFMET_10","hLogLinXPFMET_10",NMOBINS, loglinMbins);
   hLogLinXPFMET_10->Sumw2();

   TH1F * hLogLinXPFMET_JetEnUp_10= new TH1F("hLogLinXPFMET_JetEnUp_10", "hLogLinXPFMET_JetEnUp_10",NMOBINS, loglinMbins);
   hLogLinXPFMET_JetEnUp_10->Sumw2();
   TH1F * hLogLinXPFMET_JetEnDn_10= new TH1F("hLogLinXPFMET_JetEnDn_10", "hLogLinXPFMET_JetEnDn_10",NMOBINS, loglinMbins);
   hLogLinXPFMET_JetEnDn_10->Sumw2();
   TH1F * hLogLinXPFMET_ElectronEnUp_10= new TH1F("hLogLinXPFMET_ElectronEnUp_10", "hLogLinXPFMET_ElectronEnUp_10",NMOBINS, loglinMbins);
   hLogLinXPFMET_ElectronEnUp_10->Sumw2();
   TH1F * hLogLinXPFMET_ElectronEnDn_10= new TH1F("hLogLinXPFMET_ElectronEnDn_10", "hLogLinXPFMET_ElectronEnDn_10",NMOBINS, loglinMbins);
   hLogLinXPFMET_ElectronEnDn_10->Sumw2();
   TH1F * hLogLinXPFMET_MuonEnUp_10= new TH1F("hLogLinXPFMET_MuonEnUp_10", "hLogLinXPFMET_MuonEnUp_10",NMOBINS, loglinMbins);
   hLogLinXPFMET_MuonEnUp_10->Sumw2();
   TH1F * hLogLinXPFMET_MuonEnDn_10= new TH1F("hLogLinXPFMET_MuonEnDn_10", "hLogLinXPFMET_MuonEnDn_10",NMOBINS, loglinMbins);
   hLogLinXPFMET_MuonEnDn_10->Sumw2();
   TH1F * hLogLinXPFMET_JetResUp_10= new TH1F("hLogLinXPFMET_JetResUp_10", "hLogLinXPFMET_JetResUp_10",NMOBINS, loglinMbins);
   hLogLinXPFMET_JetResUp_10->Sumw2();
   TH1F * hLogLinXPFMET_JetResDn_10= new TH1F("hLogLinXPFMET_JetResDn_10", "hLogLinXPFMET_JetResDn_10",NMOBINS, loglinMbins);
   hLogLinXPFMET_JetResDn_10->Sumw2();
   TH1F * hLogLinXPFMET_UnclusteredEnUp_10= new TH1F("hLogLinXPFMET_UnclusteredEnUp_10", "hLogLinXPFMET_UnclusteredEnUp_10",NMOBINS, loglinMbins);
   hLogLinXPFMET_UnclusteredEnUp_10->Sumw2();
   TH1F * hLogLinXPFMET_UnclusteredEnDn_10= new TH1F("hLogLinXPFMET_UnclusteredEnDn_10", "hLogLinXPFMET_UnclusteredEnDn_10",NMOBINS, loglinMbins);
   hLogLinXPFMET_UnclusteredEnDn_10->Sumw2();
   TH1F * hLogLinXPFMET_PhotonEnUp_10= new TH1F("hLogLinXPFMET_PhotonEnUp_10", "hLogLinXPFMET_PhotonEnUp_10",NMOBINS, loglinMbins);
   hLogLinXPFMET_PhotonEnUp_10->Sumw2();
   TH1F * hLogLinXPFMET_PhotonEnDn_10= new TH1F("hLogLinXPFMET_PhotonEnDn_10", "hLogLinXPFMET_PhotonEnDn_10",NMOBINS, loglinMbins);
   hLogLinXPFMET_PhotonEnDn_10->Sumw2();

   TH1F * hM4l_T_10 = new TH1F("hM4l_T_10", "Transverse Mass of four leptons + MET after step 10", 1000, 0., 1000. );
   hM4l_T_10->SetXTitle("m_{T} + PF MET (GeV)");
   TH1F *hLogLinXM4l_T_10          = new TH1F("hLogLinXM4l_T_10","hLogLinXM4l_T_10",NMOBINS, loglinMbins);
   hLogLinXM4l_T_10->Sumw2();
   
   TH1F * hDPHI_10 = new TH1F("DPHI_10", "polar angle between 4l and E_{T,miss} after step 10", 1000, 0., 5. );
   hDPHI_10->SetXTitle("#DELTA#phi(4l,E_{T,miss})");

   TH1F * hPtLep1_10 = new TH1F("hPtLep1_10", "Pt of Lep1 after selection step 10", 200 , -0.5 , 199.5 );
   hPtLep1_10->SetXTitle("pt_Lep1  (GeV)");  

   TH1F * hPtLep2_10 = new TH1F("hPtLep2_10", "Pt of Lep2 after selection step 10", 200 , -0.5 , 199.5 );
   hPtLep2_10->SetXTitle("pt_Lep2  (GeV)");   

   TH1F * hPtLep3_10 = new TH1F("hPtLep3_10", "Pt of Lep3 after selection step 10", 200 , -0.5 , 199.5 );
   hPtLep3_10->SetXTitle("pt_Lep3  (GeV)");  

   TH1F * hPtLep4_10 = new TH1F("hPtLep4_10", "Pt of Lep4 after selection step 10", 200 , -0.5 , 199.5 );
   hPtLep4_10->SetXTitle("pt_Lep4  (GeV)");  

   TH1D * hNgood_10 = new TH1D("hNgood_10", "Number of good leptons", 10, -0.5, 9.5);
   hNgood_10->SetXTitle("# good leptons");

   TH1D * hNbjets_10 = new TH1D("hNbjets_10", "Number of b jets", 10, -0.5, 9.5);
   hNbjets_10->SetXTitle("# b-jets"); 

   // reverted mass peak cut
   TH1F * hPFMET_12 = new TH1F("hPFMET_12", "PF MET after selection step 12", 1000 , 0., 1000.);
   hPFMET_12->SetXTitle("PF MET (GeV)");  
   TH1F *hLogLinXPFMET_12 = new TH1F("hLogLinXPFMET_12","hLogLinXPFMET_12",NMOBINS, loglinMbins);
   hLogLinXPFMET_12->Sumw2();

   //global histos (during step 2..)
   
   TH1F * hN_loose_mu = new TH1F("hN_loose_mu", "N_loose_mu", 30 , 0. , 30. );
   hN_loose_mu->SetXTitle("N_loose_mu");
   TH1F * hN_loose_e = new TH1F("hN_loose_e", "N_loose_e", 30 , 0. , 30. );
   hN_loose_e->SetXTitle("N_loose_e");

   TH1F * hIso_loose_mu = new TH1F("hIso_loose_mu", "Isolation maxima after loose selection ", 2000 , -10. , 10. );
   hIso_loose_mu->SetXTitle("Iso");
   TH1F * hSip_loose_mu = new TH1F("hSip_loose_mu", "Sip maxima after loose selection ",  1000 , -20. , 40. );
   hSip_loose_mu->SetXTitle("Sip");
   TH1F * hIp_loose_mu = new TH1F("hIp_loose_mu", "Ip maxima after loose selection ",  1000 , -20. , 40. );
   hIp_loose_mu->SetXTitle("Ip");

   TH1F * hIso_loose_e = new TH1F("hIso_loose_e", "Isolation maxima after loose selection ", 2000 , -10. , 10. );
   hIso_loose_e->SetXTitle("Iso");
   TH1F * hSip_loose_e = new TH1F("hSip_loose_e", "Sip maxima after loose selection ",  1000 , -20. , 40. );
   hSip_loose_e->SetXTitle("Sip");
   TH1F * hIp_loose_e = new TH1F("hIp_loose_e", "Ip maxima after loose selection ",  1000 , -20. , 40. );
   hIp_loose_e->SetXTitle("Ip");


   TH1F * hN_good_lep = new TH1F("hN_good_lep", "N_good_lep", 30 , 0. , 30. );
   hN_good_lep->SetXTitle("N_good_lep");

   TH1F * hN_good_mu = new TH1F("hN_good_mu", "N_good_mu", 30 , 0. , 30. );
   hN_good_mu->SetXTitle("N_good_mu");
   TH1F * hN_good_ele = new TH1F("hN_good_ele", "N_good_ele", 30 , 0. , 30. );
   hN_good_ele->SetXTitle("N_good_ele");
   TH1F * hN_good_phot = new TH1F("hN_good_phot", "N_good_phot", 30 , 0. , 30. );
   hN_good_phot->SetXTitle("N_good_phot");

   TH1F * hMELA_8 = new TH1F("hMELA_8", "MELA after selection step 8", 300,-0.00166,1.00166 );
   hMELA_8->SetXTitle("MELA discriminant (4mu)");
   TH2F * hMELA_vs_M4l_8 = new TH2F("hMELA_vs_M4l_8", "MELA after selection step 8",1200, 4.5,1204.5,300,-0.00166,1.00166 );
   hMELA_vs_M4l_8->SetXTitle("MELA discriminant (4mu)");

   TH1F * hMELA_9 = new TH1F("hMELA_9", "MELA after selection step 9", 300,-0.00166,1.00166 );
   hMELA_9->SetXTitle("MELA discriminant (4mu)");
   TH2F * hMELA_vs_M4l_9 = new TH2F("hMELA_vs_M4l_9", "MELA after selection step 9",1200, 4.5,1204.5,300,-0.00166,1.00166 );
   hMELA_vs_M4l_9->SetXTitle("MELA discriminant (4mu)");
   
   TH1I * hVBF_PUID = new TH1I("hVBF_PUID", "PUID after step 8", 6,-1,4);
   hVBF_PUID->SetXTitle("PUID VBF");



   //PFJET Plots
   /* TH1I * hN_PFJET_6 = new TH1I("hN_PFJET_6", "Number of PFJets after step 6", 102, 0, 102);
   hN_PFJET_6->SetXTitle("Number of Jets");
   TH1F * hChg_PFJET_6 = new TH1F("hChg_PFJET_6", "Charge of PFJets after step 6", 50, -10., 10.);
   hChg_PFJET_6->SetXTitle("Charge");
   TH1F * hEt_PFJET_6 = new TH1F("hEt_PFJET_6", "Et of PFJets after step 6", 201, -0.5, 400.5 );
   hEt_PFJET_6->SetXTitle("Et_PFJET (GeV)");
   TH1F * hPt_PFJET_6 = new TH1F("hPt_PFJET_6", "Pt of PFJets after step 6", 201, -0.5, 400.5 );
   hPt_PFJET_6->SetXTitle("pt_PFJET (GeV)");
   TH1F * hEta_PFJET_6 = new TH1F("hEta_PFJET_6", "Eta of PFJets after step 6", 50, -10., 10.);
   hEta_PFJET_6->SetXTitle("eta_PFJET");
   TH1F * hPhi_PFJET_6 = new TH1F("hPhi_PFJET_6", "Phi of PFJets after step 6", 50, -10., 10.);
   hPhi_PFJET_6->SetXTitle("phi_PFJET");

   TH1I * hN_PFJET_8 = new TH1I("hN_PFJET_8", "Number of PFJets after step 8", 102, 0, 102);
   hN_PFJET_8->SetXTitle("Number of Jets");
   TH1F * hChg_PFJET_8 = new TH1F("hChg_PFJET_8", "Charge of PFJets after step 8", 50, -10., 10.);
   hChg_PFJET_8->SetXTitle("Charge");
   TH1F * hEt_PFJET_8 = new TH1F("hEt_PFJET_8", "Et of PFJets after step 8", 201, -0.5, 400.5 );
   hEt_PFJET_8->SetXTitle("Et_PFJET (GeV)");
   TH1F * hPt_PFJET_8 = new TH1F("hPt_PFJET_8", "Pt of PFJets after step 8", 201, -0.5, 400.5 );
   hPt_PFJET_8->SetXTitle("pt_PFJET (GeV)");
   TH1F * hEta_PFJET_8 = new TH1F("hEta_PFJET_8", "Eta of PFJets after step 8", 50, -10., 10.);
   hEta_PFJET_8->SetXTitle("eta_PFJET");
   TH1F * hPhi_PFJET_8 = new TH1F("hPhi_PFJET_8", "Phi of PFJets after step 8", 50, -10., 10.);
   hPhi_PFJET_8->SetXTitle("phi_PFJET");

   TH1I * hN_PFJET_VBF = new TH1I("hN_PFJET_VBF", "Number of PFJets after VBF selection", 102, 0, 102);
   hN_PFJET_VBF->SetXTitle("Number of Jets");
   TH1F * hChg_PFJET_VBF = new TH1F("hChg_PFJET_VBF", "Charge of PFJets after VBF selection", 50, -10., 10.);
   hChg_PFJET_VBF->SetXTitle("Charge");
   TH1F * hEt_PFJET_VBF = new TH1F("hEt_PFJET_VBF", "Et of PFJets after step VBF selection", 201, -0.5, 400.5 );
   hEt_PFJET_VBF->SetXTitle("Et_PFJET (GeV)");
   TH1F * hPt_PFJET_VBF = new TH1F("hPt_PFJET_VBF", "Pt of PFJets after step VBF selection", 201, -0.5, 400.5 );
   hPt_PFJET_VBF->SetXTitle("pt_PFJET (GeV)");
   TH1F * hEta_PFJET_VBF = new TH1F("hEta_PFJET_VBF", "Eta of PFJets after step VBF selection", 50, -10., 10.);
   hEta_PFJET_VBF->SetXTitle("eta_PFJET");
   TH1F * hPhi_PFJET_VBF = new TH1F("hPhi_PFJET_VBF", "Phi of PFJets after step VBF selection", 50, -10., 10.);
   hPhi_PFJET_VBF->SetXTitle("phi_PFJET");
   */


   // end book histo ***
      
   TTree *newtree = new TTree("HZZ4LeptonsAnalysisReduced", "reduced ttree");
  
    // Add branches to output rootuple 
   Float_t f_weight, f_int_weight, f_pu_weight, f_eff_weight, f_lept1_pt, f_lept1_eta, f_lept1_phi, f_lept1_charge, f_lept1_pfx, f_lept1_sip, f_lept1_mvaid, f_lept2_pt, f_lept2_eta, f_lept2_phi, f_lept2_charge, f_lept2_pfx, f_lept2_sip, f_lept2_mvaid, f_lept3_pt, f_lept3_eta, f_lept3_phi, f_lept3_charge, f_lept3_pfx, f_lept3_sip, f_lept3_mvaid, f_lept4_pt, f_lept4_eta, f_lept4_phi, f_lept4_charge, f_lept4_pfx, f_lept4_sip, f_lept4_mvaid, f_iso_max, f_sip_max, f_Z1mass, f_Z2mass, f_angle_costhetastar, f_angle_costheta1, f_angle_costheta2, f_angle_phi, f_angle_phistar1, f_eta4l, f_pt4l, f_mass4l, f_mass4lErr, f_njets_pass, f_deltajj, f_massjj, f_D_jet, f_jet1_pt, f_jet1_eta, f_jet1_phi, f_jet1_e, f_jet2_pt, f_jet2_eta, f_jet2_phi, f_jet2_e, f_jet1_highpt_pt, f_jet1_highpt_eta, f_jet1_highpt_phi, f_jet1_highpt_e, f_jet2_highpt_pt, f_jet2_highpt_eta, f_jet2_highpt_phi, f_jet2_highpt_e, f_jet3_highpt_pt, f_jet3_highpt_eta, f_jet3_highpt_phi, f_jet3_highpt_e;;
   Float_t f_D_bkg_kin,f_D_bkg,f_D_gg,f_D_g4,f_Djet_VAJHU; 
   Float_t f_genmet,f_pfmet,f_mT,f_dphi;
   Int_t f_lept1_pdgid,f_lept2_pdgid,f_lept3_pdgid,f_lept4_pdgid;
   Int_t f_category,f_Ngood,f_Nbjets;

   Int_t f_run, f_lumi, f_event;

   //@
   const int maxJets = 4;
   const int maxJetsComponents = 60;
   Int_t f_jets_highpt_charged_hadron_multiplicity[maxJets], f_jets_highpt_charged_multiplicity[maxJets], f_jets_highpt_component_pdgid[maxJets][maxJetsComponents], f_jets_highpt_electron_multiplicity[maxJets];
   Int_t f_jets_highpt_hf_em_multiplicity[maxJets], f_jets_highpt_hf_hadron_multiplicity[maxJets], f_jets_highpt_muon_multiplicity[maxJets], f_jets_highpt_ncomponents[maxJets];
   Int_t f_jets_highpt_neutral_hadron_multiplicity[maxJets], f_jets_highpt_neutral_multiplicity[maxJets], f_jets_highpt_photon_multiplicity[maxJets];
   Float_t f_jets_highpt_pt[maxJets], f_jets_highpt_eta[maxJets], f_jets_highpt_phi[maxJets], f_jets_highpt_et[maxJets], f_jets_highpt_area[maxJets], f_jets_highpt_charged_hadron_energy[maxJets], f_jets_highpt_neutral_hadron_energy[maxJets];
   Float_t f_jets_highpt_photon_energy[maxJets], f_jets_highpt_electron_energy[maxJets], f_jets_highpt_muon_energy[maxJets], f_jets_highpt_hf_hadron_energy[maxJets], f_jets_highpt_hf_em_energy[maxJets];
   Float_t f_jets_highpt_charged_em_energy[maxJets], f_jets_highpt_charged_mu_energy[maxJets], f_jets_highpt_neutral_em_energy[maxJets], f_jets_highpt_ptd[maxJets];
   Float_t f_jets_highpt_component_pt[maxJets][maxJetsComponents], f_jets_highpt_component_eta[maxJets][maxJetsComponents], f_jets_highpt_component_phi[maxJets][maxJetsComponents], f_jets_highpt_component_energy[maxJets][maxJetsComponents];
   Float_t f_jets_highpt_component_charge[maxJets][maxJetsComponents], f_jets_highpt_component_mt[maxJetsComponents][maxJetsComponents], f_jets_highpt_component_xvertex[maxJets][maxJetsComponents];
   Float_t f_jets_highpt_component_yvertex[maxJets][maxJetsComponents], f_jets_highpt_component_zvertex[maxJets][maxJetsComponents], f_jets_highpt_component_vertex_chi2[maxJets][maxJetsComponents];
   //@
   const int maxPartons = 8;
   bool f_lhe_parton_clear[maxPartons];
   Int_t f_lhe_npartons, f_lhe_parton_pdgid[maxPartons];
   Float_t f_lhe_parton_pt[maxPartons], f_lhe_parton_eta[maxPartons], f_lhe_parton_phi[maxPartons], f_lhe_parton_e[maxPartons];
   Float_t f_lept1_pt_error, f_lept2_pt_error, f_lept3_pt_error, f_lept4_pt_error, f_jets_highpt_pt_error[maxJets], f_Djet_VAJHU_UncUp, f_Djet_VAJHU_UncDn;
   Float_t f_pfmet_JetEnUp, f_pfmet_JetEnDn, f_pfmet_ElectronEnUp, f_pfmet_ElectronEnDn, f_pfmet_MuonEnUp, f_pfmet_MuonEnDn, f_pfmet_JetResUp, f_pfmet_JetResDn, f_pfmet_UnclusteredEnUp, f_pfmet_UnclusteredEnDn, f_pfmet_PhotonEnUp, f_pfmet_PhotonEnDn;  

   TBranch *b_run= newtree->Branch("f_run", &f_run,"f_run/I");
   TBranch *b_lumi= newtree->Branch("f_lumi", &f_lumi,"f_lumi/I");    
   TBranch *b_event= newtree->Branch("f_event", &f_event,"f_event/I");    
   
   TBranch *b_weight= newtree->Branch("f_weight", &f_weight,"f_weight/F");
   TBranch *b_int_weight= newtree->Branch("f_int_weight", &f_int_weight,"f_int_weight/F");
   TBranch *b_pu_weight= newtree->Branch("f_pu_weight", &f_pu_weight,"f_pu_weight/F");
   TBranch *b_eff_weight= newtree->Branch("f_eff_weight", &f_eff_weight,"f_eff_weight/F");
   TBranch *b_lept1_pt= newtree->Branch("f_lept1_pt", &f_lept1_pt,"f_lept1_pt/F");
   TBranch *b_lept1_pt_error= newtree->Branch("f_lept1_pt_error", &f_lept1_pt_error,"f_lept1_pt_error/F");
   TBranch *b_lept1_eta= newtree->Branch("f_lept1_eta", &f_lept1_eta,"f_lept1_eta/F");
   TBranch *b_lept1_phi= newtree->Branch("f_lept1_phi", &f_lept1_phi,"f_lept1_phi/F");
   TBranch *b_lept1_charge= newtree->Branch("f_lept1_charge", &f_lept1_charge,"f_lept1_charge/F");
   TBranch *b_lept1_pfx= newtree->Branch("f_lept1_pfx", &f_lept1_pfx,"f_lept1_pfx/F");
   TBranch *b_lept1_sip= newtree->Branch("f_lept1_sip", &f_lept1_sip,"f_lept1_sip/F");
   TBranch *b_lept1_pdgid= newtree->Branch("f_lept1_pdgid", &f_lept1_pdgid,"f_lept1_pdgid/I");
   TBranch *b_lept2_pt= newtree->Branch("f_lept2_pt", &f_lept2_pt,"f_lept2_pt/F");
   TBranch *b_lept2_pt_error= newtree->Branch("f_lept2_pt_error", &f_lept2_pt_error,"f_lept2_pt_error/F");
   TBranch *b_lept2_eta= newtree->Branch("f_lept2_eta", &f_lept2_eta,"f_lept2_eta/F");
   TBranch *b_lept2_phi= newtree->Branch("f_lept2_phi", &f_lept2_phi,"f_lept2_phi/F");
   TBranch *b_lept2_charge= newtree->Branch("f_lept2_charge", &f_lept2_charge,"f_lept2_charge/F");
   TBranch *b_lept2_pfx= newtree->Branch("f_lept2_pfx", &f_lept2_pfx,"f_lept2_pfx/F");
   TBranch *b_lept2_sip= newtree->Branch("f_lept2_sip", &f_lept2_sip,"f_lept2_sip/F");
   TBranch *b_lept2_pdgid= newtree->Branch("f_lept2_pdgid", &f_lept2_pdgid,"f_lept2_pdgid/I");
   TBranch *b_lept3_pt= newtree->Branch("f_lept3_pt", &f_lept3_pt,"f_lept3_pt/F");
   TBranch *b_lept3_pt_error= newtree->Branch("f_lept3_pt_error", &f_lept3_pt_error,"f_lept3_pt_error/F");
   TBranch *b_lept3_eta= newtree->Branch("f_lept3_eta", &f_lept3_eta,"f_lept3_eta/F");
   TBranch *b_lept3_phi= newtree->Branch("f_lept3_phi", &f_lept3_phi,"f_lept3_phi/F");
   TBranch *b_lept3_charge= newtree->Branch("f_lept3_charge", &f_lept3_charge,"f_lept3_charge/F");
   TBranch *b_lept3_pfx= newtree->Branch("f_lept3_pfx", &f_lept3_pfx,"f_lept3_pfx/F");
   TBranch *b_lept3_sip= newtree->Branch("f_lept3_sip", &f_lept3_sip,"f_lept3_sip/F");
   TBranch *b_lept3_pdgid= newtree->Branch("f_lept3_pdgid", &f_lept3_pdgid,"f_lept3_pdgid/I");
   TBranch *b_lept4_pt= newtree->Branch("f_lept4_pt", &f_lept4_pt,"f_lept4_pt/F");
   TBranch *b_lept4_pt_error= newtree->Branch("f_lept4_pt_error", &f_lept4_pt_error,"f_lept4_pt_error/F");
   TBranch *b_lept4_eta= newtree->Branch("f_lept4_eta", &f_lept4_eta,"f_lept4_eta/F");
   TBranch *b_lept4_phi= newtree->Branch("f_lept4_phi", &f_lept4_phi,"f_lept4_phi/F");
   TBranch *b_lept4_charge= newtree->Branch("f_lept4_charge", &f_lept4_charge,"f_lept4_charge/F");
   TBranch *b_lept4_pfx= newtree->Branch("f_lept4_pfx", &f_lept4_pfx,"f_lept4_pfx/F");
   TBranch *b_lept4_sip= newtree->Branch("f_lept4_sip", &f_lept4_sip,"f_lept4_sip/F");
   TBranch *b_lept4_pdgid= newtree->Branch("f_lept4_pdgid", &f_lept4_pdgid,"f_lept4_pdgid/I");
   TBranch *b_iso_max= newtree->Branch("f_iso_max", &f_iso_max,"f_iso_max/F");
   TBranch *b_sip_max= newtree->Branch("f_sip_max", &f_sip_max,"f_sip_max/F");
   TBranch *b_Z1mass= newtree->Branch("f_Z1mass", &f_Z1mass,"f_Z1mass/F");
   TBranch *b_Z2mass= newtree->Branch("f_Z2mass", &f_Z2mass,"f_Z2mass/F");
   TBranch *b_angle_costhetastar= newtree->Branch("f_angle_costhetastar", &f_angle_costhetastar,"f_angle_costhetastar/F");
   TBranch *b_angle_costheta1= newtree->Branch("f_angle_costheta1", &f_angle_costheta1,"f_angle_costheta1/F");
   TBranch *b_angle_costheta2= newtree->Branch("f_angle_costheta2", &f_angle_costheta2,"f_angle_costheta2/F");
   TBranch *b_angle_phi= newtree->Branch("f_angle_phi", &f_angle_phi,"f_angle_phi/F");
   TBranch *b_angle_phistar1= newtree->Branch("f_angle_phistar1", &f_angle_phistar1,"f_angle_phistar1/F");
   TBranch *b_pt4l= newtree->Branch("f_pt4l", &f_pt4l,"f_pt4l/F");
   TBranch *b_eta4l= newtree->Branch("f_eta4l", &f_eta4l,"f_eta4l/F");
   TBranch *b_mass4l= newtree->Branch("f_mass4l", &f_mass4l,"f_mass4l/F");
   TBranch *b_mass4lErr= newtree->Branch("f_mass4lErr", &f_mass4lErr,"f_mass4lErr/F");
   TBranch *b_njets_pass= newtree->Branch("f_njets_pass", &f_njets_pass,"f_njets_pass/F");
   TBranch *b_deltajj= newtree->Branch("f_deltajj", &f_deltajj,"f_deltajj/F");
   TBranch *b_massjj= newtree->Branch("f_massjj", &f_massjj,"f_massjj/F");
   TBranch *b_D_jet= newtree->Branch("f_D_jet", &f_D_jet,"f_D_jet/F");
   TBranch *b_jet1_pt= newtree->Branch("f_jet1_pt", &f_jet1_pt,"f_jet1_pt/F");
   TBranch *b_jet1_eta= newtree->Branch("f_jet1_eta", &f_jet1_eta,"f_jet1_eta/F");
   TBranch *b_jet1_phi= newtree->Branch("f_jet1_phi", &f_jet1_phi,"f_jet1_phi/F");
   TBranch *b_jet1_e= newtree->Branch("f_jet1_e", &f_jet1_e,"f_jet1_e/F");
   TBranch *b_jet2_pt= newtree->Branch("f_jet2_pt", &f_jet2_pt,"f_jet2_pt/F");
   TBranch *b_jet2_eta= newtree->Branch("f_jet2_eta", &f_jet2_eta,"f_jet2_eta/F");
   TBranch *b_jet2_phi= newtree->Branch("f_jet2_phi", &f_jet2_phi,"f_jet2_phi/F");
   TBranch *b_jet2_e= newtree->Branch("f_jet2_e", &f_jet2_e,"f_jet2_e/F");

   //
   TBranch *b_jet1_highpt_pt= newtree->Branch("f_jet1_highpt_pt", &f_jet1_highpt_pt,"f_jet1_highpt_pt/F");
   TBranch *b_jet1_highpt_eta= newtree->Branch("f_jet1_highpt_eta", &f_jet1_highpt_eta,"f_jet1_highpt_eta/F");
   TBranch *b_jet1_highpt_phi= newtree->Branch("f_jet1_highpt_phi", &f_jet1_highpt_phi,"f_jet1_highpt_phi/F");
   TBranch *b_jet1_highpt_e= newtree->Branch("f_jet1_highpt_e", &f_jet1_highpt_e,"f_jet1_highpt_e/F");

   TBranch *b_jet2_highpt_pt= newtree->Branch("f_jet2_highpt_pt", &f_jet2_highpt_pt,"f_jet2_highpt_pt/F");
   TBranch *b_jet2_highpt_eta= newtree->Branch("f_jet2_highpt_eta", &f_jet2_highpt_eta,"f_jet2_highpt_eta/F");
   TBranch *b_jet2_highpt_phi= newtree->Branch("f_jet2_highpt_phi", &f_jet2_highpt_phi,"f_jet2_highpt_phi/F");
   TBranch *b_jet2_highpt_e= newtree->Branch("f_jet2_highpt_e", &f_jet2_highpt_e,"f_jet2_highpt_e/F");

   TBranch *b_jet3_highpt_pt= newtree->Branch("f_jet3_highpt_pt", &f_jet3_highpt_pt,"f_jet3_highpt_pt/F");
   TBranch *b_jet3_highpt_eta= newtree->Branch("f_jet3_highpt_eta", &f_jet3_highpt_eta,"f_jet3_highpt_eta/F");
   TBranch *b_jet3_highpt_phi= newtree->Branch("f_jet3_highpt_phi", &f_jet3_highpt_phi,"f_jet3_highpt_phi/F");
   TBranch *b_jet3_highpt_e= newtree->Branch("f_jet3_highpt_e", &f_jet3_highpt_e,"f_jet3_highpt_e/F");
   //
   //@
   TBranch *b_jets_highpt_pt = newtree->Branch("f_jets_highpt_pt", &f_jets_highpt_pt,Form("f_jets_highpt_pt[%i]/F",maxJets));
   TBranch *b_jets_highpt_pt_error = newtree->Branch("f_jets_highpt_pt_error", &f_jets_highpt_pt_error,Form("f_jets_highpt_pt_error[%i]/F",maxJets));
   TBranch *b_jets_highpt_eta = newtree->Branch("f_jets_highpt_eta", &f_jets_highpt_eta,Form("f_jets_highpt_eta[%i]/F",maxJets));
   TBranch *b_jets_highpt_phi = newtree->Branch("f_jets_highpt_phi", &f_jets_highpt_phi,Form("f_jets_highpt_phi[%i]/F",maxJets));
   TBranch *b_jets_highpt_et = newtree->Branch("f_jets_highpt_et", &f_jets_highpt_et,Form("f_jets_highpt_et[%i]/F",maxJets));
   TBranch *b_jets_highpt_area = newtree->Branch("f_jets_highpt_area", &f_jets_highpt_area, Form("f_jets_highpt_area[%i]/F",maxJets));
   TBranch *b_jets_highpt_ptd = newtree->Branch("f_jets_highpt_ptd", &f_jets_highpt_ptd, Form("f_jets_highpt_ptd[%i]/F",maxJets));
   TBranch *b_jets_highpt_charged_hadron_energy = newtree->Branch("f_jets_highpt_charged_hadron_energy", &f_jets_highpt_charged_hadron_energy, Form("f_jets_highpt_charged_hadron_energy[%i]/F",maxJets));
   TBranch *b_jets_highpt_neutral_hadron_energy = newtree->Branch("f_jets_highpt_neutral_hadron_energy", &f_jets_highpt_neutral_hadron_energy, Form("f_jets_highpt_neutral_hadron_energy[%i]/F",maxJets));
   TBranch *b_jets_highpt_photon_energy = newtree->Branch("f_jets_highpt_photon_energy", &f_jets_highpt_photon_energy, Form("f_jets_highpt_photon_energy[%i]/F",maxJets));
   TBranch *b_jets_highpt_electron_energy = newtree->Branch("f_jets_highpt_electron_energy", &f_jets_highpt_electron_energy, Form("f_jets_highpt_electron_energy[%i]/F",maxJets));
   TBranch *b_jets_highpt_muon_energy = newtree->Branch("f_jets_highpt_muon_energy", &f_jets_highpt_muon_energy, Form("f_jets_highpt_muon_energy[%i]/F",maxJets));
   TBranch *b_jets_highpt_hf_hadron_energy = newtree->Branch("f_jets_highpt_hf_hadron_energy", &f_jets_highpt_hf_hadron_energy, Form("f_jets_highpt_hf_hadron_energy[%i]/F",maxJets));
   TBranch *b_jets_highpt_hf_em_energy = newtree->Branch("f_jets_highpt_hf_em_energy", &f_jets_highpt_hf_em_energy, Form("f_jets_highpt_hf_em_energy[%i]/F",maxJets));
   TBranch *b_jets_highpt_charged_em_energy = newtree->Branch("f_jets_highpt_charged_em_energy", &f_jets_highpt_charged_em_energy, Form("f_jets_highpt_charged_em_energy[%i]/F",maxJets));
   TBranch *b_jets_highpt_charged_mu_energy = newtree->Branch("f_jets_highpt_charged_mu_energy", &f_jets_highpt_charged_mu_energy, Form("f_jets_highpt_charged_mu_energy[%i]/F",maxJets));
   TBranch *b_jets_highpt_neutral_em_energy = newtree->Branch("f_jets_highpt_neutral_em_energy", &f_jets_highpt_neutral_em_energy, Form("f_jets_highpt_neutral_em_energy[%i]/F",maxJets));
   TBranch *b_jets_highpt_charged_hadron_multiplicity = newtree->Branch("f_jets_highpt_charged_hadron_multiplicity", &f_jets_highpt_charged_hadron_multiplicity, Form("f_jets_highpt_charged_hadron_multiplicity[%i]/I",maxJets));
   TBranch *b_jets_highpt_neutral_hadron_multiplicity = newtree->Branch("f_jets_highpt_neutral_hadron_multiplicity", &f_jets_highpt_neutral_hadron_multiplicity, Form("f_jets_highpt_neutral_hadron_multiplicity[%i]/I",maxJets));
   TBranch *b_jets_highpt_photon_multiplicity = newtree->Branch("f_jets_highpt_photon_multiplicity", &f_jets_highpt_photon_multiplicity, Form("f_jets_highpt_photon_multiplicity[%i]/I",maxJets));
   TBranch *b_jets_highpt_electron_multiplicity = newtree->Branch("f_jets_highpt_electron_multiplicity", &f_jets_highpt_electron_multiplicity, Form("f_jets_highpt_electron_multiplicity[%i]/I",maxJets));
   TBranch *b_jets_highpt_muon_multiplicity = newtree->Branch("f_jets_highpt_muon_multiplicity", &f_jets_highpt_muon_multiplicity, Form("f_jets_highpt_muon_multiplicity[%i]/I",maxJets));
   TBranch *b_jets_highpt_hf_hadron_multiplicity = newtree->Branch("f_jets_highpt_hf_hadron_multiplicity", &f_jets_highpt_hf_hadron_multiplicity, Form("f_jets_highpt_hf_hadron_multiplicity[%i]/I",maxJets));
   TBranch *b_jets_highpt_hf_em_multiplicity = newtree->Branch("f_jets_highpt_hf_em_multiplicity", &f_jets_highpt_hf_em_multiplicity, Form("f_jets_highpt_hf_em_multiplicity[%i]/I",maxJets));
   TBranch *b_jets_highpt_charged_multiplicity = newtree->Branch("f_jets_highpt_charged_multiplicity", &f_jets_highpt_charged_multiplicity, Form("f_jets_highpt_charged_multiplicity[%i]/I",maxJets));
   TBranch *b_jets_highpt_neutral_multiplicity = newtree->Branch("f_jets_highpt_neutral_multiplicity", &f_jets_highpt_neutral_multiplicity, Form("f_jets_highpt_neutral_multiplicity[%i]/I",maxJets));
   TBranch *b_jets_highpt_ncomponents = newtree->Branch("f_jets_highpt_ncomponents", &f_jets_highpt_ncomponents, Form("f_jets_highpt_ncomponents[%i]/I",maxJets));
   TBranch *b_jets_highpt_component_pdgid = newtree->Branch("f_jets_highpt_component_pdgid", &f_jets_highpt_component_pdgid, Form("f_jets_highpt_component_pdgid[%i][%i]/I",maxJets,maxJetsComponents));
   TBranch *b_jets_highpt_component_pt = newtree->Branch("f_jets_highpt_component_pt", &f_jets_highpt_component_pt, Form("f_jets_highpt_component_pt[%i][%i]/F",maxJets,maxJetsComponents));
   TBranch *b_jets_highpt_component_eta = newtree->Branch("f_jets_highpt_component_eta", &f_jets_highpt_component_eta, Form("f_jets_highpt_component_eta[%i][%i]/F",maxJets,maxJetsComponents));
   TBranch *b_jets_highpt_component_phi = newtree->Branch("f_jets_highpt_component_phi", &f_jets_highpt_component_phi, Form("f_jets_highpt_component_phi[%i][%i]/F",maxJets,maxJetsComponents));
   TBranch *b_jets_highpt_component_energy = newtree->Branch("f_jets_highpt_component_energy", &f_jets_highpt_component_energy, Form("f_jets_highpt_component_energy[%i][%i]/F",maxJets,maxJetsComponents));
   TBranch *b_jets_highpt_component_charge = newtree->Branch("f_jets_highpt_component_charge", &f_jets_highpt_component_charge, Form("f_jets_highpt_component_charge[%i][%i]/F",maxJets,maxJetsComponents));
   TBranch *b_jets_highpt_component_mt = newtree->Branch("f_jets_highpt_component_mt", &f_jets_highpt_component_mt, Form("f_jets_highpt_component_mt[%i][%i]/F",maxJets,maxJetsComponents));
   TBranch *b_jets_highpt_component_xvertex = newtree->Branch("f_jets_highpt_component_xvertex", &f_jets_highpt_component_xvertex, Form("f_jets_highpt_component_xvertex[%i][%i]/F",maxJets,maxJetsComponents));
   TBranch *b_jets_highpt_component_yvertex = newtree->Branch("f_jets_highpt_component_yvertex", &f_jets_highpt_component_yvertex, Form("f_jets_highpt_component_yvertex[%i][%i]/F",maxJets,maxJetsComponents));
   TBranch *b_jets_highpt_component_zvertex = newtree->Branch("f_jets_highpt_component_zvertex", &f_jets_highpt_component_zvertex, Form("f_jets_highpt_component_zvertex[%i][%i]/F",maxJets,maxJetsComponents));
   TBranch *b_jets_highpt_component_vertex_chi2 = newtree->Branch("f_jets_highpt_component_vertex_chi2", &f_jets_highpt_component_vertex_chi2, Form("f_jets_highpt_component_vertex_chi2[%i][%i]/F",maxJets,maxJetsComponents));
   //@
   TBranch *b_lhe_npartons = newtree->Branch("f_lhe_npartons", &f_lhe_npartons, "f_lhe_npartons/I");
   TBranch *b_lhe_parton_clear = newtree->Branch("f_lhe_parton_clear", &f_lhe_parton_clear, Form("f_lhe_parton_clear[%i]/b",maxPartons));
   TBranch *b_lhe_parton_pdgid = newtree->Branch("f_lhe_parton_pdgid", &f_lhe_parton_pdgid, Form("f_lhe_parton_pdgid[%i]/I",maxPartons));
   TBranch *b_lhe_parton_pt = newtree->Branch("f_lhe_parton_pt", &f_lhe_parton_pt, Form("f_lhe_parton_pt[%i]/F",maxPartons));
   TBranch *b_lhe_parton_eta = newtree->Branch("f_lhe_parton_eta", &f_lhe_parton_eta, Form("f_lhe_parton_eta[%i]/F",maxPartons));
   TBranch *b_lhe_parton_phi = newtree->Branch("f_lhe_parton_phi", &f_lhe_parton_phi, Form("f_lhe_parton_phi[%i]/F",maxPartons));
   TBranch *b_lhe_parton_e = newtree->Branch("f_lhe_parton_e", &f_lhe_parton_e, Form("f_lhe_parton_e[%i]/F",maxPartons));
   

   TBranch *b_D_bkg_kin= newtree->Branch("f_D_bkg_kin", &f_D_bkg_kin,"f_D_bkg_kin/F");
   TBranch *b_D_bkg= newtree->Branch("f_D_bkg", &f_D_bkg,"f_D_bkg/F");
   TBranch *b_D_gg= newtree->Branch("f_D_gg", &f_D_gg,"f_D_gg/F");
   TBranch *b_D_g4= newtree->Branch("f_D_g4", &f_D_g4,"f_D_g4/F");
   TBranch *b_Djet_VAJHU= newtree->Branch("f_Djet_VAJHU", &f_Djet_VAJHU,"f_Djet_VAJHU/F");
   TBranch *b_Djet_VAJHU_UncUp = newtree->Branch("f_Djet_VAJHU_UncUp", &f_Djet_VAJHU_UncUp,"f_Djet_VAJHU_UncUp/F");
   TBranch *b_Djet_VAJHU_UncDn = newtree->Branch("f_Djet_VAJHU_UncDn", &f_Djet_VAJHU_UncDn,"f_Djet_VAJHU_UncDn/F");

   TBranch *b_pfmet= newtree->Branch("f_pfmet", &f_pfmet,"f_pfmet/F");
    TBranch *b_pfmet_JetEnUp = newtree->Branch("f_pfmet_JetEnUp", &f_pfmet_JetEnUp, "f_pfmet_JetEnUp/F");
   TBranch *b_pfmet_JetEnDn = newtree->Branch("f_pfmet_JetEnDn", &f_pfmet_JetEnDn, "f_pfmet_JetEnDn/F");
   TBranch *b_pfmet_ElectronEnUp = newtree->Branch("f_pfmet_ElectronEnUp", &f_pfmet_ElectronEnUp, "f_pfmet_ElectronEnUp/F");
   TBranch *b_pfmet_ElectronEnDn = newtree->Branch("f_pfmet_ElectronEnDn", &f_pfmet_ElectronEnDn, "f_pfmet_ElectronEnDn/F");
   TBranch *b_pfmet_MuonEnUp = newtree->Branch("f_pfmet_MuonEnUp", &f_pfmet_MuonEnUp, "f_pfmet_MuonEnUp/F");
   TBranch *b_pfmet_MuonEnDn = newtree->Branch("f_pfmet_MuonEnDn", &f_pfmet_MuonEnDn, "f_pfmet_MuonEnDn/F");
   TBranch *b_pfmet_JetResUp = newtree->Branch("f_pfmet_JetResUp", &f_pfmet_JetResUp, "f_pfmet_JetResUp/F");
   TBranch *b_pfmet_JetResDn = newtree->Branch("f_pfmet_JetResDn", &f_pfmet_JetResDn, "f_pfmet_JetResDn/F");
   TBranch *b_pfmet_UnclusteredEnUp = newtree->Branch("f_pfmet_UnclusteredEnUp", &f_pfmet_UnclusteredEnUp, "f_pfmet_UnclusteredEnUp/F");
   TBranch *b_pfmet_UnclusteredEnDn = newtree->Branch("f_pfmet_UnclusteredEnDn", &f_pfmet_UnclusteredEnDn, "f_pfmet_UnclusteredEnDn/F");
   TBranch *b_pfmet_PhotonEnUp = newtree->Branch("f_pfmet_PhotonEnUp", &f_pfmet_PhotonEnUp, "f_pfmet_PhotonEnUp/F");
   TBranch *b_pfmet_PhotonEnDn = newtree->Branch("f_pfmet_PhotonEnDn", &f_pfmet_PhotonEnDn, "f_pfmet_PhotonEnDn/F");
   
   TBranch *b_genmet= newtree->Branch("f_genmet", &f_genmet,"f_genmet/F");
   TBranch *b_f_mT= newtree->Branch("f_mT", &f_mT,"f_mT/F");
   TBranch *b_f_dphi= newtree->Branch("f_dphi", &f_dphi,"f_dphi/F");
   TBranch *b_f_category= newtree->Branch("f_category", &f_category,"f_category/I");
   TBranch *b_f_Ngood= newtree->Branch("f_Ngood", &f_Ngood,"f_Ngood/I");
   TBranch *b_f_Nbjets=newtree->Branch("f_Nbjets", &f_Nbjets, "f_Nbjets/I");

   float newweight=1.;
   
   // New tree with clone of events passing the final selection
   TFile *skimfile = new TFile(eventsOUT,"recreate");
   // Clone tree for final events
   TTree *finaltree = fChain->CloneTree(0);

   // loop on entries
   
   Long64_t nentries = fChain->GetEntries();

   cout << "\n****************************"  <<endl;
   cout << "Analyzing " << nentries << " entries"  <<endl;     

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // Initialize reduced tree variables
      f_weight=1., f_int_weight=1., f_pu_weight=1., f_eff_weight=1., f_lept1_pt=-999., f_lept1_eta=-999., f_lept1_phi=-999., f_lept1_charge=-999., f_lept1_pfx=-999., f_lept1_sip=-999., f_lept1_mvaid=-999., f_lept2_pt=-999., f_lept2_eta=-999., f_lept2_phi=-999., f_lept2_charge=-999., f_lept2_pfx=-999., f_lept2_sip=-999., f_lept2_mvaid=-999., f_lept3_pt=-999., f_lept3_eta=-999., f_lept3_phi=-999., f_lept3_charge=-999., f_lept3_pfx=-999., f_lept3_sip=-999., f_lept3_mvaid=-999., f_lept4_pt=-999., f_lept4_eta=-999., f_lept4_phi=-999., f_lept4_charge=-999., f_lept4_pfx=-999., f_lept4_sip=-999., f_lept4_mvaid=-999., f_iso_max=-999., f_sip_max=-999., f_Z1mass=-999., f_Z2mass=-999., f_angle_costhetastar=-999., f_angle_costheta1=-999., f_angle_costheta2=-999., f_angle_phi=-999., f_angle_phistar1=-999., f_eta4l=-999., f_pt4l=-999., f_mass4l=-999., f_mass4lErr=-999., f_njets_pass=-999., f_deltajj=-999., f_massjj=-999., f_D_jet=-999., f_jet1_pt=-999., f_jet1_eta=-999., f_jet1_phi=-999., f_jet1_e=-999., f_jet2_pt=-999., f_jet2_eta=-999., f_jet2_phi=-999., f_jet2_e=-999.,f_D_bkg_kin=-999.,f_D_bkg=-999.,f_D_gg=-999.,f_D_g4=-999.,f_Djet_VAJHU=-999.,f_genmet=-999.,f_pfmet=-999.,f_mT=-999.,f_dphi=-999.,f_lept1_pdgid=-999,f_lept2_pdgid=-999,f_lept3_pdgid=-999,f_lept4_pdgid=-999,f_run=-999, f_lumi=-999, f_event=-999, f_category=-999, f_Ngood=-999, f_Nbjets=-999,f_jet1_highpt_pt=-999., f_jet1_highpt_eta=-999., f_jet1_highpt_phi=-999., f_jet1_highpt_e=-999., f_jet2_highpt_pt=-999., f_jet2_highpt_eta=-999., f_jet2_highpt_phi=-999., f_jet2_highpt_e=-999.,f_jet3_highpt_pt=-999., f_jet3_highpt_eta=-999., f_jet3_highpt_phi=-999., f_jet3_highpt_e=-999.;

       //@
      f_lept1_pt_error=-999., f_lept2_pt_error=-999., f_lept3_pt_error=-999., f_lept4_pt_error=-999.;
      f_Djet_VAJHU_UncUp = -999, f_Djet_VAJHU_UncDn = -999;
      f_pfmet_JetEnUp = -999., f_pfmet_JetEnDn = -999., f_pfmet_ElectronEnUp = -999., f_pfmet_ElectronEnDn = -999., f_pfmet_MuonEnUp = -999., f_pfmet_MuonEnDn = -999., f_pfmet_JetResUp = -999., f_pfmet_JetResDn = -999., f_pfmet_UnclusteredEnUp = -999., f_pfmet_UnclusteredEnDn = -999., f_pfmet_PhotonEnUp = -999., f_pfmet_PhotonEnDn = -999.;
      for(unsigned int ijet=0; ijet<maxJets; ++ijet){
        f_jets_highpt_charged_hadron_multiplicity[ijet] = -999;
        f_jets_highpt_charged_multiplicity[ijet] = -999;
        f_jets_highpt_electron_multiplicity[ijet] = -999;
        f_jets_highpt_hf_em_multiplicity[ijet] = -999;
        f_jets_highpt_hf_hadron_multiplicity[ijet] = -999;
        f_jets_highpt_muon_multiplicity[ijet] = -999;
        f_jets_highpt_ncomponents[ijet] = -999;
        f_jets_highpt_neutral_hadron_multiplicity[ijet] = -999;
        f_jets_highpt_neutral_multiplicity[ijet] = -999;
        f_jets_highpt_photon_multiplicity[ijet] = -999;
        f_jets_highpt_pt[ijet] = -999;
        f_jets_highpt_eta[ijet] = -999;
        f_jets_highpt_phi[ijet] = -999;
        f_jets_highpt_et[ijet] = -999;
        f_jets_highpt_pt_error[ijet] = -999;
        f_jets_highpt_area[ijet] = -999;
        f_jets_highpt_charged_hadron_energy[ijet] = -999;
        f_jets_highpt_neutral_hadron_energy[ijet] = -999;
        f_jets_highpt_photon_energy[ijet] = -999;
        f_jets_highpt_electron_energy[ijet] = -999;
        f_jets_highpt_muon_energy[ijet] = -999;
        f_jets_highpt_hf_hadron_energy[ijet] = -999;
        f_jets_highpt_hf_em_energy[ijet] = -999;
        f_jets_highpt_charged_em_energy[ijet] = -999;
        f_jets_highpt_charged_mu_energy[ijet] = -999;
        f_jets_highpt_neutral_em_energy[ijet] = -999;
        f_jets_highpt_ptd[ijet] = -999;
        for(unsigned int ijetc=0; ijetc<maxJetsComponents; ++ijetc){
          f_jets_highpt_component_pdgid[ijet][ijetc] = -999;
          f_jets_highpt_component_pt[ijet][ijetc] = -999;
          f_jets_highpt_component_eta[ijet][ijetc] = -999;
          f_jets_highpt_component_phi[ijet][ijetc] = -999;
          f_jets_highpt_component_energy[ijet][ijetc] = -999;
          f_jets_highpt_component_charge[ijet][ijetc] = -999;
          f_jets_highpt_component_mt[ijet][ijetc] = -999;
          f_jets_highpt_component_xvertex[ijet][ijetc] = -999;
          f_jets_highpt_component_yvertex[ijet][ijetc] = -999;
          f_jets_highpt_component_zvertex[ijet][ijetc] = -999;
          f_jets_highpt_component_vertex_chi2[ijet][ijetc] = -999;
        }
      }
      //@
      f_lhe_npartons = LHE_PARTON_N;
      for(unsigned int iparton=0; iparton<maxPartons; ++iparton){
	if(iparton < LHE_PARTON_N){
	  unsigned int index = -1;
	  float max_pt = -1;
	  //To sort the partons by decreasing pt
	  for(int jparton=0; jparton<LHE_PARTON_N; ++jparton){
	    if(LHE_PARTON_PT[jparton] > max_pt){
	      max_pt = LHE_PARTON_PT[jparton];
	      index = jparton;
	    }
	  }
	  
	  f_lhe_parton_clear[iparton] = LHE_PARTON_CLEAR[index];
	  f_lhe_parton_pdgid[iparton] = LHE_PARTON_PDGID[index];
	  f_lhe_parton_pt[iparton]    = LHE_PARTON_PT[index];
	  f_lhe_parton_eta[iparton]   = LHE_PARTON_ETA[index];
	  f_lhe_parton_phi[iparton]   = LHE_PARTON_PHI[index];
	  f_lhe_parton_e[iparton]     = LHE_PARTON_E[index];
	  LHE_PARTON_PT[index] = -2;
	}
	else{
	  f_lhe_parton_clear[iparton] = -999;
	  f_lhe_parton_pdgid[iparton] = -999;
	  f_lhe_parton_pt[iparton]    = -999;
	  f_lhe_parton_eta[iparton]   = -999;
	  f_lhe_parton_phi[iparton]   = -999;
	  f_lhe_parton_e[iparton]     = -999;
	}
      }

      //Reham 1:2920:259836
      //	if (!(Run==1 && LumiSection==2254 && Event==184363)) continue;	      
        cout<<"This is run= "<<Run<<"Lumi= "<<LumiSection<<"Event= "<<Event<<endl;

      //if(!(Run==257531 && LumiSection==121 && Event==178979541)) continue;
 
      if(jentry%1 == 5000) cout << "Analyzing entry: " << jentry << endl;
      

      if( RECO_NMU > 100 ) RECO_NMU = 100;
      if( RECO_NELE > 100 ) RECO_NELE = 100;
      if( RECO_NPFPHOT > 20 ) RECO_NPFPHOT = 20;
      
      bool debug=true;  //debug flag  -- default false

      newweight=weight;
      cout << "Starting weight= " << newweight << endl;

      // pT gen reweighting                                                                                                                                                                           
      double pT_weight=1.;
      if ( ( datasetBase.Contains("2HDM") || datasetBase.Contains("Baryonic") ) && datasetBase.Contains("Target") ){
        if (MC_PT[0] > 0.) {
          Int_t binx = pTweighthisto->GetXaxis()->FindBin(MC_PT[0]);
          cout << " bin x= " << binx << " " << pTweighthisto->GetBinContent(binx) << endl;
          pT_weight=double(pTweighthisto->GetBinContent(binx));
        }
        else {
          cout << "No pT reweighting of 2HDM samples" << endl;
        }
      }
      // Changing the weight for gen reweighting                                                                                                                                        
      newweight=weight*pT_weight;
      cout << "Starting weight + pT gen = " << newweight << endl;

      // pileup reweighting 2012 and 2011
      if (DATA_type=="NO" && num_PU_vertices < 0) continue; 

      hPUvertices->Fill(num_PU_vertices,weight);

      // pileup reweighting 2016
      hPUvertices->Fill(num_PU_vertices,weight);

      double pu_weight=1.;
      if (MC_type == "Fall17"){
	Int_t binx = puweight->GetXaxis()->FindBin(num_PU_vertices);
	cout << " bin x= " << binx << " " << puweight->GetBinContent(binx) << endl;	
	pu_weight=double(puweight->GetBinContent(binx));
	
      }      
       
      hPUvertices_ReWeighted->Fill(num_PU_vertices,weight*pu_weight);
      cout << "Pileup interations and weight is= " << num_PU_vertices << " " << " and weight= " << pu_weight << endl;  
      
      //if (num_PU_vertices < 0) continue;

      // Changing the weight for pileup
      newweight=weight*pT_weight*pu_weight;
      cout << "Starting weight + pTgen + pileup = " << newweight << endl;

      // ggZZ kfactor
      double ggzz_kf_wgt[9];
      float  weight_kfactor=1.;
      if (DATA_type=="NO"){
        //if( datasetName.Contains("GluGluHToZZ")){
        //  for(int f=0;f<9;f++) ggzz_kf_wgt[f] = ggZZ_kf[f]->Eval(MC_MASS[0]); // Evaluate at the true m4l 
        //  weight_kfactor=ggzz_kf_wgt[0]; // Using the nominal one                                       
        //  newweight=weight*pT_weight*pu_weight*weight_kfactor;
        //}
        if ( datasetName.Contains("GluGluToZZ") ||  datasetName.Contains("GluGluToContinToZZ")){
          for(int f=0;f<9;f++) ggzz_kf_wgt[f] = ggZZ_kf[f]->Eval(MC_ZZ_MASS[0][0]); // Evaluate at the true m4l    
          weight_kfactor=ggzz_kf_wgt[0]; // Using the nominal one                                                  
          newweight=weight*pT_weight*pu_weight*weight_kfactor;
        }
        //weight_kfactor=2.3;                                                                                                   
      }

      // qqZZ kfactor
      double qqzz_kf_wgt;
      //weight_kfactor=1.;
      int finalState=-999;
      if (DATA_type=="NO"){
	if( datasetName.Contains("ZZTo4L") )  {	
	  for (int l=0;l<4;l++){
	    if (MC_ZZ_MASS[l][0]>0. &&
		fabs(MC_ZZ_PDGID[l][3])==fabs(MC_ZZ_PDGID[l][4]) && 
		fabs(MC_ZZ_PDGID[l][3])==fabs(MC_ZZ_PDGID[l][5]) &&
		fabs(MC_ZZ_PDGID[l][3])==fabs(MC_ZZ_PDGID[l][6])) finalState=1; // 4e, 4mu, 4tau
	    else finalState=2;
	    weight_kfactor=HZZ4LeptonsAnalysis::kfactor_qqZZ_qcd_M(MC_ZZ_MASS[l][0],finalState);
	    newweight=weight*pT_weight*pu_weight*weight_kfactor;
	  }	
	}
      }

      // Weight for MCNLO samples                                                                                      
      if( datasetName.Contains("amcatnlo")) {
        cout << "Reweighting sample of amcatnlo with weight= " << MC_weighting << endl;
        if (MC_weighting!=0) newweight=weight*pT_weight*pu_weight*weight_kfactor*MC_weighting;
      }

                 
      float pFill[11];for(int pf=0;pf<11;pf++)pFill[11]=-999.;

      // ** Step 0:
      // simply number of entries...
      if( debug ) cout << "\n** Step 0: \nAnalyzing entry: " << jentry << " Run: " << Run << " Event: " << Event << " LumiSection: " << LumiSection << endl ;
      ++N_0 ;  // fill counter
      N_0_w=N_0_w+newweight;
      hweight->Fill(newweight);

      // ** Step 0.1:
      // number of 4L (4mu)...
      
      //MC_PDGID is the golden variable
      //MC_PDGID[0]  the Higgs25
      //MC_PDGID[1]  the Z1   23
      //MC_PDGID[2]  the Z2   23
      //MC_PDGID[3]  the l1+  mu- = 13  e- = 11
      //MC_PDGID[4]  the l2-  mu+ = -13  e+ = -11
      //MC_PDGID[5]  the l3+  mu- = 13  e- = 11
      //MC_PDGID[6]  the l4-  mu+ = -13  e+ = -11
     
      bool is4e=false;

      if( isSignal && debug ) cout << "\n ** Step 0.1 (MC truth): "
				   << "\nMC_PDGID[0] " << MC_PDGID[0]
				   << "\nMC_PDGID[1] " << MC_PDGID[1]	
				   << "\nMC_PDGID[2] " << MC_PDGID[2]	
				   << "\nMC_PDGID[3] " << MC_PDGID[3]	
				   << "\nMC_PDGID[4] " << MC_PDGID[4]	
				   << "\nMC_PDGID[5] " << MC_PDGID[5]	
				   << "\nMC_PDGID[6] " << MC_PDGID[6]	
				   << endl ;
      
      
      //
      if (isSignal || isggH){
	if( MC_PDGID[3] == 11 && MC_PDGID[4] == -11  
	    && MC_PDGID[5] == 11 && MC_PDGID[6] == -11 ){
	  is4e=true;
	  ++N_01 ;  // fill counter 
	  N_01_w=N_01_w+newweight;
	}
	//else continue;
	
	
	// ** Step 0.2:
	// events in the ACCEPTANCE:
	// variables contributing: eta & pt
	// NEW only eta cut
	bool isInAcceptance = 1;
	
	for( int i = 3; i  < 7; ++i){
	  
	  if( debug ) cout << "\n ** Step 0.2 (acceptance): "
			   << "\n fabs(MC_ETA[i])" << fabs(MC_ETA[i])	
			   << "\n MC_PT[i]" << MC_PT[i]	
			   << endl ;
	  
	  //if( (isSignal || isvbf) && (fabs(MC_ETA[i]) >= 2.4 || MC_PT[i] <= 5.) ){ 
	  if( (fabs(MC_ETA[i]) >= 2.4) ){ // only eta cut NEW
	    isInAcceptance = 0;
	    break;        
	  }
	  
	}// end for leptons
	
	if( is4e && isInAcceptance )       {
	  ++N_02 ;  // fill counter
	  N_02_w=N_02_w+newweight;
	}
	//else continue;	
      }
      
      ++N_1 ;  // fill counter
      N_1_w=N_1_w+newweight;
      

      // Effective AREA
      bool tag_2011=false;
      if (DATA_type=="2010" || DATA_type=="2011" || MC_type=="Fall11"){
	tag_2011=true;
      }
      
      // Loose lepton identification
      
      int N_loose_mu = 0;
      int N_loose_e = 0;
      double max_Iso_loose_mu = -1 ;
      double max_Sip_loose_mu = -1 ;
      double max_Ip_loose_mu = -1 ;
      double max_Iso_loose_e = -1 ;
      double max_Sip_loose_e = -1 ;
      double max_Ip_loose_e = -1 ;
      
      int* arraysize_mu = new int[1];
      arraysize_mu[0] = RECO_NMU;
      int iL_loose_mu[arraysize_mu[0]];
      delete [] arraysize_mu;

      for( int i = 0; i < RECO_NMU; ++i ){
	iL_loose_mu[i]=-999.;
      }
 
      for( int i = 0; i < RECO_NMU; ++i ){

        if( debug ) cout << "\n Lepton i="<< i <<" properties: "
			 << "\nRECOMU_isGlobalMu[i] " << int(RECOMU_isGlobalMu[i])
			 << "\nRECOMU_isTrackerMu[i] " << int(RECOMU_isTrackerMu[i])
			 << "\nRECOMU_PT[i] " << RECOMU_PT[i]
			 << "\nfabs(RECOMU_ETA[i]) " << fabs(RECOMU_ETA[i])
			 << "\nfabs( RECOMU_mubestrkDxy[i] ) " << fabs( RECOMU_mubesttrkDxy[i] )
			 << "\nfabs( RECOMU_mubesttrkDz[i] ) " << fabs( RECOMU_mubesttrkDz[i] )
			 << endl ;
       	
 	if( ( RECOMU_isGlobalMu[i] || (RECOMU_isTrackerMu[i] && RECOMU_numberOfMatches[i]>0) )
	    && RECOMU_mubesttrkType[i]!=2
	    && RECOMU_PT[i] > 5. 
	    && fabs(RECOMU_ETA[i]) < 2.4 
	    && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1.
	    ){ 
	  iL_loose_mu[N_loose_mu]=i;
	  ++N_loose_mu ;
	  if( RECOMU_PFX_dB[i] > max_Iso_loose_mu ) max_Iso_loose_mu = RECOMU_PFX_dB[i] ;
	  if( fabs( RECOMU_SIP[i] ) > max_Sip_loose_mu ) max_Sip_loose_mu = fabs( RECOMU_SIP[i] ) ;
	  if( fabs( RECOMU_IP[i] ) > max_Ip_loose_mu ) max_Ip_loose_mu = fabs( RECOMU_IP[i] ) ;
	}
	
      } // end loop on muons

      
      int* arraysize_e = new int[1];
      arraysize_e[0] = RECO_NELE;
      int iL_loose_e[arraysize_e[0]];
      delete [] arraysize_e;

      for( int i = 0; i < RECO_NELE; ++i ){
	iL_loose_e[i]=-999.;
      }

      for( int i = 0; i < RECO_NELE; ++i ){
	
        if( debug ) cout << "\n Lepton i="<< i <<" properties: "
			 << "\nRECOELE_PT[i] " << RECOELE_PT[i]
			 << "\nfabs(RECOELE_ETA[i]) " << fabs(RECOELE_ETA[i])
			 << "\nfabs( RECOELE_gsftrack_dxy[i] ) " << fabs( RECOELE_gsftrack_dxy[i] )
			 << "\nfabs( RECOELE_gsftrack_dz[i] ) " << fabs( RECOELE_gsftrack_dz[i] )
			 << endl ;
       	
 	if( RECOELE_PT[i] > 7. 
	    && fabs(RECOELE_ETA[i]) < 2.5 
	    //&& RECOELE_gsftrack_expected_inner_hits[i]<=1  // nto used anymore
	    && fabs(RECOELE_gsftrack_dxy[i]) < .5 
	    && fabs(RECOELE_gsftrack_dz[i]) < 1. 
	    ) {	  
	  iL_loose_e[N_loose_e]=i;
	  ++N_loose_e ;
	  if( RECOELE_PFX_rho[i] > max_Iso_loose_e ) max_Iso_loose_e = RECOELE_PFX_rho[i] ;
	  if( fabs( RECOELE_SIP[i] ) > max_Sip_loose_e ) max_Sip_loose_e = fabs( RECOELE_SIP[i] ) ;
	  if( fabs( RECOELE_IP[i] ) > max_Ip_loose_e ) max_Ip_loose_e = fabs( RECOELE_IP[i] ) ;
	}
	
      }// end loop on electrons
      
      hN_loose_mu->Fill( N_loose_mu,newweight );
      hN_loose_e->Fill( N_loose_e,newweight );
      hIso_loose_mu->Fill( max_Iso_loose_mu,newweight);
      hSip_loose_mu->Fill( max_Sip_loose_mu,newweight );
      hIp_loose_mu->Fill( max_Ip_loose_mu,newweight );
      hIso_loose_e->Fill( max_Iso_loose_e,newweight );
      hSip_loose_e->Fill( max_Sip_loose_e,newweight );
      hIp_loose_e->Fill( max_Ip_loose_e,newweight );
      
      // Electron Cross Cleaning  -- eles separated from muons (deltaR > 0.05)
      
      for(int e = 0; e < RECO_NELE; ++e)
      	for(int mu = 0; mu < RECO_NMU; ++mu){
	  	 
	  if( (RECOMU_isPFMu[mu] || (RECOMU_isTrackerHighPtMu[mu] && RECOMU_PT[mu] > 200.))
	      && (RECOMU_isGlobalMu[mu] || (RECOMU_isTrackerMu[mu] && RECOMU_numberOfMatches[mu]>0))
	      && RECOMU_mubesttrkType[mu]!=2
	      && RECOMU_PT[mu] > 5. 
	      && fabs(RECOMU_ETA[mu]) < 2.4 
	      && fabs(RECOMU_mubesttrkDxy[mu]) < .5 && fabs(RECOMU_mubesttrkDz[mu]) < 1. 
	      && fabs(RECOMU_SIP[mu])<4.  // TightID + SIP cut
	      );
	  else continue;
	  
	  double deltaR = sqrt( pow( DELTAPHI( RECOMU_PHI[mu] , RECOELE_PHI[e] ),2) + pow(RECOMU_ETA[mu] - RECOELE_ETA[e],2) );
	  
	  if( deltaR <= 0.05 ){
	    
	    if( debug )cout << "Electron not passing the cross cleaning" << endl;
	    
	    RECOELE_PT[e]  = -0.01;
	    RECOELE_ETA[e] = -99.;
	    RECOELE_PHI[e] = -99.;
	    RECOELE_SIP[e] = -99.;
	  }
	}
      
                  
      // Lepton identification -- no iso
      
      int iL[8]= { -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1};
      
      int N_good = 0 ;
      
      for( int i = 0; i < RECO_NMU; ++i ){
	
        if( debug ) cout << "\n Lepton i="<< i <<" properties: "
			 << "\nRECOMU_isPFMu[i] " << int(RECOMU_isPFMu[i])
			 << "\nRECOMU_isGlobalMu[i] " << int(RECOMU_isGlobalMu[i])
			 << "\nRECOMU_isTrackerMu[i] " << int(RECOMU_isTrackerMu[i])
			 << "\nRECOMU_PT[i] " << RECOMU_PT[i]
			 << "\nfabs(RECOMU_ETA[i]) " << fabs(RECOMU_ETA[i])
			 << "\nRECOMU_PFX_dB[i] " << RECOMU_PFX_dB[i]
			 << "\nfabs( RECOMU_SIP[i] ) " << fabs( RECOMU_SIP[i] )
			 << "\nfabs( RECOMU_mubesttrkDxy[i] ) " << fabs( RECOMU_mubesttrkDxy[i] )
			 << "\nfabs( RECOMU_mubesttrkDz[i] ) " << fabs( RECOMU_mubesttrkDz[i] )
			 << endl ;
	
       	// Tight muons
	if( (RECOMU_isPFMu[i] || (RECOMU_isTrackerHighPtMu[i] && RECOMU_PT[i] > 200.))
	    && ( RECOMU_isGlobalMu[i] || (RECOMU_isTrackerMu[i] && RECOMU_numberOfMatches[i]>0))
	    && RECOMU_mubesttrkType[i]!=2	 
	    && RECOMU_PT[i] > 5. 
	    && fabs(RECOMU_ETA[i]) < 2.4 
	    && fabs(RECOMU_mubesttrkDxy[i]) < .5 && fabs(RECOMU_mubesttrkDz[i]) < 1.
	    ){
	  
	  iL[ N_good ] = i ;
	  ++N_good ;	  
	}
      } // end loop on muons
      
      hN_good_mu->Fill( N_good,newweight );
      
      if( debug ) cout << "\nLeptons' indeces: "
		       << "\niL[0]: " << iL[0]
		       << "\niL[1]: " << iL[1]
		       << "\niL[2]: " << iL[2]
		       << "\niL[3]: " << iL[3]
		       << "\niL[4]: " << iL[4]
		       << "\niL[5]: " << iL[5]
		       << "\niL[6]: " << iL[6]
		       << "\niL[7]: " << iL[7]
		       << "\nNumber of good muons: " << N_good
		       << endl ;
      
      ++N_2 ;  // fill counter
      N_2_w=N_2_w+newweight;

   
      
      /// *** FSR
      // Photon identification & cleaning
      // ele identification is also needed
      
      //FSR photon identifications, will be used with MELA later
      int FSR_Z1_photid=-1;
      int FSR_Z2_photid=-1;
      int FSR_Z1_lepid=-1;
      int FSR_Z2_lepid=-1;

      //electrons:
      
      int Ne_good = 0 ;
      int iLe[8]= { -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1}; //electrons

      for( int i = 0; i < RECO_NELE; ++i ){

        if( debug ) cout << "\n Electron i="<< i <<" properties: "
			 << "\nRECOELE_PT[i] " << RECOELE_PT[i]
			 << "\nfabs(RECOELE_ETA[i]) " << fabs(RECOELE_ETA[i])
			 << "\nfabs(RECOELE_scl_Eta[i]) " << fabs(RECOELE_scl_Eta[i])
			 << "\nRECOELE_PFX_rho[i] " << RECOELE_PFX_rho[i]
			 << "\n chHad= " <<  RECOELE_PFchHad[i]                       
			 << "\n neuHad= " << RECOELE_PFneuHad[i]
			 << "\n photon= " <<  RECOELE_PFphoton[i]
                         << "\n PU= " << RECOELE_PFPUchAllPart[i]
			 << "\nfabs( RECOELE_SIP[i] ) " << fabs( RECOELE_SIP[i] )
			 << "\nRECOELE_mvaNonTrigV0[i] " << RECOELE_mvaNonTrigV0[i]
			 << "\nfabs( RECOELE_gsftrack_dxy[i] ) " << fabs( RECOELE_gsftrack_dxy[i] )
			 << "\nfabs( RECOELE_gsftrack_dz[i] ) " << fabs( RECOELE_gsftrack_dz[i] )
			 <<"\n charge = "<<RECOELE_CHARGE[i]
			 <<"\n Energy = "<<RECOELE_E[i]
			 << endl ;
       	
 	if( RECOELE_PT[i] > 7. && fabs(RECOELE_ETA[i]) < 2.5  ) ;
	  //&& RECOELE_gsftrack_expected_inner_hits[i]<=1 ) /* ok */ ;
	else continue ;
	cout<<"step0"<<endl;
	bool BDT_ok = 0; // Fall17
	// if( RECOELE_PT[i] > 7. &&  RECOELE_PT[i] <= 10. ){
	// 	if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > 0.5739521065342641 ) BDT_ok = 1 ;
	// 	if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) < 1.479 )
	// 					 && RECOELE_mvaNonTrigV0[i] > 0.5504628790992929 ) BDT_ok = 1 ;
	// 	if( fabs(RECOELE_scl_Eta[i]) >= 1.479 && RECOELE_mvaNonTrigV0[i] > 0.5924627534389098 ) BDT_ok = 1 ;
	// }
	// else { 
	// 	if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > -0.03391387993354392 ) BDT_ok = 1 ;
	// 	if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) <= 1.479 )
	// 					 && RECOELE_mvaNonTrigV0[i] > -0.018451958064666783 ) BDT_ok = 1 ;
	// 	if( fabs(RECOELE_scl_Eta[i]) > 1.479 && RECOELE_mvaNonTrigV0[i] > -0.38565459150737535 ) BDT_ok = 1 ;
	// }

	if( RECOELE_PT[i] > 7. && RECOELE_ID[i]==1 ){BDT_ok =1;};
	if( !BDT_ok ) continue ;
	cout<<"step1"<<endl;
	
	if( fabs(RECOELE_gsftrack_dxy[i]) < .5 
	 && fabs(RECOELE_gsftrack_dz[i])  < 1. ) /* ok */ ;
	else continue ; 
	cout<<"step2"<<endl;
	iLe[ Ne_good ] = i ;
	++Ne_good ;

      }// end loop on electrons


      hN_good_ele->Fill( Ne_good,newweight );
      hN_good_lep->Fill( N_good + Ne_good,newweight );


      if( debug ) cout << "\n Electrons' indeces: "
 		  << "\niLe[0]: " << iLe[0]
  		  << "\niLe[1]: " << iLe[1]
 		  << "\niLe[2]: " << iLe[2]
		  << "\niLe[3]: " << iLe[3]
 		  << "\niLe[4]: " << iLe[4]
  		  << "\niLe[5]: " << iLe[5]
 		  << "\niLe[6]: " << iLe[6]
		  << "\niLe[7]: " << iLe[7]
		  << "\nNumber of good eles: " << Ne_good
		  << endl ;


      // Define a new isolation array to allocate the contribution of photons
      //float RECOMU_PFX_dB_new[100],RECOELE_PFX_rho_new[100];
      for (int i=0;i<100;i++){
	RECOMU_PFX_dB_new[i]=RECOMU_PFX_dB[i];
	RECOELE_PFX_rho_new[i]=RECOELE_PFX_rho[i];
      }
      //
      
      // photon definition & cleaning:
      int iLp[30];
      for( int i = 0 ; i < 30 ; ++i )iLp[i] = -1;
      
      int Nphotons = 0;

      for( int i = 0; i < RECO_NPFPHOT; ++i ){
	
        if( debug ) cout << "\n Photon i="<< i <<" properties: "
			 << "\n RECOPFPHOT_PT[i] " << RECOPFPHOT_PT[i]
			 << "\n fabs(RECOPFPHOT_ETA[i]) " << fabs(RECOPFPHOT_ETA[i])
			 << "\n RECOPFPHOT_PHI[i] " << RECOPFPHOT_PHI[i]
			 << "\n RECOPFPHOT_PFX_rho[i] " << RECOPFPHOT_PFX_rho[i]
			 << endl ;
	
	if ( RECOPFPHOT_PT[i] > 2. && fabs(RECOPFPHOT_ETA[i]) < 2.4 && RECOPFPHOT_PFX_rho[i]<1.8 ) {
	  
	  bool is_clean = 1;
	  
	  // cleaning
	  for(int e = 0; e < N_loose_e; ++e){
  	    if (fabs( RECOELE_SIP[iL_loose_e[e]])>=4.) continue;  // loose ID + SIP cut
	    double deltaPhi = DELTAPHI( RECOPFPHOT_PHI[i] , RECOELE_scl_Phi[iL_loose_e[e]] ) ;
	    double deltaEta = fabs( RECOPFPHOT_ETA[i] - RECOELE_scl_Eta[iL_loose_e[e]] );
	    double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[i] , RECOELE_scl_Phi[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[i] - RECOELE_scl_Eta[iL_loose_e[e]],2) );
	    
	    if( ( fabs(deltaPhi) < 2 && fabs(deltaEta) < 0.05 ) || deltaR <= 0.15 ){		  
	      if( debug )cout << "Photon not passing the electron cleaning" << endl;	
	      is_clean = 0;	  
	      
	    }
	  } // end loop on eles		             	
	  
	  if( !is_clean ) continue ;
	  
	  
	  iLp[ Nphotons ] = i ;
	  ++Nphotons ;
	  
	}
      }// end loop on photons
      
      hN_good_phot->Fill( Nphotons,newweight );
      
      if( debug ) cout << "Photons' indeces: "
		       << "\niLp[0]: " << iLp[0]
		       << "\niLp[1]: " << iLp[1]
		       << "\niLp[2]: " << iLp[2]
		       << "\niLp[3]: " << iLp[3]
		       << "\niLp[4]: " << iLp[4]
		       << "\niLp[5]: " << iLp[5]
		       << "\niLp[6]: " << iLp[6]
		       << "\niLp[7]: " << iLp[7]
		       << "\nNumber of good photons: " << Nphotons
		       << endl ;
      
      
      // assign to each photon the closest lepton
      int iLp_l[30];
      for( int i = 0 ; i < 30 ; ++i )iLp_l[i] = -1;
      int iLp_tagEM[30];
      for( int i = 0 ; i < 30 ; ++i )iLp_tagEM[i] = -1;  // tag  0: mu  1: ele
      
      float RECOPFPHOT_DR[30];
      for( int i = 0 ; i < 30 ; ++i ) RECOPFPHOT_DR[i] = -999; 

      for( int i = 0; i < Nphotons; ++i ){
	
	double min_deltaR = 1000;
	int  l_min_deltaR = -1;
	int  tag_min_deltaR = -1;   // 0: mu  1: ele
	
	for(int l = 0; l < N_loose_mu; ++l){ // loop on muons
	  if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue;  //loose ID + SIP cut
	  double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOMU_PHI[iL_loose_mu[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOMU_ETA[iL_loose_mu[l]],2) );
	  if(!(deltaR < 0.5 && deltaR/pow(RECOPFPHOT_PT[iLp[i]],2)<0.012) ) continue;
	  if( deltaR<min_deltaR) { // the closest lepton
	    cout << "Possible candidate of photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " associated to a muon with pT= " << RECOMU_PT[iL_loose_mu[l]]<< endl;
	    min_deltaR = deltaR;
	    l_min_deltaR = l;
	    tag_min_deltaR = 0;
	  }
	  
	}//end loop on muons  
	
	for(int l = 0; l < N_loose_e; ++l){ // loop on electrons
	  if (fabs(RECOELE_SIP[iL_loose_e[l]])>=4.) continue;  //loose ID + SIP cut	  
	  double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_PHI[iL_loose_e[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_ETA[iL_loose_e[l]],2) );
	  //cout << "DeltaR= " << deltaR << " " << deltaR/pow(RECOPFPHOT_PT[iLp[i]],2) << endl;
	  if(!(deltaR < 0.5 && deltaR/pow(RECOPFPHOT_PT[iLp[i]],2)<0.012) ) continue;
	  if( deltaR<min_deltaR) { // the closest lepton
	    cout << "Possible candidate of photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " associated to an electron with pT= " << RECOELE_PT[iL_loose_e[l]]<< endl;
	    min_deltaR = deltaR;
	    l_min_deltaR = l;
	    tag_min_deltaR = 1;
	  }
	  
	}//end loop on electrons  
	
       
	if( min_deltaR < 0.5 ){
	    if (tag_min_deltaR==0) iLp_l[ i ] = iL_loose_mu[l_min_deltaR];
	    if (tag_min_deltaR==1) iLp_l[ i ] = iL_loose_e[l_min_deltaR];
	    iLp_tagEM[ i ] = tag_min_deltaR;
	    RECOPFPHOT_DR[iLp[i]]=min_deltaR; 	  
	}	
      }

      
      
      if( debug ) cout << "Indeces of Tight leptons associated to photons: "
		       << "\niLp_l[0]: " << iLp_l[0]
		       << "\niLp_l[1]: " << iLp_l[1]
		       << "\niLp_l[2]: " << iLp_l[2]
		       << "\niLp_l[3]: " << iLp_l[3]
		       << "\niLp_l[4]: " << iLp_l[4]
		       << "\niLp_l[5]: " << iLp_l[5]
		       << "\niLp_l[6]: " << iLp_l[6]
		       << "\niLp_l[7]: " << iLp_l[7]
		       << endl ;
      
      if( debug ) cout << "Tag of leptons associated to photons: (0: mu , 1:ele)"
 		  << "\niLp_tagEM[0]: " << iLp_tagEM[0]
  		  << "\niLp_tagEM[1]: " << iLp_tagEM[1]
 		  << "\niLp_tagEM[2]: " << iLp_tagEM[2]
		  << "\niLp_tagEM[3]: " << iLp_tagEM[3]
 		  << "\niLp_tagEM[4]: " << iLp_tagEM[4]
  		  << "\niLp_tagEM[5]: " << iLp_tagEM[5]
 		  << "\niLp_tagEM[6]: " << iLp_tagEM[6]
		  << "\niLp_tagEM[7]: " << iLp_tagEM[7]
		  << endl ;

      // for(int i=0.;i<Nphotons;i++) {
      // 	if (iLp_l[i]!=-1 && iLp_tagEM[i]==1) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " attached to an electron with pT= " << RECOELE_PT[iLp_l[i]] << endl;
      // };

      // Multiple photons associated to the same lepton: the lowest-ΔR(γ,l)/ETγ2 has to be selected.
      double min_deltaR_ET2=1000;
      int p_min_deltaR_ET2=-1;

      for(int l = 0; l < N_loose_mu; ++l){ // loop on muons
	if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue; //loose ID + SIP cut
	min_deltaR_ET2=1000;
	p_min_deltaR_ET2=-1;
	
	for( int p = 0; p < Nphotons; ++p ){
	  if( iLp_l[ p ] == iL_loose_mu[l] && iLp_tagEM[ p ] == 0 )  {
	    cout <<  "index muon" << iL_loose_mu[l] << endl;
	    double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[p]] , RECOMU_PHI[iL_loose_mu[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[p]] - RECOMU_ETA[iL_loose_mu[l]],2) );
	    double deltaR_ET2 = deltaR/pow(RECOPFPHOT_PT[iLp[p]],2);
	    if (deltaR_ET2<min_deltaR_ET2) {
	      min_deltaR_ET2=deltaR_ET2;
	      RECOPFPHOT_DR[iLp[p]]=deltaR;
	      p_min_deltaR_ET2=p;
	    }
	  }
	}
	
	if (p_min_deltaR_ET2!=-1){
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == iL_loose_mu[l] && iLp_tagEM[ p ] == 0 )  {
	      if (p!=p_min_deltaR_ET2){
		iLp_l[ p ] = -1;
		iLp_tagEM[ p ] = -1;
	      }
	    }
	  }
	}
	
      }

   
      
      //
      min_deltaR_ET2=1000;
      p_min_deltaR_ET2=-1;
      
      for(int l = 0; l < N_loose_e; ++l){ // loop on electrons
	if (fabs(RECOELE_SIP[iL_loose_e[l]])>=4.) continue; //loose ID + SIP cut
	min_deltaR_ET2=1000;
	p_min_deltaR_ET2=-1;
	
	for( int p = 0; p < Nphotons; ++p ){
	  if( iLp_l[ p ] == iL_loose_e[l] && iLp_tagEM[ p ] == 1 )  {
	    cout <<  "index electron" << iL_loose_e[l] << endl;
	    double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[p]] , RECOELE_PHI[iL_loose_e[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[p]] - RECOELE_ETA[iL_loose_e[l]],2));
	    double deltaR_ET2 = deltaR/pow(RECOPFPHOT_PT[iLp[p]],2);
	    cout << " deltaR_ET2= " << deltaR_ET2 <<endl;
	    if (deltaR_ET2<min_deltaR_ET2){
	      min_deltaR_ET2=deltaR_ET2;
	      RECOPFPHOT_DR[iLp[p]]=deltaR;
	      p_min_deltaR_ET2=p;
	      cout << " p_min_deltaR_ET2= " << p_min_deltaR_ET2 <<endl;
	    }
	  }	  
	}
	
	if (p_min_deltaR_ET2!=-1){
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == iL_loose_e[l] && iLp_tagEM[ p ] == 1 )  {
	      if (p!=p_min_deltaR_ET2){
		iLp_l[ p ] = -1;
		iLp_tagEM[ p ] = -1;
	      }
	    }
	  }	  
	}
	
      }	
     

      
      if( debug ) cout << "Indeces of loose leptons associated to the photon with lowest DeltaR/ET2: "
		       << "\niLp_l[0]: " << iLp_l[0]
		       << "\niLp_l[1]: " << iLp_l[1]
		       << "\niLp_l[2]: " << iLp_l[2]
		       << "\niLp_l[3]: " << iLp_l[3]
		       << "\niLp_l[4]: " << iLp_l[4]
		       << "\niLp_l[5]: " << iLp_l[5]
		       << "\niLp_l[6]: " << iLp_l[6]
		       << "\niLp_l[7]: " << iLp_l[7]
		       << endl ;
      
      if( debug ) cout << "Tag of leptons associated to the photon with lowest DetaR/ET2: (0: mu , 1:ele)"
		       << "\niLp_tagEM[0]: " << iLp_tagEM[0]
		       << "\niLp_tagEM[1]: " << iLp_tagEM[1]
		       << "\niLp_tagEM[2]: " << iLp_tagEM[2]
		       << "\niLp_tagEM[3]: " << iLp_tagEM[3]
		       << "\niLp_tagEM[4]: " << iLp_tagEM[4]
		       << "\niLp_tagEM[5]: " << iLp_tagEM[5]
		       << "\niLp_tagEM[6]: " << iLp_tagEM[6]
		       << "\niLp_tagEM[7]: " << iLp_tagEM[7]
		       << endl ;

      for(int i=0.;i<Nphotons;i++) {
       	if (iLp_l[i]!=-1 && iLp_tagEM[i]==1) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[i]] << " attached to an electron with pT= " << RECOELE_PT[iLp_l[i]] << endl;
      };           
      

      // Exclude that photon from the isolation cone all leptons in the event passing loose ID + SIP cut if it was in the isolation cone and outside the isolation veto (ΔR>0.01 for muons and (ele->supercluster()->eta() < 1.479 || dR > 0.08) for electrons

      double EffectiveArea=-9999.;
      cout << "Rho for electron pileup isolation correction is= " << RHO_ele << endl;
      
      for(int i=0.;i<Nphotons;i++) {
	if (iLp_l[i]==-1) continue;
	
	for(int e = 0; e < N_loose_e; ++e){
	  //if(!( iLp_l[i] == iL_loose_e[e] && iLp_tagEM[i] == 1 ) ) continue;
	  if (fabs( RECOELE_SIP[iL_loose_e[e]])>=4.) continue;
	  //double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_scl_Phi[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_scl_Eta[iL_loose_e[e]],2) );
	  double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOELE_PHI[iL_loose_e[e]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOELE_ETA[iL_loose_e[e]],2) );
	  cout << "deltaR for photon subtraction= " << deltaR << endl;
	  if( deltaR<=0.3 && (RECOELE_scl_Eta[iL_loose_e[e]]< 1.479 || deltaR>0.08) ){ // 0.4 is the isolation cone for electrons in 74x -> 0.3 in 76x              
	    if( debug )cout << "Subtracting the photon isolation from the electron isolation value " << endl;
	    
	    EffectiveArea=EAele(iL_loose_e[e],tag_2011);
	    RECOELE_PFX_rho_new[iL_loose_e[e]]=
              (RECOELE_PFchHad[iL_loose_e[e]]+
               max(0.,RECOELE_PFneuHad[iL_loose_e[e]]+
                   (RECOELE_PFphoton[iL_loose_e[e]]-RECOPFPHOT_PT[iLp[i]] )-
                   max(RHO_ele,0.0)*(EffectiveArea)))/RECOELE_PT[iL_loose_e[e]];	    
	  }
	} // end loop on ele
	
	for(int l = 0; l < N_loose_mu; ++l){ // loop on muons
	  //if(!( iLp_l[i] == iL_loose_mu[l] && iLp_tagEM[i] == 0 ) ) continue;
          if (fabs(RECOMU_SIP[iL_loose_mu[l]])>=4.) continue;
          double deltaR = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[i]] , RECOMU_PHI[iL_loose_mu[l]] ),2) + pow(RECOPFPHOT_ETA[iLp[i]] - RECOMU_ETA[iL_loose_mu[l]],2) );
	  cout << "deltaR for photon subtraction= " << deltaR << endl;
	  if( deltaR<=0.3 && deltaR>0.01){ // 0.3 is the isolation cone for muons in 76x
	    if( debug )cout << "Subtracting the photon isolation from the muon isolation value " << endl;
		    
	    RECOMU_PFX_dB_new[iL_loose_mu[l]]=
              (RECOMU_PFchHad[iL_loose_mu[l]]+
               max(0.,RECOMU_PFneuHad[iL_loose_mu[l]]+
                   (RECOMU_PFphoton[iL_loose_mu[l]]-RECOPFPHOT_PT[iLp[i]] )-
                   0.5*RECOMU_PFPUchAllPart[iL_loose_mu[l]]))/RECOMU_PT[iL_loose_mu[l]];
	    
	  }
	} // end loop on mu
	
	
      }	
      
     
      // *** end FSR
      

      // **** Step 3:
      struct candidateZ {
	float massvalue;
	int ilept1;
	float pt1;
	float isol1;
	bool ilept1_FSR;
	float eta1;
	float phi1;
	int charge1;
	int charge2;
	int ilept2;
	float pt2;
	float isol2;
	bool ilept2_FSR;
	float eta2;
	float phi2;
	float pxZ;
	float pyZ;
	float pzZ;
	float EZ;
	bool withFSR;
	float ptFSR;
      };

      vector<candidateZ> Zcandvector;
      Zcandvector.clear();
      vector<candidateZ> Zcandisolvector;
      Zcandisolvector.clear();

      // a) pair #1: mass closest to Z1
      // b) mLL in ] 40,120 [
      if( debug ) cout  << "\nStep 3: Number of good leptons: " << Ne_good << endl;
 
      if( Ne_good < 2 ) continue ; 	

      int i1 = -1; //index of the first lepton (from Z1)
      int j1 = -1; //index of the second lepton (from Z1)
      int pi1 = -1; 
      int pj1 = -1;
      
      bool has_FSR_Z1 = 0;
      TLorentzVector Lepton1,Lepton2,DiLepton,LeptonCorrection;

      for (int i=0; i<4; i++){
	cout<<"Electron properties PT= "<<RECOELE_PT[ iLe[i] ]
	    <<"\n SIP = "<<RECOELE_SIP[iLe[i]]
	    <<"\n ISO = "<<RECOELE_PFX_rho_new[iLe[i]]<<endl; }
      
      for(int i = 0; i < Ne_good; ++i){
        for(int j = i + 1; j < Ne_good; ++j){

	  if (fabs(RECOELE_SIP[iLe[i]])>=4.) continue;
	  if (fabs(RECOELE_SIP[iLe[j]])>=4.) continue;
	  // if (fabs(RECOELE_PFX_rho_new[iLe[i]])>=0.35) continue; // Isolation cut //reham
	  // if (fabs(RECOELE_PFX_rho_new[iLe[j]])>=0.35) continue; //reham
	  
	  if(RECOELE_CHARGE[ iLe[j] ] == RECOELE_CHARGE[ iLe[i] ]) continue; // opposite charge

	  cout << "Pairing electrons with pT= " << RECOELE_PT[ iLe[i] ] << " and " <<  RECOELE_PT[ iLe[j] ] << endl;
	  
	  // evaluate the mass &
	  double pxZ, pyZ, pzZ;
	  double EZ;
	  double massZ;
	  double massZ_noFSR = 0;
	  
	  int tempphotid=-1;
	  int templepid=-1;

	  float pTphot=-999.;
	  Lepton1.SetPtEtaPhiM(RECOELE_PT[iLe[i]], RECOELE_ETA[iLe[i]], RECOELE_PHI[iLe[i]], 0.000511);
	  Lepton2.SetPtEtaPhiM(RECOELE_PT[iLe[j]], RECOELE_ETA[iLe[j]], RECOELE_PHI[iLe[j]], 0.000511);
	  DiLepton=Lepton1+Lepton2;	  
	  massZ = DiLepton.M();	  
	  massZ_noFSR = massZ;
	  if (debug) cout << "Mass Z= " << massZ << endl;
	  pxZ=DiLepton.Px();
	  pyZ=DiLepton.Py();
	  pzZ=DiLepton.Pz();
	  EZ=DiLepton.E();
	  	  	 
	  // ** Association of FSR to Z
	  if( debug ) cout  << "Step Z+FSR  " << endl;
	  
	  bool has_FSR_Z = 0;
	  int N_FSR_Z = 0;
	  double max_pt_FSR_Z = -1.;
	  int pi = -1; 
	  int pj = -1;
	  
	  
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == iLe[i] && iLp_tagEM[ p ] == 1 )  {  // exists a photon associated to a lepton electron
	      
	      // evaluate the mass
	      LeptonCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]],RECOPFPHOT_PHI[iLp[p]],0);
	      Lepton1=Lepton1+LeptonCorrection;
	      DiLepton=Lepton1+Lepton2;
	      double mllp=DiLepton.M();
	      pxZ=DiLepton.Px();
	      pyZ=DiLepton.Py();
	      pzZ=DiLepton.Pz();
	      EZ=DiLepton.E();	      	      	   
	      
	      has_FSR_Z = 1; 
	      pi = p; 
	      ++N_FSR_Z;
	      if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	      massZ=mllp;

	      cout << "Mass Z with FSR= "<< massZ << endl;
	      
	    }
	    
	    if( iLp_l[ p ] == iLe[j] && iLp_tagEM[ p ] == 1 )  { 
	      
	      // evaluate the mass
	      LeptonCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]],RECOPFPHOT_PHI[iLp[p]],0);
	      Lepton2=Lepton2+LeptonCorrection;
	      DiLepton=Lepton1+Lepton2;
	      double mllp=DiLepton.M();
	      pxZ=DiLepton.Px();
	      pyZ=DiLepton.Py();
	      pzZ=DiLepton.Pz();
	      EZ=DiLepton.E();
	      	    
	      pj = p;
	      has_FSR_Z = 1;
	      ++N_FSR_Z; 
	      if( RECOPFPHOT_PT[iLp[p]] > max_pt_FSR_Z ) max_pt_FSR_Z = RECOPFPHOT_PT[iLp[p]];
	      massZ=mllp;

	      cout << "Mass Z with FSR= "<< massZ << endl;
	      
	    }
	  } // end loop on FSR photons

	 
	  
	  //if( has_FSR_Z ) debug = 1;
	  
	  if( debug && has_FSR_Z) {
	    cout  << " Z has FSR! " << endl;
	    cout  << "  N_FSR_Z " << N_FSR_Z << endl;
	    cout  << "  max_pt of photon FSR_Z " << max_pt_FSR_Z << endl;
	    if( pi > -1 ) cout  << "  pi " << pi << " --> index photon: " << iLp[pi] << " associated lepton: " << iLp_l[pi] << " (= "<< iLe[i]<<" ? )  tag: " << iLp_tagEM[pi] << endl;
	    if( pj > -1 ) cout  << "  pj " << pj << " --> index photon: " << iLp[pj] << " associated lepton: " << iLp_l[pj] << " (= "<< iLe[j]<<" ? )  tag: " << iLp_tagEM[pj] << endl;
	  }
	  else {
	    cout << "No FSR photon attached" << endl;
	  }
	  
	  
	  if( has_FSR_Z ){ // if Z has FSR
	    
	    ++N_3_FSR; // fill the counter
	    N_3_FSR_w=N_3_FSR_w+newweight;
	    	    
	    // do not recompute isolation here
	    if( debug ) cout  << "Z Isolation (not corrected for photon): "
                              << "\n RECOELE_PFX_rho[ iLe[i] ] " << RECOELE_PFX_rho[ iLe[i] ]
                              << "\n RECOELE_PFX_rho[ iLe[j] ] " << RECOELE_PFX_rho[ iLe[j] ]
                              << endl;
	    	   
	    if( pi != -1 ){
	      //double deltaR_i = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[pi]] , RECOELE_PHI[iLe[i]] ),2) + pow(RECOPFPHOT_ETA[iLp[pi]] - RECOELE_ETA[iLe[i]],2) );
              //if( deltaR_i < 0.4 && deltaR_i > 0.01 )
	      pTphot=RECOPFPHOT_PT[iLp[pi]];	      
	    }
	    else if( pj != -1 ){	      
              //double deltaR_j = sqrt( pow( DELTAPHI( RECOPFPHOT_PHI[iLp[pj]] , RECOELE_PHI[iLe[j]] ),2) + pow(RECOPFPHOT_ETA[iLp[pj]] - RECOELE_ETA[iLe[j]],2) );
              //if( deltaR_j < 0.4 && deltaR_j > 0.01 )
	      pTphot=RECOPFPHOT_PT[iLp[pj]];	      
	    }

 
	  } // end if has FSR
	  else{
	    
	    if( debug ) cout  << "Z Isolation: "  
			      << "\n RECOELE_PFX_rho_new[ iLe[i] ] " << RECOELE_PFX_rho_new[ iLe[i] ]
			      << "\n RECOELE_PFX_rho_new[ iLe[j] ] " << RECOELE_PFX_rho_new[ iLe[j] ]
			      << endl;
	    
	  }
	  // ** end association of FSR to Z
	  
	  //if( massZ == 0 || i1 == -1 || j1 == -1) continue;
	  
	  cout << "Filling a struct for Z" << endl; 
	  candidateZ *Z = new candidateZ;
	  Z->massvalue=massZ;
	  Z->ilept1=iLe[i];
	  Z->ilept2=iLe[j];
	  Z->pt1=RECOELE_PT[iLe[i]];
	  Z->pt2=RECOELE_PT[iLe[j]];
	  Z->eta1=RECOELE_ETA[iLe[i]];
	  Z->eta2=RECOELE_ETA[iLe[j]];
	  Z->phi1=RECOELE_PHI[iLe[i]];
	  Z->phi2=RECOELE_PHI[iLe[j]];
	  Z->charge1=RECOELE_CHARGE[iLe[i]];
	  Z->charge2=RECOELE_CHARGE[iLe[j]];
	  Z->isol1=RECOELE_PFX_rho_new[ iLe[i] ];
	  Z->isol2=RECOELE_PFX_rho_new[ iLe[j] ];
	  if( pi != -1 ) Z->ilept1_FSR=true;
	  if( pj != -1 ) Z->ilept2_FSR=true;
	  Z->pxZ=pxZ;
	  Z->pyZ=pyZ;
	  Z->pzZ=pzZ;
	  Z->EZ=EZ;
	  if( has_FSR_Z ) {
	    Z->withFSR=1;
	    Z->ptFSR=pTphot;
	  }	      
	  else {
	    Z->withFSR=0;
	    Z->ptFSR=0.;
	  }
	  	 	  
	  Zcandvector.push_back(*Z);	  
	  
	}
      } // end loop on couples
      
     
      if (Zcandvector.size()<2) {
	cout << "Less than two Z pairs with isolated leptons...exiting" << endl;
	continue; 
      }

      ++N_3a ;  // fill counter
      N_3a_w=N_3a_w+newweight;
      

      // Mass cut on Z
      vector<candidateZ> Zcandisolmassvector;
      Zcandisolmassvector.clear();

      for (int index=0; index<Zcandvector.size();index++){
	if (!(Zcandvector.at(index).massvalue > 12. && Zcandvector.at(index).massvalue < 120.)) continue;
	cout << "Z passing the 12 < mll < 120 cut with mass= " << Zcandvector.at(index).massvalue<< endl;
	Zcandisolmassvector.push_back(Zcandvector.at(index));
      };
      
      if (Zcandisolmassvector.size()<2) {
	cout << "No ZZ passing the mass cut"<< endl;
	continue;
      }

      cout << "Number of Z passing the isolation and the 12 << mll < 120 cut is= " << Zcandisolmassvector.size() << endl;
      
      ++N_3b ;  // fill counter
      N_3b_w=N_3b_w+newweight;
      
      Double_t eff_weight_3;

      // MZ and PFMET
      for (int index=0; index<Zcandisolmassvector.size();index++){
	int ipt1_3;
        int ipt2_3;
        ipt1_3 = Zcandisolmassvector.at(index).ilept1; 
        ipt2_3  = Zcandisolmassvector.at(index).ilept2; 
	
        // Execute Efficiency Reweighting
        int z1lept[2]={ipt1_3,ipt2_3};
	
        for(int i = 0; i < 2; ++i){
          Double_t Pt = RECOELE_PT[ z1lept[i] ]; 
          Double_t Eta = RECOELE_ETA[ z1lept[i] ]; 
	  Double_t SC_Eta = RECOELE_scl_Eta[ z1lept[i] ]; //Reham

	  Double_t Tmp_Pt = Pt;
	  if (Pt>500)Tmp_Pt = 500; //overflow
          
	  // if( (MC_type == "Spring16" || MC_type=="Moriond17") && DATA_type == "NO"){
	  if(  MC_type == "Fall17"  && DATA_type == "NO"){  
            if(RECOELE_isGap[ z1lept[i] ]==0){
	      // Int_t binx = ele_scale_2017->GetXaxis()->FindBin(Eta);
	      Int_t binx = ele_scale_2017->GetXaxis()->FindBin(SC_Eta);
              Int_t biny = ele_scale_2017->GetYaxis()->FindBin(Tmp_Pt);
              if (ele_scale_2017->GetBinContent(binx,biny)>0.) eff_weight_3*=ele_scale_2017->GetBinContent(binx,biny); 
            }
            else if(RECOELE_isGap[ z1lept[i] ]==1){
              //Int_t binx = ele_scale_2017_gap->GetXaxis()->FindBin(Eta);
	      Int_t binx = ele_scale_2017_gap->GetXaxis()->FindBin(SC_Eta);
              Int_t biny = ele_scale_2017_gap->GetYaxis()->FindBin(Tmp_Pt);
              if (ele_scale_2017_gap->GetBinContent(binx,biny)>0.) eff_weight_3*=ele_scale_2017_gap->GetBinContent(binx,biny); 
            }
	    
          }
        }//end scale

	//add Reco efficiency

        for(int i = 0; i < 2; ++i){
          Double_t Pt = RECOELE_PT[ z1lept[i] ]; 
          Double_t Eta = RECOELE_ETA[ z1lept[i] ]; 
	  Double_t SC_Eta = RECOELE_scl_Eta[ z1lept[i] ]; //Reham
          
	  // if( (MC_type == "Spring16" || MC_type=="Moriond17") && DATA_type == "NO"){
	  if(  MC_type == "Fall17"  && DATA_type == "NO"){  
            if(Pt <=20.){
	      // Int_t binx = ele_Reco_eff_2017_lowEt->GetXaxis()->FindBin(Eta);
	      Int_t binx = ele_Reco_eff_2017_lowEt->GetXaxis()->FindBin(SC_Eta);
              Int_t biny = ele_Reco_eff_2017_lowEt->GetYaxis()->FindBin(Pt);
              if (ele_Reco_eff_2017_lowEt->GetBinContent(binx,biny)>0.) eff_weight_3*=ele_Reco_eff_2017_lowEt->GetBinContent(binx,biny); 
            }
            else if(Pt>20.){
	      Double_t Tmp_Pt = Pt;
	      if (Pt>500.)Tmp_Pt = 500.;
              //Int_t binx = ele_Reco_eff_2017_highEt->GetXaxis()->FindBin(Eta);
	      Int_t binx = ele_Reco_eff_2017_highEt->GetXaxis()->FindBin(SC_Eta);
              Int_t biny = ele_Reco_eff_2017_highEt->GetYaxis()->FindBin(Tmp_Pt);
              if (ele_Reco_eff_2017_highEt->GetBinContent(binx,biny)>0.) eff_weight_3*=ele_Reco_eff_2017_highEt->GetBinContent(binx,biny); 
            }
	    
          }
        }//end Reco scale
	
        // Changing the weight for pileup and efficiency
        if (DATA_type == "2017") eff_weight_3=1.; 
        if (eff_weight_3>0.) newweight=weight*pT_weight*pu_weight*weight_kfactor*eff_weight_3;
	else newweight=weight*pT_weight*pu_weight*weight_kfactor;
	
	hMZ_3->Fill(Zcandisolmassvector.at(index).massvalue,newweight );
      }

      hPFMET_3->Fill(RECO_PFMET,newweight);
      hLogXPFMET_3->Fill(RECO_PFMET,newweight); 
      hLogLinXPFMET_3->Fill(RECO_PFMET,newweight);

      cout << "Starting weight + pileup + efficiency= " << newweight << endl;
      if(debug) cout << "Efficiency Weight for the Z1: " << eff_weight_3 << " Final weight for Z1= " << newweight << endl;
     
      /*
      //z1tree->Fill();
      
   
      int ipt1_3;
      int ipt2_3;
      if( RECOELE_PT[ iLe[i1] ] > RECOELE_PT[ iLe[j1] ] ){ ipt1_3 = iLe[i1]; ipt2_3  = iLe[j1]; }
      else { ipt1_3 = iLe[j1]; ipt2_3  = iLe[i1]; }
      
      hPtLep1_3->Fill( RECOELE_PT[ ipt1_3 ],newweight );
      hYLep1_3->Fill( RECOELE_ETA[ ipt1_3 ],newweight );
      hIsoLep1_3->Fill( RECOELE_PFX_rho[ ipt1_3 ],newweight );
      hSipLep1_3->Fill( RECOELE_SIP[ ipt1_3 ],newweight );
      hIpLep1_3->Fill( RECOELE_IP[ ipt1_3 ],newweight );
      hIpErLep1_3->Fill( RECOELE_IPERROR[ ipt1_3 ],newweight );
      
      
      hPtLep2_3->Fill( RECOELE_PT[ ipt2_3 ],newweight );
      hYLep2_3->Fill( RECOELE_ETA[ ipt2_3 ],newweight );
      hIsoLep2_3->Fill( RECOELE_PFX_rho[ ipt2_3 ],newweight );
      hSipLep2_3->Fill( RECOELE_SIP[ ipt2_3 ],newweight );
      hIpLep2_3->Fill( RECOELE_IP[ ipt2_3 ],newweight );
      hIpErLep2_3->Fill( RECOELE_IPERROR[ ipt2_3 ],newweight );

      
      float max_Iso_3 = -1 ;
      float max_Sip_3 = -1 ;
      float max_Ip_3 = -1 ;

      for (int i=0;i<RECO_NELE;i++) {
	if (i==iLe[i1] || i==iLe[j1]) continue;
	 if( RECOELE_PT[i] > 7. 
            && fabs(RECOELE_ETA[i]) < 2.5 
//            && RECOELE_gsftrack_expected_inner_hits[i]<=1 
            && fabs(RECOELE_gsftrack_dxy[i]) < .5 
            && fabs(RECOELE_gsftrack_dz[i]) < 1. 
            ) {   
	  if( RECOELE_PFX_rho[i] > max_Iso_3 ) max_Iso_3 = RECOELE_PFX_rho[i] ;
	  if( fabs( RECOELE_SIP[i] ) > max_Sip_3 ) max_Sip_3 = fabs( RECOELE_SIP[i] ) ;
	  if( fabs( RECOELE_IP[i] ) > max_Ip_3 ) max_Ip_3 = fabs( RECOELE_IP[i] ) ;
	}
      }

      hIso_3->Fill(max_Iso_3,newweight);
      hSip_3->Fill(max_Sip_3,newweight);
      hIp_3->Fill(max_Ip_3,newweight);


      // Filling plots for VBF analysis
      int njets_pass_3=0;
      float deltajj_3=-999.,massjj_3=-999.;
      float VD_3=-999.; //fisher discriminant
      int jet1_3,jet2_3;      
      int jetfail_3[100];
      double deltaR_3=0.;

      for(int i=0;i<100;i++) jetfail_3[i]=0;

      
      for(int i=0;i<RECO_PFJET_N;i++){
	//cout<<i<<" Jet with pt= "<<RECO_PFJET_PT[i]<<" ETA "<<RECO_PFJET_ETA[i]<<" PUID "<<RECO_PFJET_PUID[i]<<endl;
	if(RECO_PFJET_PUID[i]==1 && RECO_PFJET_PT[i]>30. && fabs(RECO_PFJET_ETA[i])<4.7 ){
	  //Check that jet has deltaR>0.4 away from any tight lepton 
	  
	   for(int mu = 0; mu < N_good; ++mu){
            if (RECOMU_SIP[iL[mu]]>=4.) continue;
            if (RECOMU_PFX_dB[iL[mu]]>=0.35) continue;
            deltaR_3 = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOMU_PHI[iL[mu]]),2) + pow(RECO_PFJET_ETA[i] - RECOMU_ETA[iL[mu]],2));
            cout << "1st lepton: " <<" deltaR "<< deltaR_3 <<endl;
            if (deltaR_3<0.4){
              jetfail_3[i]=1;
              cout << " jetfail " << jetfail_3[i] <<endl;
            }
          }

	  for(int ele = 0; ele < Ne_good; ++ele){
	    if (RECOELE_SIP[iLe[ele]]>=4.) continue;
	    if (RECOELE_PFX_rho_new[iLe[[ele]]>=0.5) continue;
	    deltaR_3 = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOELE_PHI[iLe[ele]]),2) + pow(RECO_PFJET_ETA[i] - RECOELE_ETA[iLe[ele]],2));
	    cout << "1st lepton: " <<" deltaR "<< deltaR_3 <<endl;
	    if (deltaR_3<0.4){
	      jetfail_3[i]=1;
	      cout << " jetfail " << jetfail_3[i] <<endl;
	    }
	  }

	  
	  
	  //Check that jet has deltaR>0.4 away from any FSR photon
	  //if (FSR_Z1_lepid != -1){
	  //  //cout<<"photon phi "<<RECOPFPHOT_PHI[FSR_Z1_photid]<<" photon eta "<<RECOPFPHOT_ETA[FSR_Z1_photid]<<endl;
	  //  deltaR_3 = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOPFPHOT_PHI[FSR_Z1_photid]),2) + pow(RECO_PFJET_ETA[i] - RECOPFPHOT_ETA[FSR_Z1_photid],2));
	  //  //cout<<FSR_Z1_photid<<" deltaR "<<deltaR_3<<endl;
	  //  if (deltaR_3<0.4){
	  //    jetfail_3[i]=1;
	  //  }
	  //}
	    
	  if (jetfail_3[i]==0){
	    cout<< " PASS jet " <<i<<" PT= "<<RECO_PFJET_PT[i]<<" ETA= "<<RECO_PFJET_ETA[i]<<" PUID= "<<RECO_PFJET_PUID[i]<<endl;
	    njets_pass_3++;
	  }
	}
	else{
	  jetfail_3[i]=1;
	}
	//cout<<" JETFAIL "<<jetfail_3[i]<<endl;
      }

      // VBF-tagged category - category 2
     int njets_vbf_3=0;
     float deltajj_3=-999.,massjj_3=-999.;
     float D_jet_3=-999.; //fisher discriminant
     TLorentzVector JET1_VBF_3, JET2_VBF_3,mJJ_3;
     int index_VBFjets_3[2]={-999,-999};
       
     for (int i=0;i<RECO_PFJET_N;i++){
       if (jetfail_3[i]!=0 ) continue; 
       for (int j=i+1;j<RECO_PFJET_N;j++){
         if (jetfail_3[j]!=0 ) continue;
         JET1_VBF_3.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
         JET2_VBF_3.SetPtEtaPhiE(RECO_PFJET_PT[j],RECO_PFJET_ETA[j],RECO_PFJET_PHI[j],RECO_PFJET_ET[j]*TMath::CosH(RECO_PFJET_ETA[j]));
         mJJ_3=JET1_VBF_3+JET2_VBF_3;
         //cout<<mJJ.M()<<endl;              
         
         deltajj_3=fabs(JET1_VBF_3.Eta()-JET2_VBF_3.Eta());
         massjj_3=mJJ_3.M();
         D_jet_£=0.18*fabs(JET1_VBF_3.Eta()-JET2_VBF_3.Eta())+1.92E-4*mJJ_3.M();
         if (debug) cout << "Mass jet jet= " << massjj_3 << " Delta jet jet= " << deltajj_3 << " Fisher D_jet" << D_jet_3 << endl;
         if (D_jet_3>0.5) {
           njets_vbf_3++;
           if (njets_vbf==1) {
             index_VBFjets_3[0]=i;
             index_VBFjets_3[1]=j;    
             JET1_VBF_3.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
             JET2_VBF_3.SetPtEtaPhiE(RECO_PFJET_PT[j],RECO_PFJET_ETA[j],RECO_PFJET_PHI[j],RECO_PFJET_ET[j]*TMath::CosH(RECO_PFJET_ETA[j]));
           }
         }
	
	hDjj_3->Fill( deltajj_3,newweight );
	hMjj_3->Fill( massjj_3,newweight );
	hVD_3->Fill( D_jet_3,newweight );
	
      }
      */

      
      // **** Step 4:
       // a) 4 leptons
      // b) pair #2
      // c) highest pt
      // d) mZ2 in ] 4,120 [

      int issamesign = 0;

      //if( debug ) cout  << "\nStep 4: Number of good leptons: " << N_good << endl;

      int N_Z2_pairs = 0;

      int i2 = -1; //index of the first lepton (from Z1)
      int j2 = -1; //index of the second lepton (from Z1)
      int pi2 = -1; 
      int pj2 = -1; 
      
      /*
      double pxZ2 = 0;  //Z2 kinematics
      double pyZ2 = 0;
      double pzZ2 = 0;
      double ptZ2 = 0;
      double EZ2 = 0;
      double Y_Z2 = -9;
      double massZ2 = 0;
      double massZ2_noFSR = 0;
      double sum_ptZ2 = 0.;
      */
      bool has_FSR_Z2 = 0;
      
      // ZZ object pairs
      array<int,2> isolgoodZs;
      vector<std::array<int, 2> > isolgoodZsv;
        
      
      // Ghost removal delta R > 0.02 
      vector<candidateZ> goodZ;

      for (int i=0;i<Zcandisolmassvector.size();i++){
	cout << "Checking mass i= " << Zcandisolmassvector.at(i).massvalue << endl;
	cout << "Ghost removal check 0: deltaR= " << sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi1, Zcandisolmassvector.at(i).phi2 ),2) 
		  + pow(Zcandisolmassvector.at(i).eta1-Zcandisolmassvector.at(i).eta2,2) ) << endl;
		if( sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi1, Zcandisolmassvector.at(i).phi2 ),2) 
		  + pow(Zcandisolmassvector.at(i).eta1-Zcandisolmassvector.at(i).eta2,2) ) <= 0.02) continue;
	
	
	for (int j=i+1;j<Zcandisolmassvector.size();j++){
	  cout << "Checking mass j= " << Zcandisolmassvector.at(j).massvalue << endl;
	  cout <<Zcandisolmassvector.at(i).pt1 << " " << Zcandisolmassvector.at(j).pt1 << " " << Zcandisolmassvector.at(j).pt2 << endl;
	  cout <<Zcandisolmassvector.at(i).pt2 << " " << Zcandisolmassvector.at(j).pt1 << " " << Zcandisolmassvector.at(j).pt2 << endl;
	  if (Zcandisolmassvector.at(i).pt1==Zcandisolmassvector.at(j).pt1 || Zcandisolmassvector.at(i).pt1==Zcandisolmassvector.at(j).pt2) continue;
	  if (Zcandisolmassvector.at(i).pt2==Zcandisolmassvector.at(j).pt1 || Zcandisolmassvector.at(i).pt2==Zcandisolmassvector.at(j).pt2) continue;
	 	  
	  cout << "Ghost removal check 1: deltaR= " << sqrt( pow( DELTAPHI(Zcandisolmassvector.at(j).phi1, Zcandisolmassvector.at(j).phi2 ),2) 
		    + pow(Zcandisolmassvector.at(j).eta1-Zcandisolmassvector.at(j).eta2,2) )<< endl;
	  if( sqrt( pow( DELTAPHI(Zcandisolmassvector.at(j).phi1, Zcandisolmassvector.at(j).phi2 ),2) 
		    + pow(Zcandisolmassvector.at(j).eta1-Zcandisolmassvector.at(j).eta2,2) ) <= 0.02) continue;
	  
	  cout << "Ghost removal check 2: deltaR= " << sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi1, Zcandisolmassvector.at(j).phi1 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta1-Zcandisolmassvector.at(j).eta1,2) )<< endl;
	  if( sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi1, Zcandisolmassvector.at(j).phi1 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta1-Zcandisolmassvector.at(j).eta1,2) ) <= 0.02) continue;

	  cout << "Ghost removal check 3: deltaR= " << sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi1, Zcandisolmassvector.at(j).phi2 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta1-Zcandisolmassvector.at(j).eta2,2) ) << endl;
	  if( sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi1, Zcandisolmassvector.at(j).phi2 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta1-Zcandisolmassvector.at(j).eta2,2) ) <= 0.02) continue;

	  cout << "Ghost removal check 4: deltaR= " << sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi2, Zcandisolmassvector.at(j).phi1 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta2-Zcandisolmassvector.at(j).eta1,2) ) << endl;
	  if( sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi2, Zcandisolmassvector.at(j).phi1 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta2-Zcandisolmassvector.at(j).eta1,2) ) <= 0.02) continue;

	  cout << "Ghost removal check 5: deltaR= " << sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi2, Zcandisolmassvector.at(j).phi2 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta2-Zcandisolmassvector.at(j).eta2,2) ) << endl;
	  if( sqrt( pow( DELTAPHI(Zcandisolmassvector.at(i).phi2, Zcandisolmassvector.at(j).phi2 ),2) 
		    + pow(Zcandisolmassvector.at(i).eta2-Zcandisolmassvector.at(j).eta2,2) ) <= 0.02) continue;

	  cout << "There is a set of 4 leptons passing the ghost removal (deltaR > 0.02)" << endl;
	  goodZ.push_back(Zcandisolmassvector.at(i));
	  goodZ.push_back(Zcandisolmassvector.at(j));
	  cout << "Filling isolgoodZsv vector " << endl;
	  isolgoodZs={i,j};
	  isolgoodZsv.push_back(isolgoodZs);
	}
      }

      //cout << "Debug ZZ= " << Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(0) ).massvalue << " " << Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(1)).massvalue << endl;
      if (isolgoodZsv.size()==0) {
	cout << "No ZZ combination passing the cuts  ...exiting " << endl;
	continue;
      }

      ++N_4a ;  // fill counter
      N_4a_w=N_4a_w+newweight;
 
      // PT,20/10 for any di-lepton
      vector<candidateZ> firstpTcleanedgoodZ;    
      vector<float> leptonspTcleaned;
      
      array<int,2> ileptonspTcleanedgoodZs;
      vector<std::array<int, 2> > ileptonspTcleanedgoodZsv;
      
      // PT,20/10 for any di-lepton      
      for (int l=0;l<isolgoodZsv.size();l++){
	leptonspTcleaned.clear();	
	leptonspTcleaned.push_back(Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(0) ).pt1);
	leptonspTcleaned.push_back(Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(0) ).pt2);
	leptonspTcleaned.push_back(Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(1) ).pt1);
	leptonspTcleaned.push_back(Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(1) ).pt2);
	std::sort(leptonspTcleaned.rbegin(),leptonspTcleaned.rend());
	
	if (leptonspTcleaned.at(0)>20. && leptonspTcleaned.at(1)>10.) {
	  firstpTcleanedgoodZ.push_back(Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(0) ));
	  firstpTcleanedgoodZ.push_back(Zcandisolmassvector.at( (isolgoodZsv.at(0)).at(1) ));
	  ileptonspTcleanedgoodZs={ (isolgoodZsv.at(l)).at(0), (isolgoodZsv.at(l)).at(1)};
	  ileptonspTcleanedgoodZsv.push_back(ileptonspTcleanedgoodZs);
	}	 	
      }
      
      cout << "Cleaned Good ZZ passing pT cuts are " << ileptonspTcleanedgoodZsv.size() << endl;            
      if (ileptonspTcleanedgoodZsv.size()==0) continue;


      // **** Step 6:
      // QCD suppression: mll>4 GeV cut on all OS-SF pairs (4/4)
      
      double min_mass_2L = 10000.;
      TLorentzVector Lepton1qcd,Lepton2qcd,DiLeptonQCD;
      vector<int> ileptonsee;
      
      array<int,2> iqcdcleanedgoodZs;
      vector<std::array<int, 2> > iqcdcleanedgoodZsv;
      vector<candidateZ> pTcleanedgoodZ; 
      
      for (int l=0;l<ileptonspTcleanedgoodZsv.size();l++){
	cout << "checking masses for QCD suppression= " 
	     << Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(0)).massvalue << " and " 
	     << Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(1)).massvalue
	     << endl;
	min_mass_2L = 10000.;
	ileptonsee.clear();
       
	ileptonsee.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(0)).ilept1);
	ileptonsee.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(0)).ilept2);
	ileptonsee.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(1)).ilept1);
	ileptonsee.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(1)).ilept2);
	
	double mass;
	
	// ee
	cout << "Checking ee pairs for QCD rejection" << endl;

	for (int i=0;i<4;i++){ 
	  for (int j=i+1;j<4;j++){
	    
	    if ( RECOELE_CHARGE[ileptonsee.at(i)] == RECOELE_CHARGE[ileptonsee.at(j)] ) continue; // shoud be OS	  
	    // evaluate the mass
	    double mass;
	    
	    Lepton1qcd.SetPtEtaPhiM(RECOELE_PT[ileptonsee.at(i)], RECOELE_ETA[ileptonsee.at(i)], RECOELE_PHI[ileptonsee.at(i)], 0.000511);
	    Lepton2qcd.SetPtEtaPhiM(RECOELE_PT[ileptonsee.at(j)], RECOELE_ETA[ileptonsee.at(j)], RECOELE_PHI[ileptonsee.at(j)], 0.000511);

	    DiLeptonQCD=Lepton1qcd+Lepton2qcd;       
	    mass = DiLeptonQCD.M();
	    cout << "Mass of ee for QCD rejection= " << mass << endl;
	    if( mass < min_mass_2L ) min_mass_2L = mass ;
	  }
	}

	hminMll_6->Fill( min_mass_2L,newweight );
	if( min_mass_2L <= 4 ) { 
	  cout << "Not passing the mll>4 cut" << endl;
	  continue ;
	}
	else {
	  pTcleanedgoodZ.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(0)));
	  pTcleanedgoodZ.push_back(Zcandisolmassvector.at( (ileptonspTcleanedgoodZsv.at(l)).at(1)));
	  
	  iqcdcleanedgoodZs={ileptonspTcleanedgoodZsv.at(l).at(0),ileptonspTcleanedgoodZsv.at(l).at(1)};
	  iqcdcleanedgoodZsv.push_back(iqcdcleanedgoodZs);
	}
	
      }
            
     
      
      // Z1 selection
      double pxZ1 = 0;  //Z1 kinematics
      double pyZ1 = 0;
      double pzZ1 = 0;
      double ptZ1 = 0;
      double EZ1 = 0;
      double Y_Z1 = -9;
      double massZ1 = 0;
      double massZ1_noFSR = 0;
      double sum_ptZ1 = 0.;
      int indexlep1Z1 = -1;
      int indexlep2Z1 = -1;
      int indexZ1= -1;
      
      // Choice of Z1 as the closest to the Z mass
      for (int i=0;i<pTcleanedgoodZ.size();++i){
	
	if( fabs(pTcleanedgoodZ.at(i).massvalue - Zmass) < fabs(massZ1 - Zmass) ){
	  
	  massZ1 = pTcleanedgoodZ.at(i).massvalue;
	  indexZ1=i;
	  
	  pxZ1 = pTcleanedgoodZ.at(i).pxZ;
	  pyZ1 = pTcleanedgoodZ.at(i).pyZ;
	  pzZ1 = pTcleanedgoodZ.at(i).pzZ;
	  EZ1  = pTcleanedgoodZ.at(i).EZ;
	  
	  ptZ1 = sqrt( pxZ1*pxZ1 + pyZ1*pyZ1 );
	  sum_ptZ1 = pTcleanedgoodZ.at(i).pt1+pTcleanedgoodZ.at(i).pt2;
	  
	  Y_Z1 = 0.5 * log ( (EZ1 + pzZ1)/(EZ1 - pzZ1) );
	  indexlep1Z1=pTcleanedgoodZ.at(i).ilept1;
	  indexlep2Z1=pTcleanedgoodZ.at(i).ilept2;
	}
      } 
      
      if (massZ1 < 40.) {
	cout << "The mass of Z1 is < 40 GeV...exiting" << endl;
	continue;
      } 
      
      if( debug ) cout  << "\n Final Z1 properties: "
			<< "\n pxZ1 " << pxZ1
			<< "\n pyZ1 " << pyZ1
			<< "\n pzZ1 " << pzZ1
			<< "\n ptZ1 " << ptZ1
			<< "\n EZ1 "  << EZ1
			<< "\n Y_Z1 " << Y_Z1
			<< "\n massZ1 " << massZ1
			<< "\n indexlep1 " << indexlep1Z1
			<< "\n indexlep2 " << indexlep2Z1
		    
			<< endl;
      
      
      // vector<candidateZ> pTcleanedgoodZ_nosharedlept;
      // pTcleanedgoodZ_nosharedlept.push_back(pTcleanedgoodZ.at(indexZ1));
      
      // for (int i=0;i<pTcleanedgoodZ.size();i++){
      //   if ( pTcleanedgoodZ.at(i).ilept1==indexlep1Z1 || pTcleanedgoodZ.at(i).ilept1==indexlep2Z1 ||  
      // 	    pTcleanedgoodZ.at(i).ilept2==indexlep1Z1 || pTcleanedgoodZ.at(i).ilept2==indexlep2Z1
      // 	    ) continue;
      //   pTcleanedgoodZ_nosharedlept.push_back(pTcleanedgoodZ.at(i));
      // }
      
      // for (int i=0;i<pTcleanedgoodZ_nosharedlept.size();i++){
      //   cout << "Z Mass surviving= " << pTcleanedgoodZ_nosharedlept.at(i).massvalue << endl;
      // }
      
      // vector<candidateZ> pTcleanedgoodZ_noZ1;
      // for (int i=0;i<pTcleanedgoodZ.size();i++){
      //   if (pTcleanedgoodZ.at(i).massvalue == massZ1 ) continue; 
      //   pTcleanedgoodZ_noZ1.push_back(pTcleanedgoodZ.at(i));
      //   cout << "Z Mass surviving _ no Z1 " << pTcleanedgoodZ.at(i).massvalue << endl;                                                                              
      // }
      
      // ZZ objects with alternative pairing
      array<int,2> icleanedgoodZs;
      vector<std::array<int, 2> > icleanedgoodZsv;
      
      for (int i=0;i<pTcleanedgoodZ.size();i++){
	
	
	float lepton1ch=-999., lepton2ch=-999.;
	lepton1ch=pTcleanedgoodZ.at(i).charge1;
	lepton2ch=pTcleanedgoodZ.at(i).charge2;
	
	TLorentzVector Lepton1,Lepton2,Lepton1Correction,Lepton2Correction;
	
	Lepton1.SetPtEtaPhiM(pTcleanedgoodZ.at(i).pt1, pTcleanedgoodZ.at(i).eta1, pTcleanedgoodZ.at(i).phi1, 0.000511);
	
	if (pTcleanedgoodZ.at(i).ilept1_FSR) {
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == pTcleanedgoodZ.at(i).ilept1 && iLp_tagEM[ p ] == 1 ){
	      Lepton1Correction.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]],RECOPFPHOT_PHI[iLp[p]],0);
	      Lepton1+=Lepton1Correction;
	    }
	  }
	}	 
	Lepton2.SetPtEtaPhiM(pTcleanedgoodZ.at(i).pt2, pTcleanedgoodZ.at(i).eta2, pTcleanedgoodZ.at(i).phi2, 0.000511);
	if (pTcleanedgoodZ.at(i).ilept2_FSR) {
	  for( int p = 0; p < Nphotons; ++p ){
	    if( iLp_l[ p ] == pTcleanedgoodZ.at(i).ilept2 && iLp_tagEM[ p ] == 1 ){
	      Lepton2Correction.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]],RECOPFPHOT_PHI[iLp[p]],0);
	      Lepton2+=Lepton2Correction;
	    }
	  }
	}
	
	
	for (int j=i+1;j<pTcleanedgoodZ.size();j++){
	  float massZa=-999.,massZb=-999.;	 	 
	  cout << "Masses= " << pTcleanedgoodZ.at(i).massvalue << " " << pTcleanedgoodZ.at(j).massvalue << endl;
	  
	  float lepton3ch=-999., lepton4ch=-999.;
	  lepton3ch=pTcleanedgoodZ.at(j).charge1;
	  lepton4ch=pTcleanedgoodZ.at(j).charge2;
	  
	  TLorentzVector Lepton3,Lepton4,Lepton3Correction,Lepton4Correction,DiLeptonZa,DiLeptonZb;
	  Lepton3.SetPtEtaPhiM(pTcleanedgoodZ.at(j).pt1, pTcleanedgoodZ.at(j).eta1, pTcleanedgoodZ.at(j).phi1, 0.000511);
	  if (pTcleanedgoodZ.at(j).ilept1_FSR) {
	    for( int p = 0; p < Nphotons; ++p ){
	      if( iLp_l[ p ] == pTcleanedgoodZ.at(j).ilept1 && iLp_tagEM[ p ] == 1 ){
		Lepton3Correction.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]],RECOPFPHOT_PHI[iLp[p]],0);
		Lepton3+=Lepton3Correction;
	      }
	    }
	  }	 
	  Lepton4.SetPtEtaPhiM(pTcleanedgoodZ.at(j).pt2, pTcleanedgoodZ.at(j).eta2, pTcleanedgoodZ.at(j).phi2, 0.000511);
	  if (pTcleanedgoodZ.at(j).ilept2_FSR) {
	    for( int p = 0; p < Nphotons; ++p ){
	      if( iLp_l[ p ] == pTcleanedgoodZ.at(j).ilept2 && iLp_tagEM[ p ] == 1 ){
		Lepton4Correction.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[p]],RECOPFPHOT_ETA[iLp[p]],RECOPFPHOT_PHI[iLp[p]],0);
		Lepton4+=Lepton4Correction;
	      }
	    }
	  }
	  
	  if (lepton1ch*lepton3ch<0){
	    DiLeptonZa=Lepton1+Lepton3;
	    DiLeptonZb=Lepton2+Lepton4;
	  }  
	  if (lepton1ch*lepton4ch<0){
	    DiLeptonZa=Lepton1+Lepton4;
	    DiLeptonZb=Lepton2+Lepton3;
	  }
	  
	  if (fabs(DiLeptonZa.M()-Zmass) < fabs(DiLeptonZb.M()-Zmass)) {
	    massZa=DiLeptonZa.M();
	    massZb=DiLeptonZb.M();
	  }
	  else {
	    massZa=DiLeptonZb.M();
	    massZb=DiLeptonZa.M();
	  }	  
	  
	  cout << "mass Za and b= " << massZa << " " << massZb << endl;
	  if ( fabs(massZa-Zmass) < fabs(massZ1-Zmass) && massZb<12) continue; // exclude the original pairs
	  
	  //	 cout << "i j " << i << " " << j << endl;
	  icleanedgoodZs={i,j};
	  icleanedgoodZsv.push_back(icleanedgoodZs);
	}
      }
      
      cout << "How many ZZ objects are present? " << icleanedgoodZsv.size() << endl;
      if (icleanedgoodZsv.size()==0) continue;
      
      // // sort Z by mass value
      // struct SortCandByClosestToZ {
      //   bool operator()( candidateZ c1, candidateZ c2) {
      //          return (fabs(c1.massvalue - Zmass) < fabs(c2.massvalue - Zmass));
      //   }
      // };
      // std::sort(vcandZ2.begin(),vcandZ2.end(),SortCandByClosestToZ());                                                                                                       
      
      double pxZ2 = 0;  //Z2 kinematics                                                                                                                                         
      double pyZ2 = 0;
      double pzZ2 = 0;
      double ptZ2 = 0;
      double EZ2 = 0;
      double Y_Z2 = -9;
      double massZ2 = 0;
      double massZ2_noFSR = 0;
      double sum_ptZ2 = 0.;
      int indexlep1Z2=-1;
      int indexlep2Z2=-1;
      int indexZ2=-1;
      
      // identify pairs  Z1+Z2s
      vector<int> vindexZZ;
      for (int l=0;l<icleanedgoodZsv.size();l++){
	cout << "Z+Z pairs with masses: " << pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(0)).massvalue << " " <<  pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(1)).massvalue << endl;
	if ( (icleanedgoodZsv.at(l).at(0)==indexZ1 &&  
	      !(pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(1)).ilept1==indexlep1Z1 || pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(1)).ilept2==indexlep1Z1 || 
		pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(1)).ilept1==indexlep2Z1 || pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(1)).ilept2==indexlep2Z1)) 
	     || 
	     (icleanedgoodZsv.at(l).at(1)==indexZ1 && 
	      !(pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(0)).ilept1==indexlep1Z1 || pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(0)).ilept2==indexlep1Z1 ||
		pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(0)).ilept1==indexlep2Z1 || pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(0)).ilept2==indexlep2Z1))
	     ) 
	  {
	    cout << "Found a Z1 + X  pair with no lepton shared with masses= " << pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(0)).massvalue << " " <<  pTcleanedgoodZ.at(icleanedgoodZsv.at(l).at(1)).massvalue << endl;
	    vindexZZ.push_back(l);
	  }
      }
      
      cout << "How many Z1+Z2s objects are present? " << vindexZZ.size() << endl;
      
      // If more than one Z2 (couple to the same Z1), choose the one with the highest-pT Z2 leptons
      if (vindexZZ.size()==1){
	cout << "Just one Z1+Z2 pair" << endl;
	if (icleanedgoodZsv.at(vindexZZ.at(0)).at(0)==indexZ1) indexZ2=icleanedgoodZsv.at(vindexZZ.at(0)).at(1); 
	if (icleanedgoodZsv.at(vindexZZ.at(0)).at(1)==indexZ1) indexZ2=icleanedgoodZsv.at(vindexZZ.at(0)).at(0);
      }
      else {
	cout << "more than one Z2 (couple to the same Z1), choose the one with the highest-pT Z2 leptons" << endl;
	int indexZ2tmp=-999;
	float sumpT=-999.,tmpsumpT=-999.;
	
	for (int ll=0;ll<vindexZZ.size();ll++){
	  if (icleanedgoodZsv.at(vindexZZ.at(ll)).at(0)==indexZ1) indexZ2tmp=icleanedgoodZsv.at(vindexZZ.at(ll)).at(1);
	  if (icleanedgoodZsv.at(vindexZZ.at(ll)).at(1)==indexZ1) indexZ2tmp=icleanedgoodZsv.at(vindexZZ.at(ll)).at(0);
	  if (indexZ2tmp>=0){
	    
	    TLorentzVector LorentzZ1,LorentzZ2,LorentzZZ;
	    LorentzZ1.SetPxPyPzE(pTcleanedgoodZ.at(indexZ1).pxZ, pTcleanedgoodZ.at(indexZ1).pyZ, pTcleanedgoodZ.at(indexZ1).pzZ, pTcleanedgoodZ.at(indexZ1).EZ);
	    LorentzZ2.SetPxPyPzE(pTcleanedgoodZ.at(indexZ2tmp).pxZ, pTcleanedgoodZ.at(indexZ2tmp).pyZ, pTcleanedgoodZ.at(indexZ2tmp).pzZ, pTcleanedgoodZ.at(indexZ2tmp).EZ);
	    LorentzZZ=LorentzZ1+LorentzZ2;
	    if (LorentzZZ.M()<70.) continue; // cut m4l>70
	    cout << "Passed m4l>70. cut"<< endl;
	    
	    tmpsumpT=
	      pTcleanedgoodZ.at(indexZ2tmp).pt1+
	      pTcleanedgoodZ.at(indexZ2tmp).pt2;
	    if (tmpsumpT>sumpT) {
	      sumpT=tmpsumpT;
	      indexZ2=indexZ2tmp;	     
	    }
	  }       
	}
      }
      
      // Z2
      if (indexZ2<0) continue;
      cout << "The highest pT leptons Z2 has mass= " <<  pTcleanedgoodZ.at(indexZ2).massvalue << endl;
      
      massZ2 = pTcleanedgoodZ.at(indexZ2).massvalue;	 
      pxZ2 = pTcleanedgoodZ.at(indexZ2).pxZ;
      pyZ2 = pTcleanedgoodZ.at(indexZ2).pyZ;
      pzZ2 = pTcleanedgoodZ.at(indexZ2).pzZ;
      EZ2  = pTcleanedgoodZ.at(indexZ2).EZ;	 
      ptZ2 = sqrt( pxZ2*pxZ2 + pyZ2*pyZ2 );
      sum_ptZ2 = pTcleanedgoodZ.at(indexZ2).pt1+pTcleanedgoodZ.at(indexZ2).pt2;	 
      Y_Z2 = 0.5 * log ( (EZ2 + pzZ2)/(EZ2 - pzZ2) );
      indexlep1Z2=pTcleanedgoodZ.at(indexZ2).ilept1;
      indexlep2Z2=pTcleanedgoodZ.at(indexZ2).ilept2;
      
      // Z1 and Z2 final 
      cout << "Z1 has index= " << indexZ1 << "  Z2 has index= " << indexZ2 << endl;
      
      cout << "PTs= " << RECOELE_PT[indexlep1Z1] << " " << RECOELE_PT[indexlep2Z1] << " " <<  RECOELE_PT[indexlep1Z2] << " " << RECOELE_PT[indexlep2Z2]<< endl; 
      
      if (std::isnan(massZ2)) {
	cout << "No Z2 found" << endl;
	continue; 
      }
      else {
	if( debug ) cout  << "\n Final Z2 properties: "  
			  << "\n pxZ2 " << pxZ2
			  << "\n pyZ2 " << pyZ2
			  << "\n pzZ2 " << pzZ2
			  << "\n ptZ2 " << ptZ2
			  << "\n EZ2 "  << EZ2
			  << "\n Y_Z2 " << Y_Z2
			  << "\n massZ2 " << massZ2
			  << "\n indexlep1_Z2 " << indexlep1Z2
			  << "\n indexlep2_Z2 " << indexlep2Z2
			  << endl;
      }
      
      if( debug && has_FSR_Z2) {
      	cout  << " Z2 has FSR! " << endl;
     	cout  << "  pi " << pi2 << " --> index: " << iLp[pi2] << " associated lepton: " << iLp_l[pi2] << " (= "<< iLe[i2]<<" ? )  tag: " << iLp_tagEM[pi2] << endl;
     	cout  << "  pj " << pj2 << " --> index: " << iLp[pj2] << " associated lepton: " << iLp_l[pj2] << " (= "<< iLe[j2]<<" ? )  tag: " << iLp_tagEM[pj2] << endl;
      }
      
      ++N_4b ;  // fill counter
      N_4b_w=N_4b_w+newweight;

      if( N_Z2_pairs > 1 ) {
	++N_4c ;  // fill counter
	N_4c_w=N_4c_w+newweight;
      }
      
      
      ++N_4d ;  // fill counter
      N_4d_w=N_4d_w+newweight;
      
      // **** Step 5:

      // // Execute Efficiency Reweighting
      int iL_4L[4] = {indexlep1Z1,indexlep2Z1,indexlep1Z2,indexlep2Z2};

      Double_t eff_weight = 1.;
      
      for(int i = 0; i < 4; ++i){
        Double_t Pt = RECOELE_PT[ iL_4L[i] ]; 
        Double_t Eta = RECOELE_ETA[ iL_4L[i] ]; 
	Double_t SC_Eta = RECOELE_scl_Eta[ iL_4L[i] ]; //Reham

	Double_t Tmp_Pt = Pt;
	if (Pt>500)Tmp_Pt = 500; //overflow
        
	// if( (MC_type == "Spring16" || MC_type=="Moriond17") && DATA_type == "NO"){
	if(  MC_type == "Fall17"  && DATA_type == "NO"){ 
          if(RECOELE_isGap[ iL_4L[i] ]==0){
            //Int_t binx = ele_scale_2017->GetXaxis()->FindBin(Eta);
	    Int_t binx = ele_scale_2017->GetXaxis()->FindBin(SC_Eta);
            Int_t biny = ele_scale_2017->GetYaxis()->FindBin(Tmp_Pt);
            if (ele_scale_2017->GetBinContent(binx,biny)>0.) eff_weight*=ele_scale_2017->GetBinContent(binx,biny); 
          }
          else if(RECOELE_isGap[ iL_4L[i] ]==1){
	    // Int_t binx = ele_scale_2017_gap->GetXaxis()->FindBin(Eta);
	    Int_t binx = ele_scale_2017_gap->GetXaxis()->FindBin(SC_Eta);
            Int_t biny = ele_scale_2017_gap->GetYaxis()->FindBin(Tmp_Pt);
            if (ele_scale_2017_gap->GetBinContent(binx,biny)>0.) eff_weight*=ele_scale_2017_gap->GetBinContent(binx,biny); 
          }
          
        }
      }//end for loop 

      //add Reco efficency

      for(int i = 0; i < 4; ++i){
        Double_t Pt = RECOELE_PT[ iL_4L[i] ]; 
        Double_t Eta = RECOELE_ETA[ iL_4L[i] ]; 
	Double_t SC_Eta = RECOELE_scl_Eta[ iL_4L[i] ]; //Reham

        
	// if( (MC_type == "Spring16" || MC_type=="Moriond17") && DATA_type == "NO"){
	if(  MC_type == "Fall17"  && DATA_type == "NO"){ 
          if(Pt<=20.){

	    
            //Int_t binx = ele_Reco_eff_2017_lowEt->GetXaxis()->FindBin(Eta);
	    Int_t binx = ele_Reco_eff_2017_lowEt->GetXaxis()->FindBin(SC_Eta);
            Int_t biny = ele_Reco_eff_2017_lowEt->GetYaxis()->FindBin(Pt);
            if (ele_Reco_eff_2017_lowEt->GetBinContent(binx,biny)>0.) eff_weight*=ele_Reco_eff_2017_lowEt->GetBinContent(binx,biny); 
          }
          else if(Pt>20.){

	    Double_t Tmp_Pt =Pt;
	    if(Pt>500.)Tmp_Pt =500.;

	    // Int_t binx = ele_Reco_eff_2017_highEt->GetXaxis()->FindBin(Eta);
	    Int_t binx = ele_Reco_eff_2017_highEt->GetXaxis()->FindBin(SC_Eta);
            Int_t biny = ele_Reco_eff_2017_highEt->GetYaxis()->FindBin(Tmp_Pt);
            if (ele_Reco_eff_2017_highEt->GetBinContent(binx,biny)>0.) eff_weight*=ele_Reco_eff_2017_highEt->GetBinContent(binx,biny); 
          }
          
        }
      }//end Reco eff.

      
      if (DATA_type == "2017") eff_weight=1.; 
      // // Changing the weight for pileup and efficiency
      if (eff_weight>0.) newweight=weight*pT_weight*pu_weight*weight_kfactor*eff_weight;
      else newweight=weight*pT_weight*pu_weight*weight_kfactor;
      
      cout << "Starting weight + pileup + efficiency= " << newweight << endl;
      if(debug) cout << "Efficiency Weight for the 4l: " << eff_weight << " Final weight= " << newweight << endl;

 
      TLorentzVector hP4,Z1P4,Z2P4;
      Z1P4.SetPxPyPzE(pxZ1,pyZ1,pzZ1,EZ1);
      Z2P4.SetPxPyPzE(pxZ2,pyZ2,pzZ2,EZ2);
      hP4 = Z1P4 + Z2P4;
      cout << "Mass of 4l best candidate= " << hP4.M() << endl;
      double mass4l=hP4.M();      

      hM4l_5->Fill( mass4l,newweight );
      
      hMZ1_5->Fill( massZ1,newweight );
      hPtZ1_5->Fill( ptZ1,newweight );
      hYZ1_5->Fill( Y_Z1,newweight );

      hMZ2_5->Fill( massZ2,newweight );
      hPtZ2_5->Fill( ptZ2,newweight );
      hYZ2_5->Fill( Y_Z2,newweight );

      
      // sort index by pt (kinematics not corrected for FSR)
      int ipt[4] ;
      double tmp_pt[4];
      cout << "PTs= " << RECOELE_PT[indexlep1Z1] << " " << RECOELE_PT[indexlep2Z1] << " " <<  RECOELE_PT[indexlep1Z2] << " " << RECOELE_PT[indexlep2Z2]<< endl;
      int indexleptonfinal[4]={indexlep1Z1,indexlep2Z1,indexlep1Z2,indexlep2Z2};
      cout << "PTs= " << RECOELE_PT[indexleptonfinal[0]] << " " << RECOELE_PT[indexleptonfinal[1]] << " " <<  RECOELE_PT[indexleptonfinal[2]] << " " << RECOELE_PT[indexleptonfinal[3]]<< endl;

      for(int i = 0; i < 4; ++i) tmp_pt[i] =  RECOELE_PT[indexleptonfinal[i]];
      for(int i = 0; i < 4; ++i){
        double tmp_max_pt = 0;
      	int jj = i;
        for(int j = 0; j < 4; ++j){
      		if( tmp_pt[j] > tmp_max_pt ){
      			tmp_max_pt = tmp_pt[j];
      			jj  = j;
      		}
      	}
      	ipt[i] = indexleptonfinal[jj];
      	tmp_pt[jj] = 0;
      }
      if(debug)
	for(int i = 0; i < 4; ++i)
	  cout << "\n ipt[" << i << "] = " << ipt[i] << "\t pt = " << RECOELE_PT[ ipt[i] ] << endl; 
      //end sorting
     

      // Format lepton syncronization                                                                                                                                      
      // {run}:{lumi}:{event}:{pdgId}:{pT:.2f}:{eta:.2f}:{phi:.2f}{SIP:.2f}:{PFChargedHadIso:.2f}:{PFNeutralHadIso:.2f}:{PFPhotonIso:.2f}:{PUCorr:.2f}:{combRelIsoPF:.3f}:{eleBDT:.3f}:{photpT:.2f}:{photDR:.2f}:{photRelIso:.2f}          
      
      Char_t leptformat[20000];

      for(int i = 0; i < 4; ++i){      

	float dummy=0.;
	int flagFSR=0;
	int pfsr=-999;

	for( int p = 0; p < Nphotons; ++p ){
	  if( iLp_l[ p ] == ipt[i] && iLp_tagEM[ p ] == 1 )  {
	    cout << "Electron with pT= " << RECOELE_PT[ipt[i]] << " has associated a photon with pT= " << RECOPFPHOT_PT[iLp[p]] <<  endl;
	    flagFSR=1;
	    pfsr=p;
	  }
	}
	if (flagFSR==1){
	  sprintf (leptformat,"FormatLept=%d:%d:%d:%d:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.3f:%.3f:%.2f:%.2f:%.2f",
		   Run,LumiSection,Event,
		   int(-11*RECOELE_CHARGE[ipt[i]]),
		   RECOELE_PT[ipt[i]],RECOELE_ETA[ipt[i]],RECOELE_PHI[ipt[i]],RECOELE_SIP[ipt[i]],
		   RECOELE_PFchHad[ipt[i]],RECOELE_PFneuHad[ipt[i]],RECOELE_PFphoton[ipt[i]],RHO_ele,RECOELE_PFX_rho[ipt[i]],RECOELE_mvaNonTrigV0[ipt[i]],
		   RECOPFPHOT_PT[iLp[pfsr]],RECOPFPHOT_DR[iLp[pfsr]],RECOPFPHOT_PFX_rho[iLp[pfsr]]
		   );
	}
	else {
	  sprintf (leptformat,"FormatLept=%d:%d:%d:%d:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.2f:%.3f:%.3f",
		   Run,LumiSection,Event,
		   int(-11*RECOELE_CHARGE[ipt[i]]),
		   RECOELE_PT[ipt[i]],RECOELE_ETA[ipt[i]],RECOELE_PHI[ipt[i]],RECOELE_SIP[ipt[i]],
		   RECOELE_PFchHad[ipt[i]],RECOELE_PFneuHad[ipt[i]],RECOELE_PFphoton[ipt[i]],RHO_ele,RECOELE_PFX_rho[ipt[i]],RECOELE_mvaNonTrigV0[ipt[i]]
		   );
	}
	  
	output_txt  << leptformat << endl;

       if( i == 0){
	 hPtLep1_5->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
	 hYLep1_5->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
	 hIsoLep1_5->Fill( RECOELE_PFX_dB[ ipt[i] ],newweight ) ;
	 hSipLep1_5->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
	 hIpLep1_5->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
	 hIpErLep1_5->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 1){
	 hPtLep2_5->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
	 hYLep2_5->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
	 hIsoLep2_5->Fill( RECOELE_PFX_dB[ ipt[i] ],newweight ) ;
	 hSipLep2_5->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
	 hIpLep2_5->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
	 hIpErLep2_5->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 2){
	 hPtLep3_5->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
	 hYLep3_5->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
	 hIsoLep3_5->Fill( RECOELE_PFX_dB[ ipt[i] ],newweight ) ;
	 hSipLep3_5->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
	 hIpLep3_5->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
	 hIpErLep3_5->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 3){
	 hPtLep4_5->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
	 hYLep4_5->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
	 hIsoLep4_5->Fill( RECOELE_PFX_dB[ ipt[i] ],newweight ) ;
	 hSipLep4_5->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
	 hIpLep4_5->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
	 hIpErLep4_5->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       
      }//end fill leptons
      
      
      
      //hIso_5->Fill( Iso_max,newweight ) ;
      //hSip_5->Fill( Sip_max,newweight ) ;
      //hIp_5->Fill( Ip_max,newweight ) ;

      // N.B. Update the Isolation values and correct the 4 momenta of leptons for FSR
      for(int i = 0; i < Ne_good; ++i){
	int flagFSR=0;
	int pfsr=-999;
	
	for( int p = 0; p < Nphotons; ++p ){
	  cout << "Index of lepton with photon ISR= " << iLp_l[ p ] << " and final lepton index= " << iLe[i] << endl;
	  if( iLp_l[ p ] == iLe[i] && iLp_tagEM[ p ] == 1 )  {
	    cout << "Electron with pT= " << RECOELE_PT[iLe[i]] << " has associated a photon with pT= " << RECOPFPHOT_PT[iLp[p]] <<  endl;
	    // RECOELE_PFX_rho_new[iLe[i]]=
	    //   (RECOELE_PFchHad[iLe[i]]+
	    //    max(0.,RECOELE_PFneuHad[iLe[i]]+
	    // 	   (RECOELE_PFphoton[iLe[i]]-RECOPFPHOT_PT[iLp[p]] )-
	    // 	   max(RHO_ele,0.0)*(EffectiveArea)))/RECOELE_PT[iLe[i]];
	   
	    flagFSR=1;
	    pfsr=p;
	  }
	}
	
	if (flagFSR==1){
	  cout << "Before correcting for FSR; electron pT= " << RECOELE_PT[iLe[i]] << " Eta= " << RECOELE_ETA[iLe[i]] << " Phi= " << RECOELE_PHI[iLe[i]] << endl;
	  TLorentzVector Lept,LeptCorrection;
	  Lept.SetPtEtaPhiM(RECOELE_PT[iLe[i]], RECOELE_ETA[iLe[i]], RECOELE_PHI[iLe[i]],0.000511);
	  LeptCorrection.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[pfsr]],RECOPFPHOT_ETA[iLp[pfsr]],RECOPFPHOT_PHI[iLp[pfsr]],0);
	  Lept+=LeptCorrection;
	  RECOELE_PT[iLe[i]]=Lept.Pt();
	  RECOELE_ETA[iLe[i]]=Lept.Eta();
	  RECOELE_PHI[iLe[i]]=Lept.Phi();
	  cout << "After correcting for FSR; muon pT= " << RECOELE_PT[iLe[i]] << " Eta= " << RECOELE_ETA[iLe[i]] << " Phi= " << RECOELE_PHI[iLe[i]] << endl;
	}
      }
      cout << "Kinematics of leptons corrected for FSR photons (if existing)" << endl;

      // sort index by pt (kinematics is now corrected for FSR)
      int iptcorr[4] ;
      double tmp_ptcorr[4];
      cout << "PTs= " << RECOELE_PT[indexlep1Z1] << " " << RECOELE_PT[indexlep2Z1] << " " <<  RECOELE_PT[indexlep1Z2] << " " << RECOELE_PT[indexlep2Z2]<< endl;
      int indexleptonfinalcorr[4]={indexlep1Z1,indexlep2Z1,indexlep1Z2,indexlep2Z2};
      cout << "PTs= " << RECOELE_PT[indexleptonfinalcorr[0]] << " " << RECOELE_PT[indexleptonfinalcorr[1]] << " " <<  RECOELE_PT[indexleptonfinalcorr[2]] << " " << RECOELE_PT[indexleptonfinalcorr[3]]<< endl;

      for(int i = 0; i < 4; ++i) tmp_ptcorr[i] =  RECOELE_PT[indexleptonfinalcorr[i]];
      for(int i = 0; i < 4; ++i){
        double tmp_max_ptcorr = 0;
      	int jj = i;
        for(int j = 0; j < 4; ++j){
      		if( tmp_ptcorr[j] > tmp_max_ptcorr ){
      			tmp_max_ptcorr = tmp_ptcorr[j];
      			jj  = j;
      		}
      	}
      	iptcorr[i] = indexleptonfinalcorr[jj];
      	tmp_ptcorr[jj] = 0;
      }
      if(debug)
      	  for(int i = 0; i < 4; ++i)
      	  	cout << "\n iptcorr[" << i << "] = " << iptcorr[i] << "\t pt = " << RECOELE_PT[ iptcorr[i] ] << endl; 
      //end sorting with FSR included
      
      
      // // **** Step 6:
      //  // QCD suppression: mll>4 GeV cut on all OS-SF pairs (4/4)           
     //if( min_mass_2L <= 4 ) continue ;
     
     ++N_6 ;  // fill counter
     N_6_w=N_6_w+newweight;

     // **** Step 7:
     // mass4l > 70 


     if( debug ) cout  << "\n Step 7: m4l > 70" 
      		 << "\n mass4l " << mass4l
		 << endl; 

     if( mass4l <= 70 ) continue ;
     
     ++N_7 ;  // fill counter
     N_7_w=N_7_w+newweight;
      

     // Calculate per-event-error 
     // We use the "simpler" method of adding the deviation to the m4l error from individual
     // muon pT error in quadrature to get the per-event-error
     
     vector<TLorentzVector> vtlv; vtlv.clear();
     vector<double> vpterr; vpterr.clear(); 
     vector<double> vpterrcorr; vpterrcorr.clear();
     
     /*

      for(unsigned int i=0; i< 4; i++){
	TLorentzVector tlv; tlv.SetPtEtaPhiM(RECOELE_PT[indexleptonfinal[i]], RECOELE_ETA[indexleptonfinal[i]], RECOELE_PHI[indexleptonfinal[i]], 0.000511); 
	vtlv.push_back(tlv);
        double pterr = ((RECOELE_PT[]/RECOELE_P[indexleptonfinal[i]])*RECOELE_PTError[indexleptonfinal[i]]);
	vpterr.push_back(pterr);

	Double_t Pt = RECOELE_PT[indexleptonfinal[i]]; 
        Double_t Eta = RECOELE_ETA[indexleptonfinal[i]]; 
        int PDGID = 11;
	double corr=1.;

	if(fabs(PDGID) == 11){    
 	  if( MC_type == "Fall11" && DATA_type == "NO"){
 	    Int_t binx = ebe_ele_mc42x->GetXaxis()->FindBin(Pt);
 	    Int_t biny = ebe_ele_mc42x->GetYaxis()->FindBin(Eta);
 	    if (ebe_ele_mc42x->GetBinContent(binx,biny)>0.) corr=ebe_ele_mc42x->GetBinContent(binx,biny);
            vpterrcorr.push_back(double((RECOELE_PT[indexleptonfinal[i]]/RECOELE_P[iL_4L[i]])*RECOELE_PTError[indexleptonfinal[i]]*corr));        	    
 	  }
 	  else if( MC_type == "Summer12" && DATA_type == "NO"){
	    Int_t binx = ebe_ele_mc53x->GetXaxis()->FindBin(Pt);
	    Int_t biny = ebe_ele_mc53x->GetYaxis()->FindBin(Eta);
	    if (ebe_ele_mc53x->GetBinContent(binx,biny)>0.) corr=ebe_ele_mc53x->GetBinContent(binx,biny);
            vpterrcorr.push_back(double((RECOELE_PT[indexleptonfinal[i]]/RECOELE_P[indexleptonfinal[i]])*RECOELE_PTError[indexleptonfinal[i]]*corr));
	  }
	  else if( MC_type == "Phys14" && DATA_type == "NO"){
	    vpterrcorr.push_back(pterr);
	  }
	  else if( MC_type == "NO" && DATA_type == "2011"){
	    Int_t binx = ebe_ele_reco42x->GetXaxis()->FindBin(Pt);
	    Int_t biny = ebe_ele_reco42x->GetYaxis()->FindBin(Eta);
	    if (ebe_ele_reco42x->GetBinContent(binx,biny)>0.) corr=ebe_ele_reco42x->GetBinContent(binx,biny);
            vpterrcorr.push_back(double((RECOELE_PT[indexleptonfinal[i]]/RECOELE_P[indexleptonfinal[i]])*RECOELE_PTError[indexleptonfinal[i]]*corr));
	  }
	  else if( MC_type == "NO" && DATA_type == "2012"){
	    Int_t binx = ebe_ele_reco53x->GetXaxis()->FindBin(Pt);
	    Int_t biny = ebe_ele_reco53x->GetYaxis()->FindBin(Eta);
	    if (ebe_ele_reco53x->GetBinContent(binx,biny)>0.) corr=ebe_ele_reco53x->GetBinContent(binx,biny);
	    vpterrcorr.push_back(double((RECOELE_PT[indexleptonfinal[i]]/RECOELE_P[indexleptonfinal[i]])*RECOELE_PTError[indexleptonfinal[i]]*corr));
	  }	  
	  else if( MC_type == "NO" && DATA_type == "2013"){
	    vpterrcorr.push_back(pterr);
	  }
 	}
      }   
     */

      double massErr = masserror(vtlv, vpterr);
      if(debug) cout << "Mass Error: " << massErr << endl;
      
      double massErrCorr = masserror(vtlv, vpterrcorr);
      if(debug) cout << "MassErrCorr: " << massErrCorr << endl;

//       // if has_FSR_Z1 sumErr=RECOPFPHOT_PTError[iphot]**2+RECOPFPHOT_PTError[iphot]**2
//       for(int i = 0; i < 4; ++i){
//         Double_t Pt = RECOMU_PT[ iL_4L[i] ]; 
//         Double_t Eta = RECOMU_ETA[ iL_4L[i] ]; 
//         int PDGID = 13;

// 	if(fabs(PDGID) == 13){    
// 	  if( MC_type == "Fall11" && DATA_type == "NO"){
// 	    Int_t binx = ebe_mu_mc42x->GetXaxis()->FindBin(Pt);
// 	    Int_t biny = ebe_mu_mc42x->GetYaxis()->FindBin(Eta);
// 	    sumErr+=pow( double(RECOMU_muInnertrkPTError[iL_4L[i]] * ebe_mu_mc42x->GetBinContent(binx,biny)),2.); 
// 	    vpterrcorr.push_back(double(RECOMU_muInnertrkPTError[iL_4L[i]]* ebe_mu_mc42x->GetBinContent(binx,biny)));
// 	  }
// 	  else if( MC_type == "Summer12" && DATA_type == "NO"){
// 	    Int_t binx = ebe_mu_mc53x->GetXaxis()->FindBin(Pt);
// 	    Int_t biny = ebe_mu_mc53x->GetYaxis()->FindBin(Eta);
// 	    sumErr+=pow (double(RECOMU_muInnertrkPTError[ iL_4L[i]] * ebe_mu_mc53x->GetBinContent(binx,biny)),2.);
// 	    vpterrcorr.push_back(double(RECOMU_muInnertrkPTError[iL_4L[i]] * ebe_mu_mc53x->GetBinContent(binx,biny)));
// 	  }
// 	  else if( MC_type == "NO" && DATA_type == "2011"){
// 	    Int_t binx = ebe_mu_reco42x->GetXaxis()->FindBin(Pt);
// 	    Int_t biny = ebe_mu_reco42x->GetYaxis()->FindBin(Eta);
// 	    sumErr+=pow(double(RECOMU_muInnertrkPTError[ iL_4L[i] ] * ebe_mu_reco42x->GetBinContent(binx,biny)),2.);
// 	    vpterrcorr.push_back(double(RECOMU_muInnertrkPTError[iL_4L[i]] * ebe_mu_reco42x->GetBinContent(binx,biny)));
// 	  }
// 	  else if( MC_type == "NO" && DATA_type == "2012"){
// 	    Int_t binx = ebe_mu_reco53x->GetXaxis()->FindBin(Pt);
// 	    Int_t biny = ebe_mu_reco53x->GetYaxis()->FindBin(Eta);
// 	    sumErr+=pow(double(RECOMU_muInnertrkPTError[ iL_4L[i] ] * ebe_mu_reco53x->GetBinContent(binx,biny)),2.);
// 	    vpterrcorr.push_back(double(RECOMU_muInnertrkPTError[iL_4L[i]] * ebe_mu_reco53x->GetBinContent(binx,biny)));
// 	  }
	  
// 	}

//       }

      ++N_7_PFMET;
      N_7_PFMET_w=N_7_PFMET_w+newweight;
      
      
     hM4l_7->Fill( mass4l,newweight );
     hPFMET_7->Fill(RECO_PFMET,newweight);
     hMZ1_7->Fill( massZ1,newweight );
     hPtZ1_7->Fill( ptZ1,newweight );
     hYZ1_7->Fill( Y_Z1,newweight );
     
     hMZ2_7->Fill( massZ2,newweight );
     hPtZ2_7->Fill( ptZ2,newweight );
     hYZ2_7->Fill( Y_Z2,newweight );
     
     // Leptons PT, ETA, Phi, Isol corrected for FSR
     for(int i = 0; i < 4; ++i){
       
       if( i == 0){
	 hPtLep1_7->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
	 hYLep1_7->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
	 hIsoLep1_7->Fill( RECOELE_PFX_rho_new[ ipt[i] ],newweight ) ;
	 hSipLep1_7->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
	 hIpLep1_7->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
	 hIpErLep1_7->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 1){
	 hPtLep2_7->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
	 hYLep2_7->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
	 hIsoLep2_7->Fill( RECOELE_PFX_rho_new[ ipt[i] ],newweight ) ;
	 hSipLep2_7->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
	 hIpLep2_7->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
	 hIpErLep2_7->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 2){
	 hPtLep3_7->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
	 hYLep3_7->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
	 hIsoLep3_7->Fill( RECOELE_PFX_rho_new[ ipt[i] ],newweight ) ;
	 hSipLep3_7->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
	 hIpLep3_7->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
	 hIpErLep3_7->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 3){
	 hPtLep4_7->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
	 hYLep4_7->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
	 hIsoLep4_7->Fill( RECOELE_PFX_rho_new[ ipt[i] ],newweight ) ;
	 hSipLep4_7->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
	 hIpLep4_7->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
	 hIpErLep4_7->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       
     }//end fill leptons
     
     //hIso_7->Fill( Iso_max,newweight ) ;
     //hSip_7->Fill( Sip_max,newweight ) ;
     //hIp_7->Fill( Ip_max,newweight ) ;
     

     //Reham to synchronize

     //output_txt << "EVENT_7:"<< Run<< ":" << LumiSection<< ":" << Event	<< ":" << mass4l<< ":" << massZ1<< ":" << massZ2 << ":" <<RECOELE_PT[ ipt[0] ]<< ":" <<RECOELE_PT[ ipt[1] ]<< ":" <<RECOELE_PT[ ipt[2] ]<<":"<<RECOELE_PT[ ipt[3] ]<<":"<< f_D_bkg_kin<<endl;
     output_txt << "EVENT_7:"<< Run<< ":" << LumiSection<< ":" << Event << ":" << mass4l<< ":" << massZ1<< ":" << massZ2 << ":" <<f_D_bkg_kin<< ":" <<RECOELE_PT[ ipt[0] ]<< ":" <<RECOELE_PT[ ipt[1] ]<< ":" <<RECOELE_PT[ ipt[2] ]<< ":" <<RECOELE_PT[ ipt[3] ]<< ":" <<RECOELE_ETA[ ipt[0] ]<< ":" <<RECOELE_ETA[ ipt[1] ]<< ":" <<RECOELE_ETA[ ipt[2] ]<< ":" <<RECOELE_ETA[ ipt[3] ]<< ":" <<RECOELE_PHI[ ipt[0] ]<< ":" <<RECOELE_PHI[ ipt[1] ]<< ":" <<RECOELE_PHI[ ipt[2] ]<< ":" <<RECOELE_PHI[ ipt[3] ]<<endl;

     // **** Step 8:
     // mass4l > 100 

     if( debug ) cout  << "\n Step 8: m4l > 100" 
      		 << "\n mass4l " << mass4l
		 << endl; 

     if( mass4l <= 70.) continue ;     
     ++N_8 ;  // fill counter     
     N_8_w=N_8_w+newweight;
     cout << " N_8_w= " << N_8_w << " " << newweight << endl;

     hM4l_8->Fill( mass4l,newweight );
     if (mass4l>=100. && mass4l<=800.) hM4l_8_100_800->Fill( mass4l,newweight );
     
     hMZ1_8->Fill( massZ1,newweight );
     hPtZ1_8->Fill( ptZ1,newweight );
     hYZ1_8->Fill( Y_Z1,newweight );
     
     hMZ2_8->Fill( massZ2,newweight );
     hPtZ2_8->Fill( ptZ2,newweight );
     hYZ2_8->Fill( Y_Z2,newweight );
     
     hMZ1_noFSR_8->Fill( massZ1_noFSR,newweight );
     hMZ2_noFSR_8->Fill( massZ2_noFSR,newweight );
     
     
     for(int i = 0; i < 4; ++i){
       
       if( i == 0){
	 hPtLep1_8->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
	 hYLep1_8->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
	 hIsoLep1_8->Fill( RECOELE_PFX_rho_new[ ipt[i] ],newweight ) ;
	 hSipLep1_8->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
	 hIpLep1_8->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
	 hIpErLep1_8->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 1){
	 hPtLep2_8->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
	 hYLep2_8->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
	 hIsoLep2_8->Fill( RECOELE_PFX_rho_new[ ipt[i] ],newweight ) ;
	 hSipLep2_8->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
	 hIpLep2_8->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
	 hIpErLep2_8->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 2){
	 hPtLep3_8->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
	 hYLep3_8->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
	 hIsoLep3_8->Fill( RECOELE_PFX_rho_new[ ipt[i] ],newweight ) ;
	 hSipLep3_8->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
	 hIpLep3_8->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
	 hIpErLep3_8->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       if( i == 3){
	 hPtLep4_8->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;
	 hYLep4_8->Fill(  RECOELE_ETA[ ipt[i] ],newweight ) ;
	 hIsoLep4_8->Fill( RECOELE_PFX_rho_new[ ipt[i] ],newweight ) ;
	 hSipLep4_8->Fill( RECOELE_SIP[ ipt[i] ],newweight ) ;
	 hIpLep4_8->Fill( RECOELE_IP[ ipt[i] ],newweight ) ;
	 hIpErLep4_8->Fill( RECOELE_IPERROR[ ipt[i] ],newweight ) ;
       }
       
     }//end fill leptons
     
     //hIso_8->Fill( Iso_max,newweight ) ;
     //hSip_8->Fill( Sip_max,newweight ) ;
     //hIp_8->Fill( Ip_max,newweight ) ;
     
     hPFMET_8->Fill(RECO_PFMET,newweight);
     hLogXPFMET_8->Fill(RECO_PFMET,newweight); 
     hLogLinXPFMET_8->Fill(RECO_PFMET,newweight);

     // //VBF
     cout << "Starting VBF analysis " << endl;
     
     //Basic cuts to jets AND delta R section
     int njets_pass=0;
     TLorentzVector JET1,JET2;
     int jet1=-999,jet2=-999;      
     int jetfail[100];

     for(int i=0;i<100;i++) jetfail[i]=0;
     
     for(int i=0;i<RECO_PFJET_N;i++){
       cout<<i<<" Jet with pt= "<<RECO_PFJET_PT[i]<<" ETA "<<RECO_PFJET_ETA[i]<<" PUID "<<RECO_PFJET_PUID[i] << " PUID_MVA "<< RECO_PFJET_PUID_MVA[i]<<endl;
       //if(RECO_PFJET_PUID[i]==1 && RECO_PFJET_PT[i]>30. && fabs(RECO_PFJET_ETA[i])<4.7 ){ // NO PU ID temporary
       if(RECO_PFJET_PT[i]>30. && fabs(RECO_PFJET_ETA[i])<4.7 ){       

      	 //Check that jet has deltaR>0.4 away from any tight lepton corrected for FSR
	 for(int mu = 0; mu < N_good; ++mu){
	   if (fabs(RECOMU_SIP[iL[mu]])>=4.) continue;
      	   if (RECOMU_PFX_dB_new[iL[mu]]>=0.35) continue;
	   double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOMU_PHI[iL[mu]]),2) + pow(RECO_PFJET_ETA[i] - RECOMU_ETA[iL[mu]],2));
	   cout << "1st lepton muon: " << " pT=" << RECOMU_PT[iL[mu]] <<" deltaR "<< deltaR <<endl;	   
	   if (deltaR<0.4){
	     jetfail[i]=1;
     	     cout << " jetfail " << jetfail[i] <<endl;
	     break;
     	   }
     	 }
	 
      	 for(int ele = 0; ele < Ne_good; ++ele){
      	   if (fabs(RECOELE_SIP[iLe[ele]])>=4.) continue;
	   // if (RECOELE_PFX_rho_new[iLe[ele]]>=0.35) continue;
      	   double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOELE_PHI[iLe[ele]]),2) + pow(RECO_PFJET_ETA[i] - RECOELE_ETA[iLe[ele]],2));
     	   cout << "1st lepton electron: " << " pT=" << RECOELE_PT[iLe[ele]] <<" deltaR "<< deltaR <<endl;
	   if (deltaR<0.4){
     	     jetfail[i]=1;
     	     cout << " jetfail " << jetfail[i] <<endl;
	     break;
     	   }
     	 }

	 // cleaning w.r.t FSR photons attached to leptons
	 for(int j=0.;j<Nphotons;j++) {
           if (iLp_l[j]!=-1 && (iLp_tagEM[j]==0 || iLp_tagEM[j]==1) ) {
	     if (iLp_tagEM[j]==0) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[j]] << " attached to a muon with pT= " << RECOMU_PT[iLp_l[j]] << endl;
	     if (iLp_tagEM[j]==1) cout << "There is photon with pT= " << RECOPFPHOT_PT[iLp[j]] << " attached to a electron with pT= " << RECOELE_PT[iLp_l[j]] << endl;
	     double deltaR = sqrt( pow(DELTAPHI(RECO_PFJET_PHI[i],RECOPFPHOT_PHI[iLp[j]]),2) + pow(RECO_PFJET_ETA[i] - RECOPFPHOT_ETA[iLp[j]],2));
	     if (deltaR<0.4){
	       jetfail[i]=1;
	       cout << " jetfail " << jetfail[i] <<endl;
	       break;
	     }
	   }
         }
	 // 
	 
	 if (jetfail[i]==0){
	   cout<< " PASS jet " <<i<<" PT= "<<RECO_PFJET_PT[i]<<" ETA= "<<RECO_PFJET_ETA[i]<<" PUID= "<<RECO_PFJET_PUID[i]<<endl;
	   njets_pass++;
	   if (njets_pass==1){
	     jet1=i;
	     JET1.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
	     
	     f_jet1_highpt_pt = RECO_PFJET_PT[i];
	     f_jet1_highpt_eta = RECO_PFJET_ETA[i];
	     f_jet1_highpt_phi = RECO_PFJET_PHI[i];
	     f_jet1_highpt_e = RECO_PFJET_ET[i];

	   }
	   if (njets_pass==2){
	     jet2=i;
	     JET2.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));

	     f_jet2_highpt_pt = RECO_PFJET_PT[i];
	     f_jet2_highpt_eta = RECO_PFJET_ETA[i];
	     f_jet2_highpt_phi = RECO_PFJET_PHI[i];
	     f_jet2_highpt_e = RECO_PFJET_ET[i];
	   }
	   if (njets_pass==3){
	     f_jet3_highpt_pt = RECO_PFJET_PT[i];
	     f_jet3_highpt_eta = RECO_PFJET_ETA[i];
	     f_jet3_highpt_phi = RECO_PFJET_PHI[i];
	     f_jet3_highpt_e = RECO_PFJET_ET[i];
	   }

	   //@ fills all jets
	  if(njets_pass <= maxJets){
            f_jets_highpt_pt[njets_pass-1]                           = RECO_PFJET_PT[i];
	    f_jets_highpt_pt_error[njets_pass-1]                     = RECO_PFJET_PT[i]-RECO_PFJET_PT_UncDn[i];//The error is symmetric
            f_jets_highpt_eta[njets_pass-1]                          = RECO_PFJET_ETA[i];
            f_jets_highpt_phi[njets_pass-1]                          = RECO_PFJET_PHI[i];
            f_jets_highpt_et[njets_pass-1]                           = RECO_PFJET_ET[i];
            f_jets_highpt_area[njets_pass-1]                         = RECO_PFJET_AREA[i];
            f_jets_highpt_charged_hadron_multiplicity[njets_pass-1]  = RECO_PFJET_CHARGED_HADRON_MULTIPLICITY[i];
            f_jets_highpt_charged_multiplicity[njets_pass-1]         = RECO_PFJET_CHARGED_MULTIPLICITY[i];
            f_jets_highpt_electron_multiplicity[njets_pass-1]        = RECO_PFJET_ELECTRON_MULTIPLICITY[i];
            f_jets_highpt_hf_em_multiplicity[njets_pass-1]           = RECO_PFJET_HF_EM_MULTIPLICITY[i];
            f_jets_highpt_hf_hadron_multiplicity[njets_pass-1]       = RECO_PFJET_HF_HADRON_MULTIPLICTY[i];
            f_jets_highpt_muon_multiplicity[njets_pass-1]            = RECO_PFJET_MUON_MULTIPLICITY[i];
            f_jets_highpt_neutral_hadron_multiplicity[njets_pass-1]  = RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY[i];
            f_jets_highpt_neutral_multiplicity[njets_pass-1]         = RECO_PFJET_NEUTRAL_MULTIPLICITY[i];
            f_jets_highpt_photon_multiplicity[njets_pass-1]          = RECO_PFJET_PHOTON_MULTIPLICITY[i];
            f_jets_highpt_charged_hadron_energy[njets_pass-1]        = RECO_PFJET_CHARGED_HADRON_ENERGY[i];
            f_jets_highpt_neutral_hadron_energy[njets_pass-1]        = RECO_PFJET_NEUTRAL_HADRON_ENERGY[i];
            f_jets_highpt_photon_energy[njets_pass-1]                = RECO_PFJET_PHOTON_ENERGY[i];
            f_jets_highpt_electron_energy[njets_pass-1]              = RECO_PFJET_ELECTRON_ENERGY[i];
            f_jets_highpt_muon_energy[njets_pass-1]                  = RECO_PFJET_MUON_ENERGY[i];
            f_jets_highpt_hf_hadron_energy[njets_pass-1]             = RECO_PFJET_HF_HADRON_ENERGY[i];
            f_jets_highpt_hf_em_energy[njets_pass-1]                 = RECO_PFJET_HF_EM_ENERGY[i];
            f_jets_highpt_charged_em_energy[njets_pass-1]            = RECO_PFJET_CHARGED_EM_ENERGY[i];
            f_jets_highpt_charged_mu_energy[njets_pass-1]            = RECO_PFJET_CHARGED_MU_ENERGY[i];
            f_jets_highpt_neutral_em_energy[njets_pass-1]            = RECO_PFJET_NEUTRAL_EM_ENERGY[i];
            f_jets_highpt_ptd[njets_pass-1]                          = RECO_PFJET_PTD[i];
            f_jets_highpt_ncomponents[njets_pass-1]                  = RECO_PFJET_NCOMPONENTS[i];
            for(int ic=0; ic<maxJetsComponents; ++ic){
              f_jets_highpt_component_pdgid[njets_pass-1][ic]        = RECO_PFJET_COMPONENT_PDGID[i][ic];
              f_jets_highpt_component_pt[njets_pass-1][ic]           = RECO_PFJET_COMPONENT_PT[i][ic];
              f_jets_highpt_component_eta[njets_pass-1][ic]          = RECO_PFJET_COMPONENT_ETA[i][ic];
              f_jets_highpt_component_phi[njets_pass-1][ic]          = RECO_PFJET_COMPONENT_PHI[i][ic];
              f_jets_highpt_component_energy[njets_pass-1][ic]       = RECO_PFJET_COMPONENT_E[i][ic];
              f_jets_highpt_component_charge[njets_pass-1][ic]       = RECO_PFJET_COMPONENT_CHARGE[i][ic];
              f_jets_highpt_component_mt[njets_pass-1][ic]           = RECO_PFJET_COMPONENT_TRANSVERSE_MASS[i][ic];
              f_jets_highpt_component_xvertex[njets_pass-1][ic]      = RECO_PFJET_COMPONENT_XVERTEX[i][ic];
              f_jets_highpt_component_yvertex[njets_pass-1][ic]      = RECO_PFJET_COMPONENT_YVERTEX[i][ic];
              f_jets_highpt_component_zvertex[njets_pass-1][ic]      = RECO_PFJET_COMPONENT_ZVERTEX[i][ic];
              f_jets_highpt_component_vertex_chi2[njets_pass-1][ic]  = RECO_PFJET_COMPONENT_VERTEX_CHI2[i][ic];
            }
          }

	 }

       }
       else{
      	 jetfail[i]=1;
       }
       //cout<<" JETFAIL "<<jetfail[i]<<endl;
     }
     
     //Number of jets and mJJ,delta eta cuts // categories
     cout << "Number of jets passing the VBF selection is = " << njets_pass << endl;
     
     // b-tagged jets - ordered in pT
     int n_bjets=0;
     int index_bjets[2]={-999,-999};
     
     for (int i=0;i<50;i++){
       //if (cSV_BTagJet_DISCR[i] > 0.460){ // Loose
       if (cSV_BTagJet_DISCR[i] > 0.8484){ // Medium     
	 if(cSV_BTagJet_PT[i]>30. && fabs(cSV_BTagJet_ETA[i])<4.7 ) cout << "Found a bjet (pT>30 and |eta|<2.4) with pT= " << cSV_BTagJet_PT[i] << endl;	 
	 n_bjets++;
	 if (n_bjets==1) index_bjets[0]=i; 
	 if (n_bjets==2) index_bjets[1]=i;	   
       }
     }
     cout << "Number of b-jets passing the selection= " << n_bjets << endl;
     
     // VBF-tagged category - category 2
     int njets_vbf=0;
     float deltajj=-999.,massjj=-999.;
     float D_jet=-999.; //fisher discriminant
     TLorentzVector JET1_VBF, JET2_VBF,mJJ;
     int index_VBFjets[2]={-999,-999};
       
     for (int i=0;i<RECO_PFJET_N;i++){
       if (jetfail[i]!=0 ) continue; 
       for (int j=i+1;j<RECO_PFJET_N;j++){
	 if (jetfail[j]!=0 ) continue;
	 JET1_VBF.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
	 JET2_VBF.SetPtEtaPhiE(RECO_PFJET_PT[j],RECO_PFJET_ETA[j],RECO_PFJET_PHI[j],RECO_PFJET_ET[j]*TMath::CosH(RECO_PFJET_ETA[j]));
	 mJJ=JET1_VBF+JET2_VBF;
	 //cout<<mJJ.M()<<endl;              
	 
	 deltajj=fabs(JET1_VBF.Eta()-JET2_VBF.Eta());
	 massjj=mJJ.M();
	 D_jet=0.18*fabs(JET1_VBF.Eta()-JET2_VBF.Eta())+1.92E-4*mJJ.M();
	 if (debug) cout << "JetpT / Jet 2 pT= " << JET1_VBF.Pt() << " " << JET2_VBF.Pt()<< " Mass jet jet= " << massjj << " Delta jet jet= " << deltajj << " Fisher D_jet= " << D_jet << endl;
	 if (D_jet>0.5) {
	   njets_vbf++;
	   if (njets_vbf==1) {
	     index_VBFjets[0]=i;
	     index_VBFjets[1]=j;    
	     JET1_VBF.SetPtEtaPhiE(RECO_PFJET_PT[i],RECO_PFJET_ETA[i],RECO_PFJET_PHI[i],RECO_PFJET_ET[i]*TMath::CosH(RECO_PFJET_ETA[i]));
	     JET2_VBF.SetPtEtaPhiE(RECO_PFJET_PT[j],RECO_PFJET_ETA[j],RECO_PFJET_PHI[j],RECO_PFJET_ET[j]*TMath::CosH(RECO_PFJET_ETA[j]));
	   }
	 }
	 hDjj_8->Fill( deltajj,newweight );
	 hMjj_8->Fill( massjj,newweight );
	 hVD_8->Fill( D_jet,newweight );
	 hMjj_Djj_8->Fill( massjj,deltajj,newweight);
       }
     }	 
     
     
     //VH-hadronic-tagged category - category 4
     int n_match_bjets=0;
     for(int i=0;i<RECO_PFJET_N;i++){
       if (jetfail[i]!=0) continue;
       for (int j=0;j<50;j++){
	 if (cSV_BTagJet_PT[j]==RECO_PFJET_PT[i] && cSV_BTagJet_DISCR[j]>0.8484 && cSV_BTagJet_PT[j]>30. && fabs(cSV_BTagJet_ETA[j])<4.7) {
	   //if (cSV_BTagJet_PT[j]==RECO_PFJET_PT[i] && cSV_BTagJet_DISCR[j]>0.460 && cSV_BTagJet_PT[j]>30. && fabs(cSV_BTagJet_ETA[j])<4.7) {
	   //if (cSV_BTagJet_DISCR[j]>0.8484 && cSV_BTagJet_PT[j]>30. && fabs(cSV_BTagJet_ETA[j])<4.7) {
	   n_match_bjets++;
	   //if (n_match_bjets==1) jet1=i;
	   //if (n_match_bjets==2) jet2=i;
	 }
       }
     }
     cout << "Number of matched b-jets to jets= " << n_match_bjets << endl;
     
     // VH jets
     int n_jets_mjj_VH=0;
     TLorentzVector JET1_VH,JET2_VH,mJJ_VH;
     int index_VHjets[2]={-999,-999};
     
     for(int j=0;j<RECO_PFJET_N;j++){
        if (jetfail[j]!=0) continue;
        if(RECO_PFJET_PT[j]>40. && fabs(RECO_PFJET_ETA[j])<2.4){
	  //cout << "Possible jet for VH categories with pT= " << RECO_PFJET_PT[j] << endl;
	 for (int k=j+1;k<RECO_PFJET_N; k++){
	   if (jetfail[k]!=0) continue;
	   if(RECO_PFJET_PT[k]>40. && fabs(RECO_PFJET_ETA[k])<2.4){	     
	     JET1_VH.SetPtEtaPhiE(RECO_PFJET_PT[j],RECO_PFJET_ETA[j],RECO_PFJET_PHI[j],RECO_PFJET_ET[j]*TMath::CosH(RECO_PFJET_ETA[j]));
	     JET2_VH.SetPtEtaPhiE(RECO_PFJET_PT[k],RECO_PFJET_ETA[k],RECO_PFJET_PHI[k],RECO_PFJET_ET[k]*TMath::CosH(RECO_PFJET_ETA[k]));
	     mJJ_VH=JET1_VH+JET2_VH;
	     
	     if (mJJ_VH.M()>60. && mJJ_VH.M()<120.) {
	       cout << "Jets for VH categories passing, pT>40, |eta|<2.4 and 60<mjj<120 with pT= " << RECO_PFJET_PT[j] << " " << RECO_PFJET_PT[k] << endl;
	       n_jets_mjj_VH++ ;
	       if (n_jets_mjj_VH==1) { // always the pair with highest pT jets
		 index_VHjets[0]=j;
		 index_VHjets[1]=k;		 
	       }	     
	     }
	   }
	 }
	}
     }	
     
     // Categorization
     int category=-999;

     //  exactly 4 leptons + at least 2 jets with Djet>0.5 + at most 1 b-tag jet in the event  - category 2
     if (Ne_good==4 && njets_vbf>=1 && n_bjets<=1){
       
       jet1=index_VBFjets[0];
       jet2=index_VBFjets[1];
       JET1.SetPtEtaPhiE(RECO_PFJET_PT[jet1],RECO_PFJET_ETA[jet1],RECO_PFJET_PHI[jet1],RECO_PFJET_ET[jet1]*TMath::CosH(RECO_PFJET_ETA[jet1]));
       JET2.SetPtEtaPhiE(RECO_PFJET_PT[jet2],RECO_PFJET_ETA[jet2],RECO_PFJET_PHI[jet2],RECO_PFJET_ET[jet2]*TMath::CosH(RECO_PFJET_ETA[jet2]));
       
       category=2;
       cout<<"PASS VBF category= "<<Event<<endl;
       ++N_VBF;
       N_VBF_w=N_VBF_w+newweight;
       
       cout << "EVENT CANDIDATE VBF: \n" 
	    << " N " << jentry 
	    << " RUN " << Run
	    << " EVENT " << Event	
	    << " LumiSection " << LumiSection
	    << " massZ1 " << massZ1
	    << " massZ2 " << massZ2
	    << " mass4l " << mass4l; 
       if( has_FSR_Z1 || has_FSR_Z2 ) cout << " FSR " << endl;	
       else cout << " noFSR " << endl;
       
       
       output_txt_vbf << "EVENT CANDIDATE VBF: \n" 
		      << " N " << jentry 
		      << " RUN " << Run
		      << " EVENT " << Event	
		      << " LumiSection " << LumiSection
		      << " massZ1 " << massZ1
		      << " massZ2 " << massZ2
		      << " mass4l " << mass4l; 
       if( has_FSR_Z1 || has_FSR_Z2 ) output_txt_vbf << " FSR " << endl;	
       else output_txt_vbf << " noFSR " << endl;	   
       
     }
     
     // exactly 4 leptons + at least one pair of jets passing { |η|<2.4, pT>40, 60<mjj<120 } + pT,4l>m4l - category 4
     else if ((n_bjets==2 && Ne_good==4) || (n_jets_mjj_VH>=1 && Ne_good==4 && hP4.Pt()>hP4.M())) {
       if (n_bjets==2 && N_good==4) {	 
	 JET1.SetPtEtaPhiE(cSV_BTagJet_PT[index_bjets[0]],cSV_BTagJet_ETA[index_bjets[0]],
			   cSV_BTagJet_PHI[index_bjets[0]],cSV_BTagJet_ET[index_bjets[0]]*TMath::CosH(cSV_BTagJet_ETA[index_bjets[0]]));
	 JET2.SetPtEtaPhiE(cSV_BTagJet_PT[index_bjets[1]],cSV_BTagJet_ETA[index_bjets[1]],
			   cSV_BTagJet_PHI[index_bjets[1]],cSV_BTagJet_ET[index_bjets[1]]*TMath::CosH(cSV_BTagJet_ETA[index_bjets[1]]));
	 cout << "category 4 with two b-jets with pT" << JET1.Pt() << " and " << JET2.Pt() << endl;
       }
       else if (n_jets_mjj_VH>=1) {	 
	 jet1=index_VHjets[0];
	 jet2=index_VHjets[1];
	 JET1.SetPtEtaPhiE(RECO_PFJET_PT[jet1],RECO_PFJET_ETA[jet1],RECO_PFJET_PHI[jet1],RECO_PFJET_ET[jet1]*TMath::CosH(RECO_PFJET_ETA[jet1]));
	 JET2.SetPtEtaPhiE(RECO_PFJET_PT[jet2],RECO_PFJET_ETA[jet2],RECO_PFJET_PHI[jet2],RECO_PFJET_ET[jet2]*TMath::CosH(RECO_PFJET_ETA[jet2]));
       }
       
       category=4;
     }

     // VH-leptonic-tagged category - category 3
     else if (njets_pass<=2 && n_bjets==0 && Ne_good>=5) category=3;

     // ttH-tagged category - category 5
     else if ((njets_pass>=3 && n_bjets>=1) || Ne_good>=5) category=5;

     // 1-jet-tagged category - category 1
     else if (njets_pass>=1) category=1;
     
     // untagged category (#0) - category 0
     else if (category<0 ) category =0;

     // MonoHiggs category
     if (abs(mass4l-125.)<=10. && Ne_good==4 && n_bjets<=1) category=6;

     cout << "Categorization completed - event category is= " << category << endl;

     // filling branches in the reduced tree
     f_weight = newweight;
     
     f_int_weight = -1;
     
     f_pu_weight = pu_weight;
     f_eff_weight = eff_weight;

     f_run = Run;
     f_event = Event;
     f_lumi = LumiSection;

     f_lept1_pt = RECOELE_PT[iptcorr[0]] ;
     f_lept1_pt_error = RECOELE_PTError[iptcorr[0]] ;
     f_lept1_eta = RECOELE_ETA[iptcorr[0]] ;
     f_lept1_phi = RECOELE_PHI[iptcorr[0]];
     f_lept1_charge = RECOELE_CHARGE[iptcorr[0]];
     f_lept1_pfx = RECOELE_PFX_rho_new[iptcorr[0]];
     f_lept1_sip = RECOELE_SIP[iptcorr[0]];
     f_lept1_mvaid = RECOELE_mvaNonTrigV0[iptcorr[0]];
     f_lept2_pt = RECOELE_PT[iptcorr[1]] ;
     f_lept2_pt_error = RECOELE_PTError[iptcorr[1]] ;
     f_lept1_pt_error = RECOMU_mubesttrkPTError[iptcorr[0]] ;
     f_lept2_eta = RECOELE_ETA[iptcorr[1]] ;
     f_lept2_phi = RECOELE_PHI[iptcorr[1]];
     f_lept2_charge = RECOELE_CHARGE[iptcorr[1]];
     f_lept2_pfx = RECOELE_PFX_rho_new[iptcorr[1]];
     f_lept2_sip = RECOELE_SIP[iptcorr[1]];
     f_lept2_mvaid = RECOELE_mvaNonTrigV0[iptcorr[1]];
     f_lept3_pt = RECOELE_PT[iptcorr[2]] ;
     f_lept3_pt_error = RECOELE_PTError[iptcorr[2]] ;
     f_lept3_eta = RECOELE_ETA[iptcorr[2]] ;
     f_lept3_phi = RECOELE_PHI[iptcorr[2]];
     f_lept3_charge = RECOELE_CHARGE[iptcorr[2]];
     f_lept3_pfx = RECOELE_PFX_rho_new[iptcorr[2]];
     f_lept3_sip = RECOELE_SIP[iptcorr[2]];
     f_lept3_mvaid = RECOELE_mvaNonTrigV0[iptcorr[2]];
     f_lept4_pt = RECOELE_PT[iptcorr[3]] ;
     f_lept4_pt_error = RECOELE_PTError[iptcorr[3]] ;
     f_lept4_eta = RECOELE_ETA[iptcorr[3]] ;
     f_lept4_phi = RECOELE_PHI[iptcorr[3]];
     f_lept4_charge = RECOELE_CHARGE[iptcorr[3]];
     f_lept4_pfx = RECOELE_PFX_rho_new[iptcorr[3]];
     f_lept4_sip = RECOELE_SIP[iptcorr[3]];
     f_lept4_mvaid = RECOELE_mvaNonTrigV0[iptcorr[3]];
     //     f_iso_max = Iso_max;
     //f_sip_max = Sip_max;
     f_Z1mass = massZ1;
     f_Z2mass = massZ2;
     f_angle_costhetastar = -99.;
     f_angle_costheta1 = -99.;
     f_angle_costheta2 = -99.;
     f_angle_phi = -99.;
     f_angle_phistar1 = -99.;
     f_pt4l = -99.;
     f_eta4l = -99;
     f_mass4l = mass4l;
     f_mass4lErr = massErr;
     f_njets_pass = njets_pass;
     f_deltajj = -999.;
     f_massjj = -999.;
     f_D_jet = -999.;
     
     if (njets_pass<1) {
       f_jet1_pt = -999;
       f_jet1_eta = -999;
       f_jet1_phi = -999;
       f_jet1_e = -999;
       f_jet2_pt = -999;
       f_jet2_eta = -999;
       f_jet2_phi = -999;
       f_jet2_e = -999;
     }
     else if (njets_pass==1) { // store the highest pT jets passing the criteria
       f_jet1_pt = RECO_PFJET_PT[jet1];
       f_jet1_eta = RECO_PFJET_ETA[jet1];
       f_jet1_phi = RECO_PFJET_PHI[jet1];
       f_jet1_e = RECO_PFJET_ET[jet1];
     }
     else if (njets_pass>=2) { // store the second highest pT jets passing the criteria
       f_jet1_pt = RECO_PFJET_PT[jet1];
       f_jet1_eta = RECO_PFJET_ETA[jet1];
       f_jet1_phi = RECO_PFJET_PHI[jet1];
       f_jet1_e = RECO_PFJET_ET[jet1];
       f_jet2_pt = RECO_PFJET_PT[jet2];
       f_jet2_eta = RECO_PFJET_ETA[jet2];
       f_jet2_phi = RECO_PFJET_PHI[jet2];
       f_jet2_e = RECO_PFJET_ET[jet2];

       TLorentzVector JET1new,JET2new,mJJnew;
       JET1new.SetPtEtaPhiE(RECO_PFJET_PT[jet1],RECO_PFJET_ETA[jet1],RECO_PFJET_PHI[jet1],RECO_PFJET_ET[jet1]*TMath::CosH(RECO_PFJET_ETA[jet1]));
       JET2new.SetPtEtaPhiE(RECO_PFJET_PT[jet2],RECO_PFJET_ETA[jet2],RECO_PFJET_PHI[jet2],RECO_PFJET_ET[jet2]*TMath::CosH(RECO_PFJET_ETA[jet2]));
       mJJnew=JET1new+JET2new;
       f_deltajj = fabs(JET1new.Eta()-JET2new.Eta());
       f_massjj = mJJnew.M();
       f_D_jet =0.18*fabs(JET1new.Eta()-JET2new.Eta())+1.92E-4*mJJnew.M();
     }

     f_pfmet=RECO_PFMET;
     f_pfmet_JetEnUp = RECO_PFMET_JetEnUp;
     f_pfmet_JetEnDn = RECO_PFMET_JetEnDn;
     f_pfmet_ElectronEnUp = RECO_PFMET_ElectronEnUp;
     f_pfmet_ElectronEnDn = RECO_PFMET_ElectronEnDn;
     f_pfmet_MuonEnUp = RECO_PFMET_MuonEnUp;
     f_pfmet_MuonEnDn = RECO_PFMET_MuonEnDn;
     f_pfmet_JetResUp = RECO_PFMET_JetResUp;
     f_pfmet_JetResDn = RECO_PFMET_JetResDn;
     f_pfmet_UnclusteredEnUp = RECO_PFMET_UnclusteredEnUp;
     f_pfmet_UnclusteredEnDn = RECO_PFMET_UnclusteredEnDn;
     f_pfmet_PhotonEnUp = RECO_PFMET_PhotonEnUp;
     f_pfmet_PhotonEnDn = RECO_PFMET_PhotonEnDn;

     f_genmet=MC_GENMET;
     f_category=category;
     f_Ngood=Ne_good;
     f_Nbjets=n_bjets;

     //MELA discriminant - Leoton kinematics already corrected for FSR
     TLorentzVector L11P4,L12P4,L21P4,L22P4;
     float angle_costheta1,angle_costheta2,angle_Phi,angle_costhetastar,angle_Phi1;

     cout << "Final: pT of final leptons from Z1 and Z2 ="
	  << RECOELE_PT[indexleptonfinal[0]] << " " << RECOELE_PT[indexleptonfinal[1]] << " " << RECOELE_PT[indexleptonfinal[2]] << " " << RECOELE_PT[indexleptonfinal[3]] <<endl;
     
     L11P4.SetPtEtaPhiM(RECOELE_PT[indexleptonfinal[0]], RECOELE_ETA[indexleptonfinal[0]], RECOELE_PHI[indexleptonfinal[0]], 0.000511);
     L12P4.SetPtEtaPhiM(RECOELE_PT[indexleptonfinal[1]], RECOELE_ETA[indexleptonfinal[1]], RECOELE_PHI[indexleptonfinal[1]], 0.000511);
     L21P4.SetPtEtaPhiM(RECOELE_PT[indexleptonfinal[2]], RECOELE_ETA[indexleptonfinal[2]], RECOELE_PHI[indexleptonfinal[2]], 0.000511);
     L22P4.SetPtEtaPhiM(RECOELE_PT[indexleptonfinal[3]], RECOELE_ETA[indexleptonfinal[3]], RECOELE_PHI[indexleptonfinal[3]], 0.000511);
     
     //Correcting 4-momenta for FSR photons
     // cout << "Correcting leptons for FSR..... "<< endl;

     // TLorentzVector L1correctionP4;
     // for(int i = 0; i < 4; ++i){

     //   int flagFSR=0;
     //   int pfsr=-999;
       
     //   for( int p = 0; p < Nphotons; ++p ){
     // 	 cout << "Index of lepton with photon ISR= " << iLp_l[ p ] << " and final lepton index= " << indexleptonfinal[i] << endl;
     // 	 if( iLp_l[ p ] == indexleptonfinal[i] && iLp_tagEM[ p ] == 1 )  {
     // 	   cout << "Muon with pT= " << RECOELE_PT[indexleptonfinal[i]] << " has associated a photon with pT= " << RECOPFPHOT_PT[iLp[p]] <<  endl;
     // 	   flagFSR=1;
     // 	   pfsr=p;
     // 	 }
     //   }

     //   if (flagFSR==1) L1correctionP4.SetPtEtaPhiM(RECOPFPHOT_PT[iLp[pfsr]],RECOPFPHOT_ETA[iLp[pfsr]],RECOPFPHOT_PHI[iLp[pfsr]],0);
     //   if (i==0) L11P4+=L1correctionP4;
     //   if (i==1) L12P4+=L1correctionP4;
     //   if (i==2) L21P4+=L1correctionP4;
     //   if (i==3) L22P4+=L1correctionP4;
     // }
     
     
     //Assigning correct PID numbers
     int L11PID,L12PID,L21PID,L22PID;
     if (RECOELE_CHARGE[indexleptonfinal[0]] == 1){
       L11PID=-11;
       L12PID=+11;
       f_lept1_pdgid=-11;
       f_lept2_pdgid=11;
     }
     else if (RECOELE_CHARGE[indexleptonfinal[0]] == -1){
       L11PID=11;
       L12PID=-11;
       f_lept1_pdgid=11;
       f_lept2_pdgid=-11;
     }
     else cout << "What the hell?!?!"<<endl;
     if (RECOELE_CHARGE[indexleptonfinal[2]] == 1){
       L21PID=-11;
       L22PID=11;
       f_lept3_pdgid=-11;
       f_lept4_pdgid=11;
     }
     else if (RECOELE_CHARGE[indexleptonfinal[2]] == -1){
       L21PID=11;
       L22PID=-11;
       f_lept3_pdgid=11;
       f_lept4_pdgid=-11;
     }
     else cout << "What the hell?!?!"<<endl;
     
          
     double massofhiggs=hP4.M();
     double massofZ1=Z1P4.M();
     double massofZ2=Z2P4.M();
     float Psig,Pbkg;
     float psPsig,psPbkg,psD,gravPsig,gravPbkg,gravD;     
     double pt4l=0.,eta4l=0.;

     // 4l transverse mass + met
     float DPHI = 0.;
     DPHI = DELTAPHI(hP4.Phi(),RECO_PFMET_PHI);

     float m4l_T=sqrt(2*hP4.Pt()*RECO_PFMET*(1-cos(DPHI)));
     hM4l_T_8->Fill(m4l_T,newweight);       
     hLogLinXM4l_T_8->Fill(m4l_T,newweight);
     f_mT=m4l_T;
     hDPHI_8->Fill(fabs(DPHI),newweight);
     f_dphi=DPHI;

     hNgood_8->Fill(f_Ngood,newweight);
     hNbjets_8->Fill(f_Nbjets,newweight);

     cout << "Mass of Higgs passed to MELA= " << massofhiggs << endl;

     //if(massofhiggs>100. && massofhiggs<1000. && massofZ2>4.){
     if(massofhiggs>0.){
       
       // MEM       
       vector<TLorentzVector> partP;
       partP.push_back(L11P4);
       partP.push_back(L12P4);
       partP.push_back(L21P4);
       partP.push_back(L22P4);
       
       vector<int> partId;
       partId.push_back(L11PID);
       partId.push_back(L12PID);
       partId.push_back(L21PID);
       partId.push_back(L22PID);
       
       
       vector<TLorentzVector> partPprod;
       vector<int> partIdprod;
       partPprod.push_back(L11P4);
       partPprod.push_back(L12P4);
       partPprod.push_back(L21P4);
       partPprod.push_back(L22P4);      
       if (JET1.Pt()>0.) partPprod.push_back(JET1);
       if (JET2.Pt()>0.) partPprod.push_back(JET2); // Can also use partPprod.push_back(nullFourVector) instead for integrated VBF MEs
       
       partIdprod.push_back(L11PID);
       partIdprod.push_back(L12PID);
       partIdprod.push_back(L21PID);
       partIdprod.push_back(L22PID);
       if (JET1.Pt()>0.) partIdprod.push_back(0);
       if (JET2.Pt()>0.) partIdprod.push_back(0); // For leptonic ZH in the future, this could actually be the opposite lepton flavor
       

       double p0plus_VAJHU, bkg_VAMCFM, p0plus_m4l, bkg_m4l;
       double p0minus_VAJHU, Dgg10_VAMCFM;
       double phjj_VAJHU, pvbf_VAJHU;
       
       int a=combinedMEM.computeME(MEMNames::kSMHiggs, MEMNames::kJHUGen, partP, partId, p0plus_VAJHU); // Calculation of SM gg->H->4l JHUGen ME      
       //cout << "a= "  << p0plus_VAJHU << endl;
       int b=combinedMEM.computeME(MEMNames::k0minus, MEMNames::kJHUGen, partP, partId, p0minus_VAJHU); // Calculation of PS (0-, fa3=1) gg->H->4l JHUGen ME 
       int c=combinedMEM.computeME(MEMNames::kggHZZ_10, MEMNames::kMCFM, partP, partId, Dgg10_VAMCFM); // Direct calculation of Dgg (D^kin for off-shell) from MCFM MEs
       int d=combinedMEM.computeME(MEMNames::kqqZZ, MEMNames::kMCFM, partP, partId, bkg_VAMCFM); // qq->4l background calculation from MCFM
       combinedMEM.computePm4l(partP,partId, MEMNames::kNone, p0plus_m4l, bkg_m4l); // m4l probabilities for signal and background, nominal resolution
       if (njets_pass>=2){
	 int f=combinedMEM.computeME(MEMNames::kJJ_SMHiggs_GG, MEMNames::kJHUGen, partPprod, partIdprod, phjj_VAJHU); // SM gg->H+2j
	 int g=combinedMEM.computeME(MEMNames::kJJ_SMHiggs_VBF, MEMNames::kJHUGen, partPprod, partIdprod, pvbf_VAJHU);  // SM VBF->H
       }

       f_D_bkg_kin = p0plus_VAJHU / ( p0plus_VAJHU + bkg_VAMCFM ); // D^kin_bkg
       f_D_bkg = p0plus_VAJHU * p0plus_m4l / ( p0plus_VAJHU * p0plus_m4l + bkg_VAMCFM * bkg_m4l ); // D^kin including superMELA
       f_D_gg = Dgg10_VAMCFM; // D_gg
       f_D_g4 = p0plus_VAJHU / ( p0plus_VAJHU + p0minus_VAJHU ); // D_0-
       if (njets_pass>=2) f_Djet_VAJHU = pvbf_VAJHU / ( pvbf_VAJHU + phjj_VAJHU ); // D^VBF_HJJ
       else f_Djet_VAJHU=-1;
       
       // TUtil::computeAngles(partP[0], partId[0], partP[1], partId[1],
       // partP[2], partId[2], partP[3], partId[3],
       //angle_costhetastar,angle_costheta1,angle_costheta2,angle_Phi,angle_Phi1);
      
       TUtil::computeAngles(angle_costhetastar,angle_costheta1,angle_costheta2,angle_Phi,angle_Phi1,
			    partP[0], partId[0], partP[1], partId[1], partP[2], partId[2], partP[3], partId[3]);
       
       
       
       TLorentzVector p4l = L11P4 + L12P4 + L21P4 + L22P4;
       pt4l = p4l.Pt();
       eta4l = p4l.Eta();
       float w; 
       //= mMEM->getMELAWeight();
       
              
       f_int_weight = w;
       f_Z1mass = massofZ1;
       f_Z2mass = massofZ2;
       f_angle_costhetastar = angle_costhetastar;
       f_angle_costheta1 = angle_costheta1;
       f_angle_costheta2 = angle_costheta2;
       f_angle_phi = angle_Phi;
       f_angle_phistar1 = angle_Phi1;
       f_pt4l = pt4l;
       f_eta4l = eta4l;


        //@ Computes the variations of MELA due to systematic uncertainties
       if(njets_pass >= 2){
	 TLorentzVector partP;
	 vector<TLorentzVector> vpartP;
	 vector<vector<TLorentzVector>> vpartsP;
	 vector<int> partsId;
       
	 partsId.push_back(L11PID);
         partsId.push_back(L12PID);
         partsId.push_back(L21PID);
         partsId.push_back(L22PID);
         partsId.push_back(0);
         partsId.push_back(0);
	 
	 vpartP.clear();
	 partP.SetPtEtaPhiE(L11P4.Pt(),L11P4.Eta(),L11P4.Phi(),L11P4.E());						vpartP.push_back( partP );//Nominal
	 partP.SetPtEtaPhiE(L11P4.Pt()+RECOELE_PTError[indexleptonfinal[0]],L11P4.Eta(),L11P4.Phi(),L11P4.E());		vpartP.push_back( partP );//+ Pt Uncertainty
	 partP.SetPtEtaPhiE(L11P4.Pt()-RECOELE_PTError[indexleptonfinal[0]],L11P4.Eta(),L11P4.Phi(),L11P4.E());		vpartP.push_back( partP );//- Pt Uncertainty
	 vpartsP.push_back( vpartP );

	 vpartP.clear();
	 partP.SetPtEtaPhiE(L12P4.Pt(),L12P4.Eta(),L12P4.Phi(),L12P4.E());						vpartP.push_back( partP );//Nominal
	 partP.SetPtEtaPhiE(L12P4.Pt()+RECOELE_PTError[indexleptonfinal[1]],L12P4.Eta(),L12P4.Phi(),L12P4.E());		vpartP.push_back( partP );//+Pt Uncertainty
	 partP.SetPtEtaPhiE(L12P4.Pt()-RECOELE_PTError[indexleptonfinal[1]],L12P4.Eta(),L12P4.Phi(),L12P4.E());		vpartP.push_back( partP );//-Pt Uncertainty
	 vpartsP.push_back( vpartP );

	 vpartP.clear();
	 partP.SetPtEtaPhiE(L21P4.Pt(),L21P4.Eta(),L21P4.Phi(),L21P4.E());						vpartP.push_back( partP );//Nominal
	 partP.SetPtEtaPhiE(L21P4.Pt()+RECOELE_PTError[indexleptonfinal[2]],L21P4.Eta(),L21P4.Phi(),L21P4.E());		vpartP.push_back( partP );//+Pt Uncertainty
	 partP.SetPtEtaPhiE(L21P4.Pt()-RECOELE_PTError[indexleptonfinal[2]],L21P4.Eta(),L21P4.Phi(),L21P4.E());		vpartP.push_back( partP );//-Pt Uncertainty
	 vpartsP.push_back( vpartP );
	 
	 vpartP.clear();
	 partP.SetPtEtaPhiE(L22P4.Pt(),L22P4.Eta(),L22P4.Phi(),L22P4.E());						vpartP.push_back( partP );//Nominal
	 partP.SetPtEtaPhiE(L22P4.Pt()+RECOELE_PTError[indexleptonfinal[3]],L22P4.Eta(),L22P4.Phi(),L22P4.E());		vpartP.push_back( partP );//+Pt Uncertainty
	 partP.SetPtEtaPhiE(L22P4.Pt()-RECOELE_PTError[indexleptonfinal[3]],L22P4.Eta(),L22P4.Phi(),L22P4.E());		vpartP.push_back( partP );//-Pt Uncertainty
	 vpartsP.push_back( vpartP );
	 
	 vpartP.clear();
	 TLorentzVector Jet1, Jet2;
	 Jet1.SetPtEtaPhiE(f_jets_highpt_pt[0],f_jets_highpt_eta[0],f_jets_highpt_phi[0],f_jets_highpt_et[0]*TMath::CosH(f_jets_highpt_eta[0]));
	 Jet2.SetPtEtaPhiE(f_jets_highpt_pt[1],f_jets_highpt_eta[1],f_jets_highpt_phi[1],f_jets_highpt_et[1]*TMath::CosH(f_jets_highpt_eta[1]));
	 partP.SetPtEtaPhiE(Jet1.Pt(),Jet1.Eta(),Jet1.Phi(),Jet1.E());								vpartP.push_back( partP );//Nominal
	 partP.SetPtEtaPhiE(Jet1.Pt()+f_jets_highpt_pt_error[0],Jet1.Eta(),Jet1.Phi(),Jet1.E());				vpartP.push_back( partP );//+Pt Uncertainty
	 partP.SetPtEtaPhiE(Jet1.Pt()-f_jets_highpt_pt_error[0],Jet1.Eta(),Jet1.Phi(),Jet1.E());				vpartP.push_back( partP );//-Pt Uncertainty
	 vpartsP.push_back( vpartP );
	 
	 vpartP.clear();
	 partP.SetPtEtaPhiE(Jet2.Pt(),Jet2.Eta(),Jet2.Phi(),Jet2.E());								vpartP.push_back( partP );//Nominal
	 partP.SetPtEtaPhiE(Jet2.Pt()+f_jets_highpt_pt_error[1],Jet2.Eta(),Jet2.Phi(),Jet2.E());				vpartP.push_back( partP );//+Pt Uncertainty
	 partP.SetPtEtaPhiE(Jet2.Pt()-f_jets_highpt_pt_error[1],Jet2.Eta(),Jet2.Phi(),Jet2.E());				vpartP.push_back( partP );//-Pt Uncertainty
	 vpartsP.push_back( vpartP );
       
	 float Djet_VAJHU_max = f_Djet_VAJHU, Djet_VAJHU_min = f_Djet_VAJHU;
         for(int il1=0; il1<3; ++il1){
	   for(int il2=0; il2<3; ++il2){
	     for(int il3=0; il3<3; ++il3){
	       for(int il4=0; il4<3; ++il4){
		 for(int ij1=0; ij1<3; ++ij1){
		   for(int ij2=0; ij2<3; ++ij2){
		     vector<TLorentzVector> partsP = {vpartsP[0][il1],vpartsP[1][il2],vpartsP[2][il3],vpartsP[3][il4],vpartsP[4][ij1],vpartsP[5][ij2]};
		     int f=combinedMEM.computeME(MEMNames::kJJ_SMHiggs_GG, MEMNames::kJHUGen, partsP, partsId, phjj_VAJHU); // SM gg->H+2j
		     int g=combinedMEM.computeME(MEMNames::kJJ_SMHiggs_VBF, MEMNames::kJHUGen, partsP, partsId, pvbf_VAJHU);  // SM VBF->H	   
		     float Djet_VAJHU = pvbf_VAJHU / ( pvbf_VAJHU + phjj_VAJHU ); // D^VBF_HJJ
		     if(Djet_VAJHU > Djet_VAJHU_max) Djet_VAJHU_max = Djet_VAJHU;
		     if(Djet_VAJHU < Djet_VAJHU_min) Djet_VAJHU_min = Djet_VAJHU;
		   }
		 }
	       }
	     }
	   }
         }
         f_Djet_VAJHU_UncUp = Djet_VAJHU_max;
	 f_Djet_VAJHU_UncDn = Djet_VAJHU_min;
	 //To assure that the f_Djet_VAJHU is computed with the 2 highest pt jets
         vector<TLorentzVector> partsP = {vpartsP[0][0],vpartsP[1][0],vpartsP[2][0],vpartsP[3][0],vpartsP[4][0],vpartsP[5][0]};
         int f=combinedMEM.computeME(MEMNames::kJJ_SMHiggs_GG, MEMNames::kJHUGen, partsP, partsId, phjj_VAJHU); // SM gg->H+2j
         int g=combinedMEM.computeME(MEMNames::kJJ_SMHiggs_VBF, MEMNames::kJHUGen, partsP, partsId, pvbf_VAJHU);  // SM VBF->H	   
         f_Djet_VAJHU = pvbf_VAJHU / ( pvbf_VAJHU + phjj_VAJHU ); // D^VBF_HJJ	 
       }else{
	 f_Djet_VAJHU_UncUp =-1;
	 f_Djet_VAJHU_UncDn =-1;
       }

       newtree->Fill();              
       
       // Filling BNN input       
       bnn_file << newweight << " " << w << " ";
       printelebnn(indexleptonfinal[0]);
       printelebnn(indexleptonfinal[1]);
       printelebnn(indexleptonfinal[2]);
       printelebnn(indexleptonfinal[3]);
       //bnn_file << Iso_max << " " << Sip_max << " " << massofZ1 << " " << massofZ2 << " " << angle_costhetastar << " " << angle_costheta1 << " " << angle_costheta2 << " " << angle_Phi << " " << angle_Phi1<< " " << D << " " << psD << " " << gravD << " " << pt4l << " " << massofhiggs << " " << RECO_PFMET << " ";
       
       // cout << "Njets=" << njets_pass << " VD= " << VD << endl;
       // if (njets_pass>=2) {
       // 	 bnn_file << RECO_PFJET_PT[jet1] << " " << RECO_PFJET_ETA[jet1] << " " << RECO_PFJET_PHI[jet1] << " " << RECO_PFJET_ET[jet1] << " ";
       // 	 bnn_file << RECO_PFJET_PT[jet2] << " " << RECO_PFJET_ETA[jet2] << " " << RECO_PFJET_PHI[jet2] << " " << RECO_PFJET_ET[jet2] << " " 
       // 		  << deltajj << " " << massjj << " " << VD  << " " << njets_pass << endl;
       // }
       // else {
       // 	 bnn_file << RECO_PFJET_PT[jet1] << " " << RECO_PFJET_ETA[jet1] << " " << RECO_PFJET_PHI[jet1] << " " << RECO_PFJET_ET[jet1] << " ";
       // 	 bnn_file << "0. 0. 0. 0. 0. 0. 0. " << njets_pass << endl;
       // }
       
       if(massofhiggs>100. && massofhiggs<180.){
	 hMELA_8->Fill(f_D_bkg_kin,newweight);
	 hMELA_vs_M4l_8->Fill(massofhiggs,f_D_bkg_kin,newweight);
       }
       if (massofhiggs>100. && f_D_bkg_kin > 0.1){       
	 hMELA_9->Fill(f_D_bkg_kin,newweight);
	 hMELA_vs_M4l_9->Fill(massofhiggs,f_D_bkg_kin,newweight);
	 ++N_9;
	 N_9_w=N_9_w+newweight;
       }
       if (f_D_g4 > 0.3){
	 ++N_9PS;
	 N_9PS_w=N_9PS_w+newweight;
       }
     }
     else {
       cout<<endl;
       continue;
     }
     
     // end of KD
     
     // Step 8a m4l window cut
     if ( abs(mass4l-125.)<=10.){ 
       ++N_8_a;
       N_8_a_w=N_8_a_w+newweight;     
       hPFMET_8a->Fill(RECO_PFMET,newweight);
       hLogLinXPFMET_8a->Fill(RECO_PFMET,newweight);
       hM4l_T_8a->Fill(m4l_T,newweight);  
       hDPHI_8a->Fill(fabs(DPHI),newweight);
       hNgood_8a->Fill(f_Ngood,newweight);
       hNbjets_8a->Fill(f_Nbjets,newweight);
     }

     // Step 8b m4l window cut
     if ( abs(mass4l-Zmass)<=10.){ 
       ++N_8_b;
       N_8_b_w=N_8_b_w+newweight;     
       hPFMET_8b->Fill(RECO_PFMET,newweight);
       hLogLinXPFMET_8b->Fill(RECO_PFMET,newweight);
     }

     // Step 9: signal enriched variables
     if ( abs(mass4l-125.)<=10. && RECO_PFMET>=50.){ 
       ++N_9_PFMET;
       N_9_PFMET_w=N_9_PFMET_w+newweight;     
       hPFMET_9->Fill(RECO_PFMET,newweight);
       hM4l_9->Fill( mass4l,newweight );
       hM4l_T_9->Fill(m4l_T,newweight);  
       hDPHI_9->Fill(fabs(DPHI),newweight);
       hNgood_9->Fill(f_Ngood,newweight);
       hNbjets_9->Fill(f_Nbjets,newweight);
     }

     // **** Step 10:
     // |mass4l-125| < 10 + Ngood=4 + Nbjets<=1 .

     if (abs(mass4l-125.)<=10. && Ne_good==4 && n_bjets<=1) {
       ++N_10;
       N_10_w=N_10_w+newweight;     
       hPFMET_10->Fill(RECO_PFMET,newweight);
       hLogLinXPFMET_10->Fill(RECO_PFMET,newweight);
       hM4l_T_10->Fill(m4l_T,newweight);
       hLogLinXM4l_T_10->Fill(m4l_T,newweight);
       hDPHI_10->Fill(fabs(DPHI),newweight);
       hNgood_10->Fill(f_Ngood,newweight);
       hNbjets_10->Fill(f_Nbjets,newweight);

       // shifts up and down
       hPFMET_JetEnUp_10->Fill(RECO_PFMET_JetEnUp,newweight);
       hPFMET_JetEnDn_10->Fill(RECO_PFMET_JetEnDn,newweight);
       hPFMET_ElectronEnUp_10->Fill(RECO_PFMET_ElectronEnUp,newweight);
       hPFMET_ElectronEnDn_10->Fill(RECO_PFMET_ElectronEnDn,newweight);
       hPFMET_MuonEnUp_10->Fill(RECO_PFMET_MuonEnUp,newweight);
       hPFMET_MuonEnDn_10->Fill(RECO_PFMET_MuonEnDn,newweight);
       hPFMET_JetResUp_10->Fill(RECO_PFMET_JetResUp,newweight);
       hPFMET_JetResDn_10->Fill(RECO_PFMET_JetResDn,newweight);
       hPFMET_UnclusteredEnUp_10->Fill(RECO_PFMET_UnclusteredEnUp,newweight);
       hPFMET_UnclusteredEnDn_10->Fill(RECO_PFMET_UnclusteredEnDn,newweight);
       hPFMET_PhotonEnUp_10->Fill(RECO_PFMET_PhotonEnUp,newweight);
       hPFMET_PhotonEnDn_10->Fill(RECO_PFMET_PhotonEnDn,newweight);

       hLogLinXPFMET_JetEnUp_10->Fill(RECO_PFMET_JetEnUp,newweight);
       hLogLinXPFMET_JetEnDn_10->Fill(RECO_PFMET_JetEnDn,newweight);
       hLogLinXPFMET_ElectronEnUp_10->Fill(RECO_PFMET_ElectronEnUp,newweight);
       hLogLinXPFMET_ElectronEnDn_10->Fill(RECO_PFMET_ElectronEnDn,newweight);
       hLogLinXPFMET_MuonEnUp_10->Fill(RECO_PFMET_MuonEnUp,newweight);
       hLogLinXPFMET_MuonEnDn_10->Fill(RECO_PFMET_MuonEnDn,newweight);
       hLogLinXPFMET_JetResUp_10->Fill(RECO_PFMET_JetResUp,newweight);
       hLogLinXPFMET_JetResDn_10->Fill(RECO_PFMET_JetResDn,newweight);
       hLogLinXPFMET_UnclusteredEnUp_10->Fill(RECO_PFMET_UnclusteredEnUp,newweight);
       hLogLinXPFMET_UnclusteredEnDn_10->Fill(RECO_PFMET_UnclusteredEnDn,newweight);
       hLogLinXPFMET_PhotonEnUp_10->Fill(RECO_PFMET_PhotonEnUp,newweight);
       hLogLinXPFMET_PhotonEnDn_10->Fill(RECO_PFMET_PhotonEnDn,newweight);

       for(int i = 0; i < 4; ++i){
	 
	 if( i == 0){
	   hPtLep1_10->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;	   
	 }
	 if( i == 1){
	   hPtLep2_10->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;	  
	 }
	 if( i == 2){
	   hPtLep3_10->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;	   
	 }
	 if( i == 3){
	   hPtLep4_10->Fill( RECOELE_PT[ ipt[i] ],newweight ) ;	   
	 }
       
       }//end fill leptons
     
     }

     if (abs(mass4l-125.)<=10. && Ne_good==4 && n_bjets<=1 && RECO_PFMET>=50.) {
       ++N_11;
       N_11_w=N_11_w+newweight;  
     }  
     
     // reverte the mass peak
     if (abs(mass4l-125.)>10. && Ne_good==4 && n_bjets<=1){
       hPFMET_12->Fill(RECO_PFMET,newweight);
       hLogLinXPFMET_12->Fill(RECO_PFMET,newweight);
     }


     //if( debug )
     cout << "EVENT CANDIDATE: \n" 
	  << " N " << jentry 
	  << " RUN " << Run
	  << " EVENT " << Event	
	  << " LumiSection " << LumiSection
	  << " massZ1 " << massZ1
	  << " massZ2 " << massZ2
	  << " mass4l " << mass4l 
	  << " massError " << massErr
          << " massErrorCorr " << massErrCorr
       //	  << " Iso_max " << Iso_max 
       //  << " Sip_max " << Sip_max
	  << " MELA LD " << f_D_bkg_kin;
     if( has_FSR_Z1 || has_FSR_Z2 ) cout << " FSR " << endl;	
     else cout << " noFSR " << endl;

     output_txt << "EVENT CANDIDATE: \n" 
		<< " N " << jentry 
		<< " RUN " << Run
		<< " EVENT " << Event	
		<< " LumiSection " << LumiSection
		<< " massZ1 " << massZ1
		<< " massZ2 " << massZ2
		<< " mass4l " << mass4l 
		<< " massError " << massErr
		<< " massErrorCorr " << massErrCorr
       //	<< " Iso_max " << Iso_max 
       //	<< " Sip_max " << Sip_max 
		<< " MELA LD " << f_D_bkg_kin;
     //<<endl;
     if( has_FSR_Z1 || has_FSR_Z2 ) output_txt << " FSR " << endl;	
     else output_txt << " noFSR " << endl;
     if (debug) for(unsigned int i=0; i< 4; i++) output_txt << "\t mu[" << i << "] pT: " << RECOELE_PT[indexleptonfinal[i]] << "\t eta: " <<  RECOELE_ETA[indexleptonfinal[i]] << "\t phi: " <<  RECOELE_PHI[indexleptonfinal[i]] << endl;

     // Format
     //{run}:{lumi}:{event}:{mass4l:.2f}:{mZ1:.2f}:{mZ2:.2f}:{massErrRaw:.2f}:{massErrCorr:.2f}:{KD:.3f}:{pt4l:.2f}:{njets30:.0}:{jet1pt:.2f}:{jet2pt:.2f}:{mjj:.2f}:{detajj:.3f}:{fishjj:.3f}:{KD(0-):.3f}:{KD(0+h):.3f}:{KD(1+):.3f}:{KD(1-):.3f}:{KD(gg2+m):.3f}:{KD(qq2+m):.3f}
     //{run}:{lumi}:{event}:{mass4l:.2f}:{mZ1:.2f}:{mZ2:.2f}::{D_bkg^kin:.3f}:{D_bkg:.3f}:{D_gg:.3f}:{D_HJJ^VBF:.3f}:{D_0-:.3f}:{njets30:d}:{jet1pt:.2f}:{jet2pt:.2f}:{category} 

     Char_t outformat[20000];
     float dummy=-1.;
     
     // if (JET1.Pt()> && jet2<0) 
     //   sprintf (outformat,"Format=%d:%d:%d:%.2f:%.2f:%.2f:%.3f:%.3f:%.3f:%.3f:%.3f:%d:%.2f:%.2f:%d",
     // 		Run,LumiSection,Event,massofhiggs,massZ1,massZ2,f_D_bkg_kin,f_D_bkg,f_D_gg,f_Djet_VAJHU,f_D_g4,njets_pass,dummy,dummy,category);
     // else if (jet1>=0 && jet2<0)
     //   sprintf (outformat,"Format=%d:%d:%d:%.2f:%.2f:%.2f:%.3f:%.3f:%.3f:%.3f:%.3f:%d:%.2f:%.2f:%d",
     //            Run,LumiSection,Event,massofhiggs,massZ1,massZ2,f_D_bkg_kin,f_D_bkg,f_D_gg,f_Djet_VAJHU,f_D_g4,njets_pass,JET1.Pt(),dummy,category);
     
     if (JET1.Pt()>0. && JET2.Pt()>0.)
       sprintf (outformat,"Format=%d:%d:%d:%.2f:%.2f:%.2f:%.3f:%.3f:%.3f:%.3f:%.3f:%d:%.2f:%.2f:%d",
		Run,LumiSection,Event,massofhiggs,massZ1,massZ2,f_D_bkg_kin,f_D_bkg,f_D_gg,f_Djet_VAJHU,f_D_g4,njets_pass,JET1.Pt(),JET2.Pt(),category);
     else if (JET1.Pt()>0.)
       sprintf (outformat,"Format=%d:%d:%d:%.2f:%.2f:%.2f:%.3f:%.3f:%.3f:%.3f:%.3f:%d:%.2f:%.2f:%d",
                Run,LumiSection,Event,massofhiggs,massZ1,massZ2,f_D_bkg_kin,f_D_bkg,f_D_gg,f_Djet_VAJHU,f_D_g4,njets_pass,JET1.Pt(),dummy,category);
     else sprintf (outformat,"Format=%d:%d:%d:%.2f:%.2f:%.2f:%.3f:%.3f:%.3f:%.3f:%.3f:%d:%.2f:%.2f:%d",
		Run,LumiSection,Event,massofhiggs,massZ1,massZ2,f_D_bkg_kin,f_D_bkg,f_D_gg,f_Djet_VAJHU,f_D_g4,njets_pass,dummy,dummy,category);
     
     output_txt  << outformat << endl;

	
     // fill FSR counters
     if( has_FSR_Z1 || has_FSR_Z2 ){
       ++N_9_1FSR; // fill the counter
       N_9_1FSR_w=N_9_1FSR_w+newweight;
     }
     if( has_FSR_Z1 && has_FSR_Z2 ){
       ++N_9_2FSR; // fill the counter
       N_9_2FSR_w=N_9_2FSR_w+newweight;
     }
     
     // fill final tree
     finaltree->Fill();

     

   } // end loop on entries

   // write on output txt file:

   
   output_txt << "N_0 "  << N_0  << " \n" 
	      << "N_01 " << N_01 << " \n"	
	      << "N_02 " << N_02 << " \n"	
	      << "N_1 "  << N_1  << " \n"	
	      << "N_2 "  << N_2  << " \n"	
	      << "N_3a " << N_3a << " \n"	
	      << "N_3_FSR " << N_3_FSR << " \n"	
	      << "N_3b " << N_3b << " \n"	
	      << "N_4a " << N_4a << " \n"	
	      << "N_4b " << N_4b << " \n"	
	      << "N_4c " << N_4c << " \n"	
	      << "N_4d " << N_4d << " \n"	
	      << "N_5 "  << N_5  << " \n"	
	      << "N_6 "  << N_6  << " \n"	
	      << "N_7 "  << N_7  << " \n"	
	      << "N_8 "  << N_8  << " \n"
	      << "N_8_a "<< N_8_a<< " \n"
      	      << "N_9 "  << N_9  << " \n"
      	      << "N_9_1FSR " << N_9_1FSR  << " \n" 
	      << "N_9_2FSR " << N_9_2FSR  << " \n" 
	      << "N_9PS "    << N_9PS << " \n"
              << "N_9GRAV"   << N_9GRAV << "\n"
	      << "N_9a_VBF " << N_9a_VBF << " \n"
              << "N_9b_VBF " << N_9b_VBF << "\n"
	      << "N_VBF "    << N_VBF << " \n"
     	      << "N_9_PFMET "<< N_9_PFMET << " \n"
	      << "N_10 " << N_10 << " \n";
   
   output_txt_vbf
              << "N_0 "  << N_0  << " \n" 
	      << "N_01 " << N_01 << " \n"	
	      << "N_02 " << N_02 << " \n"	
	      << "N_1 "  << N_1  << " \n"	
	      << "N_2 "  << N_2  << " \n"	
	      << "N_3a " << N_3a << " \n"	
	      << "N_3_FSR " << N_3_FSR << " \n"	
	      << "N_3b " << N_3b << " \n"	
	      << "N_4a " << N_4a << " \n"	
	      << "N_4b " << N_4b << " \n"	
	      << "N_4c " << N_4c << " \n"	
	      << "N_4d " << N_4d << " \n"	
	      << "N_5 "  << N_5  << " \n"	
	      << "N_6 "  << N_6  << " \n"	
	      << "N_7 "  << N_7  << " \n"	
	      << "N_8 "  << N_8  << " \n"
	      << "N_8_a "<< N_8_a<< " \n"
      	      << "N_9 "  << N_9  << " \n"
      	      << "N_9_1FSR " << N_9_1FSR  << " \n" 
	      << "N_9_2FSR " << N_9_2FSR  << " \n" 
	      << "N_9PS "    << N_9PS << " \n"
              << "N_9GRAV"   << N_9GRAV << "\n"
	      << "N_9a_VBF " << N_9a_VBF << " \n"
              << "N_9b_VBF " << N_9b_VBF << "\n"
	      << "N_VBF "    << N_VBF << " \n"
     	      << "N_9_PFMET "<< N_9_PFMET << " \n"
	      << "N_10 " << N_10 << " \n";

   bnn_file.close();
   output_txt.close();
   output_txt_vbf.close();
   


   cout << "N_0 "  << N_0  << " \n" 
	      << "N_01 " << N_01 << " \n"	
	      << "N_02 " << N_02 << " \n"	
	      << "N_1 "  << N_1  << " \n"	
	      << "N_2 "  << N_2  << " \n"	
	      << "N_3a " << N_3a << " \n"	
	      << "N_3_FSR " << N_3_FSR << " \n"	
	      << "N_3b " << N_3b << " \n"	
	      << "N_4a " << N_4a << " \n"	
	      << "N_4b " << N_4b << " \n"	
	      << "N_4c " << N_4c << " \n"	
	      << "N_4d " << N_4d << " \n"	
	      << "N_5 "  << N_5  << " \n"	
	      << "N_6 "  << N_6  << " \n"	
	      << "N_7 "  << N_7  << " \n"	
	      << "N_8 "  << N_8  << " \n"
	      << "N_8_a "<< N_8_a<< " \n"
      	      << "N_9 "  << N_9  << " \n"
      	      << "N_9_1FSR "  << N_9_1FSR  << " \n" 
	      << "N_9_2FSR "  << N_9_2FSR  << " \n" 
	      << "N_9PS "     << N_9PS << " \n"
	      << "N_9GRAV"    << N_9GRAV << "\n"
	      << "N_9a_VBF "  << N_9a_VBF << " \n"
              << "N_9b_VBF "  << N_9b_VBF << "\n"
	      << "N_VBF "     << N_VBF << " \n"
	      << "N_9_PFMET " << N_9_PFMET << " \n"
	      << "N_10 " << N_10 << " \n";

   nEvent_4l->GetXaxis()->SetBinLabel(1,"Init.");
   nEvent_4l->GetXaxis()->SetBinLabel(2,"MCTruth: 4e");
   nEvent_4l->GetXaxis()->SetBinLabel(3,"MCTruth: Acc");
   nEvent_4l->GetXaxis()->SetBinLabel(4,"Init");
   nEvent_4l->GetXaxis()->SetBinLabel(5,"HLT");
   nEvent_4l->GetXaxis()->SetBinLabel(6,"Z1 lept. cuts");
   nEvent_4l->GetXaxis()->SetBinLabel(7,"Z1+#gamma");
   nEvent_4l->GetXaxis()->SetBinLabel(8,"m_{Z1}");
   nEvent_4l->GetXaxis()->SetBinLabel(9,"4e");
   nEvent_4l->GetXaxis()->SetBinLabel(10,"at least one Z2");
   nEvent_4l->GetXaxis()->SetBinLabel(11,"Z2 lept. cuts");
   nEvent_4l->GetXaxis()->SetBinLabel(12,"m_{Z2}");
   nEvent_4l->GetXaxis()->SetBinLabel(13,"pT cuts");
   nEvent_4l->GetXaxis()->SetBinLabel(14,"mll>4 for OS-SF");
   nEvent_4l->GetXaxis()->SetBinLabel(15,"m4l > 70");
   nEvent_4l->GetXaxis()->SetBinLabel(16,"m_{Z2} > 12");
   nEvent_4l->GetXaxis()->SetBinLabel(17,"m4l > 100");
   nEvent_4l->GetXaxis()->SetBinLabel(18,"MELA KD > 0.1");
   nEvent_4l->GetXaxis()->SetBinLabel(19,"one Z+#gamma");
   nEvent_4l->GetXaxis()->SetBinLabel(20,"two Z+#gamma");
   nEvent_4l->GetXaxis()->SetBinLabel(21,"MET>100");
   nEvent_4l->GetXaxis()->SetBinLabel(22,"$|m_{4l}-125.|<=10$, $N_{good \\, l}=4$, $n_{b-jets}<=1$");
   nEvent_4l->GetXaxis()->SetBinLabel(23,"$|m_{4l}-125.|<=10$, $N_{good \\, l}=4$, $n_{b-jets}<=1$, $MET>50$");

   nEvent_4l_w->GetXaxis()->SetBinLabel(1,"Init."); 
   nEvent_4l_w->GetXaxis()->SetBinLabel(2,"MCTruth: 4e");  
   nEvent_4l_w->GetXaxis()->SetBinLabel(3,"MCTruth: Acc");
   nEvent_4l_w->GetXaxis()->SetBinLabel(4,"Init");
   nEvent_4l_w->GetXaxis()->SetBinLabel(5,"HLT");
   nEvent_4l_w->GetXaxis()->SetBinLabel(6,"Z1 lept. cuts");
   nEvent_4l_w->GetXaxis()->SetBinLabel(7,"Z1+#gamma");
   nEvent_4l_w->GetXaxis()->SetBinLabel(8,"m_{Z1}");
   nEvent_4l_w->GetXaxis()->SetBinLabel(9,"4e");
   nEvent_4l_w->GetXaxis()->SetBinLabel(10,"at least one Z2");
   nEvent_4l_w->GetXaxis()->SetBinLabel(11,"Z2 lept. cuts");
   nEvent_4l_w->GetXaxis()->SetBinLabel(12,"m_{Z2}");
   nEvent_4l_w->GetXaxis()->SetBinLabel(13,"pT cuts");
   nEvent_4l_w->GetXaxis()->SetBinLabel(14,"mll>4 for OS-SF");
   nEvent_4l_w->GetXaxis()->SetBinLabel(15,"m4l > 70");
   nEvent_4l_w->GetXaxis()->SetBinLabel(16,"m_{Z2} > 12");
   nEvent_4l_w->GetXaxis()->SetBinLabel(17,"m4l > 100");
   nEvent_4l_w->GetXaxis()->SetBinLabel(18,"MELA KD > 0.1");
   nEvent_4l_w->GetXaxis()->SetBinLabel(19,"one Z+#gamma");
   nEvent_4l_w->GetXaxis()->SetBinLabel(20,"two Z+#gamma");
   nEvent_4l_w->GetXaxis()->SetBinLabel(21,"MET>100");
   nEvent_4l_w->GetXaxis()->SetBinLabel(22,"$|m_{4l}-125.|<=10$, $N_{good \\, l}=4$, $n_{b-jets}<=1$");
   nEvent_4l_w->GetXaxis()->SetBinLabel(23,"$|m_{4l}-125.|<=10$, $N_{good \\, l}=4$, $n_{b-jets}<=1$, $MET>50$");

   nEvent_4l->SetBinContent(1,N_0);
   nEvent_4l->SetBinContent(2,N_01);
   nEvent_4l->SetBinContent(3,N_02);
   nEvent_4l->SetBinContent(4,N_1);
   nEvent_4l->SetBinContent(5,N_2);
   nEvent_4l->SetBinContent(6,N_3a);
   nEvent_4l->SetBinContent(7,N_3_FSR);
   nEvent_4l->SetBinContent(8,N_3b);
   nEvent_4l->SetBinContent(9,N_4a);
   nEvent_4l->SetBinContent(10,N_4b);
   nEvent_4l->SetBinContent(11,N_4c);
   nEvent_4l->SetBinContent(12,N_4d);
   nEvent_4l->SetBinContent(13,N_5);
   nEvent_4l->SetBinContent(14,N_6);
   nEvent_4l->SetBinContent(15,N_7);
   nEvent_4l->SetBinContent(16,N_8);
   nEvent_4l->SetBinContent(17,N_8_a);
   nEvent_4l->SetBinContent(18,N_9);
   nEvent_4l->SetBinContent(19,N_9_1FSR);
   nEvent_4l->SetBinContent(20,N_9_2FSR);
   nEvent_4l->SetBinContent(21,N_9_PFMET);
   nEvent_4l->SetBinContent(22,N_10);
   nEvent_4l->SetBinContent(23,N_11);

   nEvent_4l_w->SetBinContent(1,N_0_w);        nEvent_4l_w->SetBinError(1,sqrt(N_0)*N_0_w/N_0);
   nEvent_4l_w->SetBinContent(2,N_01_w);       nEvent_4l_w->SetBinError(2,sqrt(N_01)*N_01_w/N_01);
   nEvent_4l_w->SetBinContent(3,N_02_w);       nEvent_4l_w->SetBinError(3,sqrt(N_02)*N_02_w/N_02);
   nEvent_4l_w->SetBinContent(4,N_1_w);        nEvent_4l_w->SetBinError(4,sqrt(N_1)*N_1_w/N_1);
   nEvent_4l_w->SetBinContent(5,N_2_w);        nEvent_4l_w->SetBinError(5,sqrt(N_2)*N_2_w/N_2);
   nEvent_4l_w->SetBinContent(6,N_3a_w);       nEvent_4l_w->SetBinError(6,sqrt(N_3a)*N_3a_w/N_3a);
   nEvent_4l_w->SetBinContent(7,N_3_FSR_w);    nEvent_4l_w->SetBinError(7,sqrt(N_3_FSR)*N_3_FSR_w/N_3_FSR) ;
   nEvent_4l_w->SetBinContent(8,N_3b_w);       nEvent_4l_w->SetBinError(8,sqrt(N_3b)*N_3b_w/N_3b);
   nEvent_4l_w->SetBinContent(9,N_4a_w);       nEvent_4l_w->SetBinError(9,sqrt(N_4a)*N_4a_w/N_4a);
   nEvent_4l_w->SetBinContent(10,N_4b_w);      nEvent_4l_w->SetBinError(10,sqrt(N_4b)*N_4b_w/N_4b);
   nEvent_4l_w->SetBinContent(11,N_4c_w);      nEvent_4l_w->SetBinError(11,sqrt(N_4c)*N_4c_w/N_4c);
   nEvent_4l_w->SetBinContent(12,N_4d_w);      nEvent_4l_w->SetBinError(12,sqrt(N_4d)*N_4d_w/N_4d);
   nEvent_4l_w->SetBinContent(13,N_5_w);       nEvent_4l_w->SetBinError(13,sqrt(N_5)*N_5_w/N_5);
   nEvent_4l_w->SetBinContent(14,N_6_w);       nEvent_4l_w->SetBinError(14,sqrt(N_6)*N_6_w/N_6);
   nEvent_4l_w->SetBinContent(15,N_7_w);       nEvent_4l_w->SetBinError(15,sqrt(N_7)*N_7_w/N_7);
   nEvent_4l_w->SetBinContent(16,N_8_w);       nEvent_4l_w->SetBinError(16,sqrt(N_8)*N_8_w/N_8);
   nEvent_4l_w->SetBinContent(17,N_8_a_w);     nEvent_4l_w->SetBinError(17,sqrt(N_8_a)*N_8_a_w/N_8_a);
   nEvent_4l_w->SetBinContent(18,N_9_w);       nEvent_4l_w->SetBinError(18,sqrt(N_9)*N_9_w/N_9);
   nEvent_4l_w->SetBinContent(19,N_9_1FSR_w);  nEvent_4l_w->SetBinError(19,sqrt(N_9_1FSR)*N_9_1FSR_w/N_9_1FSR);
   nEvent_4l_w->SetBinContent(20,N_9_2FSR_w);  nEvent_4l_w->SetBinError(20,sqrt(N_9_2FSR)*N_9_2FSR_w/N_9_2FSR);
   nEvent_4l_w->SetBinContent(21,N_9_PFMET_w); nEvent_4l_w->SetBinError(21,sqrt(N_9_PFMET)*N_9_PFMET_w/N_9_PFMET);
   nEvent_4l_w->SetBinContent(22,N_10_w);      nEvent_4l_w->SetBinError(22,sqrt(N_10)*N_10_w/N_10);
   nEvent_4l_w->SetBinContent(23,N_11_w);      nEvent_4l_w->SetBinError(22,sqrt(N_11)*N_11_w/N_11);

   nEvent_red->GetXaxis()->SetBinLabel(1,"MCTruth: 4e in Geom. Acc");
   nEvent_red->GetXaxis()->SetBinLabel(2,"One Z with lept. cuts");
   nEvent_red->GetXaxis()->SetBinLabel(3,"Two Z with 12 < m_{Z} < 120");
   nEvent_red->GetXaxis()->SetBinLabel(4,"One ZZ (no lept. overlap, #Delta R > 0.02)");
   nEvent_red->GetXaxis()->SetBinLabel(5,"One ZZ with lept. pT > 20/10");
   nEvent_red->GetXaxis()->SetBinLabel(6,"Z_{1}(m_{Z1}>40) and Z_{2}");
   nEvent_red->GetXaxis()->SetBinLabel(7,"m4l > 100");

   nEvent_red->SetBinContent(1,N_02_w);
   nEvent_red->SetBinContent(2,N_3a_w);
   nEvent_red->SetBinContent(3,N_3b_w);
   nEvent_red->SetBinContent(4,N_4b_w);
   nEvent_red->SetBinContent(5,N_4c_w);
   nEvent_red->SetBinContent(6,N_4d_w);
   nEvent_red->SetBinContent(7,N_8_w);

   nEvent_ZZ->GetXaxis()->SetBinLabel(1,"At two Z with 12 < m_{Z} < 120");
   nEvent_ZZ->GetXaxis()->SetBinLabel(2,"One ZZ (no lept. overlap, #Delta R > 0.02)");
   nEvent_ZZ->GetXaxis()->SetBinLabel(3,"One ZZ with lept. pT > 20/10");
   nEvent_ZZ->GetXaxis()->SetBinLabel(4,"Z_{1}(m_{Z1}>40) and Z_{2}");
   nEvent_ZZ->GetXaxis()->SetBinLabel(5,"m4l > 100");

   nEvent_ZZ->SetBinContent(1,N_3b_w);
   nEvent_ZZ->SetBinContent(2,N_4b_w);
   nEvent_ZZ->SetBinContent(3,N_4c_w);
   nEvent_ZZ->SetBinContent(4,N_4d_w);
   nEvent_ZZ->SetBinContent(5,N_8_w);

   // write on output root file:
   _filePU->Close();
   theFile->cd();
   //z1tree->Write();
   newtree->Write();
   theFile->Write();
   theFile->Close();
   skimfile->cd();
   finaltree->Write();
   skimfile->Write();
   skimfile->Close();
} // end main


double DELTAPHI( double phi1, double phi2 ){

	if( phi1 > mPI || phi1 < -mPI || phi2 > mPI || phi2 < -mPI) {
	  // cout << "Angles out of range!!! " << endl;
	  // cout << " phi1 " << phi1 << endl;
	  // cout << " phi2 " << phi2 << endl;
	  return -999;
	}
	float dp=std::abs(phi1-phi2);
	if (dp>mPI) dp-=float(2*mPI);
	return dp;
	//return  min( fabs( phi1 - phi2 ) , 2*mPI - fabs( phi1 - phi2 ) ) ;

}

double invmass (float M1, float PT1, float ETA1, float PHI1, float M2, float PT2, float ETA2, float PHI2 ){ 
 float phi1=PHI1; 
 float eta1=ETA1; 
 float pt1=PT1; 
 float m1=M1; 

 float px1=pt1*cos(phi1); 
 float py1=pt1*sin(phi1); 
 float pz1=pt1/(2.*(exp(-1*eta1))/(1.0-exp(-2.*eta1))); 

 float phi2=PHI2; 
 float eta2=ETA2; 
 float pt2=PT2; 
 float m2=M2;

 float px2=pt2*cos(phi2); 
 float py2=pt2*sin(phi2); 
 float pz2=pt2/(2.*(exp(-1*eta2))/(1.0-exp(-2.*eta2))); 

 float e1sqr=pz1*pz1+pt1*pt1+m1*m1; 
 float e2sqr=pz2*pz2+pt2*pt2+m2*m2; 
 float e1e2=sqrt(e1sqr*e2sqr); 
 float p1dotp2=px1*px2+py1*py2+pz1*pz2; 

 float m=sqrt(m1*m1+m2*m2+2.*(e1e2-p1dotp2)); 
 //cout << "Invariant mass= " << m << endl;  
 return m;
} // float invmass closed

double HZZ4LeptonsAnalysis::masserror( std::vector<TLorentzVector> Lep, std::vector<double> pterr){ 

  //	if(Lep.size()!=4 or pterr.size()!=4) {std::cout<<" Lepsize="<<Lep.size()<<", "<<pterr.size()<<std::endl;}
  int debug_ = 0;
  TLorentzVector compositeParticle ;
  for(unsigned int i=0; i<Lep.size(); i++){
    compositeParticle+=Lep[i];
    if(debug_) std::cout<<" in mass error :  add lep  "<<i<<endl;
  }
  double mass  =  compositeParticle.M();
  
  if(debug_) std::cout<<" in mass error :  mass "<<mass<<endl;
  double masserr = 0;
  
  for(unsigned int i=0; i<Lep.size(); i++){
    if(debug_) std::cout<<" in mass error :  varying lep "<<i<<endl;
    TLorentzVector variedLep; // = Lep[i];
    
    if(debug_) std::cout<<" in mass error : pterr = "<<pterr[i]<<endl;
    variedLep.SetPtEtaPhiM(Lep[i].Pt()+ pterr[i], Lep[i].Eta(), Lep[i].Phi(), Lep[i].M());
    TLorentzVector compositeParticleVariation ;
    for(unsigned int j=0; j<Lep.size(); j++){
      if(i!=j)compositeParticleVariation+=Lep[j];
      else compositeParticleVariation+=variedLep;
    }
    
    masserr += (compositeParticleVariation.M()-mass)*(compositeParticleVariation.M()-mass);
    if(debug_) std::cout<<" in mass error :  intermediate masserr "<<masserr<<endl;
  }
  return sqrt(masserr);
}

double HZZ4LeptonsAnalysis::EAele(int index,bool use2011EA){
  
  double EffectiveArea=0.;
  if (use2011EA){
    if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.18;
    if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.20;
    if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.15;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.19;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.21;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.22;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.4 )                                       EffectiveArea = 0.29;
  }
  //else { // 7_4_x use eta
  // if (fabs(RECOELE_ETA[index]) >= 0.0   && fabs(RECOELE_ETA[index]) < 0.8 )   EffectiveArea = 0.1830;
  // if (fabs(RECOELE_ETA[index]) >= 0.8   && fabs(RECOELE_ETA[index]) < 1.3 )   EffectiveArea = 0.1734;
  // if (fabs(RECOELE_ETA[index]) >= 1.3   && fabs(RECOELE_ETA[index]) < 2.0 )   EffectiveArea = 0.1077;
  // if (fabs(RECOELE_ETA[index]) >= 2.0   && fabs(RECOELE_ETA[index]) < 2.2 )   EffectiveArea = 0.1565;
  // if (fabs(RECOELE_ETA[index]) >= 2.2 )                                       EffectiveArea = 0.2680;
  //}
  else { // 7_6_X use eta supercluster   
    if (fabs(RECOELE_scl_Eta[index]) >= 0.0   && fabs(RECOELE_scl_Eta[index]) < 1.0 )   EffectiveArea = 0.1752;
    if (fabs(RECOELE_scl_Eta[index]) >= 1.0   && fabs(RECOELE_scl_Eta[index]) < 1.479 ) EffectiveArea = 0.1862;
    if (fabs(RECOELE_scl_Eta[index]) >= 1.479 && fabs(RECOELE_scl_Eta[index]) < 2.0 )   EffectiveArea = 0.1411;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.0   && fabs(RECOELE_scl_Eta[index]) < 2.2 )   EffectiveArea = 0.1534;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.2   && fabs(RECOELE_scl_Eta[index]) < 2.3 )   EffectiveArea = 0.1903;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.3   && fabs(RECOELE_scl_Eta[index]) < 2.4 )   EffectiveArea = 0.2243;
    if (fabs(RECOELE_scl_Eta[index]) >= 2.4   && fabs(RECOELE_scl_Eta[index]) < 5.0  )  EffectiveArea = 0.2687;
  }
  
  return EffectiveArea;

}



float HZZ4LeptonsAnalysis::kfactor_qqZZ_qcd_dPhi(float GENdPhiZZ, int finalState)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {        
        k+=1.515838921760*(abs(GENdPhiZZ)>0.0&&abs(GENdPhiZZ)<=0.1);
        k+=1.496256665410*(abs(GENdPhiZZ)>0.1&&abs(GENdPhiZZ)<=0.2);
        k+=1.495522061910*(abs(GENdPhiZZ)>0.2&&abs(GENdPhiZZ)<=0.3);
        k+=1.483273154250*(abs(GENdPhiZZ)>0.3&&abs(GENdPhiZZ)<=0.4);
        k+=1.465589701130*(abs(GENdPhiZZ)>0.4&&abs(GENdPhiZZ)<=0.5);
        k+=1.491500887510*(abs(GENdPhiZZ)>0.5&&abs(GENdPhiZZ)<=0.6);
        k+=1.441183580450*(abs(GENdPhiZZ)>0.6&&abs(GENdPhiZZ)<=0.7);
        k+=1.440830603990*(abs(GENdPhiZZ)>0.7&&abs(GENdPhiZZ)<=0.8);
        k+=1.414339019120*(abs(GENdPhiZZ)>0.8&&abs(GENdPhiZZ)<=0.9);
        k+=1.422534218560*(abs(GENdPhiZZ)>0.9&&abs(GENdPhiZZ)<=1.0);
        k+=1.401037066000*(abs(GENdPhiZZ)>1.0&&abs(GENdPhiZZ)<=1.1);
        k+=1.408539428810*(abs(GENdPhiZZ)>1.1&&abs(GENdPhiZZ)<=1.2);
        k+=1.381247744080*(abs(GENdPhiZZ)>1.2&&abs(GENdPhiZZ)<=1.3);
        k+=1.370553357430*(abs(GENdPhiZZ)>1.3&&abs(GENdPhiZZ)<=1.4);
        k+=1.347323316000*(abs(GENdPhiZZ)>1.4&&abs(GENdPhiZZ)<=1.5);
        k+=1.340113437450*(abs(GENdPhiZZ)>1.5&&abs(GENdPhiZZ)<=1.6);
        k+=1.312661036510*(abs(GENdPhiZZ)>1.6&&abs(GENdPhiZZ)<=1.7);
        k+=1.290055062010*(abs(GENdPhiZZ)>1.7&&abs(GENdPhiZZ)<=1.8);
        k+=1.255322614790*(abs(GENdPhiZZ)>1.8&&abs(GENdPhiZZ)<=1.9);
        k+=1.254455642450*(abs(GENdPhiZZ)>1.9&&abs(GENdPhiZZ)<=2.0);
        k+=1.224047664420*(abs(GENdPhiZZ)>2.0&&abs(GENdPhiZZ)<=2.1);
        k+=1.178816782670*(abs(GENdPhiZZ)>2.1&&abs(GENdPhiZZ)<=2.2);
        k+=1.162624827140*(abs(GENdPhiZZ)>2.2&&abs(GENdPhiZZ)<=2.3);
        k+=1.105401140940*(abs(GENdPhiZZ)>2.3&&abs(GENdPhiZZ)<=2.4);
        k+=1.074749265690*(abs(GENdPhiZZ)>2.4&&abs(GENdPhiZZ)<=2.5);
        k+=1.021864599380*(abs(GENdPhiZZ)>2.5&&abs(GENdPhiZZ)<=2.6);
        k+=0.946334793286*(abs(GENdPhiZZ)>2.6&&abs(GENdPhiZZ)<=2.7);
        k+=0.857458082628*(abs(GENdPhiZZ)>2.7&&abs(GENdPhiZZ)<=2.8);
        k+=0.716607670482*(abs(GENdPhiZZ)>2.8&&abs(GENdPhiZZ)<=2.9);
        k+=1.132841784840*(abs(GENdPhiZZ)>2.9&&abs(GENdPhiZZ)<=3.1416);
    }

    if (finalState==2) {
       k+=1.513834489150*(abs(GENdPhiZZ)>0.0&&abs(GENdPhiZZ)<=0.1);
       k+=1.541738780180*(abs(GENdPhiZZ)>0.1&&abs(GENdPhiZZ)<=0.2);
       k+=1.497829632510*(abs(GENdPhiZZ)>0.2&&abs(GENdPhiZZ)<=0.3);
       k+=1.534956782920*(abs(GENdPhiZZ)>0.3&&abs(GENdPhiZZ)<=0.4);
       k+=1.478217033060*(abs(GENdPhiZZ)>0.4&&abs(GENdPhiZZ)<=0.5);
       k+=1.504330859290*(abs(GENdPhiZZ)>0.5&&abs(GENdPhiZZ)<=0.6);
       k+=1.520626246850*(abs(GENdPhiZZ)>0.6&&abs(GENdPhiZZ)<=0.7);
       k+=1.507013090030*(abs(GENdPhiZZ)>0.7&&abs(GENdPhiZZ)<=0.8);
       k+=1.494243156250*(abs(GENdPhiZZ)>0.8&&abs(GENdPhiZZ)<=0.9);
       k+=1.450536096150*(abs(GENdPhiZZ)>0.9&&abs(GENdPhiZZ)<=1.0);
       k+=1.460812521660*(abs(GENdPhiZZ)>1.0&&abs(GENdPhiZZ)<=1.1);
       k+=1.471603622200*(abs(GENdPhiZZ)>1.1&&abs(GENdPhiZZ)<=1.2);
       k+=1.467700038200*(abs(GENdPhiZZ)>1.2&&abs(GENdPhiZZ)<=1.3);
       k+=1.422408690640*(abs(GENdPhiZZ)>1.3&&abs(GENdPhiZZ)<=1.4);
       k+=1.397184022730*(abs(GENdPhiZZ)>1.4&&abs(GENdPhiZZ)<=1.5);
       k+=1.375593447520*(abs(GENdPhiZZ)>1.5&&abs(GENdPhiZZ)<=1.6);
       k+=1.391901318370*(abs(GENdPhiZZ)>1.6&&abs(GENdPhiZZ)<=1.7);
       k+=1.368564350560*(abs(GENdPhiZZ)>1.7&&abs(GENdPhiZZ)<=1.8);
       k+=1.317884804290*(abs(GENdPhiZZ)>1.8&&abs(GENdPhiZZ)<=1.9);
       k+=1.314019950800*(abs(GENdPhiZZ)>1.9&&abs(GENdPhiZZ)<=2.0);
       k+=1.274641749910*(abs(GENdPhiZZ)>2.0&&abs(GENdPhiZZ)<=2.1);
       k+=1.242346606820*(abs(GENdPhiZZ)>2.1&&abs(GENdPhiZZ)<=2.2);
       k+=1.244727403840*(abs(GENdPhiZZ)>2.2&&abs(GENdPhiZZ)<=2.3);
       k+=1.146259351670*(abs(GENdPhiZZ)>2.3&&abs(GENdPhiZZ)<=2.4);
       k+=1.107804993520*(abs(GENdPhiZZ)>2.4&&abs(GENdPhiZZ)<=2.5);
       k+=1.042053646740*(abs(GENdPhiZZ)>2.5&&abs(GENdPhiZZ)<=2.6);
       k+=0.973608545141*(abs(GENdPhiZZ)>2.6&&abs(GENdPhiZZ)<=2.7);
       k+=0.872169942668*(abs(GENdPhiZZ)>2.7&&abs(GENdPhiZZ)<=2.8);
       k+=0.734505279177*(abs(GENdPhiZZ)>2.8&&abs(GENdPhiZZ)<=2.9);
       k+=1.163152837230*(abs(GENdPhiZZ)>2.9&&abs(GENdPhiZZ)<=3.1416);       
    }
    if (k==0.0) return 1.1; // if something goes wrong return inclusive k-factor
    else return k;

}

float HZZ4LeptonsAnalysis::kfactor_qqZZ_qcd_M(float GENmassZZ, int finalState)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {
        k+=1.23613311013*(abs(GENmassZZ)>0.0&&abs(GENmassZZ)<=25.0);
        k+=1.17550314639*(abs(GENmassZZ)>25.0&&abs(GENmassZZ)<=50.0);
        k+=1.17044565911*(abs(GENmassZZ)>50.0&&abs(GENmassZZ)<=75.0);
        k+=1.03141209689*(abs(GENmassZZ)>75.0&&abs(GENmassZZ)<=100.0);
        k+=1.05285574912*(abs(GENmassZZ)>100.0&&abs(GENmassZZ)<=125.0);
        k+=1.11287217794*(abs(GENmassZZ)>125.0&&abs(GENmassZZ)<=150.0);
        k+=1.13361441158*(abs(GENmassZZ)>150.0&&abs(GENmassZZ)<=175.0);
        k+=1.10355603327*(abs(GENmassZZ)>175.0&&abs(GENmassZZ)<=200.0);
        k+=1.10053981637*(abs(GENmassZZ)>200.0&&abs(GENmassZZ)<=225.0);
        k+=1.10972676811*(abs(GENmassZZ)>225.0&&abs(GENmassZZ)<=250.0);
        k+=1.12069120525*(abs(GENmassZZ)>250.0&&abs(GENmassZZ)<=275.0);
        k+=1.11589101635*(abs(GENmassZZ)>275.0&&abs(GENmassZZ)<=300.0);
        k+=1.13906170314*(abs(GENmassZZ)>300.0&&abs(GENmassZZ)<=325.0);
        k+=1.14854594271*(abs(GENmassZZ)>325.0&&abs(GENmassZZ)<=350.0);
        k+=1.14616229031*(abs(GENmassZZ)>350.0&&abs(GENmassZZ)<=375.0);
        k+=1.14573157789*(abs(GENmassZZ)>375.0&&abs(GENmassZZ)<=400.0);
        k+=1.13829430515*(abs(GENmassZZ)>400.0&&abs(GENmassZZ)<=425.0);
        k+=1.15521193686*(abs(GENmassZZ)>425.0&&abs(GENmassZZ)<=450.0);
        k+=1.13679822698*(abs(GENmassZZ)>450.0&&abs(GENmassZZ)<=475.0);
        k+=1.13223956942*(abs(GENmassZZ)>475.0);
    }

    if (finalState==2) {
        k+=1.25094466582*(abs(GENmassZZ)>0.0&&abs(GENmassZZ)<=25.0);
        k+=1.22459455362*(abs(GENmassZZ)>25.0&&abs(GENmassZZ)<=50.0);
        k+=1.19287368979*(abs(GENmassZZ)>50.0&&abs(GENmassZZ)<=75.0);
        k+=1.04597506451*(abs(GENmassZZ)>75.0&&abs(GENmassZZ)<=100.0);
        k+=1.08323413771*(abs(GENmassZZ)>100.0&&abs(GENmassZZ)<=125.0);
        k+=1.09994968030*(abs(GENmassZZ)>125.0&&abs(GENmassZZ)<=150.0);
        k+=1.16698455800*(abs(GENmassZZ)>150.0&&abs(GENmassZZ)<=175.0);
        k+=1.10399053155*(abs(GENmassZZ)>175.0&&abs(GENmassZZ)<=200.0);
        k+=1.10592664340*(abs(GENmassZZ)>200.0&&abs(GENmassZZ)<=225.0);
        k+=1.10690381480*(abs(GENmassZZ)>225.0&&abs(GENmassZZ)<=250.0);
        k+=1.11194928918*(abs(GENmassZZ)>250.0&&abs(GENmassZZ)<=275.0);
        k+=1.13522586553*(abs(GENmassZZ)>275.0&&abs(GENmassZZ)<=300.0);
        k+=1.11895090244*(abs(GENmassZZ)>300.0&&abs(GENmassZZ)<=325.0);
        k+=1.13898508615*(abs(GENmassZZ)>325.0&&abs(GENmassZZ)<=350.0);
        k+=1.15463977506*(abs(GENmassZZ)>350.0&&abs(GENmassZZ)<=375.0);
        k+=1.17341664594*(abs(GENmassZZ)>375.0&&abs(GENmassZZ)<=400.0);
        k+=1.20093349763*(abs(GENmassZZ)>400.0&&abs(GENmassZZ)<=425.0);
        k+=1.18915554919*(abs(GENmassZZ)>425.0&&abs(GENmassZZ)<=450.0);
        k+=1.18546007375*(abs(GENmassZZ)>450.0&&abs(GENmassZZ)<=475.0);
        k+=1.12864505708*(abs(GENmassZZ)>475.0);
    }

    if (k==0.0) return 1.1;
    else return k; // if something goes wrong return inclusive k-factor

}

float HZZ4LeptonsAnalysis::kfactor_qqZZ_qcd_Pt(float GENpTZZ, int finalState)
{

    // finalState=1 : 4e/4mu/4tau
    // finalState=2 : 2e2mu/2mutau/2e2tau

    float k=0.0;

    if (finalState==1) {
        k+=0.64155491983*(abs(GENpTZZ)>0.0&&abs(GENpTZZ)<=5.0);
        k+=1.09985240531*(abs(GENpTZZ)>5.0&&abs(GENpTZZ)<=10.0);
        k+=1.29390628654*(abs(GENpTZZ)>10.0&&abs(GENpTZZ)<=15.0);
        k+=1.37859998571*(abs(GENpTZZ)>15.0&&abs(GENpTZZ)<=20.0);
        k+=1.42430263312*(abs(GENpTZZ)>20.0&&abs(GENpTZZ)<=25.0);
        k+=1.45038493266*(abs(GENpTZZ)>25.0&&abs(GENpTZZ)<=30.0);
        k+=1.47015377651*(abs(GENpTZZ)>30.0&&abs(GENpTZZ)<=35.0);
        k+=1.48828685748*(abs(GENpTZZ)>35.0&&abs(GENpTZZ)<=40.0);
        k+=1.50573440448*(abs(GENpTZZ)>40.0&&abs(GENpTZZ)<=45.0);
        k+=1.50211655928*(abs(GENpTZZ)>45.0&&abs(GENpTZZ)<=50.0);
        k+=1.50918720827*(abs(GENpTZZ)>50.0&&abs(GENpTZZ)<=55.0);
        k+=1.52463089491*(abs(GENpTZZ)>55.0&&abs(GENpTZZ)<=60.0);
        k+=1.52400838378*(abs(GENpTZZ)>60.0&&abs(GENpTZZ)<=65.0);
        k+=1.52418067701*(abs(GENpTZZ)>65.0&&abs(GENpTZZ)<=70.0);
        k+=1.55424382578*(abs(GENpTZZ)>70.0&&abs(GENpTZZ)<=75.0);
        k+=1.52544284222*(abs(GENpTZZ)>75.0&&abs(GENpTZZ)<=80.0);
        k+=1.57896384602*(abs(GENpTZZ)>80.0&&abs(GENpTZZ)<=85.0);
        k+=1.53034682567*(abs(GENpTZZ)>85.0&&abs(GENpTZZ)<=90.0);
        k+=1.56147329708*(abs(GENpTZZ)>90.0&&abs(GENpTZZ)<=95.0);
        k+=1.54468169268*(abs(GENpTZZ)>95.0&&abs(GENpTZZ)<=100.0);
        k+=1.57222952415*(abs(GENpTZZ)>100.0);
    }

    if (finalState==2) {
        k+=0.743602533303*(abs(GENpTZZ)>0.0&&abs(GENpTZZ)<=5.0);
        k+=1.14789453219*(abs(GENpTZZ)>5.0&&abs(GENpTZZ)<=10.0);
        k+=1.33815867892*(abs(GENpTZZ)>10.0&&abs(GENpTZZ)<=15.0);
        k+=1.41420044104*(abs(GENpTZZ)>15.0&&abs(GENpTZZ)<=20.0);
        k+=1.45511318916*(abs(GENpTZZ)>20.0&&abs(GENpTZZ)<=25.0);
        k+=1.47569225244*(abs(GENpTZZ)>25.0&&abs(GENpTZZ)<=30.0);
        k+=1.49053003693*(abs(GENpTZZ)>30.0&&abs(GENpTZZ)<=35.0);
        k+=1.50622827695*(abs(GENpTZZ)>35.0&&abs(GENpTZZ)<=40.0);
        k+=1.50328889799*(abs(GENpTZZ)>40.0&&abs(GENpTZZ)<=45.0);
        k+=1.52186945281*(abs(GENpTZZ)>45.0&&abs(GENpTZZ)<=50.0);
        k+=1.52043468754*(abs(GENpTZZ)>50.0&&abs(GENpTZZ)<=55.0);
        k+=1.53977869986*(abs(GENpTZZ)>55.0&&abs(GENpTZZ)<=60.0);
        k+=1.53491994434*(abs(GENpTZZ)>60.0&&abs(GENpTZZ)<=65.0);
        k+=1.51772882172*(abs(GENpTZZ)>65.0&&abs(GENpTZZ)<=70.0);
        k+=1.54494489131*(abs(GENpTZZ)>70.0&&abs(GENpTZZ)<=75.0);
        k+=1.57762411697*(abs(GENpTZZ)>75.0&&abs(GENpTZZ)<=80.0);
        k+=1.55078339014*(abs(GENpTZZ)>80.0&&abs(GENpTZZ)<=85.0);
        k+=1.57078191891*(abs(GENpTZZ)>85.0&&abs(GENpTZZ)<=90.0);
        k+=1.56162666568*(abs(GENpTZZ)>90.0&&abs(GENpTZZ)<=95.0);
        k+=1.54183774627*(abs(GENpTZZ)>95.0&&abs(GENpTZZ)<=100.0);
        k+=1.58485762205*(abs(GENpTZZ)>100.0);
    }

    if (k==0.0) return 1.1;
    else return k; // if something goes wrong return inclusive k-factor

}


