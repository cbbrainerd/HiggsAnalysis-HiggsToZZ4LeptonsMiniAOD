//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun  7 12:14:04 2012 by ROOT version 5.32/00
// from TTree HZZ4LeptonsAnalysis/HZZ4Leptons Analysis Tree
// found on file: roottree_leptons_Fall11_0706.root
//////////////////////////////////////////////////////////

#ifndef HZZ4LeptonsAnalysis_h
#define HZZ4LeptonsAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <TLorentzVector.h>

using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class HZZ4LeptonsAnalysis {
public :

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          Run;
   UInt_t          Event;
   UInt_t          LumiSection;
   Float_t         Avginstlumi;
   Int_t           num_PU_vertices;
   Int_t           PU_BunchCrossing;
   Float_t         MC_weighting;
   Int_t           RECO_nMuHLTMatch;
   Float_t         RECOMU_PT_MuHLTMatch[100];  
   Char_t          HLTPathsFired[20000];
   Float_t         MC_E[7];
   Float_t         MC_PT[7];
   Float_t         MC_ETA[7];
   Float_t         MC_THETA[7];
   Float_t         MC_PHI[7];
   Float_t         MC_MASS[7];
   Float_t         MC_PDGID[7];
   Float_t         MC_LEPT_PT[4];
   Float_t         MC_LEPT_ETA[4];
   Float_t         MC_LEPT_PHI[4];
   Float_t         MC_LEPT_THETA[4];
   Float_t         MC_LEPT_PDGID[4];
   Float_t         MC_Z_PT[2][5];
   Float_t         MC_Z_ETA[2][5];
   Float_t         MC_Z_PHI[2][5];
   Float_t         MC_Z_THETA[2][5];
   Float_t         MC_Z_MASS[2][5];
   Float_t         MC_Z_PDGID[2][5];
   Float_t         MC_fourl_MASS[50][5];
   Float_t         MC_fourl_PT[50][5];
   Float_t         MC_fourl_PDGID[50][5];
   Float_t         MC_ZZ_MASS[4][7];
   Float_t         MC_ZZ_PT[4][7];
   Float_t         MC_ZZ_ETA[4][7];
   Float_t         MC_ZZ_PHI[4][7];
   Float_t         MC_ZZ_THETA[4][7];
   Float_t         MC_ZZ_PDGID[4][7];
   Float_t         MC_GENJET_PT[4];
   Float_t         MC_GENJET_ETA[4];
   Float_t         MC_GENJET_PHI[4];
   Float_t         MC_GENMET;
   Double_t        RECORF_2e2mu_cosTheta1_spin[100];
   Double_t        RECORF_2e2mu_cosTheta2_spin[100];
   Double_t        RECORF_2e2mu_cosThetaStar_spin[100];
   Double_t        RECORF_2e2mu_Phi_spin[100];
   Double_t        RECORF_2e2mu_Phi1_spin[100];
   Double_t        RECORF_2e2mu_Phi2_spin[100];
   Double_t        RECORF_2e2mu_phi1RF_spin[100];
   Double_t        RECORF_2e2mu_phi2RF_spin[100];
   Double_t        RECORF_2e2mu_MELA[100];
   Double_t        RECORF_4e_cosTheta1_spin[100];
   Double_t        RECORF_4e_cosTheta2_spin[100];
   Double_t        RECORF_4e_cosThetaStar_spin[100];
   Double_t        RECORF_4e_Phi_spin[100];
   Double_t        RECORF_4e_Phi1_spin[100];
   Double_t        RECORF_4e_Phi2_spin[100];
   Double_t        RECORF_4e_phi1RF_spin[100];
   Double_t        RECORF_4e_phi2RF_spin[100];
   Double_t        RECORF_4e_MELA[100];
   Double_t        RECORF_4mu_cosTheta1_spin[100];
   Double_t        RECORF_4mu_cosTheta2_spin[100];
   Double_t        RECORF_4mu_cosThetaStar_spin[100];
   Double_t        RECORF_4mu_Phi_spin[100];
   Double_t        RECORF_4mu_Phi1_spin[100];
   Double_t        RECORF_4mu_Phi2_spin[100];
   Double_t        RECORF_4mu_phi1RF_spin[100];
   Double_t        RECORF_4mu_phi2RF_spin[100];
   Double_t        RECORF_4mu_MELA[100];
   Float_t         RECO_ZMM_MASS[50];
   Float_t         RECO_ZEE_MASS[50];
   Float_t         RECO_DiLep_MASS[50];
   Float_t         RECO_ZMM_PT[3][50];
   Float_t         RECO_ZEE_PT[3][50];
   Float_t         RECO_DiLep_PT[3][50];
   Float_t         RECO_ZMM_ETA[3][50];
   Float_t         RECO_ZEE_ETA[3][50];
   Float_t         RECO_DiLep_ETA[3][50];
   Float_t         RECO_ZMM_PHI[3][50];
   Float_t         RECO_ZEE_PHI[3][50];
   Float_t         RECO_DiLep_PHI[3][50];
   Float_t         RECO_ZMMss_MASS[50];
   Float_t         RECO_ZEEss_MASS[50];
   Float_t         RECO_ZEM_MASS[50];
   Float_t         RECO_ZMMss_PT[3][50];
   Float_t         RECO_ZEEss_PT[3][50];
   Float_t         RECO_ZEM_PT[3][50];
   Float_t         RECO_ZMMss_ETA[3][50];
   Float_t         RECO_ZEEss_ETA[3][50];
   Float_t         RECO_ZEM_ETA[3][50];
   Float_t         RECO_ZMMss_PHI[3][50];
   Float_t         RECO_ZEEss_PHI[3][50];
   Float_t         RECO_ZEM_PHI[3][50];
   Float_t         RECO_MMMM_MASS[7][100];
   Float_t         RECO_MMMM_PT[7][100];
   Float_t         RECO_MMMM_ETA[7][100];
   Float_t         RECO_MMMM_PHI[7][100];
   Float_t         RECO_MMMM_MASS_REFIT[100];
   Float_t         RECO_EEEE_MASS[7][100];
   Float_t         RECO_EEEE_PT[7][100];
   Float_t         RECO_EEEE_ETA[7][100];
   Float_t         RECO_EEEE_PHI[7][100];
   Float_t         RECO_EEEE_MASS_REFIT[100];
   Float_t         RECO_EEMM_MASS[7][100];
   Float_t         RECO_EEMM_PT[7][100];
   Float_t         RECO_EEMM_ETA[7][100];
   Float_t         RECO_EEMM_PHI[7][100];
   Float_t         RECO_EEMM_MASS_REFIT[100];
   Float_t         RECO_LLL0_MASS[50];
   Float_t         RECO_LLL1_MASS[50];
   Float_t         RECO_LLL2_MASS[50];
   Float_t         RECO_LLL3_MASS[50];
   Float_t         RECO_LLL0_PT[4][50];
   Float_t         RECO_LLL1_PT[4][50];
   Float_t         RECO_LLL2_PT[4][50];
   Float_t         RECO_LLL3_PT[4][50];
   Float_t         RECO_LLLl0_MASS[20];
   Float_t         RECO_LLLl1_MASS[20];
   Float_t         RECO_LLLl0_PT[5][20];
   Float_t         RECO_LLLl1_PT[5][20];
   Float_t         RECO_LLLL0ss_MASS[20];
   Float_t         RECO_LLLL0ss_PT[5][20];
   Float_t         RECO_LLLL1ss_MASS[20];
   Float_t         RECO_LLLL1ss_PT[5][20];
   Float_t         RECO_LLLL2ss_MASS[20];
   Float_t         RECO_LLLL2ss_PT[5][20];
   Float_t         RECO_LLLL_MASS[7][100];
   Float_t         RECO_LLLL_PT[7][100];
   Float_t         RECO_LLLL_ETA[7][100];
   Float_t         RECO_LLLL_PHI[7][100];
   Float_t         RECOELE_E[100];
   Float_t         RECOELE_PT[100];
   Float_t         RECOELE_PTError[100];
   Float_t         RECOELE_P[100];
   Float_t         RECOELE_ETA[100];
   Float_t         RECOELE_THETA[100];
   Float_t         RECOELE_PHI[100];
   Float_t         RECOELE_MASS[100];
   Float_t         RECOELE_CHARGE[100];
   UChar_t         RECOELE_isEcalDriven[100];
   UChar_t         RECOELE_isTrackerDriven[100];
   Float_t         RECOELE_gsftrack_NPixHits[100];
   Float_t         RECOELE_gsftrack_NStripHits[100];
   Float_t         RECOELE_gsftrack_chi2[100];
   Float_t         RECOELE_gsftrack_dxyB[100];
   Float_t         RECOELE_gsftrack_dxy[100];
   Float_t         RECOELE_gsftrack_dxyError[100];
   Float_t         RECOELE_gsftrack_dzB[100];
   Float_t         RECOELE_gsftrack_dz[100];
   Float_t         RECOELE_gsftrack_dzError[100];
   Int_t           RECOELE_gsftrack_losthits[100];
   Int_t           RECOELE_gsftrack_validhits[100];
   Int_t           RECOELE_gsftrack_expected_inner_hits[100];
   Float_t         RECOELE_scl_E[100];
   Float_t         RECOELE_scl_Et[100];
   Float_t         RECOELE_scl_Eta[100];
   Float_t         RECOELE_scl_Phi[100];
   Float_t         RECOELE_ep[100];
   Float_t         RECOELE_eSeedp[100];
   Float_t         RECOELE_eSeedpout[100];
   Float_t         RECOELE_eElepout[100];
   Float_t         RECOELE_deltaEtaIn[100];
   Float_t         RECOELE_deltaEtaSeed[100];
   Float_t         RECOELE_deltaEtaEle[100];
   Float_t         RECOELE_deltaPhiIn[100];
   Float_t         RECOELE_deltaPhiSeed[100];
   Float_t         RECOELE_deltaPhiEle[100];
   Int_t           RECOELE_isbarrel[100];
   Int_t           RECOELE_isendcap[100];
   Int_t           RECOELE_isEBetaGap[100];
   Int_t           RECOELE_isEBphiGap[100];
   Int_t           RECOELE_isEEdeeGap[100];
   Int_t           RECOELE_isEEringGap[100];
   Float_t         RECOELE_sigmaIetaIeta[100];
   Float_t         RECOELE_sigmaEtaEta[100];
   Float_t         RECOELE_e15[100];
   Float_t         RECOELE_e25max[100];
   Float_t         RECOELE_e55[100];
   Float_t         RECOELE_he[100];
   Float_t         RECOELE_r9[100];
   Float_t         RECOELE_mva[100];
   Float_t         RECOELE_fbrem[100];
   Int_t           RECOELE_nbrems[100];
   Int_t           RECOELE_Class[100];
   Double_t        RECOELE_fbrem_mode[100];
   Double_t        RECOELE_fbrem_mean[100];
   Float_t         RECOELE_EGMTRACKISO[100];
   Float_t         RECOELE_EGMHCALISO[100];
   Float_t         RECOELE_EGMECALISO[100];
   Float_t         RECOELE_EGMX[100];
   Double_t        RECOELE_PFchAllPart[100];
   Double_t        RECOELE_PFchHad[100];
   Double_t        RECOELE_PFneuHad[100];
   Double_t        RECOELE_PFphoton[100];
   Double_t        RECOELE_PFPUchAllPart[100];
   Double_t        RECOELE_PFX_dB[100];
   Double_t        RECOELE_PFX_rho[100];
   Double_t        RECOELE_PFX_rho_new[100];
   Double_t        RECOELE_regEnergy[100];
   Double_t        RECOELE_regEnergyError[100];
   Float_t         RECOELE_SIP[100];
   Float_t         RECOELE_IP[100];
   Float_t         RECOELE_IPERROR[100];
   Float_t         RECOELE_SIP_KF[100];
   Float_t         RECOELE_IP_KF[100];
   Float_t         RECOELE_IPERROR_KF[100];
   Float_t         RECOELE_SIP_GD[100];
   Float_t         RECOELE_SIP_GDEEEE[100];
   Float_t         RECOELE_SIP_Std[100];
   Float_t         RECOELE_SIP_StdEEEE[100];
   Float_t         RECOELE_SIP_Kin[100];
   Float_t         RECOELE_SIP_KinEEEE[100];
   Float_t         RECOELE_STIP[100];
   Float_t         RECOELE_SLIP[100];
   Float_t         RECOELE_TIP[100];
   Float_t         RECOELE_LIP[100];
   Float_t         RECOELE_TIPERROR[100];
   Float_t         RECOELE_LIPERROR[100];
   Double_t        RECOELE_sclRawE[100];
   Double_t        RECOELE_sclX[100];
   Double_t        RECOELE_sclY[100];
   Double_t        RECOELE_sclZ[100];
   Double_t        RECOELE_seedSubdet1[100];
   Double_t        RECOELE_seedDphi1[100];
   Double_t        RECOELE_seedDrz1[100];
   Double_t        RECOELE_seedSubdet2[100];
   Double_t        RECOELE_seedDphi2[100];
   Double_t        RECOELE_seedDrz2[100];
   Double_t        RECOELE_eidVeryLoose[100];
   Double_t        RECOELE_eidLoose[100];
   Double_t        RECOELE_eidMedium[100];
   Double_t        RECOELE_eidTight[100];
   Double_t        RECOELE_eidHZZVeryLoose[100];
   Double_t        RECOELE_eidHZZLoose[100];
   Double_t        RECOELE_eidHZZMedium[100];
   Double_t        RECOELE_eidHZZTight[100];
   Double_t        RECOELE_mvaTrigV0[100];
   Double_t        RECOELE_mvaNonTrigV0[100];
   Float_t         RECOELE_ID[100];
   Double_t        RECOELE_COV[100][3][3];
   UChar_t         RECOMU_isPFMu[100];
   UChar_t         RECOMU_isGlobalMu[100];
   UChar_t         RECOMU_isStandAloneMu[100];
   UChar_t         RECOMU_isTrackerMu[100];
   UChar_t         RECOMU_isCaloMu[100];
   UChar_t         RECOMU_isTrackerHighPtMu[100];
   Float_t         RECOMU_E[100];
   Float_t         RECOMU_PT[100];
   Float_t         RECOMU_P[100];
   Float_t         RECOMU_ETA[100];
   Float_t         RECOMU_THETA[100];
   Float_t         RECOMU_PHI[100];
   Float_t         RECOMU_MASS[100];
   Float_t         RECOMU_CHARGE[100];
   Double_t        RECOMU_COV[100][3][3];
   Float_t         RECOMU_TRACKISO[100];
   Float_t         RECOMU_TRACKISO_SUMPT[100];
   Float_t         RECOMU_HCALISO[100];
   Float_t         RECOMU_ECALISO[100];
   Float_t         RECOMU_X[100];
   Double_t        RECOMU_PFchHad[100];
   Double_t        RECOMU_PFneuHad[100];
   Double_t        RECOMU_PFphoton[100];
   Double_t        RECOMU_PFPUchAllPart[100];
   Double_t        RECOMU_PFX_dB[100];
   Double_t        RECOMU_PFX_dB_new[100];
   Double_t        RECOMU_PFX_rho[100];
   Double_t        RECOPFPHOT_PFchHad[20];
   Double_t        RECOPFPHOT_PFneuHad[20];
   Double_t        RECOPFPHOT_PFphoton[20];
   Double_t        RECOPFPHOT_PFPUchAllPart[20];
   Double_t        RECOPFPHOT_PFX_rho[20];
   Float_t         RECOMU_SIP[100];
   Float_t         RECOMU_IP[100];
   Float_t         RECOMU_IPERROR[100];
   Float_t         RECOMU_SIP_KF[100];
   Float_t         RECOMU_IP_KF[100];
   Float_t         RECOMU_IPERROR_KF[100];
   Float_t         RECOMU_SIP_GD[100];
   Float_t         RECOMU_SIP_GDMMMM[100];
   Float_t         RECOMU_SIP_Std[100];
   Float_t         RECOMU_SIP_StdMMMM[100];
   Float_t         RECOMU_SIP_Kin[100];
   Float_t         RECOMU_SIP_KinMMMM[100];
   Float_t         RECOMU_STIP[100];
   Float_t         RECOMU_SLIP[100];
   Float_t         RECOMU_TIP[100];
   Float_t         RECOMU_LIP[100];
   Float_t         RECOMU_TIPERROR[100];
   Float_t         RECOMU_LIPERROR[100];
   Float_t         RECOMU_caloCompatibility[100];
   Float_t         RECOMU_segmentCompatibility[100];
   UInt_t          RECOMU_numberOfMatches[100];
   UInt_t          RECOMU_numberOfMatchedStations[100];
   UChar_t         RECOMU_glbmuPromptTight[100];
   UChar_t         RECOMU_trkmuArbitration[100];
   UChar_t         RECOMU_trkmu2DCompatibilityLoose[100];
   UChar_t         RECOMU_trkmu2DCompatibilityTight[100];
   UChar_t         RECOMU_trkmuOneStationLoose[100];
   UChar_t         RECOMU_trkmuOneStationTight[100];
   UChar_t         RECOMU_trkmuLastStationLoose[100];
   UChar_t         RECOMU_trkmuLastStationTight[100];
   UChar_t         RECOMU_trkmuOneStationAngLoose[100];
   UChar_t         RECOMU_trkmuOneStationAngTight[100];
   UChar_t         RECOMU_trkmuLastStationAngLoose[100];
   UChar_t         RECOMU_trkmuLastStationAngTight[100];
   UChar_t         RECOMU_trkmuLastStationOptimizedLowPtLoose[100];
   UChar_t         RECOMU_trkmuLastStationOptimizedLowPtTight[100];
   Float_t         RECOMU_mutrkPT[100];
   Float_t         RECOMU_mutrkPTError[100];
   Float_t         RECOMU_mutrkDxy[100];
   Float_t         RECOMU_mutrkDxyError[100];
   Float_t         RECOMU_mutrkDxyB[100];
   Float_t         RECOMU_mutrkDz[100];
   Float_t         RECOMU_mutrkDzError[100];
   Float_t         RECOMU_mutrkDzB[100];
   Float_t         RECOMU_mutrkChi2PerNdof[100];
   Float_t         RECOMU_mutrkCharge[100];
   Float_t         RECOMU_mutrkNHits[100];
   Float_t         RECOMU_mutrkNStripHits[100];
   Float_t         RECOMU_mutrkNPixHits[100];
   Float_t         RECOMU_mutrkNMuonHits[100];
   Float_t         RECOMU_mutrktrackerLayersWithMeasurement[100];
   Float_t         RECOMU_muInnertrkDxy[100];
   Float_t         RECOMU_muInnertrkDxyError[100];
   Float_t         RECOMU_muInnertrkDxyB[100];
   Float_t         RECOMU_muInnertrkDz[100];
   Float_t         RECOMU_muInnertrkDzError[100];
   Float_t         RECOMU_muInnertrkDzB[100];
   Float_t         RECOMU_muInnertrkChi2PerNdof[100];
   Float_t         RECOMU_muInnertrktrackerLayersWithMeasurement[100];
   Float_t         RECOMU_muInnertrkPT[100];
   Float_t         RECOMU_muInnertrkPTError[100];
   Float_t         RECOMU_muInnertrkCharge[100];
   Float_t         RECOMU_muInnertrkNHits[100];
   Float_t         RECOMU_muInnertrkNStripHits[100];
   Float_t         RECOMU_muInnertrkNPixHits[100];
   Int_t           RECOMU_mubesttrkType[100];
   Float_t         RECOMU_mubesttrkDxy[100];
   Float_t         RECOMU_mubesttrkDxyError[100];
   Float_t         RECOMU_mubesttrkDz[100];
   Float_t         RECOMU_mubesttrkDzError[100];
   Float_t         RECOMU_mubesttrkPTError[100];
   Double_t        ftsigma[100];
   Double_t        gdX[100];
   Double_t        gdY[100];
   Double_t        gdZ[100];
   Double_t        ftsigmalag[100];
   Double_t        gdlagX[100];
   Double_t        gdlagY[100];
   Double_t        gdlagZ[100];
   Double_t        gdlagProb[100];
   Double_t        gdlagNdof[100];
   Double_t        ftsigmaMMMM[100];
   Double_t        gdXMMMM[100];
   Double_t        gdYMMMM[100];
   Double_t        gdZMMMM[100];
   Double_t        ftsigmalagMMMM[100];
   Double_t        gdlagXMMMM[100];
   Double_t        gdlagYMMMM[100];
   Double_t        gdlagZMMMM[100];
   Double_t        gdlagProbMMMM[100];
   Double_t        gdlagNdofMMMM[100];
   Double_t        ftsigmaEEEE[100];
   Double_t        gdXEEEE[100];
   Double_t        gdYEEEE[100];
   Double_t        gdZEEEE[100];
   Double_t        ftsigmalagEEEE[100];
   Double_t        gdlagXEEEE[100];
   Double_t        gdlagYEEEE[100];
   Double_t        gdlagZEEEE[100];
   Double_t        gdlagProbEEEE[100];
   Double_t        gdlagNdofEEEE[100];
   Double_t        StdFitVertexX[100];
   Double_t        StdFitVertexY[100];
   Double_t        StdFitVertexZ[100];
   Double_t        StdFitVertexChi2r[100];
   Double_t        StdFitVertexProb[100];
   Float_t         StdFitVertexTrack_PT[4][100];
   Float_t         StdFitVertexTrack_ETA[4][100];
   Float_t         StdFitVertexTrack_PHI[4][100];
   Double_t        KinFitVertexX[100];
   Double_t        KinFitVertexY[100];
   Double_t        KinFitVertexZ[100];
   Double_t        KinFitVertexChi2r[100];
   Double_t        KinFitVertexProb[100];
   Double_t        StdFitVertexXMMMM[100];
   Double_t        StdFitVertexYMMMM[100];
   Double_t        StdFitVertexZMMMM[100];
   Double_t        StdFitVertexChi2rMMMM[100];
   Double_t        StdFitVertexProbMMMM[100];
   Float_t         StdFitVertexTrackMMMM_PT[4][100];
   Float_t         StdFitVertexTrackMMMM_ETA[4][100];
   Float_t         StdFitVertexTrackMMMM_PHI[4][100];
   Double_t        KinFitVertexXMMMM[100];
   Double_t        KinFitVertexYMMMM[100];
   Double_t        KinFitVertexZMMMM[100];
   Double_t        KinFitVertexChi2rMMMM[100];
   Double_t        KinFitVertexProbMMMM[100];
   Double_t        StdFitVertexXEEEE[100];
   Double_t        StdFitVertexYEEEE[100];
   Double_t        StdFitVertexZEEEE[100];
   Double_t        StdFitVertexChi2rEEEE[100];
   Double_t        StdFitVertexProbEEEE[100];
   Float_t         StdFitVertexTrackEEEE_PT[4][100];
   Float_t         StdFitVertexTrackEEEE_ETA[4][100];
   Float_t         StdFitVertexTrackEEEE_PHI[4][100];
   Double_t        KinFitVertexXEEEE[100];
   Double_t        KinFitVertexYEEEE[100];
   Double_t        KinFitVertexZEEEE[100];
   Double_t        KinFitVertexChi2rEEEE[100];
   Double_t        KinFitVertexProbEEEE[100];
   Double_t        StdFitVertexChi2rMMM[50];
   Double_t        StdFitVertexProbMMM[50];
   Double_t        StdFitVertexChi2rMME[50];
   Double_t        StdFitVertexProbMME[50];
   Double_t        StdFitVertexChi2rEEE[50];
   Double_t        StdFitVertexProbEEE[50];
   Double_t        StdFitVertexChi2rMEE[50];
   Double_t        StdFitVertexProbMEE[50];
   Double_t        StdFitVertexChi2rDiLep[40];
   Double_t        StdFitVertexProbDiLep[40];
   Float_t         ConvMapDist[100];
   Float_t         ConvMapDcot[100];
   UChar_t         RECOMU_MatchingMCTruth[100];
   Float_t         RECOMU_MatchingMCpT[100];
   Float_t         RECOMU_MatchingMCEta[100];
   Float_t         RECOMU_MatchingMCPhi[100];
   UChar_t         RECOELE_MatchingMCTruth[100];
   Float_t         RECOELE_MatchingMCpT[100];
   Float_t         RECOELE_MatchingMCEta[100];
   Float_t         RECOELE_MatchingMCPhi[100];
   UChar_t         RECOPHOT_MatchingMCTruth[50];
   Float_t         RECOPHOT_MatchingMCpT[50];
   Float_t         RECOPHOT_MatchingMCEta[50];
   Float_t         RECOPHOT_MatchingMCPhi[50];
   UChar_t         RECOzMuMu_MatchingMCTruth[50];
   Float_t         RECOzMuMu_MatchingMCpT[50];
   Float_t         RECOzMuMu_MatchingMCmass[50];
   Float_t         RECOzMuMu_MatchingMCEta[50];
   Float_t         RECOzMuMu_MatchingMCPhi[50];
   UChar_t         RECOzEE_MatchingMCTruth[50];
   Float_t         RECOzEE_MatchingMCpT[50];
   Float_t         RECOzEE_MatchingMCmass[50];
   Float_t         RECOzEE_MatchingMCEta[50];
   Float_t         RECOzEE_MatchingMCPhi[50];
   UChar_t         RECOHzzEEEE_MatchingMCTruth[100];
   Float_t         RECOHzzEEEE_MatchingMCpT[100];
   Float_t         RECOHzzEEEE_MatchingMCmass[100];
   Float_t         RECOHzzEEEE_MatchingMCEta[100];
   Float_t         RECOHzzEEEE_MatchingMCPhi[100];
   UChar_t         RECOHzzEEMM_MatchingMCTruth[100];
   Float_t         RECOHzzEEMM_MatchingMCpT[100];
   Float_t         RECOHzzEEMM_MatchingMCmass[100];
   Float_t         RECOHzzEEMM_MatchingMCEta[100];
   Float_t         RECOHzzEEMM_MatchingMCPhi[100];
   UChar_t         RECOHzzMMMM_MatchingMCTruth[100];
   Float_t         RECOHzzMMMM_MatchingMCpT[100];
   Float_t         RECOHzzMMMM_MatchingMCmass[100];
   Float_t         RECOHzzMMMM_MatchingMCEta[100];
   Float_t         RECOHzzMMMM_MatchingMCPhi[100];
   Int_t           RECO_NMU;
   Int_t           RECO_NELE;
   Int_t           RECO_NTRACK;
   Float_t         RECO_TRACK_PT[200];
   Float_t         RECO_TRACK_ETA[200];
   Float_t         RECO_TRACK_PHI[200];
   Float_t         RECO_TRACK_CHI2[200];
   Float_t         RECO_TRACK_CHI2RED[200];
   Float_t         RECO_TRACK_CHI2PROB[200];
   Int_t           RECO_TRACK_NHITS[200];
   Float_t         RECO_TRACK_DXY[200];
   Float_t         RECO_TRACK_DXYERR[200];
   Float_t         RECO_TRACK_DZ[200];
   Float_t         RECO_TRACK_DZERR[200];
   Int_t           RECO_NPHOT;
   Float_t         RECOPHOT_PT[20];
   Float_t         RECOPHOT_ETA[20];
   Float_t         RECOPHOT_PHI[20];
   Float_t         RECOPHOT_THETA[20];
   Int_t           RECO_NPFPHOT;
   Float_t         RECOPFPHOT_PT[20];
   Float_t         RECOPFPHOT_PTError[20];
   Float_t         RECOPFPHOT_ETA[20];
   Float_t         RECOPFPHOT_PHI[20];
   Float_t         RECOPFPHOT_THETA[20];
   Double_t        BeamSpot_X;
   Double_t        BeamSpot_Y;
   Double_t        BeamSpot_Z;
   Int_t           RECO_NVTX;
   Float_t         RECO_VERTEX_x[15];
   Float_t         RECO_VERTEX_y[15];
   Float_t         RECO_VERTEX_z[15];
   Float_t         RECO_VERTEX_ndof[15];
   Float_t         RECO_VERTEX_chi2[15];
   Int_t           RECO_VERTEX_ntracks[15];
   Float_t         RECO_VERTEXPROB[15];
   UChar_t         RECO_VERTEX_isValid[15];
   Float_t         RECO_VERTEX_TRACK_PT[15][100];
   Int_t           RECO_PFJET_N;
   Int_t           RECO_PFJET_CHARGE[100];
   Float_t         RECO_PFJET_ET[100];
   Float_t         RECO_PFJET_PT[100];
   Float_t         RECO_PFJET_ETA[100];
   Float_t         RECO_PFJET_PHI[100];
   Int_t           RECO_PFJET_PUID[100];
   Float_t         RECO_PFJET_PUID_MVA[100];
   Double_t        RHO_ele;
   Double_t        RHO_mu;

     Float_t         RECO_PFJET_PT_UncUp[200];
   Float_t         RECO_PFJET_PT_UncDn[200];
   Float_t         RECO_PFJET_AREA[200];
   Float_t         RECO_PFJET_PTD[200];
   Float_t         RECO_PFJET_CHARGED_HADRON_ENERGY[200];
   Float_t         RECO_PFJET_NEUTRAL_HADRON_ENERGY[200];
   Float_t         RECO_PFJET_PHOTON_ENERGY[200];
   Float_t         RECO_PFJET_ELECTRON_ENERGY[200];
   Float_t         RECO_PFJET_MUON_ENERGY[200];
   Float_t         RECO_PFJET_HF_HADRON_ENERGY[200];
   Float_t         RECO_PFJET_HF_EM_ENERGY[200];
   Float_t         RECO_PFJET_CHARGED_EM_ENERGY[200];
   Float_t         RECO_PFJET_CHARGED_MU_ENERGY[200];
   Float_t         RECO_PFJET_NEUTRAL_EM_ENERGY[200];
   Int_t           RECO_PFJET_CHARGED_HADRON_MULTIPLICITY[200];
   Int_t           RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY[200];
   Int_t           RECO_PFJET_PHOTON_MULTIPLICITY[200];
   Int_t           RECO_PFJET_ELECTRON_MULTIPLICITY[200];
   Int_t           RECO_PFJET_MUON_MULTIPLICITY[200];
   Int_t           RECO_PFJET_HF_HADRON_MULTIPLICTY[200];
   Int_t           RECO_PFJET_HF_EM_MULTIPLICITY[200];
   Int_t           RECO_PFJET_CHARGED_MULTIPLICITY[200];
   Int_t           RECO_PFJET_NEUTRAL_MULTIPLICITY[200];
   Int_t           RECO_PFJET_NCOMPONENTS[200];
   Int_t           RECO_PFJET_COMPONENT_PDGID[200][100];
   Float_t         RECO_PFJET_COMPONENT_PT[200][100];
   Float_t         RECO_PFJET_COMPONENT_ETA[200][100];
   Float_t         RECO_PFJET_COMPONENT_PHI[200][100];
   Float_t         RECO_PFJET_COMPONENT_E[200][100];
   Float_t         RECO_PFJET_COMPONENT_CHARGE[200][100];
   Float_t         RECO_PFJET_COMPONENT_TRANSVERSE_MASS[200][100];
   Float_t         RECO_PFJET_COMPONENT_XVERTEX[200][100];
   Float_t         RECO_PFJET_COMPONENT_YVERTEX[200][100];
   Float_t         RECO_PFJET_COMPONENT_ZVERTEX[200][100];
   Float_t         RECO_PFJET_COMPONENT_VERTEX_CHI2[200][100];
   //@
   Int_t           LHE_PARTON_N;
   UChar_t         LHE_PARTON_CLEAR[10];
   Int_t           LHE_PARTON_PDGID[10];
   Float_t         LHE_PARTON_PT[10];
   Float_t         LHE_PARTON_ETA[10];
   Float_t         LHE_PARTON_PHI[10];
   Float_t         LHE_PARTON_E[10];

   Float_t         RECO_CALOMET;
   Float_t         RECO_PFMET;
   Float_t         RECO_PFMET_X;
   Float_t         RECO_PFMET_Y;
   Float_t         RECO_PFMET_PHI;
   Float_t         RECO_PFMET_THETA;
   Float_t         RECO_PFMETT1;
   Float_t         RECO_PFMETT1_X;
   Float_t         RECO_PFMETT1_Y;
   Float_t         RECO_PFMETT1_PHI;
   Float_t         RECO_PFMETT1_THETA;
   Float_t         RECO_TCMET;
   Float_t         RECO_CORMETMUONS;

 //@Type I correction MET
   Float_t         RECO_PFMET_JetEnUp;
   Float_t         RECO_PFMET_JetEnDn;
   Float_t         RECO_PFMET_ElectronEnUp;
   Float_t         RECO_PFMET_ElectronEnDn;
   Float_t         RECO_PFMET_MuonEnUp;
   Float_t         RECO_PFMET_MuonEnDn;
   Float_t         RECO_PFMET_JetResUp;
   Float_t         RECO_PFMET_JetResDn;
   Float_t         RECO_PFMET_UnclusteredEnUp;
   Float_t         RECO_PFMET_UnclusteredEnDn;
   Float_t         RECO_PFMET_PhotonEnUp;
   Float_t         RECO_PFMET_PhotonEnDn;

   Float_t         tCHighEff_BTagJet_PT[50];
   Float_t         tCHighPur_BTagJet_PT[50];
   Float_t         cSV_BTagJet_PT[50];
   Float_t         tCHighEff_BTagJet_ETA[50];
   Float_t         tCHighPur_BTagJet_ETA[50];
   Float_t         cSV_BTagJet_ETA[50];
   Float_t         tCHighEff_BTagJet_PHI[50];
   Float_t         tCHighPur_BTagJet_PHI[50];
   Float_t         cSV_BTagJet_PHI[50];
   Float_t         tCHighEff_BTagJet_DISCR[50];
   Float_t         tCHighPur_BTagJet_DISCR[50];
   Float_t         cSV_BTagJet_DISCR[50];
   Float_t         cSV_BTagJet_ET[50];

   // List of branches
   TBranch        *b_irun;   //!
   TBranch        *b_ievt;   //!
   TBranch        *b_ils;   //!
   TBranch        *b_Avginstlumi;   //!
   TBranch        *b_num_PU_vertices;   //!
   TBranch        *b_PU_BunchCrossing;   //!
   TBranch        *b_MC_weighting;   //!
   TBranch        *b_RECO_nMuHLTMatch;   //!
   TBranch        *b_RECOMU_PT_MuHLTMatch;   //!
   TBranch        *b_HLTPathsFired;   //!
   TBranch        *b_MC_E;   //!
   TBranch        *b_MC_PT;   //!
   TBranch        *b_MC_ETA;   //!
   TBranch        *b_MC_THETA;   //!
   TBranch        *b_MC_PHI;   //!
   TBranch        *b_MC_MASS;   //!
   TBranch        *b_MC_PDGID;   //!
   TBranch        *b_MC_LEPT_PT;   //!
   TBranch        *b_MC_LEPT_ETA;   //!
   TBranch        *b_MC_LEPT_PHI;   //!
   TBranch        *b_MC_LEPT_THETA;   //!
   TBranch        *b_MC_LEPT_PDGID;   //!
   TBranch        *b_MC_Z_PT;   //!
   TBranch        *b_MC_Z_ETA;   //!
   TBranch        *b_MC_Z_PHI;   //!
   TBranch        *b_MC_Z_THETA;   //!
   TBranch        *b_MC_Z_MASS;   //!
   TBranch        *b_MC_Z_PDGID;   //!
   TBranch        *b_MC_fourl_MASS;   //!
   TBranch        *b_MC_fourl_PT;   //!
   TBranch        *b_MC_fourl_PDGID;   //!
   TBranch        *b_MC_ZZ_MASS;   //!
   TBranch        *b_MC_ZZ_PT;   //!
   TBranch        *b_MC_ZZ_ETA;   //!
   TBranch        *b_MC_ZZ_PHI;   //!
   TBranch        *b_MC_ZZ_THETA;   //!
   TBranch        *b_MC_ZZ_PDGID;   //!
   TBranch        *b_MC_GENJET_PT;  //!
   TBranch        *b_MC_GENJET_ETA; //!
   TBranch        *b_MC_GENJET_PHI;  //!
   TBranch        *b_MC_GENMET;   //!
   TBranch        *b_RECORF_2e2mu_cosTheta1_spin;   //!
   TBranch        *b_RECORF_2e2mu_cosTheta2_spin;   //!
   TBranch        *b_RECORF_2e2mu_cosThetaStar_spin;   //!
   TBranch        *b_RECORF_2e2mu_Phi_spin;   //!
   TBranch        *b_RECORF_2e2mu_Phi1_spin;   //!
   TBranch        *b_RECORF_2e2mu_Phi2_spin;   //!
   TBranch        *b_RECORF_2e2mu_phi1RF_spin;   //!
   TBranch        *b_RECORF_2e2mu_phi2RF_spin;   //!
   TBranch        *b_RECORF_2e2mu_MELA;   //!
   TBranch        *b_RECORF_4e_cosTheta1_spin;   //!
   TBranch        *b_RECORF_4e_cosTheta2_spin;   //!
   TBranch        *b_RECORF_4e_cosThetaStar_spin;   //!
   TBranch        *b_RECORF_4e_Phi_spin;   //!
   TBranch        *b_RECORF_4e_Phi1_spin;   //!
   TBranch        *b_RECORF_4e_Phi2_spin;   //!
   TBranch        *b_RECORF_4e_phi1RF_spin;   //!
   TBranch        *b_RECORF_4e_phi2RF_spin;   //!
   TBranch        *b_RECORF_4e_MELA;   //!
   TBranch        *b_RECORF_4mu_cosTheta1_spin;   //!
   TBranch        *b_RECORF_4mu_cosTheta2_spin;   //!
   TBranch        *b_RECORF_4mu_cosThetaStar_spin;   //!
   TBranch        *b_RECORF_4mu_Phi_spin;   //!
   TBranch        *b_RECORF_4mu_Phi1_spin;   //!
   TBranch        *b_RECORF_4mu_Phi2_spin;   //!
   TBranch        *b_RECORF_4mu_phi1RF_spin;   //!
   TBranch        *b_RECORF_4mu_phi2RF_spin;   //!
   TBranch        *b_RECORF_4mu_MELA;   //!
   TBranch        *b_RECO_ZMM_MASS;   //!
   TBranch        *b_RECO_ZEE_MASS;   //!
   TBranch        *b_RECO_DiLep_MASS;   //!
   TBranch        *b_RECO_ZMM_PT;   //!
   TBranch        *b_RECO_ZEE_PT;   //!
   TBranch        *b_RECO_DiLep_PT;   //!
   TBranch        *b_RECO_ZMM_ETA;   //!
   TBranch        *b_RECO_ZEE_ETA;   //!
   TBranch        *b_RECO_DiLep_ETA;   //!
   TBranch        *b_RECO_ZMM_PHI;   //!
   TBranch        *b_RECO_ZEE_PHI;   //!
   TBranch        *b_RECO_DiLep_PHI;   //!
   TBranch        *b_RECO_ZMMss_MASS;   //!
   TBranch        *b_RECO_ZEEss_MASS;   //!
   TBranch        *b_RECO_ZEM_MASS;   //!
   TBranch        *b_RECO_ZMMss_PT;   //!
   TBranch        *b_RECO_ZEEss_PT;   //!
   TBranch        *b_RECO_ZEM_PT;   //!
   TBranch        *b_RECO_ZMMss_ETA;   //!
   TBranch        *b_RECO_ZEEss_ETA;   //!
   TBranch        *b_RECO_ZEM_ETA;   //!
   TBranch        *b_RECO_ZMMss_PHI;   //!
   TBranch        *b_RECO_ZEEss_PHI;   //!
   TBranch        *b_RECO_ZEM_PHI;   //!
   TBranch        *b_RECO_MMMM_MASS;   //!
   TBranch        *b_RECO_MMMM_PT;   //!
   TBranch        *b_RECO_MMMM_ETA;   //!
   TBranch        *b_RECO_MMMM_PHI;   //!
   TBranch        *b_RECO_MMMM_MASS_REFIT;   //!
   TBranch        *b_RECO_EEEE_MASS;   //!
   TBranch        *b_RECO_EEEE_PT;   //!
   TBranch        *b_RECO_EEEE_ETA;   //!
   TBranch        *b_RECO_EEEE_PHI;   //!
   TBranch        *b_RECO_EEEE_MASS_REFIT;   //!
   TBranch        *b_RECO_EEMM_MASS;   //!
   TBranch        *b_RECO_EEMM_PT;   //!
   TBranch        *b_RECO_EEMM_ETA;   //!
   TBranch        *b_RECO_EEMM_PHI;   //!
   TBranch        *b_RECO_EEMM_MASS_REFIT;   //!
   TBranch        *b_RECO_LLL0_MASS;   //!
   TBranch        *b_RECO_LLL1_MASS;   //!
   TBranch        *b_RECO_LLL2_MASS;   //!
   TBranch        *b_RECO_LLL3_MASS;   //!
   TBranch        *b_RECO_LLL0_PT;   //!
   TBranch        *b_RECO_LLL1_PT;   //!
   TBranch        *b_RECO_LLL2_PT;   //!
   TBranch        *b_RECO_LLL3_PT;   //!
   TBranch        *b_RECO_LLLl0_MASS;   //!
   TBranch        *b_RECO_LLLl1_MASS;   //!
   TBranch        *b_RECO_LLLl0_PT;   //!
   TBranch        *b_RECO_LLLl1_PT;   //!
   TBranch        *b_RECO_LLLL0ss_MASS;   //!
   TBranch        *b_RECO_LLLL0ss_PT;   //!
   TBranch        *b_RECO_LLLL1ss_MASS;   //!
   TBranch        *b_RECO_LLLL1ss_PT;   //!
   TBranch        *b_RECO_LLLL2ss_MASS;   //!
   TBranch        *b_RECO_LLLL2ss_PT;   //!
   TBranch        *b_RECO_LLLL_MASS;   //!
   TBranch        *b_RECO_LLLL_PT;   //!
   TBranch        *b_RECO_LLLL_ETA;   //!
   TBranch        *b_RECO_LLLL_PHI;   //!
   TBranch        *b_RECOELE_E;   //!
   TBranch        *b_RECOELE_PT;   //!
   TBranch        *b_RECOELE_PTError;   //!
   TBranch        *b_RECOELE_P;   //!
   TBranch        *b_RECOELE_ETA;   //!
   TBranch        *b_RECOELE_THETA;   //!
   TBranch        *b_RECOELE_PHI;   //!
   TBranch        *b_RECOELE_MASS;   //!
   TBranch        *b_RECOELE_CHARGE;   //!
   TBranch        *b_RECOELE_isEcalDriven;   //!
   TBranch        *b_RECOELE_isTrackerDriven;   //!
   TBranch        *b_RECOELE_gsftrack_NPixHits;   //!
   TBranch        *b_RECOELE_gsftrack_NStripHits;   //!
   TBranch        *b_RECOELE_gsftrack_chi2;   //!
   TBranch        *b_RECOELE_gsftrack_dxyB;   //!
   TBranch        *b_RECOELE_gsftrack_dxy;   //!
   TBranch        *b_RECOELE_gsftrack_dxyError;   //!
   TBranch        *b_RECOELE_gsftrack_dzB;   //!
   TBranch        *b_RECOELE_gsftrack_dz;   //!
   TBranch        *b_RECOELE_gsftrack_dzError;   //!
   TBranch        *b_RECOELE_gsftrack_losthits;   //!
   TBranch        *b_RECOELE_gsftrack_validhits;   //!
   TBranch        *b_RECOELE_gsftrack_expected_inner_hits;   //!
   TBranch        *b_RECOELE_scl_E;   //!
   TBranch        *b_RECOELE_scl_Et;   //!
   TBranch        *b_RECOELE_scl_Eta;   //!
   TBranch        *b_RECOELE_scl_Phi;   //!
   TBranch        *b_RECOELE_ep;   //!
   TBranch        *b_RECOELE_eSeedp;   //!
   TBranch        *b_RECOELE_eSeedpout;   //!
   TBranch        *b_RECOELE_eElepout;   //!
   TBranch        *b_RECOELE_deltaEtaIn;   //!
   TBranch        *b_RECOELE_deltaEtaSeed;   //!
   TBranch        *b_RECOELE_deltaEtaEle;   //!
   TBranch        *b_RECOELE_deltaPhiIn;   //!
   TBranch        *b_RECOELE_deltaPhiSeed;   //!
   TBranch        *b_RECOELE_deltaPhiEle;   //!
   TBranch        *b_RECOELE_isbarrel;   //!
   TBranch        *b_RECOELE_isendcap;   //!
   TBranch        *b_RECOELE_isEBetaGap;   //!
   TBranch        *b_RECOELE_isEBphiGap;   //!
   TBranch        *b_RECOELE_isEEdeeGap;   //!
   TBranch        *b_RECOELE_isEEringGap;   //!
   TBranch        *b_RECOELE_sigmaIetaIeta;   //!
   TBranch        *b_RECOELE_sigmaEtaEta;   //!
   TBranch        *b_RECOELE_e15;   //!
   TBranch        *b_RECOELE_e25max;   //!
   TBranch        *b_RECOELE_e55;   //!
   TBranch        *b_RECOELE_he;   //!
   TBranch        *b_RECOELE_r9;   //!
   TBranch        *b_RECOELE_mva;   //!
   TBranch        *b_RECOELE_fbrem;   //!
   TBranch        *b_RECOELE_nbrems;   //!
   TBranch        *b_RECOELE_Class;   //!
   TBranch        *b_RECOELE_fbrem_mode;   //!
   TBranch        *b_RECOELE_fbrem_mean;   //!
   TBranch        *b_RECOELE_EGMTRACKISO;   //!
   TBranch        *b_RECOELE_EGMHCALISO;   //!
   TBranch        *b_RECOELE_EGMECALISO;   //!
   TBranch        *b_RECOELE_EGMX;   //!
   TBranch        *b_RECOELE_PFchAllPart;   //!
   TBranch        *b_RECOELE_PFchHad;   //!
   TBranch        *b_RECOELE_PFneuHad;   //!
   TBranch        *b_RECOELE_PFphoton;   //!
   TBranch        *b_RECOELE_PFPUchAllPart;   //!
   TBranch        *b_RECOELE_PFX_dB;   //!
   TBranch        *b_RECOELE_PFX_rho;   //!
   TBranch        *b_RECOELE_regEnergy;   //!
   TBranch        *b_RECOELE_regEnergyError;   //!
   TBranch        *b_RECOELE_SIP;   //!
   TBranch        *b_RECOELE_IP;   //!
   TBranch        *b_RECOELE_IPERROR;   //!
   TBranch        *b_RECOELE_SIP_KF;   //!
   TBranch        *b_RECOELE_IP_KF;   //!
   TBranch        *b_RECOELE_IPERROR_KF;   //!
   TBranch        *b_RECOELE_SIP_GD;   //!
   TBranch        *b_RECOELE_SIP_GDEEEE;   //!
   TBranch        *b_RECOELE_SIP_Std;   //!
   TBranch        *b_RECOELE_SIP_StdEEEE;   //!
   TBranch        *b_RECOELE_SIP_Kin;   //!
   TBranch        *b_RECOELE_SIP_KinEEEE;   //!
   TBranch        *b_RECOELE_STIP;   //!
   TBranch        *b_RECOELE_SLIP;   //!
   TBranch        *b_RECOELE_TIP;   //!
   TBranch        *b_RECOELE_LIP;   //!
   TBranch        *b_RECOELE_TIPERROR;   //!
   TBranch        *b_RECOELE_LIPERROR;   //!
   TBranch        *b_RECOELE_sclRawE;   //!
   TBranch        *b_RECOELE_sclX;   //!
   TBranch        *b_RECOELE_sclY;   //!
   TBranch        *b_RECOELE_sclZ;   //!
   TBranch        *b_RECOELE_seedSubdet1;   //!
   TBranch        *b_RECOELE_seedDphi;   //!
   TBranch        *b_RECOELE_seedDrz1;   //!
   TBranch        *b_RECOELE_seedSubdet2;   //!
   TBranch        *b_RECOELE_seedDphi2;   //!
   TBranch        *b_RECOELE_seedDrz2;   //!
   TBranch        *b_RECOELE_eidVeryLoose;   //!
   TBranch        *b_RECOELE_eidLoose;   //!
   TBranch        *b_RECOELE_eidMedium;   //!
   TBranch        *b_RECOELE_eidTight;   //!
   TBranch        *b_RECOELE_eidHZZVeryLoose;   //!
   TBranch        *b_RECOELE_eidHZZLoose;   //!
   TBranch        *b_RECOELE_eidHZZMedium;   //!
   TBranch        *b_RECOELE_eidHZZTight;   //!
   TBranch        *b_RECOELE_mvaTrigV0;   //!
   TBranch        *b_RECOELE_mvaNonTrigV0;   //!
   TBranch        *b_RECOELE_ID;   //!
   TBranch        *b_RECOELE_COV;   //!
   TBranch        *b_RECOMU_isPFMu;   //!
   TBranch        *b_RECOMU_isGlobalMu;   //!
   TBranch        *b_RECOMU_isStandAloneMu;   //!
   TBranch        *b_RECOMU_isTrackerMu;   //!
   TBranch        *b_RECOMU_isCaloMu;   //!
   TBranch        *b_RECOMU_isTrackerHighPtMu;   //!
   TBranch        *b_RECOMU_E;   //!
   TBranch        *b_RECOMU_PT;   //!
   TBranch        *b_RECOMU_P;   //!
   TBranch        *b_RECOMU_ETA;   //!
   TBranch        *b_RECOMU_THETA;   //!
   TBranch        *b_RECOMU_PHI;   //!
   TBranch        *b_RECOMU_MASS;   //!
   TBranch        *b_RECOMU_CHARGE;   //!
   TBranch        *b_RECOMU_COV;   //!
   TBranch        *b_RECOMU_TRACKISO;   //!
   TBranch        *b_RECOMU_TRACKISO_SUMPT;   //!
   TBranch        *b_RECOMU_HCALISO;   //!
   TBranch        *b_RECOMU_ECALISO;   //!
   TBranch        *b_RECOMU_X;   //!
   TBranch        *b_RECOMU_PFchHad;   //!
   TBranch        *b_RECOMU_PFneuHad;   //!
   TBranch        *b_RECOMU_PFphoton;   //!
   TBranch        *b_RECOMU_PFPUchAllPart;   //!
   TBranch        *b_RECOMU_PFX_dB;   //!
   TBranch        *b_RECOMU_PFX_rho;   //!
   TBranch        *b_RECOPFPHOT_PFchHad;   //!
   TBranch        *b_RECOPFPHOT_PFneuHad;   //!
   TBranch        *b_RECOPFPHOT_PFphoton;   //!
   TBranch        *b_RECOPFPHOT_PFPUchAllPart;   //!
   TBranch        *b_RECOPFPHOT_PFX_rho;   //!
   TBranch        *b_RECOMU_SIP;   //!
   TBranch        *b_RECOMU_IP;   //!
   TBranch        *b_RECOMU_IPERROR;   //!
   TBranch        *b_RECOMU_SIP_KF;   //!
   TBranch        *b_RECOMU_IP_KF;   //!
   TBranch        *b_RECOMU_IPERROR_KF;   //!
   TBranch        *b_RECOMU_SIP_GD;   //!
   TBranch        *b_RECOMU_SIP_GDMMMM;   //!
   TBranch        *b_RECOMU_SIP_Std;   //!
   TBranch        *b_RECOMU_SIP_StdMMMM;   //!
   TBranch        *b_RECOMU_SIP_Kin;   //!
   TBranch        *b_RECOMU_SIP_KinMMMM;   //!
   TBranch        *b_RECOMU_STIP;   //!
   TBranch        *b_RECOMU_SLIP;   //!
   TBranch        *b_RECOMU_TIP;   //!
   TBranch        *b_RECOMU_LIP;   //!
   TBranch        *b_RECOMU_TIPERROR;   //!
   TBranch        *b_RECOMU_LIPERROR;   //!
   TBranch        *b_RECOMU_caloCompatibility;   //!
   TBranch        *b_RECOMU_segmentCompatibility;   //!
   TBranch        *b_RECOMU_numberOfMatches;   //!
   TBranch        *b_RECOMU_numberOfMatchedStations;   //!
   TBranch        *b_RECOMU_glbmuPromptTight;   //!
   TBranch        *b_RECOMU_trkmuArbitration;   //!
   TBranch        *b_RECOMU_trkmu2DCompatibilityLoose;   //!
   TBranch        *b_RECOMU_trkmu2DCompatibilityTight;   //!
   TBranch        *b_RECOMU_trkmuOneStationLoose;   //!
   TBranch        *b_RECOMU_trkmuOneStationTight;   //!
   TBranch        *b_RECOMU_trkmuLastStationLoose;   //!
   TBranch        *b_RECOMU_trkmuLastStationTight;   //!
   TBranch        *b_RECOMU_trkmuOneStationAngLoose;   //!
   TBranch        *b_RECOMU_trkmuOneStationAngTight;   //!
   TBranch        *b_RECOMU_trkmuLastStationAngLoose;   //!
   TBranch        *b_RECOMU_trkmuLastStationAngTight;   //!
   TBranch        *b_RECOMU_trkmuLastStationOptimizedLowPtLoose;   //!
   TBranch        *b_RECOMU_trkmuLastStationOptimizedLowPtTight;   //!
   TBranch        *b_RECOMU_mutrkPT;   //!
   TBranch        *b_RECOMU_mutrkPTError;   //!
   TBranch        *b_RECOMU_mutrkDxy;   //!
   TBranch        *b_RECOMU_mutrkDxyError;   //!
   TBranch        *b_RECOMU_mutrkDxyB;   //!
   TBranch        *b_RECOMU_mutrkDz;   //!
   TBranch        *b_RECOMU_mutrkDzError;   //!
   TBranch        *b_RECOMU_mutrkDzB;   //!
   TBranch        *b_RECOMU_mutrkChi2PerNdof;   //!
   TBranch        *b_RECOMU_mutrkCharge;   //!
   TBranch        *b_RECOMU_mutrkNHits;   //!
   TBranch        *b_RECOMU_mutrkNStripHits;   //!
   TBranch        *b_RECOMU_mutrkNPixHits;   //!
   TBranch        *b_RECOMU_mutrkNMuonHits;   //!
   TBranch        *b_RECOMU_mutrktrackerLayersWithMeasurement;   //!
   TBranch        *b_RECOMU_muInnertrkDxy;   //!
   TBranch        *b_RECOMU_muInnertrkDxyError;   //!
   TBranch        *b_RECOMU_muInnertrkDxyB;   //!
   TBranch        *b_RECOMU_muInnertrkDz;   //!
   TBranch        *b_RECOMU_muInnertrkDzError;   //!
   TBranch        *b_RECOMU_muInnertrkDzB;   //!
   TBranch        *b_RECOMU_muInnertrkChi2PerNdof;   //!
   TBranch        *b_RECOMU_muInnertrktrackerLayersWithMeasurement;   //!
   TBranch        *b_RECOMU_muInnertrkPT;   //!
   TBranch        *b_RECOMU_muInnertrkPTError;   //!
   TBranch        *b_RECOMU_muInnertrkCharge;   //!
   TBranch        *b_RECOMU_muInnertrkNHits;   //!
   TBranch        *b_RECOMU_muInnertrkNStripHits;   //!
   TBranch        *b_RECOMU_muInnertrkNPixHits;   //!
   TBranch        *b_RECOMU_mubesttrkType;   //!
   TBranch        *b_RECOMU_mubesttrkDxy;   //!
   TBranch        *b_RECOMU_mubesttrkDxyError;   //!
   TBranch        *b_RECOMU_mubesttrkDz;   //!
   TBranch        *b_RECOMU_mubesttrkDzError;   //!
   TBranch        *b_RECOMU_mubesttrkPTError; //!
   TBranch        *b_ftsigma;   //!
   TBranch        *b_gdX;   //!
   TBranch        *b_gdY;   //!
   TBranch        *b_gdZ;   //!
   TBranch        *b_ftsigmalag;   //!
   TBranch        *b_gdlagX;   //!
   TBranch        *b_gdlagY;   //!
   TBranch        *b_gdlagZ;   //!
   TBranch        *b_gdlagProb;   //!
   TBranch        *b_gdlagNdof;   //!
   TBranch        *b_ftsigmaMMMM;   //!
   TBranch        *b_gdXMMMM;   //!
   TBranch        *b_gdYMMMM;   //!
   TBranch        *b_gdZMMMM;   //!
   TBranch        *b_ftsigmalagMMMM;   //!
   TBranch        *b_gdlagXMMMM;   //!
   TBranch        *b_gdlagYMMMM;   //!
   TBranch        *b_gdlagZMMMM;   //!
   TBranch        *b_gdlagProbMMMM;   //!
   TBranch        *b_gdlagNdofMMMM;   //!
   TBranch        *b_ftsigmaEEEE;   //!
   TBranch        *b_gdXEEEE;   //!
   TBranch        *b_gdYEEEE;   //!
   TBranch        *b_gdZEEEE;   //!
   TBranch        *b_ftsigmalagEEEE;   //!
   TBranch        *b_gdlagXEEEE;   //!
   TBranch        *b_gdlagYEEEE;   //!
   TBranch        *b_gdlagZEEEE;   //!
   TBranch        *b_gdlagProbEEEE;   //!
   TBranch        *b_gdlagNdofEEEE;   //!
   TBranch        *b_StdFitVertexX;   //!
   TBranch        *b_StdFitVertexY;   //!
   TBranch        *b_StdFitVertexZ;   //!
   TBranch        *b_StdFitVertexChi2r;   //!
   TBranch        *b_StdFitVertexProb;   //!
   TBranch        *b_StdFitVertexTrack_PT;   //!
   TBranch        *b_StdFitVertexTrack_ETA;   //!
   TBranch        *b_StdFitVertexTrack_PHI;   //!
   TBranch        *b_KinFitVertexX;   //!
   TBranch        *b_KinFitVertexY;   //!
   TBranch        *b_KinFitVertexZ;   //!
   TBranch        *b_KinFitVertexChi2r;   //!
   TBranch        *b_KinFitVertexProb;   //!
   TBranch        *b_StdFitVertexXMMMM;   //!
   TBranch        *b_StdFitVertexYMMMM;   //!
   TBranch        *b_StdFitVertexZMMMM;   //!
   TBranch        *b_StdFitVertexChi2rMMMM;   //!
   TBranch        *b_StdFitVertexProbMMMM;   //!
   TBranch        *b_StdFitVertexTrackMMMM_PT;   //!
   TBranch        *b_StdFitVertexTrackMMMM_ETA;   //!
   TBranch        *b_StdFitVertexTrackMMMM_PHI;   //!
   TBranch        *b_KinFitVertexXMMMM;   //!
   TBranch        *b_KinFitVertexYMMMM;   //!
   TBranch        *b_KinFitVertexZMMMM;   //!
   TBranch        *b_KinFitVertexChi2rMMMM;   //!
   TBranch        *b_KinFitVertexProbMMMM;   //!
   TBranch        *b_StdFitVertexXEEEE;   //!
   TBranch        *b_StdFitVertexYEEEE;   //!
   TBranch        *b_StdFitVertexZEEEE;   //!
   TBranch        *b_StdFitVertexChi2rEEEE;   //!
   TBranch        *b_StdFitVertexProbEEEE;   //!
   TBranch        *b_StdFitVertexTrackEEEE_PT;   //!
   TBranch        *b_StdFitVertexTrackEEEE_ETA;   //!
   TBranch        *b_StdFitVertexTrackEEEE_PHI;   //!
   TBranch        *b_KinFitVertexXEEEE;   //!
   TBranch        *b_KinFitVertexYEEEE;   //!
   TBranch        *b_KinFitVertexZEEEE;   //!
   TBranch        *b_KinFitVertexChi2rEEEE;   //!
   TBranch        *b_KinFitVertexProbEEEE;   //!
   TBranch        *b_StdFitVertexChi2rMMM;   //!
   TBranch        *b_StdFitVertexProbMMM;   //!
   TBranch        *b_StdFitVertexChi2rMME;   //!
   TBranch        *b_StdFitVertexProbMME;   //!
   TBranch        *b_StdFitVertexChi2rEEE;   //!
   TBranch        *b_StdFitVertexProbEEE;   //!
   TBranch        *b_StdFitVertexChi2rMEE;   //!
   TBranch        *b_StdFitVertexProbMEE;   //!
   TBranch        *b_StdFitVertexChi2rDiLep;   //!
   TBranch        *b_StdFitVertexProbDiLep;   //!
   TBranch        *b_ConvMapDist;   //!
   TBranch        *b_ConvMapDcot;   //!
   TBranch        *b_RECOMU_MatchingMCTruth;   //!
   TBranch        *b_RECOMU_MatchingMCpT;   //!
   TBranch        *b_RECOMU_MatchingMCEta;   //!
   TBranch        *b_RECOMU_MatchingMCPhi;   //!
   TBranch        *b_RECOELE_MatchingMCTruth;   //!
   TBranch        *b_RECOELE_MatchingMCpT;   //!
   TBranch        *b_RECOELE_MatchingMCEta;   //!
   TBranch        *b_RECOELE_MatchingMCPhi;   //!
   TBranch        *b_RECOPHOT_MatchingMCTruth;   //!
   TBranch        *b_RECOPHOT_MatchingMCpT;   //!
   TBranch        *b_RECOPHOT_MatchingMCEta;   //!
   TBranch        *b_RECOPHOT_MatchingMCPhi;   //!
   TBranch        *b_RECOzMuMu_MatchingMCTruth;   //!
   TBranch        *b_RECOzMuMu_MatchingMCpT;   //!
   TBranch        *b_RECOzMuMu_MatchingMCmass;   //!
   TBranch        *b_RECOzMuMu_MatchingMCEta;   //!
   TBranch        *b_RECOzMuMu_MatchingMCPhi;   //!
   TBranch        *b_RECOzEE_MatchingMCTruth;   //!
   TBranch        *b_RECOzEE_MatchingMCpT;   //!
   TBranch        *b_RECOzEE_MatchingMCmass;   //!
   TBranch        *b_RECOzEE_MatchingMCEta;   //!
   TBranch        *b_RECOzEE_MatchingMCPhi;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCTruth;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCpT;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCmass;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCEta;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCPhi;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCTruth;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCpT;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCmass;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCEta;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCPhi;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCTruth;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCpT;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCmass;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCEta;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCPhi;   //!
   TBranch        *b_RECO_NMU;   //!
   TBranch        *b_RECO_NELE;   //!
   TBranch        *b_RECO_NTRACK;   //!
   TBranch        *b_RECO_TRACK_PT;   //!
   TBranch        *b_RECO_TRACK_ETA;   //!
   TBranch        *b_RECO_TRACK_PHI;   //!
   TBranch        *b_RECO_TRACK_CHI2;   //!
   TBranch        *b_RECO_TRACK_CHI2RED;   //!
   TBranch        *b_RECO_TRACK_CHI2PROB;   //!
   TBranch        *b_RECO_TRACK_NHITS;   //!
   TBranch        *b_RECO_TRACK_DXY;   //!
   TBranch        *b_RECO_TRACK_DXYERR;   //!
   TBranch        *b_RECO_TRACK_DZ;   //!
   TBranch        *b_RECO_TRACK_DZERR;   //!
   TBranch        *b_RECO_NPHOT;   //!
   TBranch        *b_RECOPHOT_PT;   //!
   TBranch        *b_RECOPHOT_ETA;   //!
   TBranch        *b_RECOPHOT_PHI;   //!
   TBranch        *b_RECOPHOT_THETA;   //!
   TBranch        *b_RECO_NPFPHOT;   //!
   TBranch        *b_RECOPFPHOT_PT;   //!
   TBranch        *b_RECOPFPHOT_PTError;   //!
   TBranch        *b_RECOPFPHOT_ETA;   //!
   TBranch        *b_RECOPFPHOT_PHI;   //!
   TBranch        *b_RECOPFPHOT_THETA;   //!
   TBranch        *b_BeamSpot_X;   //!
   TBranch        *b_BeamSpot_Y;   //!
   TBranch        *b_BeamSpot_Z;   //!
   TBranch        *b_RECO_NVTX;   //!
   TBranch        *b_RECO_VERTEX_x;   //!
   TBranch        *b_RECO_VERTEX_y;   //!
   TBranch        *b_RECO_VERTEX_z;   //!
   TBranch        *b_RECO_VERTEX_ndof;   //!
   TBranch        *b_RECO_VERTEX_chi2;   //!
   TBranch        *b_RECO_VERTEX_ntracks;   //!
   TBranch        *b_RECO_VERTEXPROB;   //!
   TBranch        *b_RECO_VERTEX_isValid;   //!
   TBranch        *b_RECO_VERTEX_TRACK_PT;   //!
   TBranch        *b_RECO_PFJET_N;   //!
   TBranch        *b_RECO_PFJET_CHARGE;   //!
   TBranch        *b_RECO_PFJET_ET;   //!
   TBranch        *b_RECO_PFJET_PT;   //!
   TBranch        *b_RECO_PFJET_ETA;   //!
   TBranch        *b_RECO_PFJET_PHI;   //!
   TBranch        *b_RECO_PFJET_PUID;   //!
   TBranch        *b_RECO_PFJET_PUID_MVA; //!
   TBranch        *b_RHO_ele;   //!
   TBranch        *b_RHO_mu;   //!

 TBranch        *b_RECO_PFJET_PT_UncUp;//!
   TBranch        *b_RECO_PFJET_PT_UncDn;//!
   TBranch        *b_RECO_PFJET_AREA;   //!
   TBranch        *b_RECO_PFJET_PTD;   //!
   TBranch        *b_RECO_PFJET_CHARGED_HADRON_ENERGY;   //!
   TBranch        *b_RECO_PFJET_NEUTRAL_HADRON_ENERGY;   //!
   TBranch        *b_RECO_PFJET_PHOTON_ENERGY;   //!
   TBranch        *b_RECO_PFJET_ELECTRON_ENERGY;   //!
   TBranch        *b_RECO_PFJET_MUON_ENERGY;   //!
   TBranch        *b_RECO_PFJET_HF_HADRON_ENERGY;   //!
   TBranch        *b_RECO_PFJET_HF_EM_ENERGY;   //!
   TBranch        *b_RECO_PFJET_CHARGED_EM_ENERGY;   //!
   TBranch        *b_RECO_PFJET_CHARGED_MU_ENERGY;   //!
   TBranch        *b_RECO_PFJET_NEUTRAL_EM_ENERGY;   //!
   TBranch        *b_RECO_PFJET_CHARGED_HADRON_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_PHOTON_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_ELECTRON_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_MUON_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_HF_HADRON_MULTIPLICTY;   //!
   TBranch        *b_RECO_PFJET_HF_EM_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_CHARGED_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_NEUTRAL_MULTIPLICITY;   //!
   TBranch        *b_RECO_PFJET_NCOMPONENTS;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_PDGID;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_PT;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_ETA;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_PHI;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_E;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_CHARGE;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_TRANSVERSE_MASS;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_XVERTEX;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_YVERTEX;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_ZVERTEX;   //!
   TBranch        *b_RECO_PFJET_COMPONENT_VERTEX_CHI2;   //!
   //@
   TBranch        *b_LHE_PARTON_N;   //!
   TBranch        *b_LHE_PARTON_CLEAR;   //!
   TBranch        *b_LHE_PARTON_PDGID;   //!
   TBranch        *b_LHE_PARTON_PT;   //!
   TBranch        *b_LHE_PARTON_ETA;   //!
   TBranch        *b_LHE_PARTON_PHI;   //!
   TBranch        *b_LHE_PARTON_E;   //!

   TBranch        *b_RECO_CALOMET;   //!
   TBranch        *b_RECO_PFMET;   //!
   TBranch        *b_RECO_PFMET_X;   //!
   TBranch        *b_RECO_PFMET_Y;   //!
   TBranch        *b_RECO_PFMET_PHI;   //!
   TBranch        *b_RECO_PFMET_THETA;   //!
   TBranch        *b_RECO_PFMETT1;   //!
   TBranch        *b_RECO_PFMETT1_X;   //!
   TBranch        *b_RECO_PFMETT1_Y;   //!
   TBranch        *b_RECO_PFMETT1_PHI;   //!
   TBranch        *b_RECO_PFMETT1_THETA;   //!
   TBranch        *b_RECO_TCMET;   //!
   TBranch        *b_RECO_CORMETMUONS;   //!

 //@Type I correction MET
   TBranch        *b_RECO_PFMET_JetEnUp;
   TBranch        *b_RECO_PFMET_JetEnDn;
   TBranch        *b_RECO_PFMET_ElectronEnUp;
   TBranch        *b_RECO_PFMET_ElectronEnDn;
   TBranch        *b_RECO_PFMET_MuonEnUp;
   TBranch        *b_RECO_PFMET_MuonEnDn;
   TBranch        *b_RECO_PFMET_JetResUp;
   TBranch        *b_RECO_PFMET_JetResDn;
   TBranch        *b_RECO_PFMET_UnclusteredEnUp;
   TBranch        *b_RECO_PFMET_UnclusteredEnDn;
   TBranch        *b_RECO_PFMET_PhotonEnUp;
   TBranch        *b_RECO_PFMET_PhotonEnDn;

   TBranch        *b_tCHighEff_BTagJet_PT;   //!
   TBranch        *b_tCHighPur_BTagJet_PT;   //!
   TBranch        *b_cSV_BTagJet_PT;   //!
   TBranch        *b_tCHighEff_BTagJet_ETA;   //!
   TBranch        *b_tCHighPur_BTagJet_ETA;   //!
   TBranch        *b_cSV_BTagJet_ETA;   //!
   TBranch        *b_tCHighEff_BTagJet_PHI;   //!
   TBranch        *b_tCHighPur_BTagJet_PHI;   //!
   TBranch        *b_cSV_BTagJet_PHI;   //!
   TBranch        *b_tCHighEff_BTagJet_DISCR;   //!
   TBranch        *b_tCHighPur_BTagJet_DISCR;   //!
   TBranch        *b_cSV_BTagJet_DISCR;   //!
   TBranch        *b_cSV_BTagJet_ET;   //!


   HZZ4LeptonsAnalysis(TTree *tree=0,Double_t weight_=1.,string DATA_type_="DATA",string MC_type_="MC");
   virtual ~HZZ4LeptonsAnalysis();
   Double_t weight;
   string DATA_type,MC_type;
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Char_t *name);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   ofstream bnn_file;
   double EAele(int ,bool );
   double masserror( std::vector<TLorentzVector> Lep, std::vector<double> pterr );
   void printmubnn(int i);
   float kfactor_qqZZ_qcd_dPhi(float GENdPhiZZ, int finalState);
   float kfactor_qqZZ_qcd_M(float GENmassZZ, int finalState);
   float kfactor_qqZZ_qcd_Pt(float GENpTZZ, int finalState);
   float kfactor_ggZZ(float GENmassZZ, int finalState);
};

#endif

#ifdef HZZ4LeptonsAnalysis_cxx
HZZ4LeptonsAnalysis::HZZ4LeptonsAnalysis(TTree *tree,Double_t weight_, string DATA_type_, string MC_type_) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   weight = weight_;
   DATA_type = DATA_type_;
   MC_type = MC_type_;

   if (tree == 0) {
      TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
      chain->Add("/lustre/cms/store/user/ndefilip/Summer12_52X_merged/roottree_leptons_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball.root_a");
      chain->Add("/lustre/cms/store/user/ndefilip/Summer12_52X_merged/roottree_leptons_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_1.root_b");

      tree = chain;

   }
   Init(tree);
}

HZZ4LeptonsAnalysis::~HZZ4LeptonsAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HZZ4LeptonsAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HZZ4LeptonsAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HZZ4LeptonsAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);


   fChain->SetBranchAddress("Run", &Run, &b_irun);
   fChain->SetBranchAddress("Event", &Event, &b_ievt);
   fChain->SetBranchAddress("LumiSection", &LumiSection, &b_ils);
   fChain->SetBranchAddress("Avginstlumi", &Avginstlumi, &b_Avginstlumi);
   fChain->SetBranchAddress("num_PU_vertices", &num_PU_vertices, &b_num_PU_vertices);
   fChain->SetBranchAddress("PU_BunchCrossing", &PU_BunchCrossing, &b_PU_BunchCrossing);
   fChain->SetBranchAddress("MC_weighting", &MC_weighting, &b_MC_weighting);
   fChain->SetBranchAddress("RECO_nMuHLTMatch", &RECO_nMuHLTMatch, &b_RECO_nMuHLTMatch);
   fChain->SetBranchAddress("RECOMU_PT_MuHLTMatch", RECOMU_PT_MuHLTMatch, &b_RECOMU_PT_MuHLTMatch);
   fChain->SetBranchAddress("HLTPathsFired", HLTPathsFired, &b_HLTPathsFired);
   fChain->SetBranchAddress("MC_E", MC_E, &b_MC_E);
   fChain->SetBranchAddress("MC_PT", MC_PT, &b_MC_PT);
   fChain->SetBranchAddress("MC_ETA", MC_ETA, &b_MC_ETA);
   fChain->SetBranchAddress("MC_THETA", MC_THETA, &b_MC_THETA);
   fChain->SetBranchAddress("MC_PHI", MC_PHI, &b_MC_PHI);
   fChain->SetBranchAddress("MC_MASS", MC_MASS, &b_MC_MASS);
   fChain->SetBranchAddress("MC_PDGID", MC_PDGID, &b_MC_PDGID);
   fChain->SetBranchAddress("MC_LEPT_PT", MC_LEPT_PT, &b_MC_LEPT_PT);
   fChain->SetBranchAddress("MC_LEPT_ETA", MC_LEPT_ETA, &b_MC_LEPT_ETA);
   fChain->SetBranchAddress("MC_LEPT_PHI", MC_LEPT_PHI, &b_MC_LEPT_PHI);
   fChain->SetBranchAddress("MC_LEPT_THETA", MC_LEPT_THETA, &b_MC_LEPT_THETA);
   fChain->SetBranchAddress("MC_LEPT_PDGID", MC_LEPT_PDGID, &b_MC_LEPT_PDGID);
   fChain->SetBranchAddress("MC_Z_PT", MC_Z_PT, &b_MC_Z_PT);
   fChain->SetBranchAddress("MC_Z_ETA", MC_Z_ETA, &b_MC_Z_ETA);
   fChain->SetBranchAddress("MC_Z_PHI", MC_Z_PHI, &b_MC_Z_PHI);
   fChain->SetBranchAddress("MC_Z_THETA", MC_Z_THETA, &b_MC_Z_THETA);
   fChain->SetBranchAddress("MC_Z_MASS", MC_Z_MASS, &b_MC_Z_MASS);
   fChain->SetBranchAddress("MC_Z_PDGID", MC_Z_PDGID, &b_MC_Z_PDGID);
   fChain->SetBranchAddress("MC_fourl_MASS", MC_fourl_MASS, &b_MC_fourl_MASS);
   fChain->SetBranchAddress("MC_fourl_PT", MC_fourl_PT, &b_MC_fourl_PT);
   fChain->SetBranchAddress("MC_fourl_PDGID", MC_fourl_PDGID, &b_MC_fourl_PDGID);
   fChain->SetBranchAddress("MC_ZZ_MASS", MC_ZZ_MASS, &b_MC_ZZ_MASS);
   fChain->SetBranchAddress("MC_ZZ_PT", MC_ZZ_PT, &b_MC_ZZ_PT);
   fChain->SetBranchAddress("MC_ZZ_ETA", MC_ZZ_ETA, &b_MC_ZZ_ETA);
   fChain->SetBranchAddress("MC_ZZ_PHI", MC_ZZ_PHI, &b_MC_ZZ_PHI);
   fChain->SetBranchAddress("MC_ZZ_THETA", MC_ZZ_THETA, &b_MC_ZZ_THETA);
   fChain->SetBranchAddress("MC_ZZ_PDGID", MC_ZZ_PDGID, &b_MC_ZZ_PDGID);
   fChain->SetBranchAddress("MC_GENJET_PT", MC_GENJET_PT, &b_MC_GENJET_PT);
   fChain->SetBranchAddress("MC_GENJET_ETA", MC_GENJET_ETA, &b_MC_GENJET_ETA);
   fChain->SetBranchAddress("MC_GENJET_PHI", MC_GENJET_PHI, &b_MC_GENJET_PHI);       
   fChain->SetBranchAddress("MC_GENMET", &MC_GENMET, &b_MC_GENMET);
   fChain->SetBranchAddress("RECORF_2e2mu_cosTheta1_spin", RECORF_2e2mu_cosTheta1_spin, &b_RECORF_2e2mu_cosTheta1_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_cosTheta2_spin", RECORF_2e2mu_cosTheta2_spin, &b_RECORF_2e2mu_cosTheta2_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_cosThetaStar_spin", RECORF_2e2mu_cosThetaStar_spin, &b_RECORF_2e2mu_cosThetaStar_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_Phi_spin", RECORF_2e2mu_Phi_spin, &b_RECORF_2e2mu_Phi_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_Phi1_spin", RECORF_2e2mu_Phi1_spin, &b_RECORF_2e2mu_Phi1_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_Phi2_spin", RECORF_2e2mu_Phi2_spin, &b_RECORF_2e2mu_Phi2_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_phi1RF_spin", RECORF_2e2mu_phi1RF_spin, &b_RECORF_2e2mu_phi1RF_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_phi2RF_spin", RECORF_2e2mu_phi2RF_spin, &b_RECORF_2e2mu_phi2RF_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_MELA", RECORF_2e2mu_MELA, &b_RECORF_2e2mu_MELA);
   fChain->SetBranchAddress("RECORF_4e_cosTheta1_spin", RECORF_4e_cosTheta1_spin, &b_RECORF_4e_cosTheta1_spin);
   fChain->SetBranchAddress("RECORF_4e_cosTheta2_spin", RECORF_4e_cosTheta2_spin, &b_RECORF_4e_cosTheta2_spin);
   fChain->SetBranchAddress("RECORF_4e_cosThetaStar_spin", RECORF_4e_cosThetaStar_spin, &b_RECORF_4e_cosThetaStar_spin);
   fChain->SetBranchAddress("RECORF_4e_Phi_spin", RECORF_4e_Phi_spin, &b_RECORF_4e_Phi_spin);
   fChain->SetBranchAddress("RECORF_4e_Phi1_spin", RECORF_4e_Phi1_spin, &b_RECORF_4e_Phi1_spin);
   fChain->SetBranchAddress("RECORF_4e_Phi2_spin", RECORF_4e_Phi2_spin, &b_RECORF_4e_Phi2_spin);
   fChain->SetBranchAddress("RECORF_4e_phi1RF_spin", RECORF_4e_phi1RF_spin, &b_RECORF_4e_phi1RF_spin);
   fChain->SetBranchAddress("RECORF_4e_phi2RF_spin", RECORF_4e_phi2RF_spin, &b_RECORF_4e_phi2RF_spin);
   fChain->SetBranchAddress("RECORF_4e_MELA", RECORF_4e_MELA, &b_RECORF_4e_MELA);
   fChain->SetBranchAddress("RECORF_4mu_cosTheta1_spin", RECORF_4mu_cosTheta1_spin, &b_RECORF_4mu_cosTheta1_spin);
   fChain->SetBranchAddress("RECORF_4mu_cosTheta2_spin", RECORF_4mu_cosTheta2_spin, &b_RECORF_4mu_cosTheta2_spin);
   fChain->SetBranchAddress("RECORF_4mu_cosThetaStar_spin", RECORF_4mu_cosThetaStar_spin, &b_RECORF_4mu_cosThetaStar_spin);
   fChain->SetBranchAddress("RECORF_4mu_Phi_spin", RECORF_4mu_Phi_spin, &b_RECORF_4mu_Phi_spin);
   fChain->SetBranchAddress("RECORF_4mu_Phi1_spin", RECORF_4mu_Phi1_spin, &b_RECORF_4mu_Phi1_spin);
   fChain->SetBranchAddress("RECORF_4mu_Phi2_spin", RECORF_4mu_Phi2_spin, &b_RECORF_4mu_Phi2_spin);
   fChain->SetBranchAddress("RECORF_4mu_phi1RF_spin", RECORF_4mu_phi1RF_spin, &b_RECORF_4mu_phi1RF_spin);
   fChain->SetBranchAddress("RECORF_4mu_phi2RF_spin", RECORF_4mu_phi2RF_spin, &b_RECORF_4mu_phi2RF_spin);
   fChain->SetBranchAddress("RECORF_4mu_MELA", RECORF_4mu_MELA, &b_RECORF_4mu_MELA);
   fChain->SetBranchAddress("RECO_ZMM_MASS", RECO_ZMM_MASS, &b_RECO_ZMM_MASS);
   fChain->SetBranchAddress("RECO_ZEE_MASS", RECO_ZEE_MASS, &b_RECO_ZEE_MASS);
   fChain->SetBranchAddress("RECO_DiLep_MASS", RECO_DiLep_MASS, &b_RECO_DiLep_MASS);
   fChain->SetBranchAddress("RECO_ZMM_PT", RECO_ZMM_PT, &b_RECO_ZMM_PT);
   fChain->SetBranchAddress("RECO_ZEE_PT", RECO_ZEE_PT, &b_RECO_ZEE_PT);
   fChain->SetBranchAddress("RECO_DiLep_PT", RECO_DiLep_PT, &b_RECO_DiLep_PT);
   fChain->SetBranchAddress("RECO_ZMM_ETA", RECO_ZMM_ETA, &b_RECO_ZMM_ETA);
   fChain->SetBranchAddress("RECO_ZEE_ETA", RECO_ZEE_ETA, &b_RECO_ZEE_ETA);
   fChain->SetBranchAddress("RECO_DiLep_ETA", RECO_DiLep_ETA, &b_RECO_DiLep_ETA);
   fChain->SetBranchAddress("RECO_ZMM_PHI", RECO_ZMM_PHI, &b_RECO_ZMM_PHI);
   fChain->SetBranchAddress("RECO_ZEE_PHI", RECO_ZEE_PHI, &b_RECO_ZEE_PHI);
   fChain->SetBranchAddress("RECO_DiLep_PHI", RECO_DiLep_PHI, &b_RECO_DiLep_PHI);
   fChain->SetBranchAddress("RECO_ZMMss_MASS", RECO_ZMMss_MASS, &b_RECO_ZMMss_MASS);
   fChain->SetBranchAddress("RECO_ZEEss_MASS", RECO_ZEEss_MASS, &b_RECO_ZEEss_MASS);
   fChain->SetBranchAddress("RECO_ZEM_MASS", RECO_ZEM_MASS, &b_RECO_ZEM_MASS);
   fChain->SetBranchAddress("RECO_ZMMss_PT", RECO_ZMMss_PT, &b_RECO_ZMMss_PT);
   fChain->SetBranchAddress("RECO_ZEEss_PT", RECO_ZEEss_PT, &b_RECO_ZEEss_PT);
   fChain->SetBranchAddress("RECO_ZEM_PT", RECO_ZEM_PT, &b_RECO_ZEM_PT);
   fChain->SetBranchAddress("RECO_ZMMss_ETA", RECO_ZMMss_ETA, &b_RECO_ZMMss_ETA);
   fChain->SetBranchAddress("RECO_ZEEss_ETA", RECO_ZEEss_ETA, &b_RECO_ZEEss_ETA);
   fChain->SetBranchAddress("RECO_ZEM_ETA", RECO_ZEM_ETA, &b_RECO_ZEM_ETA);
   fChain->SetBranchAddress("RECO_ZMMss_PHI", RECO_ZMMss_PHI, &b_RECO_ZMMss_PHI);
   fChain->SetBranchAddress("RECO_ZEEss_PHI", RECO_ZEEss_PHI, &b_RECO_ZEEss_PHI);
   fChain->SetBranchAddress("RECO_ZEM_PHI", RECO_ZEM_PHI, &b_RECO_ZEM_PHI);
   fChain->SetBranchAddress("RECO_MMMM_MASS", RECO_MMMM_MASS, &b_RECO_MMMM_MASS);
   fChain->SetBranchAddress("RECO_MMMM_PT", RECO_MMMM_PT, &b_RECO_MMMM_PT);
   fChain->SetBranchAddress("RECO_MMMM_ETA", RECO_MMMM_ETA, &b_RECO_MMMM_ETA);
   fChain->SetBranchAddress("RECO_MMMM_PHI", RECO_MMMM_PHI, &b_RECO_MMMM_PHI);
   fChain->SetBranchAddress("RECO_MMMM_MASS_REFIT", RECO_MMMM_MASS_REFIT, &b_RECO_MMMM_MASS_REFIT);
   fChain->SetBranchAddress("RECO_EEEE_MASS", RECO_EEEE_MASS, &b_RECO_EEEE_MASS);
   fChain->SetBranchAddress("RECO_EEEE_PT", RECO_EEEE_PT, &b_RECO_EEEE_PT);
   fChain->SetBranchAddress("RECO_EEEE_ETA", RECO_EEEE_ETA, &b_RECO_EEEE_ETA);
   fChain->SetBranchAddress("RECO_EEEE_PHI", RECO_EEEE_PHI, &b_RECO_EEEE_PHI);
   fChain->SetBranchAddress("RECO_EEEE_MASS_REFIT", RECO_EEEE_MASS_REFIT, &b_RECO_EEEE_MASS_REFIT);
   fChain->SetBranchAddress("RECO_EEMM_MASS", RECO_EEMM_MASS, &b_RECO_EEMM_MASS);
   fChain->SetBranchAddress("RECO_EEMM_PT", RECO_EEMM_PT, &b_RECO_EEMM_PT);
   fChain->SetBranchAddress("RECO_EEMM_ETA", RECO_EEMM_ETA, &b_RECO_EEMM_ETA);
   fChain->SetBranchAddress("RECO_EEMM_PHI", RECO_EEMM_PHI, &b_RECO_EEMM_PHI);
   fChain->SetBranchAddress("RECO_EEMM_MASS_REFIT", RECO_EEMM_MASS_REFIT, &b_RECO_EEMM_MASS_REFIT);
   fChain->SetBranchAddress("RECO_LLL0_MASS", RECO_LLL0_MASS, &b_RECO_LLL0_MASS);
   fChain->SetBranchAddress("RECO_LLL1_MASS", RECO_LLL1_MASS, &b_RECO_LLL1_MASS);
   fChain->SetBranchAddress("RECO_LLL2_MASS", RECO_LLL2_MASS, &b_RECO_LLL2_MASS);
   fChain->SetBranchAddress("RECO_LLL3_MASS", RECO_LLL3_MASS, &b_RECO_LLL3_MASS);
   fChain->SetBranchAddress("RECO_LLL0_PT", RECO_LLL0_PT, &b_RECO_LLL0_PT);
   fChain->SetBranchAddress("RECO_LLL1_PT", RECO_LLL1_PT, &b_RECO_LLL1_PT);
   fChain->SetBranchAddress("RECO_LLL2_PT", RECO_LLL2_PT, &b_RECO_LLL2_PT);
   fChain->SetBranchAddress("RECO_LLL3_PT", RECO_LLL3_PT, &b_RECO_LLL3_PT);
   fChain->SetBranchAddress("RECO_LLLl0_MASS", RECO_LLLl0_MASS, &b_RECO_LLLl0_MASS);
   fChain->SetBranchAddress("RECO_LLLl1_MASS", RECO_LLLl1_MASS, &b_RECO_LLLl1_MASS);
   fChain->SetBranchAddress("RECO_LLLl0_PT", RECO_LLLl0_PT, &b_RECO_LLLl0_PT);
   fChain->SetBranchAddress("RECO_LLLl1_PT", RECO_LLLl1_PT, &b_RECO_LLLl1_PT);
   fChain->SetBranchAddress("RECO_LLLL0ss_MASS", RECO_LLLL0ss_MASS, &b_RECO_LLLL0ss_MASS);
   fChain->SetBranchAddress("RECO_LLLL0ss_PT", RECO_LLLL0ss_PT, &b_RECO_LLLL0ss_PT);
   fChain->SetBranchAddress("RECO_LLLL1ss_MASS", RECO_LLLL1ss_MASS, &b_RECO_LLLL1ss_MASS);
   fChain->SetBranchAddress("RECO_LLLL1ss_PT", RECO_LLLL1ss_PT, &b_RECO_LLLL1ss_PT);
   fChain->SetBranchAddress("RECO_LLLL2ss_MASS", RECO_LLLL2ss_MASS, &b_RECO_LLLL2ss_MASS);
   fChain->SetBranchAddress("RECO_LLLL2ss_PT", RECO_LLLL2ss_PT, &b_RECO_LLLL2ss_PT);
   fChain->SetBranchAddress("RECO_LLLL_MASS", RECO_LLLL_MASS, &b_RECO_LLLL_MASS);
   fChain->SetBranchAddress("RECO_LLLL_PT", RECO_LLLL_PT, &b_RECO_LLLL_PT);
   fChain->SetBranchAddress("RECO_LLLL_ETA", RECO_LLLL_ETA, &b_RECO_LLLL_ETA);
   fChain->SetBranchAddress("RECO_LLLL_PHI", RECO_LLLL_PHI, &b_RECO_LLLL_PHI);
   fChain->SetBranchAddress("RECOELE_E", RECOELE_E, &b_RECOELE_E);
   fChain->SetBranchAddress("RECOELE_PT", RECOELE_PT, &b_RECOELE_PT);
   fChain->SetBranchAddress("RECOELE_PTError", RECOELE_PTError, &b_RECOELE_PTError);
   fChain->SetBranchAddress("RECOELE_P", RECOELE_P, &b_RECOELE_P);
   fChain->SetBranchAddress("RECOELE_ETA", RECOELE_ETA, &b_RECOELE_ETA);
   fChain->SetBranchAddress("RECOELE_THETA", RECOELE_THETA, &b_RECOELE_THETA);
   fChain->SetBranchAddress("RECOELE_PHI", RECOELE_PHI, &b_RECOELE_PHI);
   fChain->SetBranchAddress("RECOELE_MASS", RECOELE_MASS, &b_RECOELE_MASS);
   fChain->SetBranchAddress("RECOELE_CHARGE", RECOELE_CHARGE, &b_RECOELE_CHARGE);
   fChain->SetBranchAddress("RECOELE_isEcalDriven", RECOELE_isEcalDriven, &b_RECOELE_isEcalDriven);
   fChain->SetBranchAddress("RECOELE_isTrackerDriven", RECOELE_isTrackerDriven, &b_RECOELE_isTrackerDriven);
   fChain->SetBranchAddress("RECOELE_gsftrack_NPixHits", RECOELE_gsftrack_NPixHits, &b_RECOELE_gsftrack_NPixHits);
   fChain->SetBranchAddress("RECOELE_gsftrack_NStripHits", RECOELE_gsftrack_NStripHits, &b_RECOELE_gsftrack_NStripHits);
   fChain->SetBranchAddress("RECOELE_gsftrack_chi2", RECOELE_gsftrack_chi2, &b_RECOELE_gsftrack_chi2);
   fChain->SetBranchAddress("RECOELE_gsftrack_dxyB", RECOELE_gsftrack_dxyB, &b_RECOELE_gsftrack_dxyB);
   fChain->SetBranchAddress("RECOELE_gsftrack_dxy", RECOELE_gsftrack_dxy, &b_RECOELE_gsftrack_dxy);
   fChain->SetBranchAddress("RECOELE_gsftrack_dxyError", RECOELE_gsftrack_dxyError, &b_RECOELE_gsftrack_dxyError);
   fChain->SetBranchAddress("RECOELE_gsftrack_dzB", RECOELE_gsftrack_dzB, &b_RECOELE_gsftrack_dzB);
   fChain->SetBranchAddress("RECOELE_gsftrack_dz", RECOELE_gsftrack_dz, &b_RECOELE_gsftrack_dz);
   fChain->SetBranchAddress("RECOELE_gsftrack_dzError", RECOELE_gsftrack_dzError, &b_RECOELE_gsftrack_dzError);
   fChain->SetBranchAddress("RECOELE_gsftrack_losthits", RECOELE_gsftrack_losthits, &b_RECOELE_gsftrack_losthits);
   fChain->SetBranchAddress("RECOELE_gsftrack_validhits", RECOELE_gsftrack_validhits, &b_RECOELE_gsftrack_validhits);
   fChain->SetBranchAddress("RECOELE_gsftrack_expected_inner_hits", RECOELE_gsftrack_expected_inner_hits, &b_RECOELE_gsftrack_expected_inner_hits);
   fChain->SetBranchAddress("RECOELE_scl_E", RECOELE_scl_E, &b_RECOELE_scl_E);
   fChain->SetBranchAddress("RECOELE_scl_Et", RECOELE_scl_Et, &b_RECOELE_scl_Et);
   fChain->SetBranchAddress("RECOELE_scl_Eta", RECOELE_scl_Eta, &b_RECOELE_scl_Eta);
   fChain->SetBranchAddress("RECOELE_scl_Phi", RECOELE_scl_Phi, &b_RECOELE_scl_Phi);
   fChain->SetBranchAddress("RECOELE_ep", RECOELE_ep, &b_RECOELE_ep);
   fChain->SetBranchAddress("RECOELE_eSeedp", RECOELE_eSeedp, &b_RECOELE_eSeedp);
   fChain->SetBranchAddress("RECOELE_eSeedpout", RECOELE_eSeedpout, &b_RECOELE_eSeedpout);
   fChain->SetBranchAddress("RECOELE_eElepout", RECOELE_eElepout, &b_RECOELE_eElepout);
   fChain->SetBranchAddress("RECOELE_deltaEtaIn", RECOELE_deltaEtaIn, &b_RECOELE_deltaEtaIn);
   fChain->SetBranchAddress("RECOELE_deltaEtaSeed", RECOELE_deltaEtaSeed, &b_RECOELE_deltaEtaSeed);
   fChain->SetBranchAddress("RECOELE_deltaEtaEle", RECOELE_deltaEtaEle, &b_RECOELE_deltaEtaEle);
   fChain->SetBranchAddress("RECOELE_deltaPhiIn", RECOELE_deltaPhiIn, &b_RECOELE_deltaPhiIn);
   fChain->SetBranchAddress("RECOELE_deltaPhiSeed", RECOELE_deltaPhiSeed, &b_RECOELE_deltaPhiSeed);
   fChain->SetBranchAddress("RECOELE_deltaPhiEle", RECOELE_deltaPhiEle, &b_RECOELE_deltaPhiEle);
   fChain->SetBranchAddress("RECOELE_isbarrel", RECOELE_isbarrel, &b_RECOELE_isbarrel);
   fChain->SetBranchAddress("RECOELE_isendcap", RECOELE_isendcap, &b_RECOELE_isendcap);
   fChain->SetBranchAddress("RECOELE_isEBetaGap", RECOELE_isEBetaGap, &b_RECOELE_isEBetaGap);
   fChain->SetBranchAddress("RECOELE_isEBphiGap", RECOELE_isEBphiGap, &b_RECOELE_isEBphiGap);
   fChain->SetBranchAddress("RECOELE_isEEdeeGap", RECOELE_isEEdeeGap, &b_RECOELE_isEEdeeGap);
   fChain->SetBranchAddress("RECOELE_isEEringGap", RECOELE_isEEringGap, &b_RECOELE_isEEringGap);
   fChain->SetBranchAddress("RECOELE_sigmaIetaIeta", RECOELE_sigmaIetaIeta, &b_RECOELE_sigmaIetaIeta);
   fChain->SetBranchAddress("RECOELE_sigmaEtaEta", RECOELE_sigmaEtaEta, &b_RECOELE_sigmaEtaEta);
   fChain->SetBranchAddress("RECOELE_e15", RECOELE_e15, &b_RECOELE_e15);
   fChain->SetBranchAddress("RECOELE_e25max", RECOELE_e25max, &b_RECOELE_e25max);
   fChain->SetBranchAddress("RECOELE_e55", RECOELE_e55, &b_RECOELE_e55);
   fChain->SetBranchAddress("RECOELE_he", RECOELE_he, &b_RECOELE_he);
   fChain->SetBranchAddress("RECOELE_r9", RECOELE_r9, &b_RECOELE_r9);
   fChain->SetBranchAddress("RECOELE_mva", RECOELE_mva, &b_RECOELE_mva);
   fChain->SetBranchAddress("RECOELE_fbrem", RECOELE_fbrem, &b_RECOELE_fbrem);
   fChain->SetBranchAddress("RECOELE_nbrems", RECOELE_nbrems, &b_RECOELE_nbrems);
   fChain->SetBranchAddress("RECOELE_Class", RECOELE_Class, &b_RECOELE_Class);
   fChain->SetBranchAddress("RECOELE_fbrem_mode", RECOELE_fbrem_mode, &b_RECOELE_fbrem_mode);
   fChain->SetBranchAddress("RECOELE_fbrem_mean", RECOELE_fbrem_mean, &b_RECOELE_fbrem_mean);
   fChain->SetBranchAddress("RECOELE_EGMTRACKISO", RECOELE_EGMTRACKISO, &b_RECOELE_EGMTRACKISO);
   fChain->SetBranchAddress("RECOELE_EGMHCALISO", RECOELE_EGMHCALISO, &b_RECOELE_EGMHCALISO);
   fChain->SetBranchAddress("RECOELE_EGMECALISO", RECOELE_EGMECALISO, &b_RECOELE_EGMECALISO);
   fChain->SetBranchAddress("RECOELE_EGMX", RECOELE_EGMX, &b_RECOELE_EGMX);
   fChain->SetBranchAddress("RECOELE_PFchAllPart", RECOELE_PFchAllPart, &b_RECOELE_PFchAllPart);
   fChain->SetBranchAddress("RECOELE_PFchHad", RECOELE_PFchHad, &b_RECOELE_PFchHad);
   fChain->SetBranchAddress("RECOELE_PFneuHad", RECOELE_PFneuHad, &b_RECOELE_PFneuHad);
   fChain->SetBranchAddress("RECOELE_PFphoton", RECOELE_PFphoton, &b_RECOELE_PFphoton);
   fChain->SetBranchAddress("RECOELE_PFPUchAllPart", RECOELE_PFPUchAllPart, &b_RECOELE_PFPUchAllPart);
   fChain->SetBranchAddress("RECOELE_PFX_dB", RECOELE_PFX_dB, &b_RECOELE_PFX_dB);
   fChain->SetBranchAddress("RECOELE_PFX_rho", RECOELE_PFX_rho, &b_RECOELE_PFX_rho);
   fChain->SetBranchAddress("RECOELE_regEnergy", RECOELE_regEnergy, &b_RECOELE_regEnergy);
   fChain->SetBranchAddress("RECOELE_regEnergyError", RECOELE_regEnergyError, &b_RECOELE_regEnergyError);
   fChain->SetBranchAddress("RECOELE_SIP", RECOELE_SIP, &b_RECOELE_SIP);
   fChain->SetBranchAddress("RECOELE_IP", RECOELE_IP, &b_RECOELE_IP);
   fChain->SetBranchAddress("RECOELE_IPERROR", RECOELE_IPERROR, &b_RECOELE_IPERROR);
   fChain->SetBranchAddress("RECOELE_SIP_KF", RECOELE_SIP_KF, &b_RECOELE_SIP_KF);
   fChain->SetBranchAddress("RECOELE_IP_KF", RECOELE_IP_KF, &b_RECOELE_IP_KF);
   fChain->SetBranchAddress("RECOELE_IPERROR_KF", RECOELE_IPERROR_KF, &b_RECOELE_IPERROR_KF);
   fChain->SetBranchAddress("RECOELE_SIP_GD", RECOELE_SIP_GD, &b_RECOELE_SIP_GD);
   fChain->SetBranchAddress("RECOELE_SIP_GDEEEE", RECOELE_SIP_GDEEEE, &b_RECOELE_SIP_GDEEEE);
   fChain->SetBranchAddress("RECOELE_SIP_Std", RECOELE_SIP_Std, &b_RECOELE_SIP_Std);
   fChain->SetBranchAddress("RECOELE_SIP_StdEEEE", RECOELE_SIP_StdEEEE, &b_RECOELE_SIP_StdEEEE);
   fChain->SetBranchAddress("RECOELE_SIP_Kin", RECOELE_SIP_Kin, &b_RECOELE_SIP_Kin);
   fChain->SetBranchAddress("RECOELE_SIP_KinEEEE", RECOELE_SIP_KinEEEE, &b_RECOELE_SIP_KinEEEE);
   fChain->SetBranchAddress("RECOELE_STIP", RECOELE_STIP, &b_RECOELE_STIP);
   fChain->SetBranchAddress("RECOELE_SLIP", RECOELE_SLIP, &b_RECOELE_SLIP);
   fChain->SetBranchAddress("RECOELE_TIP", RECOELE_TIP, &b_RECOELE_TIP);
   fChain->SetBranchAddress("RECOELE_LIP", RECOELE_LIP, &b_RECOELE_LIP);
   fChain->SetBranchAddress("RECOELE_TIPERROR", RECOELE_TIPERROR, &b_RECOELE_TIPERROR);
   fChain->SetBranchAddress("RECOELE_LIPERROR", RECOELE_LIPERROR, &b_RECOELE_LIPERROR);
   fChain->SetBranchAddress("RECOELE_sclRawE", RECOELE_sclRawE, &b_RECOELE_sclRawE);
   fChain->SetBranchAddress("RECOELE_sclX", RECOELE_sclX, &b_RECOELE_sclX);
   fChain->SetBranchAddress("RECOELE_sclY", RECOELE_sclY, &b_RECOELE_sclY);
   fChain->SetBranchAddress("RECOELE_sclZ", RECOELE_sclZ, &b_RECOELE_sclZ);
   fChain->SetBranchAddress("RECOELE_seedSubdet1", RECOELE_seedSubdet1, &b_RECOELE_seedSubdet1);
   fChain->SetBranchAddress("RECOELE_seedDphi1", RECOELE_seedDphi1, &b_RECOELE_seedDphi);
   fChain->SetBranchAddress("RECOELE_seedDrz1", RECOELE_seedDrz1, &b_RECOELE_seedDrz1);
   fChain->SetBranchAddress("RECOELE_seedSubdet2", RECOELE_seedSubdet2, &b_RECOELE_seedSubdet2);
   fChain->SetBranchAddress("RECOELE_seedDphi2", RECOELE_seedDphi2, &b_RECOELE_seedDphi2);
   fChain->SetBranchAddress("RECOELE_seedDrz2", RECOELE_seedDrz2, &b_RECOELE_seedDrz2);
   fChain->SetBranchAddress("RECOELE_eidVeryLoose", RECOELE_eidVeryLoose, &b_RECOELE_eidVeryLoose);
   fChain->SetBranchAddress("RECOELE_eidLoose", RECOELE_eidLoose, &b_RECOELE_eidLoose);
   fChain->SetBranchAddress("RECOELE_eidMedium", RECOELE_eidMedium, &b_RECOELE_eidMedium);
   fChain->SetBranchAddress("RECOELE_eidTight", RECOELE_eidTight, &b_RECOELE_eidTight);
   fChain->SetBranchAddress("RECOELE_eidHZZVeryLoose", RECOELE_eidHZZVeryLoose, &b_RECOELE_eidHZZVeryLoose);
   fChain->SetBranchAddress("RECOELE_eidHZZLoose", RECOELE_eidHZZLoose, &b_RECOELE_eidHZZLoose);
   fChain->SetBranchAddress("RECOELE_eidHZZMedium", RECOELE_eidHZZMedium, &b_RECOELE_eidHZZMedium);
   fChain->SetBranchAddress("RECOELE_eidHZZTight", RECOELE_eidHZZTight, &b_RECOELE_eidHZZTight);
   fChain->SetBranchAddress("RECOELE_mvaTrigV0", RECOELE_mvaTrigV0, &b_RECOELE_mvaTrigV0);
   fChain->SetBranchAddress("RECOELE_mvaNonTrigV0", RECOELE_mvaNonTrigV0, &b_RECOELE_mvaNonTrigV0);
   fChain->SetBranchAddress("RECOELE_ID", RECOELE_ID, &b_RECOELE_ID);
   fChain->SetBranchAddress("RECOELE_COV", RECOELE_COV, &b_RECOELE_COV);
   fChain->SetBranchAddress("RECOMU_isPFMu", RECOMU_isPFMu, &b_RECOMU_isPFMu);
   fChain->SetBranchAddress("RECOMU_isGlobalMu", RECOMU_isGlobalMu, &b_RECOMU_isGlobalMu);
   fChain->SetBranchAddress("RECOMU_isStandAloneMu", RECOMU_isStandAloneMu, &b_RECOMU_isStandAloneMu);
   fChain->SetBranchAddress("RECOMU_isTrackerMu", RECOMU_isTrackerMu, &b_RECOMU_isTrackerMu);
   fChain->SetBranchAddress("RECOMU_isCaloMu", RECOMU_isCaloMu, &b_RECOMU_isCaloMu);
   fChain->SetBranchAddress("RECOMU_isTrackerHighPtMu", RECOMU_isTrackerHighPtMu, &b_RECOMU_isTrackerHighPtMu);
   fChain->SetBranchAddress("RECOMU_E", RECOMU_E, &b_RECOMU_E);
   fChain->SetBranchAddress("RECOMU_PT", RECOMU_PT, &b_RECOMU_PT);
   fChain->SetBranchAddress("RECOMU_P", RECOMU_P, &b_RECOMU_P);
   fChain->SetBranchAddress("RECOMU_ETA", RECOMU_ETA, &b_RECOMU_ETA);
   fChain->SetBranchAddress("RECOMU_THETA", RECOMU_THETA, &b_RECOMU_THETA);
   fChain->SetBranchAddress("RECOMU_PHI", RECOMU_PHI, &b_RECOMU_PHI);
   fChain->SetBranchAddress("RECOMU_MASS", RECOMU_MASS, &b_RECOMU_MASS);
   fChain->SetBranchAddress("RECOMU_CHARGE", RECOMU_CHARGE, &b_RECOMU_CHARGE);
   fChain->SetBranchAddress("RECOMU_COV", RECOMU_COV, &b_RECOMU_COV);
   fChain->SetBranchAddress("RECOMU_TRACKISO", RECOMU_TRACKISO, &b_RECOMU_TRACKISO);
   fChain->SetBranchAddress("RECOMU_TRACKISO_SUMPT", RECOMU_TRACKISO_SUMPT, &b_RECOMU_TRACKISO_SUMPT);
   fChain->SetBranchAddress("RECOMU_HCALISO", RECOMU_HCALISO, &b_RECOMU_HCALISO);
   fChain->SetBranchAddress("RECOMU_ECALISO", RECOMU_ECALISO, &b_RECOMU_ECALISO);
   fChain->SetBranchAddress("RECOMU_X", RECOMU_X, &b_RECOMU_X);
   fChain->SetBranchAddress("RECOMU_PFchHad", RECOMU_PFchHad, &b_RECOMU_PFchHad);
   fChain->SetBranchAddress("RECOMU_PFneuHad", RECOMU_PFneuHad, &b_RECOMU_PFneuHad);
   fChain->SetBranchAddress("RECOMU_PFphoton", RECOMU_PFphoton, &b_RECOMU_PFphoton);
   fChain->SetBranchAddress("RECOMU_PFPUchAllPart", RECOMU_PFPUchAllPart, &b_RECOMU_PFPUchAllPart);
   fChain->SetBranchAddress("RECOMU_PFX_dB", RECOMU_PFX_dB, &b_RECOMU_PFX_dB);
   fChain->SetBranchAddress("RECOMU_PFX_rho", RECOMU_PFX_rho, &b_RECOMU_PFX_rho);
   fChain->SetBranchAddress("RECOPFPHOT_PFchHad", RECOPFPHOT_PFchHad, &b_RECOPFPHOT_PFchHad);
   fChain->SetBranchAddress("RECOPFPHOT_PFneuHad", RECOPFPHOT_PFneuHad, &b_RECOPFPHOT_PFneuHad);
   fChain->SetBranchAddress("RECOPFPHOT_PFphoton", RECOPFPHOT_PFphoton, &b_RECOPFPHOT_PFphoton);
   fChain->SetBranchAddress("RECOPFPHOT_PFPUchAllPart", RECOPFPHOT_PFPUchAllPart, &b_RECOPFPHOT_PFPUchAllPart);
   fChain->SetBranchAddress("RECOPFPHOT_PFX_rho", RECOPFPHOT_PFX_rho, &b_RECOPFPHOT_PFX_rho);
   fChain->SetBranchAddress("RECOMU_SIP", RECOMU_SIP, &b_RECOMU_SIP);
   fChain->SetBranchAddress("RECOMU_IP", RECOMU_IP, &b_RECOMU_IP);
   fChain->SetBranchAddress("RECOMU_IPERROR", RECOMU_IPERROR, &b_RECOMU_IPERROR);
   fChain->SetBranchAddress("RECOMU_SIP_KF", RECOMU_SIP_KF, &b_RECOMU_SIP_KF);
   fChain->SetBranchAddress("RECOMU_IP_KF", RECOMU_IP_KF, &b_RECOMU_IP_KF);
   fChain->SetBranchAddress("RECOMU_IPERROR_KF", RECOMU_IPERROR_KF, &b_RECOMU_IPERROR_KF);
   fChain->SetBranchAddress("RECOMU_SIP_GD", RECOMU_SIP_GD, &b_RECOMU_SIP_GD);
   fChain->SetBranchAddress("RECOMU_SIP_GDMMMM", RECOMU_SIP_GDMMMM, &b_RECOMU_SIP_GDMMMM);
   fChain->SetBranchAddress("RECOMU_SIP_Std", RECOMU_SIP_Std, &b_RECOMU_SIP_Std);
   fChain->SetBranchAddress("RECOMU_SIP_StdMMMM", RECOMU_SIP_StdMMMM, &b_RECOMU_SIP_StdMMMM);
   fChain->SetBranchAddress("RECOMU_SIP_Kin", RECOMU_SIP_Kin, &b_RECOMU_SIP_Kin);
   fChain->SetBranchAddress("RECOMU_SIP_KinMMMM", RECOMU_SIP_KinMMMM, &b_RECOMU_SIP_KinMMMM);
   fChain->SetBranchAddress("RECOMU_STIP", RECOMU_STIP, &b_RECOMU_STIP);
   fChain->SetBranchAddress("RECOMU_SLIP", RECOMU_SLIP, &b_RECOMU_SLIP);
   fChain->SetBranchAddress("RECOMU_TIP", RECOMU_TIP, &b_RECOMU_TIP);
   fChain->SetBranchAddress("RECOMU_LIP", RECOMU_LIP, &b_RECOMU_LIP);
   fChain->SetBranchAddress("RECOMU_TIPERROR", RECOMU_TIPERROR, &b_RECOMU_TIPERROR);
   fChain->SetBranchAddress("RECOMU_LIPERROR", RECOMU_LIPERROR, &b_RECOMU_LIPERROR);
   fChain->SetBranchAddress("RECOMU_caloCompatibility", RECOMU_caloCompatibility, &b_RECOMU_caloCompatibility);
   fChain->SetBranchAddress("RECOMU_segmentCompatibility", RECOMU_segmentCompatibility, &b_RECOMU_segmentCompatibility);
   fChain->SetBranchAddress("RECOMU_numberOfMatches", RECOMU_numberOfMatches, &b_RECOMU_numberOfMatches);
   fChain->SetBranchAddress("RECOMU_numberOfMatchedStations", RECOMU_numberOfMatchedStations, &b_RECOMU_numberOfMatchedStations);
   fChain->SetBranchAddress("RECOMU_glbmuPromptTight", RECOMU_glbmuPromptTight, &b_RECOMU_glbmuPromptTight);
   fChain->SetBranchAddress("RECOMU_trkmuArbitration", RECOMU_trkmuArbitration, &b_RECOMU_trkmuArbitration);
   fChain->SetBranchAddress("RECOMU_trkmu2DCompatibilityLoose", RECOMU_trkmu2DCompatibilityLoose, &b_RECOMU_trkmu2DCompatibilityLoose);
   fChain->SetBranchAddress("RECOMU_trkmu2DCompatibilityTight", RECOMU_trkmu2DCompatibilityTight, &b_RECOMU_trkmu2DCompatibilityTight);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationLoose", RECOMU_trkmuOneStationLoose, &b_RECOMU_trkmuOneStationLoose);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationTight", RECOMU_trkmuOneStationTight, &b_RECOMU_trkmuOneStationTight);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationLoose", RECOMU_trkmuLastStationLoose, &b_RECOMU_trkmuLastStationLoose);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationTight", RECOMU_trkmuLastStationTight, &b_RECOMU_trkmuLastStationTight);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationAngLoose", RECOMU_trkmuOneStationAngLoose, &b_RECOMU_trkmuOneStationAngLoose);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationAngTight", RECOMU_trkmuOneStationAngTight, &b_RECOMU_trkmuOneStationAngTight);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationAngLoose", RECOMU_trkmuLastStationAngLoose, &b_RECOMU_trkmuLastStationAngLoose);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationAngTight", RECOMU_trkmuLastStationAngTight, &b_RECOMU_trkmuLastStationAngTight);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationOptimizedLowPtLoose", RECOMU_trkmuLastStationOptimizedLowPtLoose, &b_RECOMU_trkmuLastStationOptimizedLowPtLoose);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationOptimizedLowPtTight", RECOMU_trkmuLastStationOptimizedLowPtTight, &b_RECOMU_trkmuLastStationOptimizedLowPtTight);
   fChain->SetBranchAddress("RECOMU_mutrkPT", RECOMU_mutrkPT, &b_RECOMU_mutrkPT);
   fChain->SetBranchAddress("RECOMU_mutrkPTError", RECOMU_mutrkPTError, &b_RECOMU_mutrkPTError);
   fChain->SetBranchAddress("RECOMU_mutrkDxy", RECOMU_mutrkDxy, &b_RECOMU_mutrkDxy);
   fChain->SetBranchAddress("RECOMU_mutrkDxyError", RECOMU_mutrkDxyError, &b_RECOMU_mutrkDxyError);
   fChain->SetBranchAddress("RECOMU_mutrkDxyB", RECOMU_mutrkDxyB, &b_RECOMU_mutrkDxyB);
   fChain->SetBranchAddress("RECOMU_mutrkDz", RECOMU_mutrkDz, &b_RECOMU_mutrkDz);
   fChain->SetBranchAddress("RECOMU_mutrkDzError", RECOMU_mutrkDzError, &b_RECOMU_mutrkDzError);
   fChain->SetBranchAddress("RECOMU_mutrkDzB", RECOMU_mutrkDzB, &b_RECOMU_mutrkDzB);
   fChain->SetBranchAddress("RECOMU_mutrkChi2PerNdof", RECOMU_mutrkChi2PerNdof, &b_RECOMU_mutrkChi2PerNdof);
   fChain->SetBranchAddress("RECOMU_mutrkCharge", RECOMU_mutrkCharge, &b_RECOMU_mutrkCharge);
   fChain->SetBranchAddress("RECOMU_mutrkNHits", RECOMU_mutrkNHits, &b_RECOMU_mutrkNHits);
   fChain->SetBranchAddress("RECOMU_mutrkNStripHits", RECOMU_mutrkNStripHits, &b_RECOMU_mutrkNStripHits);
   fChain->SetBranchAddress("RECOMU_mutrkNPixHits", RECOMU_mutrkNPixHits, &b_RECOMU_mutrkNPixHits);
   fChain->SetBranchAddress("RECOMU_mutrkNMuonHits", RECOMU_mutrkNMuonHits, &b_RECOMU_mutrkNMuonHits);
   fChain->SetBranchAddress("RECOMU_mutrktrackerLayersWithMeasurement", RECOMU_mutrktrackerLayersWithMeasurement, &b_RECOMU_mutrktrackerLayersWithMeasurement);
   fChain->SetBranchAddress("RECOMU_muInnertrkDxy", RECOMU_muInnertrkDxy, &b_RECOMU_muInnertrkDxy);
   fChain->SetBranchAddress("RECOMU_muInnertrkDxyError", RECOMU_muInnertrkDxyError, &b_RECOMU_muInnertrkDxyError);
   fChain->SetBranchAddress("RECOMU_muInnertrkDxyB", RECOMU_muInnertrkDxyB, &b_RECOMU_muInnertrkDxyB);
   fChain->SetBranchAddress("RECOMU_muInnertrkDz", RECOMU_muInnertrkDz, &b_RECOMU_muInnertrkDz);
   fChain->SetBranchAddress("RECOMU_muInnertrkDzError", RECOMU_muInnertrkDzError, &b_RECOMU_muInnertrkDzError);
   fChain->SetBranchAddress("RECOMU_muInnertrkDzB", RECOMU_muInnertrkDzB, &b_RECOMU_muInnertrkDzB);
   fChain->SetBranchAddress("RECOMU_muInnertrkChi2PerNdof", RECOMU_muInnertrkChi2PerNdof, &b_RECOMU_muInnertrkChi2PerNdof);
   fChain->SetBranchAddress("RECOMU_muInnertrktrackerLayersWithMeasurement", RECOMU_muInnertrktrackerLayersWithMeasurement, &b_RECOMU_muInnertrktrackerLayersWithMeasurement);
   fChain->SetBranchAddress("RECOMU_muInnertrkPT", RECOMU_muInnertrkPT, &b_RECOMU_muInnertrkPT);
   fChain->SetBranchAddress("RECOMU_muInnertrkPTError", RECOMU_muInnertrkPTError, &b_RECOMU_muInnertrkPTError);
   fChain->SetBranchAddress("RECOMU_muInnertrkCharge", RECOMU_muInnertrkCharge, &b_RECOMU_muInnertrkCharge);
   fChain->SetBranchAddress("RECOMU_muInnertrkNHits", RECOMU_muInnertrkNHits, &b_RECOMU_muInnertrkNHits);
   fChain->SetBranchAddress("RECOMU_muInnertrkNStripHits", RECOMU_muInnertrkNStripHits, &b_RECOMU_muInnertrkNStripHits);
   fChain->SetBranchAddress("RECOMU_muInnertrkNPixHits", RECOMU_muInnertrkNPixHits, &b_RECOMU_muInnertrkNPixHits);
   fChain->SetBranchAddress("RECOMU_mubesttrkType", RECOMU_mubesttrkType, &b_RECOMU_mubesttrkType);
   fChain->SetBranchAddress("RECOMU_mubesttrkDxy", RECOMU_mubesttrkDxy, &b_RECOMU_mubesttrkDxy);
   fChain->SetBranchAddress("RECOMU_mubesttrkDxyError", RECOMU_mubesttrkDxyError, &b_RECOMU_mubesttrkDxyError);
   fChain->SetBranchAddress("RECOMU_mubesttrkDz", RECOMU_mubesttrkDz, &b_RECOMU_mubesttrkDz);
   fChain->SetBranchAddress("RECOMU_mubesttrkDzError", RECOMU_mubesttrkDzError, &b_RECOMU_mubesttrkDzError);
   fChain->SetBranchAddress("RECOMU_mubesttrkPTError", RECOMU_mubesttrkPTError, &b_RECOMU_mubesttrkPTError);
   fChain->SetBranchAddress("ftsigma", ftsigma, &b_ftsigma);
   fChain->SetBranchAddress("gdX", gdX, &b_gdX);
   fChain->SetBranchAddress("gdY", gdY, &b_gdY);
   fChain->SetBranchAddress("gdZ", gdZ, &b_gdZ);
   fChain->SetBranchAddress("ftsigmalag", ftsigmalag, &b_ftsigmalag);
   fChain->SetBranchAddress("gdlagX", gdlagX, &b_gdlagX);
   fChain->SetBranchAddress("gdlagY", gdlagY, &b_gdlagY);
   fChain->SetBranchAddress("gdlagZ", gdlagZ, &b_gdlagZ);
   fChain->SetBranchAddress("gdlagProb", gdlagProb, &b_gdlagProb);
   fChain->SetBranchAddress("gdlagNdof", gdlagNdof, &b_gdlagNdof);
   fChain->SetBranchAddress("ftsigmaMMMM", ftsigmaMMMM, &b_ftsigmaMMMM);
   fChain->SetBranchAddress("gdXMMMM", gdXMMMM, &b_gdXMMMM);
   fChain->SetBranchAddress("gdYMMMM", gdYMMMM, &b_gdYMMMM);
   fChain->SetBranchAddress("gdZMMMM", gdZMMMM, &b_gdZMMMM);
   fChain->SetBranchAddress("ftsigmalagMMMM", ftsigmalagMMMM, &b_ftsigmalagMMMM);
   fChain->SetBranchAddress("gdlagXMMMM", gdlagXMMMM, &b_gdlagXMMMM);
   fChain->SetBranchAddress("gdlagYMMMM", gdlagYMMMM, &b_gdlagYMMMM);
   fChain->SetBranchAddress("gdlagZMMMM", gdlagZMMMM, &b_gdlagZMMMM);
   fChain->SetBranchAddress("gdlagProbMMMM", gdlagProbMMMM, &b_gdlagProbMMMM);
   fChain->SetBranchAddress("gdlagNdofMMMM", gdlagNdofMMMM, &b_gdlagNdofMMMM);
   fChain->SetBranchAddress("ftsigmaEEEE", ftsigmaEEEE, &b_ftsigmaEEEE);
   fChain->SetBranchAddress("gdXEEEE", gdXEEEE, &b_gdXEEEE);
   fChain->SetBranchAddress("gdYEEEE", gdYEEEE, &b_gdYEEEE);
   fChain->SetBranchAddress("gdZEEEE", gdZEEEE, &b_gdZEEEE);
   fChain->SetBranchAddress("ftsigmalagEEEE", ftsigmalagEEEE, &b_ftsigmalagEEEE);
   fChain->SetBranchAddress("gdlagXEEEE", gdlagXEEEE, &b_gdlagXEEEE);
   fChain->SetBranchAddress("gdlagYEEEE", gdlagYEEEE, &b_gdlagYEEEE);
   fChain->SetBranchAddress("gdlagZEEEE", gdlagZEEEE, &b_gdlagZEEEE);
   fChain->SetBranchAddress("gdlagProbEEEE", gdlagProbEEEE, &b_gdlagProbEEEE);
   fChain->SetBranchAddress("gdlagNdofEEEE", gdlagNdofEEEE, &b_gdlagNdofEEEE);
   fChain->SetBranchAddress("StdFitVertexX", StdFitVertexX, &b_StdFitVertexX);
   fChain->SetBranchAddress("StdFitVertexY", StdFitVertexY, &b_StdFitVertexY);
   fChain->SetBranchAddress("StdFitVertexZ", StdFitVertexZ, &b_StdFitVertexZ);
   fChain->SetBranchAddress("StdFitVertexChi2r", StdFitVertexChi2r, &b_StdFitVertexChi2r);
   fChain->SetBranchAddress("StdFitVertexProb", StdFitVertexProb, &b_StdFitVertexProb);
   fChain->SetBranchAddress("StdFitVertexTrack_PT", StdFitVertexTrack_PT, &b_StdFitVertexTrack_PT);
   fChain->SetBranchAddress("StdFitVertexTrack_ETA", StdFitVertexTrack_ETA, &b_StdFitVertexTrack_ETA);
   fChain->SetBranchAddress("StdFitVertexTrack_PHI", StdFitVertexTrack_PHI, &b_StdFitVertexTrack_PHI);
   fChain->SetBranchAddress("KinFitVertexX", KinFitVertexX, &b_KinFitVertexX);
   fChain->SetBranchAddress("KinFitVertexY", KinFitVertexY, &b_KinFitVertexY);
   fChain->SetBranchAddress("KinFitVertexZ", KinFitVertexZ, &b_KinFitVertexZ);
   fChain->SetBranchAddress("KinFitVertexChi2r", KinFitVertexChi2r, &b_KinFitVertexChi2r);
   fChain->SetBranchAddress("KinFitVertexProb", KinFitVertexProb, &b_KinFitVertexProb);
   fChain->SetBranchAddress("StdFitVertexXMMMM", StdFitVertexXMMMM, &b_StdFitVertexXMMMM);
   fChain->SetBranchAddress("StdFitVertexYMMMM", StdFitVertexYMMMM, &b_StdFitVertexYMMMM);
   fChain->SetBranchAddress("StdFitVertexZMMMM", StdFitVertexZMMMM, &b_StdFitVertexZMMMM);
   fChain->SetBranchAddress("StdFitVertexChi2rMMMM", StdFitVertexChi2rMMMM, &b_StdFitVertexChi2rMMMM);
   fChain->SetBranchAddress("StdFitVertexProbMMMM", StdFitVertexProbMMMM, &b_StdFitVertexProbMMMM);
   fChain->SetBranchAddress("StdFitVertexTrackMMMM_PT", StdFitVertexTrackMMMM_PT, &b_StdFitVertexTrackMMMM_PT);
   fChain->SetBranchAddress("StdFitVertexTrackMMMM_ETA", StdFitVertexTrackMMMM_ETA, &b_StdFitVertexTrackMMMM_ETA);
   fChain->SetBranchAddress("StdFitVertexTrackMMMM_PHI", StdFitVertexTrackMMMM_PHI, &b_StdFitVertexTrackMMMM_PHI);
   fChain->SetBranchAddress("KinFitVertexXMMMM", KinFitVertexXMMMM, &b_KinFitVertexXMMMM);
   fChain->SetBranchAddress("KinFitVertexYMMMM", KinFitVertexYMMMM, &b_KinFitVertexYMMMM);
   fChain->SetBranchAddress("KinFitVertexZMMMM", KinFitVertexZMMMM, &b_KinFitVertexZMMMM);
   fChain->SetBranchAddress("KinFitVertexChi2rMMMM", KinFitVertexChi2rMMMM, &b_KinFitVertexChi2rMMMM);
   fChain->SetBranchAddress("KinFitVertexProbMMMM", KinFitVertexProbMMMM, &b_KinFitVertexProbMMMM);
   fChain->SetBranchAddress("StdFitVertexXEEEE", StdFitVertexXEEEE, &b_StdFitVertexXEEEE);
   fChain->SetBranchAddress("StdFitVertexYEEEE", StdFitVertexYEEEE, &b_StdFitVertexYEEEE);
   fChain->SetBranchAddress("StdFitVertexZEEEE", StdFitVertexZEEEE, &b_StdFitVertexZEEEE);
   fChain->SetBranchAddress("StdFitVertexChi2rEEEE", StdFitVertexChi2rEEEE, &b_StdFitVertexChi2rEEEE);
   fChain->SetBranchAddress("StdFitVertexProbEEEE", StdFitVertexProbEEEE, &b_StdFitVertexProbEEEE);
   fChain->SetBranchAddress("StdFitVertexTrackEEEE_PT", StdFitVertexTrackEEEE_PT, &b_StdFitVertexTrackEEEE_PT);
   fChain->SetBranchAddress("StdFitVertexTrackEEEE_ETA", StdFitVertexTrackEEEE_ETA, &b_StdFitVertexTrackEEEE_ETA);
   fChain->SetBranchAddress("StdFitVertexTrackEEEE_PHI", StdFitVertexTrackEEEE_PHI, &b_StdFitVertexTrackEEEE_PHI);
   fChain->SetBranchAddress("KinFitVertexXEEEE", KinFitVertexXEEEE, &b_KinFitVertexXEEEE);
   fChain->SetBranchAddress("KinFitVertexYEEEE", KinFitVertexYEEEE, &b_KinFitVertexYEEEE);
   fChain->SetBranchAddress("KinFitVertexZEEEE", KinFitVertexZEEEE, &b_KinFitVertexZEEEE);
   fChain->SetBranchAddress("KinFitVertexChi2rEEEE", KinFitVertexChi2rEEEE, &b_KinFitVertexChi2rEEEE);
   fChain->SetBranchAddress("KinFitVertexProbEEEE", KinFitVertexProbEEEE, &b_KinFitVertexProbEEEE);
   fChain->SetBranchAddress("StdFitVertexChi2rMMM", StdFitVertexChi2rMMM, &b_StdFitVertexChi2rMMM);
   fChain->SetBranchAddress("StdFitVertexProbMMM", StdFitVertexProbMMM, &b_StdFitVertexProbMMM);
   fChain->SetBranchAddress("StdFitVertexChi2rMME", StdFitVertexChi2rMME, &b_StdFitVertexChi2rMME);
   fChain->SetBranchAddress("StdFitVertexProbMME", StdFitVertexProbMME, &b_StdFitVertexProbMME);
   fChain->SetBranchAddress("StdFitVertexChi2rEEE", StdFitVertexChi2rEEE, &b_StdFitVertexChi2rEEE);
   fChain->SetBranchAddress("StdFitVertexProbEEE", StdFitVertexProbEEE, &b_StdFitVertexProbEEE);
   fChain->SetBranchAddress("StdFitVertexChi2rMEE", StdFitVertexChi2rMEE, &b_StdFitVertexChi2rMEE);
   fChain->SetBranchAddress("StdFitVertexProbMEE", StdFitVertexProbMEE, &b_StdFitVertexProbMEE);
   fChain->SetBranchAddress("StdFitVertexChi2rDiLep", StdFitVertexChi2rDiLep, &b_StdFitVertexChi2rDiLep);
   fChain->SetBranchAddress("StdFitVertexProbDiLep", StdFitVertexProbDiLep, &b_StdFitVertexProbDiLep);
   fChain->SetBranchAddress("ConvMapDist", ConvMapDist, &b_ConvMapDist);
   fChain->SetBranchAddress("ConvMapDcot", ConvMapDcot, &b_ConvMapDcot);
   fChain->SetBranchAddress("RECOMU_MatchingMCTruth", RECOMU_MatchingMCTruth, &b_RECOMU_MatchingMCTruth);
   fChain->SetBranchAddress("RECOMU_MatchingMCpT", RECOMU_MatchingMCpT, &b_RECOMU_MatchingMCpT);
   fChain->SetBranchAddress("RECOMU_MatchingMCEta", RECOMU_MatchingMCEta, &b_RECOMU_MatchingMCEta);
   fChain->SetBranchAddress("RECOMU_MatchingMCPhi", RECOMU_MatchingMCPhi, &b_RECOMU_MatchingMCPhi);
   fChain->SetBranchAddress("RECOELE_MatchingMCTruth", RECOELE_MatchingMCTruth, &b_RECOELE_MatchingMCTruth);
   fChain->SetBranchAddress("RECOELE_MatchingMCpT", RECOELE_MatchingMCpT, &b_RECOELE_MatchingMCpT);
   fChain->SetBranchAddress("RECOELE_MatchingMCEta", RECOELE_MatchingMCEta, &b_RECOELE_MatchingMCEta);
   fChain->SetBranchAddress("RECOELE_MatchingMCPhi", RECOELE_MatchingMCPhi, &b_RECOELE_MatchingMCPhi);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCTruth", RECOPHOT_MatchingMCTruth, &b_RECOPHOT_MatchingMCTruth);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCpT", RECOPHOT_MatchingMCpT, &b_RECOPHOT_MatchingMCpT);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCEta", RECOPHOT_MatchingMCEta, &b_RECOPHOT_MatchingMCEta);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCPhi", RECOPHOT_MatchingMCPhi, &b_RECOPHOT_MatchingMCPhi);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCTruth", RECOzMuMu_MatchingMCTruth, &b_RECOzMuMu_MatchingMCTruth);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCpT", RECOzMuMu_MatchingMCpT, &b_RECOzMuMu_MatchingMCpT);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCmass", RECOzMuMu_MatchingMCmass, &b_RECOzMuMu_MatchingMCmass);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCEta", RECOzMuMu_MatchingMCEta, &b_RECOzMuMu_MatchingMCEta);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCPhi", RECOzMuMu_MatchingMCPhi, &b_RECOzMuMu_MatchingMCPhi);
   fChain->SetBranchAddress("RECOzEE_MatchingMCTruth", RECOzEE_MatchingMCTruth, &b_RECOzEE_MatchingMCTruth);
   fChain->SetBranchAddress("RECOzEE_MatchingMCpT", RECOzEE_MatchingMCpT, &b_RECOzEE_MatchingMCpT);
   fChain->SetBranchAddress("RECOzEE_MatchingMCmass", RECOzEE_MatchingMCmass, &b_RECOzEE_MatchingMCmass);
   fChain->SetBranchAddress("RECOzEE_MatchingMCEta", RECOzEE_MatchingMCEta, &b_RECOzEE_MatchingMCEta);
   fChain->SetBranchAddress("RECOzEE_MatchingMCPhi", RECOzEE_MatchingMCPhi, &b_RECOzEE_MatchingMCPhi);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCTruth", RECOHzzEEEE_MatchingMCTruth, &b_RECOHzzEEEE_MatchingMCTruth);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCpT", RECOHzzEEEE_MatchingMCpT, &b_RECOHzzEEEE_MatchingMCpT);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCmass", RECOHzzEEEE_MatchingMCmass, &b_RECOHzzEEEE_MatchingMCmass);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCEta", RECOHzzEEEE_MatchingMCEta, &b_RECOHzzEEEE_MatchingMCEta);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCPhi", RECOHzzEEEE_MatchingMCPhi, &b_RECOHzzEEEE_MatchingMCPhi);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCTruth", RECOHzzEEMM_MatchingMCTruth, &b_RECOHzzEEMM_MatchingMCTruth);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCpT", RECOHzzEEMM_MatchingMCpT, &b_RECOHzzEEMM_MatchingMCpT);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCmass", RECOHzzEEMM_MatchingMCmass, &b_RECOHzzEEMM_MatchingMCmass);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCEta", RECOHzzEEMM_MatchingMCEta, &b_RECOHzzEEMM_MatchingMCEta);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCPhi", RECOHzzEEMM_MatchingMCPhi, &b_RECOHzzEEMM_MatchingMCPhi);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCTruth", RECOHzzMMMM_MatchingMCTruth, &b_RECOHzzMMMM_MatchingMCTruth);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCpT", RECOHzzMMMM_MatchingMCpT, &b_RECOHzzMMMM_MatchingMCpT);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCmass", RECOHzzMMMM_MatchingMCmass, &b_RECOHzzMMMM_MatchingMCmass);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCEta", RECOHzzMMMM_MatchingMCEta, &b_RECOHzzMMMM_MatchingMCEta);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCPhi", RECOHzzMMMM_MatchingMCPhi, &b_RECOHzzMMMM_MatchingMCPhi);
   fChain->SetBranchAddress("RECO_NMU", &RECO_NMU, &b_RECO_NMU);
   fChain->SetBranchAddress("RECO_NELE", &RECO_NELE, &b_RECO_NELE);
   fChain->SetBranchAddress("RECO_NTRACK", &RECO_NTRACK, &b_RECO_NTRACK);
   fChain->SetBranchAddress("RECO_TRACK_PT", RECO_TRACK_PT, &b_RECO_TRACK_PT);
   fChain->SetBranchAddress("RECO_TRACK_ETA", RECO_TRACK_ETA, &b_RECO_TRACK_ETA);
   fChain->SetBranchAddress("RECO_TRACK_PHI", RECO_TRACK_PHI, &b_RECO_TRACK_PHI);
   fChain->SetBranchAddress("RECO_TRACK_CHI2", RECO_TRACK_CHI2, &b_RECO_TRACK_CHI2);
   fChain->SetBranchAddress("RECO_TRACK_CHI2RED", RECO_TRACK_CHI2RED, &b_RECO_TRACK_CHI2RED);
   fChain->SetBranchAddress("RECO_TRACK_CHI2PROB", RECO_TRACK_CHI2PROB, &b_RECO_TRACK_CHI2PROB);
   fChain->SetBranchAddress("RECO_TRACK_NHITS", RECO_TRACK_NHITS, &b_RECO_TRACK_NHITS);
   fChain->SetBranchAddress("RECO_TRACK_DXY", RECO_TRACK_DXY, &b_RECO_TRACK_DXY);
   fChain->SetBranchAddress("RECO_TRACK_DXYERR", RECO_TRACK_DXYERR, &b_RECO_TRACK_DXYERR);
   fChain->SetBranchAddress("RECO_TRACK_DZ", RECO_TRACK_DZ, &b_RECO_TRACK_DZ);
   fChain->SetBranchAddress("RECO_TRACK_DZERR", RECO_TRACK_DZERR, &b_RECO_TRACK_DZERR);
   fChain->SetBranchAddress("RECO_NPHOT", &RECO_NPHOT, &b_RECO_NPHOT);
   fChain->SetBranchAddress("RECOPHOT_PT", RECOPHOT_PT, &b_RECOPHOT_PT);
   fChain->SetBranchAddress("RECOPHOT_ETA", RECOPHOT_ETA, &b_RECOPHOT_ETA);
   fChain->SetBranchAddress("RECOPHOT_PHI", RECOPHOT_PHI, &b_RECOPHOT_PHI);
   fChain->SetBranchAddress("RECOPHOT_THETA", RECOPHOT_THETA, &b_RECOPHOT_THETA);
   fChain->SetBranchAddress("RECO_NPFPHOT", &RECO_NPFPHOT, &b_RECO_NPFPHOT);
   fChain->SetBranchAddress("RECOPFPHOT_PT", RECOPFPHOT_PT, &b_RECOPFPHOT_PT);
   fChain->SetBranchAddress("RECOPFPHOT_PTError", RECOPFPHOT_PTError, &b_RECOPFPHOT_PTError);
   fChain->SetBranchAddress("RECOPFPHOT_ETA", RECOPFPHOT_ETA, &b_RECOPFPHOT_ETA);
   fChain->SetBranchAddress("RECOPFPHOT_PHI", RECOPFPHOT_PHI, &b_RECOPFPHOT_PHI);
   fChain->SetBranchAddress("RECOPFPHOT_THETA", RECOPFPHOT_THETA, &b_RECOPFPHOT_THETA);
   fChain->SetBranchAddress("BeamSpot_X", &BeamSpot_X, &b_BeamSpot_X);
   fChain->SetBranchAddress("BeamSpot_Y", &BeamSpot_Y, &b_BeamSpot_Y);
   fChain->SetBranchAddress("BeamSpot_Z", &BeamSpot_Z, &b_BeamSpot_Z);
   fChain->SetBranchAddress("RECO_NVTX", &RECO_NVTX, &b_RECO_NVTX);
   fChain->SetBranchAddress("RECO_VERTEX_x", RECO_VERTEX_x, &b_RECO_VERTEX_x);
   fChain->SetBranchAddress("RECO_VERTEX_y", RECO_VERTEX_y, &b_RECO_VERTEX_y);
   fChain->SetBranchAddress("RECO_VERTEX_z", RECO_VERTEX_z, &b_RECO_VERTEX_z);
   fChain->SetBranchAddress("RECO_VERTEX_ndof", RECO_VERTEX_ndof, &b_RECO_VERTEX_ndof);
   fChain->SetBranchAddress("RECO_VERTEX_chi2", RECO_VERTEX_chi2, &b_RECO_VERTEX_chi2);
   fChain->SetBranchAddress("RECO_VERTEX_ntracks", RECO_VERTEX_ntracks, &b_RECO_VERTEX_ntracks);
   fChain->SetBranchAddress("RECO_VERTEXPROB", RECO_VERTEXPROB, &b_RECO_VERTEXPROB);
   fChain->SetBranchAddress("RECO_VERTEX_isValid", RECO_VERTEX_isValid, &b_RECO_VERTEX_isValid);
   fChain->SetBranchAddress("RECO_VERTEX_TRACK_PT", RECO_VERTEX_TRACK_PT, &b_RECO_VERTEX_TRACK_PT);
   fChain->SetBranchAddress("RECO_PFJET_N", &RECO_PFJET_N, &b_RECO_PFJET_N);
   fChain->SetBranchAddress("RECO_PFJET_CHARGE", RECO_PFJET_CHARGE, &b_RECO_PFJET_CHARGE);
   fChain->SetBranchAddress("RECO_PFJET_ET", RECO_PFJET_ET, &b_RECO_PFJET_ET);
   fChain->SetBranchAddress("RECO_PFJET_PT", RECO_PFJET_PT, &b_RECO_PFJET_PT);
   fChain->SetBranchAddress("RECO_PFJET_ETA", RECO_PFJET_ETA, &b_RECO_PFJET_ETA);
   fChain->SetBranchAddress("RECO_PFJET_PHI", RECO_PFJET_PHI, &b_RECO_PFJET_PHI);
   fChain->SetBranchAddress("RECO_PFJET_PUID", RECO_PFJET_PUID, &b_RECO_PFJET_PUID);
   fChain->SetBranchAddress("RECO_PFJET_PUID_MVA", RECO_PFJET_PUID_MVA, &b_RECO_PFJET_PUID_MVA);
   fChain->SetBranchAddress("RHO_ele", &RHO_ele, &b_RHO_ele);
   fChain->SetBranchAddress("RHO_mu", &RHO_mu, &b_RHO_mu);

  fChain->SetBranchAddress("RECO_PFJET_PT_UncUp",  RECO_PFJET_PT_UncUp,  &b_RECO_PFJET_PT_UncUp);
   fChain->SetBranchAddress("RECO_PFJET_PT_UncDn",  RECO_PFJET_PT_UncDn,  &b_RECO_PFJET_PT_UncDn);   
   fChain->SetBranchAddress("RECO_PFJET_AREA", RECO_PFJET_AREA, &b_RECO_PFJET_AREA);
   fChain->SetBranchAddress("RECO_PFJET_PTD", RECO_PFJET_PTD, &b_RECO_PFJET_PTD);
   fChain->SetBranchAddress("RECO_PFJET_CHARGED_HADRON_ENERGY", RECO_PFJET_CHARGED_HADRON_ENERGY, &b_RECO_PFJET_CHARGED_HADRON_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_NEUTRAL_HADRON_ENERGY", RECO_PFJET_NEUTRAL_HADRON_ENERGY, &b_RECO_PFJET_NEUTRAL_HADRON_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_PHOTON_ENERGY", RECO_PFJET_PHOTON_ENERGY, &b_RECO_PFJET_PHOTON_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_ELECTRON_ENERGY", RECO_PFJET_ELECTRON_ENERGY, &b_RECO_PFJET_ELECTRON_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_MUON_ENERGY", RECO_PFJET_MUON_ENERGY, &b_RECO_PFJET_MUON_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_HF_HADRON_ENERGY", RECO_PFJET_HF_HADRON_ENERGY, &b_RECO_PFJET_HF_HADRON_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_HF_EM_ENERGY", RECO_PFJET_HF_EM_ENERGY, &b_RECO_PFJET_HF_EM_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_CHARGED_EM_ENERGY", RECO_PFJET_CHARGED_EM_ENERGY, &b_RECO_PFJET_CHARGED_EM_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_CHARGED_MU_ENERGY", RECO_PFJET_CHARGED_MU_ENERGY, &b_RECO_PFJET_CHARGED_MU_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_NEUTRAL_EM_ENERGY", RECO_PFJET_NEUTRAL_EM_ENERGY, &b_RECO_PFJET_NEUTRAL_EM_ENERGY);
   fChain->SetBranchAddress("RECO_PFJET_CHARGED_HADRON_MULTIPLICITY", RECO_PFJET_CHARGED_HADRON_MULTIPLICITY, &b_RECO_PFJET_CHARGED_HADRON_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY", RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY, &b_RECO_PFJET_NEUTRAL_HADRON_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_PHOTON_MULTIPLICITY", RECO_PFJET_PHOTON_MULTIPLICITY, &b_RECO_PFJET_PHOTON_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_ELECTRON_MULTIPLICITY", RECO_PFJET_ELECTRON_MULTIPLICITY, &b_RECO_PFJET_ELECTRON_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_MUON_MULTIPLICITY", RECO_PFJET_MUON_MULTIPLICITY, &b_RECO_PFJET_MUON_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_HF_HADRON_MULTIPLICTY", RECO_PFJET_HF_HADRON_MULTIPLICTY, &b_RECO_PFJET_HF_HADRON_MULTIPLICTY);
   fChain->SetBranchAddress("RECO_PFJET_HF_EM_MULTIPLICITY", RECO_PFJET_HF_EM_MULTIPLICITY, &b_RECO_PFJET_HF_EM_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_CHARGED_MULTIPLICITY", RECO_PFJET_CHARGED_MULTIPLICITY, &b_RECO_PFJET_CHARGED_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_NEUTRAL_MULTIPLICITY", RECO_PFJET_NEUTRAL_MULTIPLICITY, &b_RECO_PFJET_NEUTRAL_MULTIPLICITY);
   fChain->SetBranchAddress("RECO_PFJET_NCOMPONENTS", RECO_PFJET_NCOMPONENTS, &b_RECO_PFJET_NCOMPONENTS);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_PDGID", RECO_PFJET_COMPONENT_PDGID, &b_RECO_PFJET_COMPONENT_PDGID);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_PT", RECO_PFJET_COMPONENT_PT, &b_RECO_PFJET_COMPONENT_PT);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_ETA", RECO_PFJET_COMPONENT_ETA, &b_RECO_PFJET_COMPONENT_ETA);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_PHI", RECO_PFJET_COMPONENT_PHI, &b_RECO_PFJET_COMPONENT_PHI);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_E", RECO_PFJET_COMPONENT_E, &b_RECO_PFJET_COMPONENT_E);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_CHARGE", RECO_PFJET_COMPONENT_CHARGE, &b_RECO_PFJET_COMPONENT_CHARGE);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_TRANSVERSE_MASS", RECO_PFJET_COMPONENT_TRANSVERSE_MASS, &b_RECO_PFJET_COMPONENT_TRANSVERSE_MASS);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_XVERTEX", RECO_PFJET_COMPONENT_XVERTEX, &b_RECO_PFJET_COMPONENT_XVERTEX);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_YVERTEX", RECO_PFJET_COMPONENT_YVERTEX, &b_RECO_PFJET_COMPONENT_YVERTEX);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_ZVERTEX", RECO_PFJET_COMPONENT_ZVERTEX, &b_RECO_PFJET_COMPONENT_ZVERTEX);
   fChain->SetBranchAddress("RECO_PFJET_COMPONENT_VERTEX_CHI2", RECO_PFJET_COMPONENT_VERTEX_CHI2, &b_RECO_PFJET_COMPONENT_VERTEX_CHI2);
   //@
   fChain->SetBranchAddress("LHE_PARTON_N", &LHE_PARTON_N, &b_LHE_PARTON_N);
   fChain->SetBranchAddress("LHE_PARTON_CLEAR", LHE_PARTON_CLEAR, &b_LHE_PARTON_CLEAR);
   fChain->SetBranchAddress("LHE_PARTON_PDGID", LHE_PARTON_PDGID, &b_LHE_PARTON_PDGID);
   fChain->SetBranchAddress("LHE_PARTON_PT", LHE_PARTON_PT, &b_LHE_PARTON_PT);
   fChain->SetBranchAddress("LHE_PARTON_ETA", LHE_PARTON_ETA, &b_LHE_PARTON_ETA);
   fChain->SetBranchAddress("LHE_PARTON_PHI", LHE_PARTON_PHI, &b_LHE_PARTON_PHI);
   fChain->SetBranchAddress("LHE_PARTON_E", LHE_PARTON_E, &b_LHE_PARTON_E);

   fChain->SetBranchAddress("RECO_CALOMET", &RECO_CALOMET, &b_RECO_CALOMET);
   fChain->SetBranchAddress("RECO_PFMET", &RECO_PFMET, &b_RECO_PFMET);
   fChain->SetBranchAddress("RECO_PFMET_X", &RECO_PFMET_X, &b_RECO_PFMET_X);
   fChain->SetBranchAddress("RECO_PFMET_Y", &RECO_PFMET_Y, &b_RECO_PFMET_Y);
   fChain->SetBranchAddress("RECO_PFMET_PHI", &RECO_PFMET_PHI, &b_RECO_PFMET_PHI);
   fChain->SetBranchAddress("RECO_PFMET_THETA", &RECO_PFMET_THETA, &b_RECO_PFMET_THETA);
   fChain->SetBranchAddress("RECO_PFMETT1", &RECO_PFMETT1, &b_RECO_PFMETT1);
   fChain->SetBranchAddress("RECO_PFMETT1_X", &RECO_PFMETT1_X, &b_RECO_PFMETT1_X);
   fChain->SetBranchAddress("RECO_PFMETT1_Y", &RECO_PFMETT1_Y, &b_RECO_PFMETT1_Y);
   fChain->SetBranchAddress("RECO_PFMETT1_PHI", &RECO_PFMETT1_PHI, &b_RECO_PFMETT1_PHI);
   fChain->SetBranchAddress("RECO_PFMETT1_THETA", &RECO_PFMETT1_THETA, &b_RECO_PFMETT1_THETA);
   fChain->SetBranchAddress("RECO_TCMET", &RECO_TCMET, &b_RECO_TCMET);
   fChain->SetBranchAddress("RECO_CORMETMUONS", &RECO_CORMETMUONS, &b_RECO_CORMETMUONS);

 //@Type I correction MET
   fChain->SetBranchAddress("RECO_PFMET_JetEnUp", &RECO_PFMET_JetEnUp, &b_RECO_PFMET_JetEnUp);
   fChain->SetBranchAddress("RECO_PFMET_JetEnDn", &RECO_PFMET_JetEnDn, &b_RECO_PFMET_JetEnDn);
   fChain->SetBranchAddress("RECO_PFMET_ElectronEnUp", &RECO_PFMET_ElectronEnUp, &b_RECO_PFMET_ElectronEnUp);
   fChain->SetBranchAddress("RECO_PFMET_ElectronEnDn", &RECO_PFMET_ElectronEnDn, &b_RECO_PFMET_ElectronEnDn);
   fChain->SetBranchAddress("RECO_PFMET_MuonEnUp", &RECO_PFMET_MuonEnUp, &b_RECO_PFMET_MuonEnUp);
   fChain->SetBranchAddress("RECO_PFMET_MuonEnDn", &RECO_PFMET_MuonEnDn, &b_RECO_PFMET_MuonEnDn);
   fChain->SetBranchAddress("RECO_PFMET_JetResUp", &RECO_PFMET_JetResUp, &b_RECO_PFMET_JetResUp);
   fChain->SetBranchAddress("RECO_PFMET_JetResDn", &RECO_PFMET_JetResDn, &b_RECO_PFMET_JetResDn);
   fChain->SetBranchAddress("RECO_PFMET_UnclusteredEnUp", &RECO_PFMET_UnclusteredEnUp, &b_RECO_PFMET_UnclusteredEnUp);
   fChain->SetBranchAddress("RECO_PFMET_UnclusteredEnDn", &RECO_PFMET_UnclusteredEnDn, &b_RECO_PFMET_UnclusteredEnDn);
   fChain->SetBranchAddress("RECO_PFMET_PhotonEnUp", &RECO_PFMET_PhotonEnUp, &b_RECO_PFMET_PhotonEnUp);
   fChain->SetBranchAddress("RECO_PFMET_PhotonEnDn", &RECO_PFMET_PhotonEnDn, &b_RECO_PFMET_PhotonEnDn);
   

   fChain->SetBranchAddress("tCHighEff_BTagJet_PT", tCHighEff_BTagJet_PT, &b_tCHighEff_BTagJet_PT);
   fChain->SetBranchAddress("tCHighPur_BTagJet_PT", tCHighPur_BTagJet_PT, &b_tCHighPur_BTagJet_PT);
   fChain->SetBranchAddress("cSV_BTagJet_PT", cSV_BTagJet_PT, &b_cSV_BTagJet_PT);
   fChain->SetBranchAddress("tCHighEff_BTagJet_ETA", tCHighEff_BTagJet_ETA, &b_tCHighEff_BTagJet_ETA);
   fChain->SetBranchAddress("tCHighPur_BTagJet_ETA", tCHighPur_BTagJet_ETA, &b_tCHighPur_BTagJet_ETA);
   fChain->SetBranchAddress("cSV_BTagJet_ETA", cSV_BTagJet_ETA, &b_cSV_BTagJet_ETA);
   fChain->SetBranchAddress("tCHighEff_BTagJet_PHI", tCHighEff_BTagJet_PHI, &b_tCHighEff_BTagJet_PHI);
   fChain->SetBranchAddress("tCHighPur_BTagJet_PHI", tCHighPur_BTagJet_PHI, &b_tCHighPur_BTagJet_PHI);
   fChain->SetBranchAddress("cSV_BTagJet_PHI", cSV_BTagJet_PHI, &b_cSV_BTagJet_PHI);
   fChain->SetBranchAddress("tCHighEff_BTagJet_DISCR", tCHighEff_BTagJet_DISCR, &b_tCHighEff_BTagJet_DISCR);
   fChain->SetBranchAddress("tCHighPur_BTagJet_DISCR", tCHighPur_BTagJet_DISCR, &b_tCHighPur_BTagJet_DISCR);
   fChain->SetBranchAddress("cSV_BTagJet_DISCR", cSV_BTagJet_DISCR, &b_cSV_BTagJet_DISCR);
   fChain->SetBranchAddress("cSV_BTagJet_ET", cSV_BTagJet_ET, &b_cSV_BTagJet_ET);
 
   Notify();
}

Bool_t HZZ4LeptonsAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HZZ4LeptonsAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HZZ4LeptonsAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void HZZ4LeptonsAnalysis::printmubnn(int i){    
  bnn_file 
                << RECOMU_PT[i]  << " "  
                << RECOMU_ETA[i] << " "  
                << RECOMU_PHI[i] << " "  
                << RECOMU_CHARGE[i] << " "
                << RECOMU_PFX_dB[i] << " "
                << RECOMU_SIP[i] << " ";
}

#endif // #ifdef HZZ4LeptonsAnalysis_cxx
