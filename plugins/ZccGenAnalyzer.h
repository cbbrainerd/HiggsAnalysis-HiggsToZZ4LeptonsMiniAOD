#ifndef ZccGenAnalyzer_h
#define ZccGenAnalyzer_h

// system inculde files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <TH1.h>

class ZccGenAnalyzer : public edm::EDAnalyzer {
public:
  explicit ZccGenAnalyzer(const edm::ParameterSet &params);
  virtual ~ZccGenAnalyzer();
   
protected:
  virtual void analyze(const edm::Event &event,
		       const edm::EventSetup &es);
  std::string getParticleName(int id) const;
 
private:
  edm::ESHandle<ParticleDataTable>        pdt_;
  edm::InputTag				sourceLabel;
  float mass,mass1,Tmass;
  bool skimm,preselptE,preseletaE,preselmB,preselmE,preselpte,preseletae,preselmb,preselme,pre2e,pre2m,pre4e,pre4m,preE12,pre2E12,preM12,pre2M12 ,pre4l,ccll;
  std::vector<float> leptonpt,elept,muonpt,muonp,eleeta,muoneta,elen,elep,mun,mup,elepair,muonpair,fourl;

  TH1 *histoPt_c,*histoPt_dc,*histoPt_dcbar,*histoeta_c,*histophi_c,*histoeta_cbar,*histophi_cbar,*histoPt_cbar,*histoeta_l_c,*histoMass_l_c,*histophi_l_c,*histoPt_l_c,*histoPt_z,*histoMass_z,*histoMass_l_z,*histophi_z,*histoeta_z,*histoMass_Zccbar,*histoPt_Zccbar,*histoeta_Zccbar,*histophi_Zccbar,*histoMass_trilepton,*histoPt_trilepton,*histophi_trilepton,*histoeta_trilepton,*histoMass_diele4l,*histoMass_dimu4l,*histoMass_dileptonsz,*histoMass_dileptonsbb,*histoMass_dileptonselepair,*histoMass_dileptonsmuonpair,*histoMass_zdileptons,*histoMass_dileptonpresel,*histoMass_fourleptons,*histoPt_g,*histoPt_l_z,*histoeta_l_z,*histoeta_muonpair,*histoeta_elepair,*histopt_muonpair,*histopt_elepair,*histophi_l_z,*histoPt_l_g,*histoPt_l_stable,*histoeta_l_stable,*histophi_l_stable,*histoPt_e_stable,*histoeta_e_stable,*histophi_e_stable,*histoPt_m_stable,*histoeta_m_stable,*histophi_m_stable,*histop_m_stable,*histoMass_elepair12,* histoPt_fourleptons,* histoPt_diele4l,*histoPt_dimu4l,*histoeta0_diele4l,*histoeta1_diele4l,*histophi0_diele4l,*histophi1_diele4l,*histoZmass,*histocpt,*histocptcm;
  TH1 *histoPt_1,*histoPt_2,*histoPt_3,*histoPt_4,*histoPt_1s,*histoPt_2s,*histoPt_3s,*histoPt_1e,*histoPt_2e,*histoPt_3e,*histoPt_4e,*histoeta_1e,*histoeta_2e,*histoeta_3e,*histoeta_4e,*histoPt_1m,*histoPt_2m,*histoPt_3m,*histoPt_4m,*histop_1m,*histop_2m,*histop_3m,*histop_4m,*histoeta_1m,*histoeta_2m,*histoeta_3m,*histoeta_4m;
  TH1 *histoMass_ccbar,*histoMass_ccbard1,*histoMass_ccbard2,*histoMass_ccbarnc,*histophi_ccbar,*histoPt_ccbar,*histoMass_c,*histoMass_cbar,*histoeta_ccbar,*histoMass_fourl,*histoPt_fourl,*histoPt_D4,*histoeta_D3,*histoeta_D2,*histoeta_D1,*histoeta_D4,*histoPt_D1,*histoPt_D2,*histoPt_D3,*histoMass_D2,*histoMass_D1,*histoeta_fourl,*histophi_fourl;
  TH1 *histonlept,*histophi_D1, *histophi_D2,*histophi_D3,*histophi_D4,*histoelept,*histoeleeta,*histomuonpt,*histomuonp,*histomuoneta,*histonlept_c;
  
};

#endif
