
/* Original author:  Nicola De Filippis - Politecnico and INFN Bari 
 * Contribution by:
 *              Sherif Elgammal - British University of Egypt
 *              Giorgia Miniello - Università di Bari
 */


#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Nicola
// code for accessing registry info
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/Registry.h"

#include "FWCore/Framework/interface/ConstProductRegistry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Provenance/interface/ProcessHistoryRegistry.h"
//nicola


// Chi2 Prob
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

// Run
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "FWCore/Framework/interface/Run.h"

// Maps
#include "DataFormats/Common/interface/ValueMap.h"

// CommonRootTree
#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsCommonRootTree.h"


// Namespaces
using namespace std;


// Constructor
HZZ4LeptonsCommonRootTree::HZZ4LeptonsCommonRootTree(const edm::ParameterSet& pset) {

  // Read parameters files from cfg
  ReadParameters(pset);

  // Create the root file
  theFile_ = new TFile(rootFileName.c_str(), "RECREATE");
  theFile_->cd();

  theTree_ = new TTree("HZZ4LeptonsAnalysis", "HZZ4Leptons Analysis Tree");

  //cout << "This is" << pset.getUntrackedParameter("fileName", std::string()) << endl;

  // Creating branches
  DefineBranches(theTree_);

  // Counter of number of analyzed events
  nevt=0;

  
}

// Destructor
HZZ4LeptonsCommonRootTree::~HZZ4LeptonsCommonRootTree() {

  // Write the histos to file
  theFile_->cd();
  theFile_->Write() ;
  theFile_->Close();

  cout << "Number of events analysed for the ROOT tree= " << nevt << std::endl;

  }


void HZZ4LeptonsCommonRootTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

// nicola
//  bool firstevent_=true;
//  if (firstevent_){
//     bool last=true;
//     // To get access to the process history registry
//     edm::ProcessHistoryRegistry* reghist = edm::ProcessHistoryRegistry::instance();
//     unsigned int i=1;
//     for (edm::ProcessHistoryRegistry::const_iterator ita=reghist->begin(); ita!=reghist->end();++ita){
//       last=false;
//       if (i==reghist->size()){
//         std::cout << "Processing History:"<<std::endl;
//         // printHistory(ita->second);
//       }
//       i++;
//     }
    
    
//     // To get access to the product registry
//     edm::Service<edm::ConstProductRegistry> regserv;
//     const edm::ProductRegistry& pReg = regserv->productRegistry();

//     typedef std::map<std::string,std::vector<edm::BranchDescription> > IdToBranches;
//     typedef std::map<std::pair<std::string,std::string>,IdToBranches> ModuleToIdBranches;
//     ModuleToIdBranches moduleToIdBranches;

//     for (edm::ProductRegistry::ProductList::const_iterator it = 
//               pReg.productList().begin(), itEnd = pReg.productList().end();
//            it != itEnd;
//            ++it) {
//          it->second.init();
//          for (std::map<edm::ProcessConfigurationID, edm::ParameterSetID>::const_iterator
//            itId = it->second.parameterSetIDs().begin(),
//            itIdEnd = it->second.parameterSetIDs().end();
//            itId != itIdEnd;
//            ++itId) {
         
//             std::stringstream s;
//             s << itId->second;
//             moduleToIdBranches[std::make_pair(it->second.processName(),it->second.moduleLabel())][s.str()].push_back(it->second);
//          }
//     }

//     for (ModuleToIdBranches::const_iterator it = moduleToIdBranches.begin(),
//           itEnd = moduleToIdBranches.end();
//           it != itEnd;
//           ++it) {
//       //std::cout <<"Module: "<<it->first.second<<" "<<it->first.first<<std::endl;

//          //Look only to modules for  monte carlo generation
//          if (it->first.second==module_to_search[0] && it->first.first==module_to_search[1]){
//            std::cout <<"Module: "<<it->first.second<<" "<<it->first.first<<std::endl;
//            const IdToBranches& idToBranches = it->second;
//            for (IdToBranches::const_iterator itIdBranch = idToBranches.begin(),
//                   itIdBranchEnd = idToBranches.end();
//                 itIdBranch != itIdBranchEnd;
//                 ++itIdBranch) {
//              std::cout <<" PSet id:"<<itIdBranch->first<<std::endl;
//              std::cout <<" products: {"<<std::endl;
//              for (std::vector<edm::BranchDescription>::const_iterator itBranch = itIdBranch->second.begin(),
//                     itBranchEnd = itIdBranch->second.end();
//                   itBranch != itBranchEnd;
//                   ++itBranch) {
//                std::cout << "  "<< itBranch->branchName()<<std::endl;
//              }
//              std::cout <<"}"<<std::endl;
//              edm::ParameterSetID psid(itIdBranch->first);
//              cout << "Psid= " << psid << endl;
             
//              // to get access to the parameter set registry
//              edm::pset::Registry* reg = edm::pset::Registry::instance();
//              bool found=false;
//              for (edm::pset::Registry::const_iterator it=reg->begin(); it!=reg->end();++it){
//                if (it->first==psid) {
//                  std::cout <<" parameters: "<< it->second << endl;
//     //             printParameters(it->second,par_to_search,originalmH);
//                  // const struct std::pair<const edm::ParameterSetID, edm::ParameterSet>
//                  found=true;
//                }             
//              }      
//              if (found==false) std::cout << "No ParameterSetID for " << psid << std::endl;
//            }
//          }
//     }
//   }

//



  // Incrementing counter of events
  nevt++;
  cout << "Dumping the information of the event in a ROOT tree: " << nevt << std::endl;
  
  //Initialize variables
  Initialize();
  
  // Dump Run, Event, LumiSection
  irun=iEvent.id().run();
  ievt=iEvent.id().event();
  ils=iEvent.luminosityBlock();

  // cout << "Dumping the information of run=" << irun << "  event=" << ievt << "  lumisection=" << ils << std::endl;
  cout << "Dumping_the_information_of_run:" << irun << ":" << ievt << ":" << ils << std::endl;

  edm::Handle<LumiSummary> l;
  iEvent.getLuminosityBlock().getByLabel("lumiProducer", l); 
  // Check that there is something
  if (l.isValid()){
    Avginstlumi=l->avgInsDelLumi();
    cout << "Instataneous luminosity= " << Avginstlumi << endl;
  }

  // file PU block
  if (fillPUinfo) fillPU(iEvent);
  EventsMCReWeighting(iEvent);

  // fill HLT block
  //fillHLTFired(iEvent);

  //trigger matching
  triggermatching(iEvent);
  
  // Get the MC Truth particles, H, ZZ and 4 leptons
  if ( fillMCTruth) {
    cout << "Filling MCtruth variables" << endl;
    fillgenparticles(iEvent,iSetup);
    fillmc(iEvent);
  }

  //Vertices
  edm::Handle<reco::VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByToken(verticesTag_,recoPrimaryVertexCollection);
  PV = *recoPrimaryVertexCollection;
  
  RECO_NVTX=recoPrimaryVertexCollection->size();
  cout << "Number of Vertices in the event= " << RECO_NVTX << endl;

  int index_vertex = 0;
  
  for (VertexCollection::const_iterator i=recoPrimaryVertexCollection->begin(); i!=recoPrimaryVertexCollection->end();i++) {
    if(index_vertex>14) break;
    RECO_VERTEX_x[index_vertex] = i->x();
    RECO_VERTEX_y[index_vertex] = i->y();
    RECO_VERTEX_z[index_vertex] = i->z();
    RECO_VERTEX_ndof[index_vertex] = i->ndof();
    RECO_VERTEX_chi2[index_vertex] = i->chi2();
    RECO_VERTEX_ntracks[index_vertex] = i->tracksSize();
    RECO_VERTEXPROB[index_vertex] = ChiSquaredProbability(i->chi2(),i->ndof());
    RECO_VERTEX_isValid[index_vertex] = i->isValid();
    cout << "Vertex made by " << i->tracksSize() << " tracks with chi2="<< i->chi2() << " and ndof=" << i->ndof() << " and prob=" << RECO_VERTEXPROB[index_vertex] << endl;
    
    int indice=0;
    for(std::vector<reco::TrackBaseRef>::const_iterator iter = i->tracks_begin();
	iter != i->tracks_end(); iter++) {
      cout << "pT of tracks building the vertex= " << (**iter).pt() << endl; 
      if (indice <100) RECO_VERTEX_TRACK_PT[index_vertex][indice]= (**iter).pt();
      indice++;
    }
   
    index_vertex++;
  } // loop on vertices
  
  

  // Fill RECO block in the rootple
  // PF Jets
  cout << "fill jet test" << endl;
  filljets(iEvent,iSetup);
  cout << "fill jet test end" << endl;

//  if (useAdditionalRECO==true) {
//    fillAdditionalRECO(iEvent);
//  cout << "RECOcollNameZ size " << RECOcollNameZ.size() << endl;
//  }


  // Filling electron and muons vectors
  cout <<"fill ele test" << endl;
  fillElectrons(iEvent,iSetup);
  cout <<"fill Muon test" << endl;
  fillMuons(iEvent,iSetup);
  cout <<"fill photon test" << endl;
  fillPhotons(iEvent);



  // ConstraintVertexFit
  //if (useAdditionalRECO==false){
  //  if (decaychannel=="2e2mu" ) fillConstraintVtx2e2mu(iEvent);
  //  if (decaychannel=="4mu" )   fillConstraintVtx4mu(iEvent);
  //  if (decaychannel=="4e" )    fillConstraintVtx4e(iEvent);
  //}
  //else if (useAdditionalRECO==true) {
  //  fillConstraintVtx2e2mu(iEvent);
  //  fillConstraintVtx4mu(iEvent);
  //  fillConstraintVtx4e(iEvent);
  //}
  // fillConstraintVtxTriLeptons(iEvent);
  // fillConstraintVtxDiLeptons(iEvent);

  
  //Tracks

  cout << "fill track test" << endl;
  fillTracks(iEvent); 

 //GENJets
  if (fillMCTruth) fillgenjets(iEvent);

    
  cout << "GENMET test" << endl;

  //GENMET
  if (fillMCTruth) {
    edm::Handle<pat::METCollection> genmetHandle;
    iEvent.getByToken(pfmetTag_,genmetHandle);
    for(unsigned int i=0; i< genmetHandle->size(); i++){
       const pat::MET &c =(*genmetHandle)[i];
       genmet = c.genMET()->pt();
    }
  }
 

 
  cout << "fill MET" << endl;
  //RECO MET
  fillMET(iEvent);

 
  cout << "fill MET Filters" << endl;
  //RECO MET
  fillMETFilters(iEvent);


   cout << " beam spot met " << endl;
  // Beam Spot
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(offlineBeamSpot_,recoBeamSpotHandle) ;
  bs = *recoBeamSpotHandle;
  
  BeamSpot_X=bs.position().x();
  BeamSpot_Y=bs.position().y();
  BeamSpot_Z=bs.position().z();
  
  cout << "BeamSpot:"
       << "  bs_X=" << BeamSpot_X
       << "  bs_Y=" << BeamSpot_Y
       << "  bs_Z=" << BeamSpot_Z
       << endl;

  cout << "Btag test " << endl;
  // btagging
  fillBTagging(iEvent);

  // fill the tree at end of loop
  theTree_->Fill();
  

}

void HZZ4LeptonsCommonRootTree::beginJob() {
  cout << "Filling a ROOT tree for offline selection" << endl;
}


void HZZ4LeptonsCommonRootTree::endJob() {
}


