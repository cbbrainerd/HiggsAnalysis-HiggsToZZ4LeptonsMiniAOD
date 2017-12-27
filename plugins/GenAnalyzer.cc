
#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/GenAnalyzer.h"

using namespace std;
using namespace edm;
using namespace reco;

GenAnalyzer::GenAnalyzer(const edm::ParameterSet &params) :
	sourceLabel(params.getParameter<edm::InputTag>("src"))
{


  // Create the root file
  theFile_ = new TFile("GenAnalyzer.root", "RECREATE");
  theFile_->cd();
  
  theTree_ = new TTree("GenAnalyzerTree", "GenAnalyzer Tree");
  
  theTree_->Branch("Run",&irun,"irun/i");
  theTree_->Branch("Event",&ievt,"ievt/i");
  theTree_->Branch("LumiSection",&ils,"ils/i");
  theTree_->Branch("weightgen",&weightgen,"weightgen/D");
  theTree_->Branch("PT_h",&PT_h,"PT_h/F");
  theTree_->Branch("PT_z",PT_z,"PT_z[2]/F");
  theTree_->Branch("PT_l",&PT_l,"PT_l/F");
  theTree_->Branch("PT_l_z",PT_l_z,"PT_l_z[4]/F");
  theTree_->Branch("PT_l_stable",PT_l_stable,"PT_l_stable[50]/F");
  theTree_->Branch("PT_l_1",&PT_l_1,"PT_l_1/F");
  theTree_->Branch("PT_l_2",&PT_l_2,"PT_l_2/F");
  theTree_->Branch("PT_l_3",&PT_l_3,"PT_l_3/F");
  theTree_->Branch("PT_l_4",&PT_l_4,"PT_l_4/F");
  theTree_->Branch("Mass_fourl",Mass_fourl,"Mass_fourl[50]/F");
  theTree_->Branch("Mass_ZZ",Mass_ZZ,"Mass_ZZ[50]/F");
  theTree_->Branch("PT_ZZ",PT_ZZ,"PT_ZZ[50]/F");

  nevt=0;

}

GenAnalyzer::~GenAnalyzer()
{

  // Write the histos to file
  theFile_->cd();
  theFile_->Write() ;
  theFile_->Close();

  cout << "Number of events analysed for the ROOT tree= " << nevt << std::endl;


}


void GenAnalyzer::analyze(const edm::Event &event, const edm::EventSetup &es)
{

  nevt++;
  //Initialize variables
  Initialize();

 // Dump Run and Event
  irun=event.id().run();
  ievt=event.id().event();
  ils=event.luminosityBlock();

  // get the weight
  edm::Handle<GenEventInfoProduct> hEventInfo;
  event.getByLabel("generator", hEventInfo);

  if (hEventInfo.isValid())
    weightgen=hEventInfo->weight();
  else 
    weightgen=1.;

  cout << "Dumping the run=" << irun << "  event=" << ievt << "  lumisection=" << ils << "  weight=" << weightgen << std::endl;


  // get gen particle candidates 
  edm::Handle<reco::GenParticleCollection> genCandidates;	    
  event.getByLabel(sourceLabel, genCandidates);
  es.getData( pdt_ );

  leptonpt.clear();
  
  for ( GenParticleCollection::const_iterator mcIter=genCandidates->begin(); mcIter!=genCandidates->end(); ++mcIter ) {
    
    // higgs pt
    if ( mcIter->pdgId()==25 ){
       PT_h=mcIter->pt();
       std::cout << "Found Higgs with mass= " << mcIter->mass() << " and PT= " << PT_h << std::endl;	    
    }
    
    // z pt
    int i=0;
    if ( abs(mcIter->pdgId()==23) && mcIter->status()==3 ) {
      PT_z[i]=mcIter->pt();
      std::cout << "Found Z with mass= " << mcIter->mass() << " and charge= " << mcIter->charge() << std::endl;	    
      i++;
    }	

    i=0;
    if ( ( abs(mcIter->pdgId())==11 || abs(mcIter->pdgId())==13 || abs(mcIter->pdgId())==15  ) && 
	 ( abs(mcIter->mother(0)->pdgId())==23) ) {	    
      std::cout << "Found lepton from z" << std::endl;
      PT_l_z[i]=mcIter->pt();
      //leptonpt.push_back(mcIter->pt());
      i++;
    }
    
    
   
    // lepton stable
    i=0;
    if ( ( abs(mcIter->pdgId())==11 || abs(mcIter->pdgId())==13 || abs(mcIter->pdgId())==15  ) && mcIter->status()==1 ){
      PT_l_stable[i]=mcIter->pt();
      leptonpt.push_back(mcIter->pt());
      i++;
      std::cout << " \n Found lepton with Id= " << mcIter->pdgId() << "  in the final state with mother= " ; 
      //if (mcIter->numberOfMothers()>0 &&  mcIter->mother(0)->status()!=3){
      if (mcIter->numberOfMothers()>0){
	std::cout << getParticleName(int(mcIter->mother(0)->pdgId())) << " << "; 
	if (mcIter->mother(0)->status()==3) continue;
	if (mcIter->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->status()!=3 ){
	  std::cout << getParticleName(int(mcIter->mother(0)->mother(0)->pdgId())) << " << ";
	  if (mcIter->mother(0)->mother(0)->status()==3) continue;
	  if (mcIter->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->status()!=3 ){
	    std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
	    if (mcIter->mother(0)->mother(0)->mother(0)->status()==3) continue;
	    if (mcIter->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
	      std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
	      if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
	      if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		std::cout<< getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << ";
		if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->numberOfMothers()>0 && mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()!=3 ){
		  std::cout << getParticleName(int(mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->pdgId())) << " << " << std::endl;
		  if (mcIter->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->mother(0)->status()==3) continue;
		}
	      }
	    }
	  }
	}
      }      
    }     
    
  }

  sort(leptonpt.rbegin(),leptonpt.rend());

  if (leptonpt.size()>=4){
    std::cout << "\n Sorted vector in decreasing order= " << leptonpt.at(0) << " " << leptonpt.at(1) << " " << leptonpt.at(2) << " " << leptonpt.at(3) << endl;
    PT_l_1=leptonpt.at(0);
    PT_l_2=leptonpt.at(1);
    PT_l_3=leptonpt.at(2);
    PT_l_4=leptonpt.at(3);
  }

 
  // // get clean 4l candidates
  int i=0;
  edm::Handle<edm::View<Candidate> >  fourlCandidates;
  event.getByLabel("fourleptons", fourlCandidates);
  for (edm::View<Candidate>::const_iterator mcIter=fourlCandidates->begin(); mcIter!=fourlCandidates->end(); ++mcIter ) {
    
    cout << "4l Mass= " << mcIter->mass()
         << " Charge= " 
	 << mcIter->daughter(0)->daughter(0)->charge() << " " 
	 << mcIter->daughter(0)->daughter(1)->charge() << " "
         << mcIter->daughter(0)->daughter(2)->charge() << " " 
         << mcIter->daughter(1)->charge() << " " 
	 << endl;
    Mass_fourl[i]=mcIter->mass();
    if (i>49 ) continue;
    i++;
  }

  // // get clean ZZ candidates
  int j=0;
  edm::Handle<edm::View<Candidate> >  ZZCandidates;
  event.getByLabel("diZ", ZZCandidates);
  for (edm::View<Candidate>::const_iterator mcIterZZ=ZZCandidates->begin(); mcIterZZ!=ZZCandidates->end(); ++mcIterZZ ) {
    
    cout << "ZZ Mass= " << mcIterZZ->mass() 
	 << " and pT= " << mcIterZZ->pt()  
         << " Charge= " 
	 << endl;
    Mass_ZZ[j]=mcIterZZ->mass();
    PT_ZZ[j]=mcIterZZ->pt();
    if (j>49 ) continue;
    j++;
  }


  // fill the tree at end of loop
  theTree_->Fill();


}

std::string GenAnalyzer::getParticleName(int id) const
{
  const ParticleData * pd = pdt_->particle( id );
  if (!pd) {
    std::ostringstream ss;
    ss << "P" << id;
    return ss.str();
  } else
    return pd->name();
}

void GenAnalyzer::Initialize(){

  irun=0;
  ievt=0;
  ils=0;
  weightgen=-99999.;

  PT_h=-999.;
  PT_l=-999.;
  for (int i=0;i<2;i++){
    PT_z[i]=-999.;
  }
  for (int i=0;i<4;i++){
    PT_l_z[i]=-999.;
  }
  for (int i=0;i<50;i++){
    PT_l_stable[i]=-999.;
    Mass_fourl[i]=-999.;
    Mass_ZZ[i]=-999.;
    PT_ZZ[i]=-999.;
  }
  PT_l_1=-999.;
  PT_l_2=-999.;
  PT_l_3=-999.;
  PT_l_4=-999.;

}
