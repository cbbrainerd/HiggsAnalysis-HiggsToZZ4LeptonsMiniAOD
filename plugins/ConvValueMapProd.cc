#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/ConvValueMapProd.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/ValueMap.h"


ConvValueMapProd::ConvValueMapProd(const edm::ParameterSet& iConfig){ 
    gsfLabel_= consumes<std::vector<reco::GsfElectron> >(iConfig.getUntrackedParameter<edm::InputTag>("gsfLabel"));
    tkLabel_ = consumes<std::vector<reco::Track> >(iConfig.getUntrackedParameter<edm::InputTag>("tkLabel"));

    produces<edm::ValueMap<float> >("dist").setBranchAlias("dist");
    produces<edm::ValueMap<float> >("dcot").setBranchAlias("dcot");

}

void ConvValueMapProd::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<reco::GsfElectronCollection> gsfH;
    iEvent.getByToken(gsfLabel_,gsfH);

    edm::Handle<reco::TrackCollection> tkH;
    iEvent.getByToken(tkLabel_,tkH);

    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
    const  MagneticField *mField = magneticField.product();
    double bfield = mField->inTesla(GlobalPoint(0.,0.,0.)).z();

    std::vector<float> distsV;
    std::auto_ptr<edm::ValueMap<float> > distsM(new edm::ValueMap<float> ());
    edm::ValueMap<float>::Filler distsF(*distsM);

    std::vector<float> dcotsV;
    std::auto_ptr<edm::ValueMap<float> > dcotsM(new edm::ValueMap<float> ());
    edm::ValueMap<float>::Filler dcotsF(*dcotsM);

    ConversionFinder convFinder;
    ConversionInfo convInfo;

    for(size_t i=0;i<gsfH->size();++i) {

        convInfo = convFinder.getConversionInfo(gsfH->at(i), tkH, bfield);
        distsV.push_back(convInfo.dist());
        dcotsV.push_back(convInfo.dcot());

    }

    distsF.insert(gsfH,distsV.begin(),distsV.end());
    dcotsF.insert(gsfH,dcotsV.begin(),dcotsV.end());

    distsF.fill();
    dcotsF.fill();
    iEvent.put(distsM,"dist");
    iEvent.put(dcotsM,"dcot");


}

ConvValueMapProd::~ConvValueMapProd() { }
void ConvValueMapProd::beginJob() { }
void ConvValueMapProd::endJob() { }
