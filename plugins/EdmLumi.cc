//
// $Id: EdmLumi.cc,v 1.3 2010/10/08 05:34:25 ndefilip Exp $
//

/**
  \class    gEdmLumi PATMuonKinematics.h "PhysicsTools/PatAlgos/interface/PATMuonKinematics.h"
  \brief    Use StandAlone track to define the 4-momentum of a PAT Muon (normally the global one is used)
            
  \author   Giovanni Petrucciani
  \version  $Id: EdmLumi.cc,v 1.3 2010/10/08 05:34:25 ndefilip Exp $
*/


#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <FWCore/ParameterSet/interface/FileInPath.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <map>
#include <vector>
#include <iostream>

class EdmLumi : public edm::EDFilter {
    public:
        explicit EdmLumi(const edm::ParameterSet & iConfig);
        virtual ~EdmLumi() { }

        virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);
        virtual bool beginLuminosityBlock(edm::LuminosityBlock &lumi, const edm::EventSetup &iSetup);
        virtual bool endLuminosityBlock(edm::LuminosityBlock &lumi, const edm::EventSetup &iSetup);
        virtual void endJob();

    private:
        bool printMissingRuns, printMissingLumis;

        enum Mode { Online, OfflineVtx, OfflineHF };
        Mode mode;

        double lumi, tlumi, scale;

        std::map<uint32_t, std::vector<float> > lumiData;
        std::map<uint32_t, uint32_t>            prescaleData;

        bool isZeroLumi;
        uint32_t nInLumi, nInZeroLumi;
        std::string label;
};

EdmLumi::EdmLumi(const edm::ParameterSet & iConfig) :
    printMissingRuns( iConfig.getUntrackedParameter<bool>("printMissingRuns",  true)),
    printMissingLumis(iConfig.getUntrackedParameter<bool>("printMissingLumis", true)),
    lumi(0), tlumi(0),
    scale(0),
    lumiData(), prescaleData(),
    nInLumi(0), nInZeroLumi(0),
    label(iConfig.getParameter<std::string>("@module_label"))
{
    if (label.find("edmLumi") == 0) label = label.substr(7);
    if (!label.empty()) label = label + ": ";
    if (iConfig.existsAs<std::string>("lumi_by_LS_all_csv",false)) {
        edm::FileInPath csv(iConfig.getUntrackedParameter<std::string>("lumi_by_LS_all_csv"));
        FILE *f = fopen(csv.fullPath().c_str(), "r");
        if (f == 0) throw cms::Exception("InputError") << "Can't open file '" << csv.fullPath() << "'\n";

        std::string method = iConfig.getParameter<std::string>("method");
        if (method == "vertex" || method == "vtx") { mode = OfflineVtx; }
        else if (method == "hf") { mode = OfflineHF; }
        else if (method == "online") { mode = Online; }
        else throw cms::Exception("Configuration") << "Parameter 'method' must be either 'vertex' or 'hf'.\n";

        char buff[1024];
        int run, lumi, c1, c2; float lvtx, lhf;
        int ret;
        int totalSize = 0;
        // skip header (5 for offilne, 1 for online)
        for (size_t i = 0; i < (mode == Online ? 1 : 5); ++i) { fgets(buff, 1024, f);  }
        while(true) {
            if (mode == Online) ret = fscanf(f, "%d,%d,%g,%g\n", &run, &lumi, &lvtx, &lhf);
            else ret = fscanf(f, " %d, %d, %d, %d, %g, %g\n", &run, &lumi, &c1, &c2, &lhf, &lvtx);
            if (ret == EOF) break;
            if (ret == (mode == Online ? 4 : 6)) {
                std::vector<float> & lvec =  lumiData[run];
                if (lvec.size() <= size_t(lumi)) lvec.resize(2*lumi,0);
                lvec[lumi] = (mode == OfflineVtx ? lvtx : lhf);
                totalSize++;
            }
        }
        fclose(f);
        scale = iConfig.getParameter<double>("scale");
        edm::LogSystem("EdmLumi") << label << "Loaded " << lumiData.size() << " runs, " << totalSize << " lumisections" <<
                                              " from lumi_by_LS_all_csv file '" << csv.fullPath() << "' .";
    } else {
        edm::LogSystem("EdmLumi") << label << "No lumi_by_LS_all_csv file: assuming all lumiblocks have unit luminosity. ";
    }

    if (iConfig.existsAs<std::string>("prescale_by_run",false)) {
        edm::FileInPath prs(iConfig.getUntrackedParameter<std::string>("prescale_by_run"));
        FILE *f = fopen(prs.fullPath().c_str(), "r");
        int ret, run, prescale;
        while((ret = fscanf(f, "%d %d\n", &run, &prescale)) != EOF) {
            if (ret == 2) {
                prescaleData[run] = prescale;
            }
        }
        edm::LogSystem("EdmLumi") << label << "Loaded prescales for " << prescaleData.size() << " runs" <<
                                              " from prescale_by_run file '" << prs.fullPath() << "' .";
    }
}


bool
EdmLumi::beginLuminosityBlock(edm::LuminosityBlock &iLumi, const edm::EventSetup &iSetup) {
    if (lumiData.empty()) { isZeroLumi = false; return true; }
    isZeroLumi = true;
    std::map<uint32_t, std::vector<float> >::iterator itrun = lumiData.find(iLumi.run());
    if (itrun != lumiData.end() &&
        iLumi.luminosityBlock() < itrun->second.size() &&
        itrun->second[iLumi.luminosityBlock()] != 0) {
        isZeroLumi = false;
    }
    return true;
}

bool
EdmLumi::endLuminosityBlock(edm::LuminosityBlock &iLumi, const edm::EventSetup &iSetup) {
    if (lumiData.empty()) { 
        lumi++; return true; 
    } 

    std::map<uint32_t, std::vector<float> >::iterator itrun = lumiData.find(iLumi.run());
    if (itrun == lumiData.end()) { 
        if (printMissingRuns) edm::LogError("MissingRun") << "Missing run " << iLumi.run();
        return true;
    } 
    if (itrun->second.size() <= iLumi.luminosityBlock()) {
        if (printMissingLumis) edm::LogError("MissingLumi") << "Missing run " << iLumi.run() << ", lumi block " << iLumi.luminosityBlock();
        return true;
    }
    float & thislumi = itrun->second[iLumi.luminosityBlock()];
    if (thislumi == 0) {
        if (printMissingLumis) edm::LogWarning("MissingLumi") << "Zero lumi for run " << iLumi.run() << ", lumi block " << iLumi.luminosityBlock();
        return true;
    }
    if (!prescaleData.empty()) {
        std::map<uint32_t, uint32_t>::const_iterator itp = prescaleData.find(iLumi.run());
        if (itp != prescaleData.end()) {
            if (itp->second != 0) {
                tlumi += thislumi/itp->second;
            }
        }
    }
    lumi += thislumi;
    thislumi = 0; // avoid double-counting in case of duplicated lumis
    return true;
}

void
EdmLumi::endJob() {
    if (!lumiData.empty()) {
        lumi  /= scale; 
        tlumi /= scale; 
        edm::LogSystem("EdmLumi") << label << "Total integrated luminosity: " << lumi/1000 << " /nb";
        if (nInLumi + nInZeroLumi > 0) {
            edm::LogSystem("EdmLumi") << label << "Fraction of events in zero-lumi blocks: " << 
                nInZeroLumi << "/" << (nInLumi + nInZeroLumi) << " = " << double(nInZeroLumi)/double(nInLumi + nInZeroLumi);
        }
        if (!prescaleData.empty()) {
            edm::LogSystem("EdmLumi") << label << "Integrated luminosity after prescale: " << tlumi/1000 << " /nb  (average prescale " << (lumi == 0 ? 1 : tlumi/lumi) << ")";
        }
    } else {
        edm::LogSystem("EdmLumi") << label << "Total integrated luminosity: " << lumi << " /au";
    }
}

bool 
EdmLumi::filter(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    if (isZeroLumi) nInZeroLumi++; else nInLumi++;
    return true;
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EdmLumi);
