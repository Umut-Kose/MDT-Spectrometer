#include "TrackAnalysis.hh"
#include <iostream>

#include "TrackAnalysis.hh"
#include <RKTrackRep.h>
#include <FieldManager.h>
#include <DAF.h>
#include <PlanarMeasurement.h>
#include <SpacepointMeasurement.h>
#include <MaterialEffects.h> 

#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoVolume.h>
#include <TGeoMaterialInterface.h>

#include <TApplication.h>


int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root> --verbose=" << std::endl;
        return 1;
    }

    std::string input_file = argv[1]; // Capture the filename FIRST

    int verbose = 0; // default
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.find("--verbose=") == 0) {
            verbose = std::stoi(arg.substr(10)); // "--verbose=" is 10 chars
        } else if (arg == "-v" && i+1 < argc) {
            verbose = std::stoi(argv[i+1]);
            ++i;
        }
    }


    TApplication app("MuonTrackingAnalysis", &argc, argv);


    const char* gdml_file = "/Users/ukose/sw/Work/Magnet_a_la_babymind/MuonSpectrometerSim/G4Sim/Geometry/detector_geometry.gdml";
    TGeoManager* geo = TGeoManager::Import(gdml_file);
    if (!geo) {
        std::cerr << "Could not load GDML geometry!" << std::endl;
        return 1;
    }
    
    //geo->GetTopVolume()->Print();
    //geo->GetListOfMaterials()->Print();

    // Setup GENFIT material effects from TGeoManager
    genfit::TGeoMaterialInterface* matInterface = new genfit::TGeoMaterialInterface();
    genfit::MaterialEffects::getInstance()->init(matInterface);
    //genfit::MaterialEffects::getInstance()->setNoEffects();
    


    TrackAnalysis analyzer(geo);
    analyzer.verbose = verbose; // Set verbosity level
    analyzer.ProcessFile(input_file.c_str());
    analyzer.WriteResults("tracking_results.root");

    //app.Run(); // This keeps canvases alive!
    return 0;
}
