#include "MyDisplay.h"
#include <TApplication.h>
#include <TSystem.h>
#include <iostream>

int main(int argc, char** argv)
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " mdt_hits.root detector_geometry.gdml" << std::endl;
    return 1;
  }
  // Print current directory and full paths
  std::cout << "Current directory: ";
  gSystem->Exec("pwd");
  
  TString rootFile = gSystem->ExpandPathName(argv[1]);
    TString gdmlFile = gSystem->ExpandPathName(argv[2]);
    
    std::cout << "Root file path: " << rootFile << std::endl;
    std::cout << "GDML file path: " << gdmlFile << std::endl;

    // Verify files exist
    if (gSystem->AccessPathName(rootFile)) {
        std::cerr << "ERROR: Root file not found: " << rootFile << std::endl;
        return 1;
    }
    if (gSystem->AccessPathName(gdmlFile)) {
        std::cerr << "ERROR: GDML file not found: " << gdmlFile << std::endl;
        return 1;
    }
    
    TApplication app("app", &argc, argv);
    
    try {
        MyDisplay display;
        display.Initialize(rootFile.Data(), gdmlFile.Data());
        display.Run();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    
    return 0;
}
