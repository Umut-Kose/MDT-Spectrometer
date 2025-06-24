#include "TrackReconstruction.h"
#include <iostream>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <input_file.root>" << std::endl;
        return 1;
    }
    
    TrackReconstructor reconstructor(argv[1]);
    
    std::cout << "Loading data..." << std::endl;
    reconstructor.LoadData();
    
    std::cout << "Reconstructing tracks..." << std::endl;
    reconstructor.ReconstructTracks();
    
    std::cout << "Analyzing performance..." << std::endl;
    reconstructor.AnalyzePerformance();
    
    std::cout << "Drawing results..." << std::endl;
    reconstructor.DrawResults();
    
    std::cout << "Analysis complete!" << std::endl;
    
    return 0;
}
