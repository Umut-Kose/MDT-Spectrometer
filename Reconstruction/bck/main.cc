#include "TrackAnalysis.hh"
#include <iostream>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root>" << std::endl;
        return 1;
    }

    TrackAnalysis analyzer;
    analyzer.ProcessFile(argv[1]);
    analyzer.WriteResults("tracking_results.root");

    return 0;
}
