#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TStyle.h"
#include <iostream>
#include <vector>
#include "TSystem.h"

void DisplayEvents(const std::string& filename) {
    // Open the ROOT file
    TFile *file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("Hits");
    if (!tree) {
        std::cerr << "TTree 'Hits' not found." << std::endl;
        return;
    }

    // Variables to read
    std::vector<double> *x = nullptr, *y = nullptr, *z = nullptr;
    int eventID;

    // Set branch addresses
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("z", &z);
    tree->SetBranchAddress("eventID", &eventID);

    TCanvas *c1 = new TCanvas("c1", "Event Display", 800, 600);
    gStyle->SetOptStat(0);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries && i < 10; ++i) {
        tree->GetEntry(i);

        std::cout << "Displaying Event ID: " << eventID << "  (Entry: " << i << ")" << std::endl;

        TGraph *gr = new TGraph(z->size());
        for (size_t j = 0; j < z->size(); ++j) {
            gr->SetPoint(j, z->at(j), y->at(j)); // Z vs Y projection
        }

        gr->SetTitle(Form("Event ID %d (z vs y);z [mm];y [mm]", eventID));
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.2);
        gr->Draw("AP");

        c1->Update();
	gSystem->ProcessEvents();
	
        std::cout << "Press Enter to show next event..." << std::endl;
        std::cin.get();

        delete gr;
    }

    //delete c1;
    file->Close();
}
int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);
    DisplayEvents("scifi_hits.root");
    app.Run();
    return 0;
}
