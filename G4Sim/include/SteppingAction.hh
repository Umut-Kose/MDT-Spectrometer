// SteppingAction.hh
#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "globals.hh"
#include <fstream>
#include "TFile.h"
#include "TTree.h"
class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction();
    virtual ~SteppingAction();
    virtual void UserSteppingAction(const G4Step* step) override;
    
    TFile* fRootFile;
    TTree* fStepTree;
    // Variables for branches:
    int eventID, trackID, pdgCode, stepNumber;
    double x, y, z, px, py, pz, theta, phi, theta0, Bx, By, Bz;

private:
    std::ofstream fOutputFile;
};