// SteppingAction.cc
#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Event.hh"          // Include for event access
#include "G4RunManager.hh"     // Needed to get current event
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagneticField.hh"

#include "TFile.h"
#include "TTree.h"


SteppingAction::SteppingAction() {

    fRootFile = new TFile("tracking_info.root", "RECREATE");
    fStepTree = new TTree("StepTree", "Step-by-step track info");
    // Book branches
    fStepTree->Branch("eventID",     &eventID,     "eventID/I");
    fStepTree->Branch("trackID",     &trackID,     "trackID/I");
    fStepTree->Branch("pdgCode",     &pdgCode,     "pdgCode/I");
    fStepTree->Branch("stepNumber",  &stepNumber,  "stepNumber/I");
    fStepTree->Branch("x",           &x,           "x/D");
    fStepTree->Branch("y",           &y,           "y/D");
    fStepTree->Branch("z",           &z,           "z/D");
    fStepTree->Branch("px",          &px,          "px/D");
    fStepTree->Branch("py",          &py,          "py/D");
    fStepTree->Branch("pz",          &pz,          "pz/D");
    fStepTree->Branch("theta",       &theta,       "theta/D");
    fStepTree->Branch("phi",         &phi,         "phi/D");
    fStepTree->Branch("theta0",      &theta0,      "theta0/D");
    fStepTree->Branch("Bx",          &Bx,          "Bx/D");
    fStepTree->Branch("By",          &By,          "By/D");
    fStepTree->Branch("Bz",          &Bz,          "Bz/D");


    //fOutputFile.open("tracking_info.dat");
    //fOutputFile << "EventID\tTrackID\tStep#\tX[mm]\tY[mm]\tZ[mm]\tPx[MeV/c]\tPy[MeV/c]\tPz[MeV/c]\tTheta[deg]\tDeltaPhi[deg]\tThetaTh[deg]\n";
    //fOutputFile << "EventID\tTrackID\tStep#\tX[mm]\tY[mm]\tZ[mm]\tPx[MeV/c]\tPy[MeV/c]\tPz[MeV/c]\tTheta[deg]\tDeltaPhi[deg]\tThetaTh[deg]\tBx[T]\tBy[T]\tBz[T]\n";
    //fOutputFile << "EventID\tTrackID\tPDG\tStep#\tX[mm]\tY[mm]\tZ[mm]\tPx[MeV/c]\tPy[MeV/c]\tPz[MeV/c]\tTheta[deg]\tDeltaPhi[deg]\tThetaTh[deg]\tBx[T]\tBy[T]\tBz[T]\n";

}

SteppingAction::~SteppingAction() {
    fRootFile->cd();
    fStepTree->Write();
    fRootFile->Close();
    delete fRootFile;}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    G4Track* track = step->GetTrack();
    // Get the event ID via G4RunManager
    auto* event = G4RunManager::GetRunManager()->GetCurrentEvent();
    if (!event) return;
    
    eventID = event->GetEventID();
    trackID = track->GetTrackID();
    stepNumber = track->GetCurrentStepNumber();

    pdgCode = track->GetDefinition()->GetPDGEncoding();

    // Get position and momentum
    G4ThreeVector position = track->GetPosition(); //where the particle is after the step is finished
    G4ThreeVector momentum = track->GetMomentum();

    // Calculate theta (angle w.r.t. Z-axis) and delta phi (bending angle)
    theta = momentum.theta() / deg;
    phi = momentum.phi() / deg;

    /////////////////////////////////////////////
    // MCS theoretical angle calculation
    const G4Material* material = step->GetPreStepPoint()->GetMaterial();
    G4double radLength = material->GetRadlen();  // Radiation length
    G4double beta = step->GetPreStepPoint()->GetBeta();
    G4double p = step->GetPreStepPoint()->GetMomentum().mag();
    G4double charge = step->GetTrack()->GetDefinition()->GetPDGCharge();
    
    // Highland formula for MCS angle (approximate)
    G4double stepLength = step->GetStepLength();
    theta0 = 0.;
    if (radLength > 0. && stepLength > 0.) {
        theta0 = (13.6 * MeV) / (beta * p) * charge *
            std::sqrt(stepLength / radLength) * 
            (1.0 + 0.038 * std::log(stepLength / radLength));
    }

    //G4cout << "MCS angle (theta0): " << theta0/deg << " deg" << G4endl;
    /////////////////////////////////////////////
    // Get the current position --> The position at the start of the step (i.e., the PreStepPoint).
    // This is where the particle was before it took the step.
    G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();
    // Get field manager for the current logical volume
    const G4FieldManager* fieldManager = step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetFieldManager();
    G4ThreeVector BfieldValue(0,0,0);
    if (fieldManager && fieldManager->GetDetectorField()) {
        G4double point[4] = {pos.x(), pos.y(), pos.z(), 0.}; // time=0
        G4double B[3] = {0,0,0};
        fieldManager->GetDetectorField()->GetFieldValue(point, B);
        BfieldValue.set(B[0], B[1], B[2]);
    }

    // Fill variables for the tree
    x = position.x()/mm;
    y = position.y()/mm;
    z = position.z()/mm;
    px = momentum.x()/MeV;
    py = momentum.y()/MeV;
    pz = momentum.z()/MeV;
    // theta, phi, theta0 already calculated above
    Bx = BfieldValue.x()/tesla;
    By = BfieldValue.y()/tesla;
    Bz = BfieldValue.z()/tesla;

    fStepTree->Fill();
}