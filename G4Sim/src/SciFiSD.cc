#include "SciFiSD.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "Randomize.hh"  // For G4RandGauss


SciFiSD::SciFiSD(const G4String& name)
  : G4VSensitiveDetector(name) {}

SciFiSD::~SciFiSD() {}

//void SciFiSD::Clear() {
//    fSmearedPositions.clear();
//}


G4bool SciFiSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
  G4StepPoint* pre = step->GetPreStepPoint();
  G4TouchableHandle touchable = pre->GetTouchableHandle();
  G4Track* track = step->GetTrack();
  
  G4ThreeVector pos = pre->GetPosition();
  G4ThreeVector momentum = pre->GetMomentum();
  
  G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0) return false;	
  
  // Smearing: 50 μm Gaussian resolution on x and y
  const double sigma = 50.0 * micrometer;
  G4double xSmear = G4RandGauss::shoot(pos.x(), sigma);
  G4double ySmear = G4RandGauss::shoot(pos.y(), sigma);
  G4double zSmear = pos.z();  // Assume perfect z for simplicity

  G4int layerID = touchable->GetCopyNumber();

  //G4ThreeVector smeared(xSmear, ySmear, zSmear);
  //fSmearedPositions.push_back(smeared);
  
  //print smeared hit
  //G4cout << "SciFi smeared hit: " << smeared << G4endl;
  
  SciFiHit hit;	
  hit.smearedPos = G4ThreeVector(xSmear, ySmear, zSmear);
  hit.trueMomentum = momentum;
  hit.trackID = track->GetTrackID();
  hit.pdgID = track->GetDefinition()->GetPDGEncoding();
  hit.layerID = layerID;
  fHits.push_back(hit);
  
  
  return true;
}

