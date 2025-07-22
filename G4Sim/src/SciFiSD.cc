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
  
  // Debug output
    G4cout << "Processing hit in volume: " << pre->GetPhysicalVolume()->GetName()
           << " at Z = " << pre->GetPosition().z()/mm << " mm"
           << " with layerID: " << touchable->GetCopyNumber() << G4endl;

  G4ThreeVector pos = pre->GetPosition();
  G4ThreeVector momentum = pre->GetMomentum();
  
  G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0) return false;	

  // Skip steps with minimal energy deposition
  if(step->GetTotalEnergyDeposit() < 0.1*keV) return false;

  // Smearing: 50 μm Gaussian resolution on x and y
  const double sigma = 50.0 * micrometer;
  G4double xSmear = G4RandGauss::shoot(pos.x(), sigma);
  G4double ySmear = G4RandGauss::shoot(pos.y(), sigma);
  G4double zSmear = pos.z();  // Assume perfect z for simplicity

  G4int layerID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(); // from 1 to 40
  G4int stationID = layerID / 4 + 1; // from 1 to 10
  
  
  G4cout << "Hit: station=" << stationID << " layer=" << layerID
       << " trackID=" << track->GetTrackID()
       << " pdgID=" << track->GetDefinition()->GetPDGEncoding()
       << " x=" << xSmear/mm << " y=" << ySmear/mm << " z=" << zSmear/mm
       << " edep=" << edep/keV << " keV"
       << G4endl;


  SciFiHit hit;	
  hit.smearedPos = G4ThreeVector(xSmear, ySmear, zSmear);
  hit.trueMomentum = momentum;
  hit.trackID = track->GetTrackID();
  hit.pdgID = track->GetDefinition()->GetPDGEncoding();
  hit.stationID = stationID;
  hit.layerID = layerID;
  fHits.push_back(hit);
  
  
  return true;
}

