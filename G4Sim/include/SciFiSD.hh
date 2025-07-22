#ifndef SCIFI_SD_HH
#define SCIFI_SD_HH

#include "G4VSensitiveDetector.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include "globals.hh"

class G4Step;
class G4TouchableHistory;


struct SciFiHit {
  G4ThreeVector smearedPos;
  G4ThreeVector trueMomentum;
  G4int trackID;
  G4int pdgID;
  G4int eventID;
  G4int stationID;
  G4int layerID;
};


class SciFiSD : public G4VSensitiveDetector {
public:
  SciFiSD(const G4String& name);
  virtual ~SciFiSD();
  
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
  //const std::vector<G4ThreeVector>& GetSmearedHits() const { return fSmearedPositions; }
  //void ClearHits() { fSmearedPositions.clear(); }
  const std::vector<SciFiHit>& GetHits() const { return fHits; }
  void ClearHits() { fHits.clear(); }
  
private:
  //std::vector<G4ThreeVector> fSmearedPositions;
  std::vector<SciFiHit> fHits;
};

#endif

