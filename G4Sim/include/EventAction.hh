#ifndef EVENT_ACTION_HH
#define EVENT_ACTION_HH

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "TTree.h"
#include "TFile.h"
#include "G4ThreeVector.hh"
#include <vector>

class EventAction : public G4UserEventAction {
public:
  EventAction();
  virtual ~EventAction();
  
  virtual void BeginOfEventAction(const G4Event*) override;
  virtual void EndOfEventAction(const G4Event*) override;

private:
  TFile* fRootFile = nullptr;
  TTree* fTree = nullptr;
  
  // MDT hit data
  std::vector<double> fTubeCenterX, fTubeCenterY, fTubeCenterZ;  // Tube center positions
  std::vector<double> fHitX;                                      // X position along tube where particle crossed
  std::vector<double> fTrueDriftRadius;                           // True drift radius (no smearing) — smear at reco
  std::vector<double> fDriftAngle;                                // True drift angle in Y-Z plane (MC truth)
  std::vector<double> fTrueX, fTrueY, fTrueZ;                    // True hit position (MC truth)
  std::vector<double> fpx, fpy, fpz;                              // Particle momentum
  std::vector<int> fpdg, ftrackID;                                // Particle info
  std::vector<int> fstationID, fplaneID, ftubeID;                // Detector IDs
  std::vector<double> fedep, fdriftTime;                         // Energy deposit and drift time
  int feventID;
};

#endif

