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
  
  std::vector<double> fx, fy, fz;
  std::vector<double> fpx, fpy, fpz;
  std::vector<int> fpdg, ftrackID, flayerID, fstationID;
  int feventID;
};

#endif

