#include "EventAction.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4AnalysisManager.hh"
#include "G4SDManager.hh"
#include "MDTSD.hh"


EventAction::EventAction() : G4UserEventAction() {
  fRootFile = new TFile("mdt_hits.root", "RECREATE");
  fTree = new TTree("Hits", "MDT Hits");
  
  // Branch for event ID
  fTree->Branch("eventID", &feventID);
  
  // Branches for particle information
  fTree->Branch("trackID", &ftrackID);
  fTree->Branch("pdg", &fpdg);
  
  // Branches for detector geometry
  fTree->Branch("stationID", &fstationID);
  fTree->Branch("planeID", &fplaneID);
  fTree->Branch("tubeID", &ftubeID);
  
  // Branches for tube center position
  fTree->Branch("tubeCenterX", &fTubeCenterX);
  fTree->Branch("tubeCenterY", &fTubeCenterY);
  fTree->Branch("tubeCenterZ", &fTubeCenterZ);
  
  // Branches for drift information
  fTree->Branch("hitX",           &fHitX);           // X position along tube
  fTree->Branch("trueDriftRadius",&fTrueDriftRadius); // True drift radius (smear at reco level)
  fTree->Branch("driftAngle",     &fDriftAngle);      // MC truth: angle in Y-Z plane
  
  // Branches for MC truth (true hit position)
  fTree->Branch("trueX", &fTrueX);
  fTree->Branch("trueY", &fTrueY);
  fTree->Branch("trueZ", &fTrueZ);
  
  // Branches for particle momentum
  fTree->Branch("px", &fpx);
  fTree->Branch("py", &fpy);
  fTree->Branch("pz", &fpz);
  
  // Branches for energy deposit and drift time
  fTree->Branch("edep", &fedep);
  fTree->Branch("driftTime", &fdriftTime);
}

EventAction::~EventAction() {
  if (fRootFile) {
    fRootFile->cd();

    if (fTree) {
      fTree->Write();
    }

    fRootFile->Close();
    delete fRootFile;
    fRootFile = nullptr;
  }
}

void EventAction::BeginOfEventAction(const G4Event* evt) {
  // Reset per-event data
  feventID = evt->GetEventID();

  ftrackID.clear();
  fpdg.clear();

  fstationID.clear();
  fplaneID.clear();
  ftubeID.clear();

  fTubeCenterX.clear();
  fTubeCenterY.clear();
  fTubeCenterZ.clear();

  fHitX.clear();
  fTrueDriftRadius.clear();
  fDriftAngle.clear();
  fdriftTime.clear();
  fedep.clear();

  fTrueX.clear();
  fTrueY.clear();
  fTrueZ.clear();

  fpx.clear();
  fpy.clear();
  fpz.clear();

  auto sdManager = G4SDManager::GetSDMpointer();
  auto mdtSD = static_cast<MDTSD*>(
      sdManager->FindSensitiveDetector("MDTSD", false));

  if (mdtSD) {
    mdtSD->ClearHits();  // Defensive clear; Initialize() handles this via Geant4 SD lifecycle
  }
}

void EventAction::EndOfEventAction(const G4Event* evt) {
  auto sdManager = G4SDManager::GetSDMpointer();
  auto mdtSD = static_cast<MDTSD*>(sdManager->FindSensitiveDetector("MDTSD", false));
  
  if (!mdtSD) return;
  
  const auto& hits = mdtSD->GetHits();
  for (const auto& hit : hits) {
    // Tube center position
    fTubeCenterX.push_back(hit.tubeCenter.x());
    fTubeCenterY.push_back(hit.tubeCenter.y());
    fTubeCenterZ.push_back(hit.tubeCenter.z());
    
    // Drift information
    fHitX.push_back(hit.hitX);                          // X position along tube
    fTrueDriftRadius.push_back(hit.trueDriftRadius);    // True radius (smear at reco)
    fDriftAngle.push_back(hit.driftAngle);              // MC truth: drift angle
    
    // True hit position (MC truth)
    fTrueX.push_back(hit.truePosition.x());
    fTrueY.push_back(hit.truePosition.y());
    fTrueZ.push_back(hit.truePosition.z());
    
    // Particle momentum
    fpx.push_back(hit.trueMomentum.x());
    fpy.push_back(hit.trueMomentum.y());
    fpz.push_back(hit.trueMomentum.z());
    
    // Particle information
    fpdg.push_back(hit.pdgID);
    ftrackID.push_back(hit.trackID);
    
    // Detector IDs
    fstationID.push_back(hit.stationID);
    fplaneID.push_back(hit.planeID);
    ftubeID.push_back(hit.tubeID);
    
    // Energy deposit and drift time
    fedep.push_back(hit.edep);
    fdriftTime.push_back(hit.driftTime);
  }
  
  fTree->Fill();
  // Note: hits are cleared at the start of each event via MDTSD::Initialize()
}

