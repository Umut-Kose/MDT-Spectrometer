#include "EventAction.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4AnalysisManager.hh"
#include "G4SDManager.hh"
#include "SciFiSD.hh"


EventAction::EventAction() : G4UserEventAction() {
  fRootFile = new TFile("scifi_hits.root", "RECREATE");
  fTree = new TTree("Hits", "SciFi Hits");
  
  fTree->Branch("eventID", &feventID);
  fTree->Branch("trackID", &ftrackID);
  fTree->Branch("pdg", &fpdg);
  fTree->Branch("layerID", &flayerID);
  fTree->Branch("x", &fx);
  fTree->Branch("y", &fy);
  fTree->Branch("z", &fz);
  fTree->Branch("px", &fpx);
  fTree->Branch("py", &fpy);
  fTree->Branch("pz", &fpz);
}

EventAction::~EventAction() {
  fRootFile->cd();
  fTree->Write();
  fRootFile->Close();
}

void EventAction::BeginOfEventAction(const G4Event* evt) {
  // reset per-event data
  fx.clear(); fy.clear(); fz.clear();
  fpx.clear(); fpy.clear(); fpz.clear();
  fpdg.clear(); ftrackID.clear(); flayerID.clear();
  feventID = evt->GetEventID();
}

void EventAction::EndOfEventAction(const G4Event* evt) {
  auto sdManager = G4SDManager::GetSDMpointer();
  auto scifiSD = static_cast<SciFiSD*>(sdManager->FindSensitiveDetector("SciFiSD"));
  
  if (!scifiSD) return;
  
  const auto& hits = scifiSD->GetHits();
  for (const auto& hit : hits) {
    fx.push_back(hit.smearedPos.x());
    fy.push_back(hit.smearedPos.y());
    fz.push_back(hit.smearedPos.z());
    fpx.push_back(hit.trueMomentum.x());
    fpy.push_back(hit.trueMomentum.y());
    fpz.push_back(hit.trueMomentum.z());
    fpdg.push_back(hit.pdgID);
    ftrackID.push_back(hit.trackID);
    feventID = evt->GetEventID();
    flayerID.push_back(hit.layerID);
  }
  fTree->Fill();
  scifiSD->ClearHits();
}


