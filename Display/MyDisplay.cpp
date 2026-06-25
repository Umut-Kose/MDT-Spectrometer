 /*
 * MyDisplay.cpp - Event Display for Muon Spectrometer Simulation
 * 
 * Detector Configuration:
 * - 4 Tracking Stations (each with 4 layers of MDTs - Muon Drift Tubes)
 * - 3 Magnet Iron pieces (40 cm thick each)
 * - Total: 16 MDT layers + 3 magnets
 * 
 * Version: MDT configuration (updated from SciFi detectors)
 * Backup of SciFi version available in: backup_scifi_version/
 */

#include "MyDisplay.h"
#include <TSystem.h>
#include <iostream>
#include "TVector3.h"
#include <map>
#include <vector>
#include <cmath>

MyDisplay::MyDisplay() : 
  fEventTree(nullptr),
  fNumberEntry(nullptr),
  fEventNumber(0),
  fHitElements(nullptr),
  fDetectorElements(nullptr),
  ApplyIsolation(kFALSE)
{
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MyDisplay::~MyDisplay()
{
  if (fHitElements) {
    delete fHitElements;
  }
  if (fDetectorElements) {
    gEve->RemoveElement(fDetectorElements, nullptr);
    delete fDetectorElements;
  }
  fInputFile->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MyDisplay::Initialize(const std::string& rootFile, const std::string& gdmlFile)
{
  // Load data
  fEventTree = ReadDataFile(rootFile);
  if (!fEventTree) {
    throw std::runtime_error("Failed to load event data");
  }
  
  // Initialize EVE
  TEveManager::Create(kTRUE, "V");

  TEveBrowser* browser = gEve->GetBrowser();
  browser->StartEmbedding(TRootBrowser::kLeft);

  CreateGUI();

  browser->StopEmbedding();
  browser->SetTabTitle("Control", 0);
  
  TEveWindowSlot* slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
  TEveWindowPack* pack = slot->MakePack();
  pack->SetElementName("3D View");
  pack->SetShowTitleBar(kFALSE);
  TEveViewer* viewer = new TEveViewer("3DView");
  pack->NewSlot()->ReplaceWindow(viewer);
  viewer->SpawnGLViewer(nullptr);

  viewer->AddScene(gEve->GetGlobalScene());
  viewer->AddScene(gEve->GetEventScene());
  
  
  // Load geometry
  LoadGeometry(gdmlFile);
  SetupViewer();
  //AddAxes();

  // Configure viewer
  TGLViewer* glViewer = viewer->GetGLViewer();
  glViewer->SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
  glViewer->SetGuideState(TGLUtil::kAxesEdge, kTRUE, kFALSE, 0);
  glViewer->SetClearColor(kBlack);
  glViewer->ResetCameras();
  
  gEve->GetBrowser()->GetTabRight()->MapWindow();

  
  gEve->Redraw3D(kTRUE);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TTree* MyDisplay::ReadDataFile(const std::string& rootFileName)
{
  std::cout << "Attempting to open: " << rootFileName << std::endl;

  fInputFile = TFile::Open(rootFileName.c_str());
  if (!fInputFile || fInputFile->IsZombie()) {
    std::cerr << "Error opening file: " << rootFileName << std::endl;
    std::cerr << "Current working directory: ";
    gSystem->Exec("pwd");
    return nullptr;
  }

  TTree* tree = dynamic_cast<TTree*>(fInputFile->Get("Hits"));
  if (!tree) {
    std::cerr << "Error: 'Hits' tree not found in file" << std::endl;
    fInputFile->ls();
    fInputFile->Close();
    return nullptr;
  }

  // Verify tree has entries
  if (tree->GetEntries() == 0) {
    std::cerr << "WARNING: Tree has 0 entries" << std::endl;
  }
  std::cerr << "WHATTTT" << std::endl;
  return tree;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MyDisplay::SetupViewer() {
  TEveViewer *viewer = gEve->GetDefaultViewer();
  TGLViewer *glviewer = viewer->GetGLViewer();
  
  // Set camera type
  glviewer->SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
  
  // Configure viewport
  glviewer->SetClearColor(kBlack);
  glviewer->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
  
  // Auto-scale to content
  glviewer->ResetCameras();
  glviewer->RefreshPadEditor();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MyDisplay::LoadGeometry(const std::string& gdmlFile) {
  std::cout << "Loading geometry from: " << gdmlFile << std::endl;
  
  if (gGeoManager) {
    delete gGeoManager;
  }
  TGeoManager::Import(gdmlFile.c_str());
  if (!gGeoManager) {
    throw std::runtime_error("Failed to load GDML geometry");
  }
  if (!gGeoManager->GetTopNode()) {
    throw std::runtime_error("No top node in geometry");
  }
  //TEveGeoTopNode* topNode = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
  //gEve->AddGlobalElement(topNode);
  std::cout << "GDML file imported successfully." << std::endl;
  //
  TGeoVolume* gdmlTop = gGeoManager->GetTopVolume();
  if (!gdmlTop) {
    std::cerr << "Failed to get top volume." << std::endl;
    return;
  }
  //std::cout << "Top volume set: " << gdmlTop->GetName() << std::endl;
  TGeoIterator nextNode(gdmlTop);
  TGeoNode* curNode;

  fDetectorElements = new TEveElementList("Detector");
  //
  while ((curNode = nextNode())) {
    TGeoVolume* vol = curNode->GetVolume();
    if (!vol) {
      std::cerr << "Volume is null for node: " << curNode->GetName() << std::endl;
      continue;
    }
    //
    TGeoShape* shape = vol->GetShape();
    if (!shape) {
      std::cerr << "Shape is null for volume: " << vol->GetName() << std::endl;
      continue;
    }
    //
    TString nodeName(curNode->GetName());
    TString nodePath;
    nextNode.GetPath(nodePath);
    const TGeoMatrix* matrix = nextNode.GetCurrentMatrix();
    if (!matrix) {
      std::cerr << "Matrix is null for node: " << curNode->GetName() << std::endl;
      continue;
    }
    //
    const Double_t* trans = matrix->GetTranslation();
    const Double_t* rotMatrix = matrix->GetRotationMatrix();
    if (!trans || !rotMatrix) {
      std::cerr << "Transformation matrix is null for node: " << curNode->GetName() << std::endl;
      continue;
    }
    //
    TGeoRotation rotation;
    rotation.SetMatrix(rotMatrix);
    TGeoCombiTrans transform(trans[0], trans[1], trans[2], &rotation);
    // Create the TEveGeoShape for visualization
    TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
    eveShape->SetShape(shape);
    eveShape->SetMainTransparency(90); // Set transparency
    eveShape->SetTransMatrix(transform);
    
    // Color-coding based on volume name
    // New detector configuration: 4 layers of MDTs per station, 4 stations total, 3 magnet iron pieces (40cm thick)
    TString volName = vol->GetName();
    
    if (volName.Contains("MDT") || volName.Contains("mdt"))
      eveShape->SetMainColor(kGreen + 1);  // MDT layers in green
    else if (volName.Contains("Magnet") || volName.Contains("magnet"))
      eveShape->SetMainColor(kGray + 1);   // Magnet iron in gray
    else if (volName.Contains("AluStrip"))
      eveShape->SetMainColor(kOrange + 7);
    else
      eveShape->SetMainColor(kBlue);
    
    
    fDetectorElements->AddElement(eveShape);
  }
  gEve->AddGlobalElement(fDetectorElements);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MyDisplay::AddAxes()
{
  TEveLine* xAxis = new TEveLine(2);
  xAxis->SetLineColor(kRed); xAxis->SetLineWidth(2);
  xAxis->SetPoint(0, 0, 0, 0); xAxis->SetPoint(1, 50, 0, 0);
  gEve->AddElement(xAxis);
  
  TEveLine* yAxis = new TEveLine(2);
  yAxis->SetLineColor(kGreen); yAxis->SetLineWidth(2);
  yAxis->SetPoint(0, 0, 0, 0); yAxis->SetPoint(1, 0, 50, 0);
  gEve->AddElement(yAxis);
  
  TEveLine* zAxis = new TEveLine(2);
  zAxis->SetLineColor(kBlue); zAxis->SetLineWidth(2);
  zAxis->SetPoint(0, 0, 0, 0); zAxis->SetPoint(1, 0, 0, 50);
  gEve->AddElement(zAxis);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MyDisplay::CreateGUI()
{
  TEveBrowser* browser = gEve->GetBrowser();
  browser->StartEmbedding(TRootBrowser::kLeft);
  
  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1200, 800);
  frmMain->SetWindowName("Event Display");
  
  TGVerticalFrame* vf = new TGVerticalFrame(frmMain);
  
  // Event control group
  TGGroupFrame* gfEvents = new TGGroupFrame(vf, "Event Control");
  {
    fNumberEntry = new TGNumberEntry(gfEvents, 0, 6, -1, 
				     TGNumberFormat::kNESInteger,
				     TGNumberFormat::kNEAAnyNumber,
				     TGNumberFormat::kNELLimitMinMax, 0, 10000);
    
    TGTextButton* showBtn = new TGTextButton(gfEvents, "Show Event");
    showBtn->Connect("Clicked()", "MyDisplay", this, "SetEventNumber()");
    
    TGTextButton* nextBtn = new TGTextButton(gfEvents, "Next");
    nextBtn->Connect("Clicked()", "MyDisplay", this, "NextEvent()");
    
    TGTextButton* prevBtn = new TGTextButton(gfEvents, "Previous");
    prevBtn->Connect("Clicked()", "MyDisplay", this, "PreviousEvent()");
    
    // Layout
    gfEvents->AddFrame(fNumberEntry, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
    gfEvents->AddFrame(showBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    gfEvents->AddFrame(nextBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    gfEvents->AddFrame(prevBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
  }
  vf->AddFrame(gfEvents, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

  fIsolate = new TGCheckButton(vf, "Hide Detector Geometry");
  fIsolate->Connect("Clicked()", "MyDisplay", this, "Isolate()");
  vf->AddFrame(fIsolate, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

  fbackground = new TGTextButton(vf, "Change Background Colour");
  fbackground->Connect("Clicked()", "MyDisplay", this, "BackgroundColor()");
  vf->AddFrame(fbackground, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

  
  // Exit button
  TGTextButton* exitBtn = new TGTextButton(vf, "Exit");
  exitBtn->Connect("Clicked()", "MyDisplay", this, "DoExit()");
  vf->AddFrame(exitBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
  
  frmMain->AddFrame(vf);
  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();
  
  browser->StopEmbedding();
  browser->SetTabTitle("Control", 0);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MyDisplay::SetEventNumber()
{
  fEventNumber = fNumberEntry->GetNumber();
  std::cout << "Loading event " << fEventNumber << std::endl;
  DrawEvent(fEventNumber);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MyDisplay::NextEvent()
{
  fEventNumber++;
  fNumberEntry->SetNumber(fEventNumber);
  std::cout << "Loading event " << fEventNumber << std::endl;
  DrawEvent(fEventNumber);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MyDisplay::PreviousEvent()
{
  fEventNumber--;
  fNumberEntry->SetNumber(fEventNumber);
  std::cout << "Loading event " << fEventNumber << std::endl;
  DrawEvent(fEventNumber);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MyDisplay::Isolate()
{
  if (!fDetectorElements) {
    std::cerr << "Error: Detector elements not initialized!" << std::endl;
    return;
  }
  bool hideGeometry = fIsolate->IsOn();
  fDetectorElements->SetRnrState(!hideGeometry);
  gEve->Redraw3D(kTRUE);
  std::cout << "Detector geometry " << (hideGeometry ? "hidden" : "shown") << std::endl;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MyDisplay::BackgroundColor()
{
  if (!gEve->GetDefaultGLViewer()) return;
  // Toggle between black and white background
  static Color_t bgColor = kBlack;
  bgColor = (bgColor == kBlack) ? kWhite : kBlack;
  
  gEve->GetDefaultGLViewer()->SetClearColor(bgColor);
  gEve->FullRedraw3D(kTRUE);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MyDisplay::DoExit()
{
  gApplication->Terminate(0);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MyDisplay::DrawEvent(Long64_t ientry)
{
  if (!gEve) {
    std::cerr << "TEveManager is not initialized!" << std::endl;
    return;
  }
  
  TEveEventManager* currEvent = gEve->GetCurrentEvent();
  if( currEvent ) currEvent->DestroyElements();    
  
  fHitElements = new TEveElementList("Hits");
  fHitElements->SetRnrSelf(kTRUE);
  // Display new event
  DisplayEventHits(ientry);
  gEve->AddElement(fHitElements);
 
  gEve->Redraw3D(kTRUE);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MyDisplay::DisplayEventHits(int targetEventID)
{
  std::cout << "Displaying hits for event " << targetEventID << std::endl;

  std::vector<double>* x = nullptr;
  std::vector<double>* y = nullptr;
  std::vector<double>* z = nullptr;
  std::vector<int>* pdg = nullptr;
  Int_t eventID = -1;
  std::vector<int>* trackID = nullptr;
  std::vector<double>* px = nullptr;
  std::vector<double>* py = nullptr;
  std::vector<double>* pz = nullptr;
  
  // MDT simulation uses trueX, trueY, trueZ for actual hit positions
  fEventTree->SetBranchAddress("trueX", &x);
  fEventTree->SetBranchAddress("trueY", &y);
  fEventTree->SetBranchAddress("trueZ", &z);
  fEventTree->SetBranchAddress("pdg", &pdg);
  fEventTree->SetBranchAddress("eventID", &eventID);
  fEventTree->SetBranchAddress("trackID", &trackID);
  fEventTree->SetBranchAddress("px", &px);
  fEventTree->SetBranchAddress("py", &py);
  fEventTree->SetBranchAddress("pz", &pz);

  //TEveBoxSet* hits = new TEveBoxSet("Hits");
  //hits->Reset(TEveBoxSet::kBT_AABox, true, 64);
  //hits->SetMainColor(kRed);  // More visible color
  //hits->SetMainTransparency(0);  // Slightly transparent
  //hits->SetPickable(kTRUE);

  TGeoMedium *air = gGeoManager->GetMedium("AIR");
  // Larger hit markers for better visibility (1 cm cube)
  TGeoShape *box = new TGeoBBox("box", 1/2.0, 1/2.0, 1/2.0);

  TEveElementList* hitList = new TEveElementList("Hits");
  
  // Store momentum info per track (trackID -> momentum vector)
  std::map<int, TVector3> trackMomentum;
  std::map<int, int> trackPDG;

  for (Long64_t i = 0; i < fEventTree->GetEntries(); ++i)
    {
      fEventTree->GetEntry(i);
      if (eventID != targetEventID) continue;
      
      std::cout << "\n=== Event " << targetEventID << " Hits ===" << std::endl;
      
      for (size_t j = 0; j < x->size(); ++j)
      {
        int trkid = trackID->at(j);
        int pdgCode = pdg->at(j);
        double p_x = px->at(j);
        double p_y = py->at(j);
        double p_z = pz->at(j);
        double p_tot = std::sqrt(p_x*p_x + p_y*p_y + p_z*p_z);
        
        // Store momentum for each track (first hit)
        if (trackMomentum.find(trkid) == trackMomentum.end()) {
          trackMomentum[trkid] = TVector3(p_x, p_y, p_z);
          trackPDG[trkid] = pdgCode;
          
          // Print momentum info for muons (PDG = ±13)
          if (std::abs(pdgCode) == 13) {
            std::cout << "Track " << trkid << " (Muon, PDG=" << pdgCode << "):" << std::endl;
            std::cout << "  Momentum: Ptot = " << p_tot << " MeV/c" << std::endl;
            std::cout << "            Px = " << p_x << " MeV/c" << std::endl;
            std::cout << "            Py = " << p_y << " MeV/c" << std::endl;
            std::cout << "            Pz = " << p_z << " MeV/c" << std::endl;
          }
        }
        
        // Create hit marker
        TGeoTranslation *trans = new TGeoTranslation(x->at(j) / 10.0, y->at(j) / 10.0, z->at(j) / 10.0);
        TGeoVolume* hitVolume = new TGeoVolume("HitVolume", box, air);
        hitVolume->SetLineColor(kRed); 
        TEveGeoShape* eveShape = new TEveGeoShape(hitVolume->GetName());
        eveShape->SetShape(hitVolume->GetShape());
        eveShape->SetMainColor(hitVolume->GetLineColor());
        eveShape->SetTransMatrix(*trans);
        hitList->AddElement(eveShape);
      }
    }
  
  // Display momentum vectors for muon tracks
  for (const auto& [trkid, pVec] : trackMomentum) {
    if (std::abs(trackPDG[trkid]) == 13) {
      // Create momentum vector arrow (scaled for visibility)
      TEveLine* momArrow = new TEveLine(2);
      momArrow->SetLineColor(kBlue);
      momArrow->SetLineWidth(3);
      momArrow->SetLineStyle(2); // Dashed line
      
      // Scale momentum vector for display (adjust scale factor as needed)
      double scale = 0.01; // 1 cm per 100 MeV/c
      momArrow->SetPoint(0, 0, 0, 0);
      momArrow->SetPoint(1, pVec.X()*scale, pVec.Y()*scale, pVec.Z()*scale);
      momArrow->SetTitle(Form("Track %d: P=%.1f MeV/c", trkid, pVec.Mag()));
      hitList->AddElement(momArrow);
    }
  }

  fHitElements->AddElement(hitList);
  //gEve->AddElement(fHitElements);
  gEve->Redraw3D(kTRUE);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MyDisplay::Run()
{
  if (gApplication) {
    gApplication->Run();
    
  }
}
