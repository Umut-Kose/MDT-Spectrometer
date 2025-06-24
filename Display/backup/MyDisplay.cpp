#include "MyDisplay.h"
#include <TSystem.h>
#include <iostream>

MyDisplay::MyDisplay() : 
  fEventTree(nullptr),
  fNumberEntry(nullptr),
  fEventNumber(0),
  fHitElements(nullptr) 
{
}

MyDisplay::~MyDisplay()
{
  if (fHitElements) {
    delete fHitElements;
  }
}

void MyDisplay::Initialize(const std::string& rootFile, const std::string& gdmlFile)
{
  // Load data
  fEventTree = ReadDataFile(rootFile);
  if (!fEventTree) {
    throw std::runtime_error("Failed to load event data");
  }
  
  // Initialize EVE
  TEveManager::Create(kTRUE, "V");
  
  TEveWindowSlot* slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
  TEveWindowPack* pack = slot->MakePack();
  pack->SetElementName("Main View");
  pack->SetShowTitleBar(kFALSE);
  pack->NewSlot()->MakeCurrent();
  gEve->GetDefaultViewer()->GetGLViewer()->ResetCameras();
  if (gEve->GetDefaultGLViewer())
    {
      gEve->GetDefaultGLViewer()->SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
      gEve->GetDefaultGLViewer()->ResetCurrentCamera();
      gEve->GetDefaultGLViewer()->RefreshPadEditor();
    }
  
  // Load geometry
  LoadGeometry(gdmlFile);
  SetupViewer();
  AddAxes();
  
  // Create GUI
  CreateGUI();
  gEve->GetBrowser()->GetTabRight()->MapWindow();
  gEve->GetDefaultGLViewer()->ResetCameras();
  gEve->Redraw3D(kTRUE);
}

TTree* MyDisplay::ReadDataFile(const std::string& rootFileName) {
  std::cout << "Attempting to open: " << rootFileName << std::endl;
  TFile* file = TFile::Open(rootFileName.c_str());
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file: " << rootFileName << std::endl;
    std::cerr << "Current working directory: ";
    gSystem->Exec("pwd");
    return nullptr;
  }
  TTree* tree = dynamic_cast<TTree*>(file->Get("Hits"));
  if (!tree) {
    std::cerr << "Error: 'Hits' tree not found in file" << std::endl;
    file->ls();
    file->Close();
    return nullptr;
  }
  // Verify tree has entries
  if (tree->GetEntries() == 0) {
    std::cerr << "WARNING: Tree has 0 entries" << std::endl;
  }
  
  tree->SetDirectory(0); 
  file->Close();
  return tree;
}

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
  TEveElementList* fDetectorElements = new TEveElementList("Detector");
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
    TString volName = vol->GetName();
    
    if (volName.Contains("SciFi"))
      eveShape->SetMainColor(kGreen + 1);
    else if (volName.Contains("Magnet"))
      eveShape->SetMainColor(kGray + 1);
    else if (volName.Contains("AluStrip"))
      eveShape->SetMainColor(kOrange + 7);
    else
      eveShape->SetMainColor(kBlue);
    
    
    fDetectorElements->AddElement(eveShape);
  }
  gEve->AddGlobalElement(fDetectorElements);
}

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

void MyDisplay::SetEventNumber()
{
  fEventNumber = fNumberEntry->GetNumber();
  std::cout << "Loading event " << fEventNumber << std::endl;
  DrawEvent(fEventNumber);
}

void MyDisplay::NextEvent()
{
  fEventNumber++;
  fNumberEntry->SetNumber(fEventNumber);
  std::cout << "Loading event " << fEventNumber << std::endl;
  DrawEvent(fEventNumber);
}

void MyDisplay::PreviousEvent()
{
  fEventNumber--;
  fNumberEntry->SetNumber(fEventNumber);
  std::cout << "Loading event " << fEventNumber << std::endl;
  DrawEvent(fEventNumber);
}

void MyDisplay::DoExit()
{
  gApplication->Terminate(0);
}

void MyDisplay::DrawEvent(Long64_t ientry)
{
  if (!gEve) {
    std::cerr << "TEveManager is not initialized!" << std::endl;
    return;
  }
  
  TEveEventManager* currEvent = gEve->GetCurrentEvent();
  if( currEvent ) currEvent->DestroyElements();    
  
  
  fHitElements = new TEveElementList("Hits");
  
  // Display new event
  DisplayEventHits(ientry);
  gEve->Redraw3D(kTRUE);
}

void MyDisplay::DisplayEventHits(int targetEventID)
{
  std::vector<double>* x = nullptr;
  std::vector<double>* y = nullptr;
  std::vector<double>* z = nullptr;
  Int_t eventID = -1;
  
  fEventTree->SetBranchAddress("x", &x);
  fEventTree->SetBranchAddress("y", &y);
  fEventTree->SetBranchAddress("z", &z);
  fEventTree->SetBranchAddress("eventID", &eventID);
  
  TEveBoxSet* hits = new TEveBoxSet("Hits");
  hits->Reset(TEveBoxSet::kBT_AABox, true, 64);
  hits->SetMainColor(kBlue);
  
  const float halfSize = 5; // mm
  for (Long64_t i = 0; i < fEventTree->GetEntries(); ++i) {
    fEventTree->GetEntry(i);
    if (eventID != targetEventID) continue;
    
    for (size_t j = 0; j < x->size(); ++j) {
      hits->AddBox(x->at(j), y->at(j), z->at(j), halfSize, halfSize, halfSize);
    }
  }
  
  fHitElements->AddElement(hits);
  gEve->AddElement(fHitElements);
}

void MyDisplay::Run()
{
  if (gApplication) {
    gApplication->Run();
    
  }
}
