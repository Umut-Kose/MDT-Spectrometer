
#include "EventDisplayTEve.h"


  
Int_t EventNumber =0;
Int_t CurrentEventNumber;
TTree* gEventTree = nullptr;
TGNumberEntry* fNumberEntry = nullptr;


//*********************************************//
//*********************************************//
void LoadGeometry(const std::string& gdmlFile)
{
  // Load GDML into ROOT TGeoManager
  TGeoManager::Import(gdmlFile.c_str());
  if (!gGeoManager) {
    std::cerr << "Failed to load GDML geometry: " << gdmlFile << std::endl;
    return;
  }
  TEveGeoTopNode* topNode = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
  gEve->AddGlobalElement(topNode);
}
//*********************************************//
//*********************************************//
void GetDetector(const std::string& gdmlFile)
{
  std::cout << "Starting GetDetector()" << std::endl;
  TGeoManager::Import(gdmlFile.c_str());
  //
  if (!gGeoManager) {
    std::cerr << "Failed to import GDML file." << std::endl;
    return;
  }
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
  gEve->Redraw3D(kTRUE);
  std::cout << "GetDetector() completed." << std::endl;
}
//*********************************************//
//*********************************************//
void AddAxes()
{
  TEveLine* xAxis = new TEveLine(2);
  xAxis->SetLineColor(kRed); xAxis->SetLineWidth(2);
  xAxis->SetPoint(0, 0, 0, 0);
  xAxis->SetPoint(1, 50, 0, 0);
  xAxis->SetTitle("X Axis");
  gEve->AddElement(xAxis);
  
  TEveLine* yAxis = new TEveLine(2);
  yAxis->SetLineColor(kGreen); yAxis->SetLineWidth(2);
  yAxis->SetPoint(0, 0, 0, 0);
  yAxis->SetPoint(1, 0, 50, 0);
  yAxis->SetTitle("Y Axis");
  gEve->AddElement(yAxis);
  
  TEveLine* zAxis = new TEveLine(2);
  zAxis->SetLineColor(kBlue); zAxis->SetLineWidth(2);
  zAxis->SetPoint(0, 0, 0, 0);
  zAxis->SetPoint(1, 0, 0, 50); 
  zAxis->SetTitle("Z Axis");
  gEve->AddElement(zAxis);
  gEve->Redraw3D(kTRUE);
}
//*********************************************//
//*********************************************//
void DisplayEventHits(int targetEventID, TTree *tree)
{
  std::vector<double>* x = nullptr;
  std::vector<double>* y = nullptr;
  std::vector<double>* z = nullptr;
  Int_t eventID = -1;
  
  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("z", &z);
  tree->SetBranchAddress("eventID", &eventID);
  
    
  TEvePointSet* points = new TEvePointSet("Hits");
  points->SetMarkerColor(kBlue);
  points->SetMarkerStyle(4);
  points->SetMarkerSize(1.0);
  // Convert cm --> mm 
  
  TEveBoxSet* hitsBoxSet = new TEveBoxSet("Hits");
  hitsBoxSet->SetPickable(kTRUE);
  hitsBoxSet->Reset(TEveBoxSet::kBT_AABox, true, 64);  // Axis-aligned boxes
  hitsBoxSet->SetMainColor(kBlue); // or kRed, etc.
  //hitsBoxSet->SetMainTransparency(20); // 0 = opaque, 100 = invisible
  hitsBoxSet->SetOwnIds(kTRUE); 
  
  // Size of each hit box (e.g., 5 mm cube)
  const float halfSize = 5; // in mm
  for (Long64_t i = 0; i < tree->GetEntries(); ++i)
    {
      tree->GetEntry(i);
      
      if (eventID != targetEventID) continue;
      for (size_t j = 0; j < x->size(); ++j)
	{
	  std::cout << j << " "
		    << x->at(j) << " "
		    << y->at(j) << " "
		    << z->at(j) << std::endl;
	  float cx = x->at(j);
	  float cy = y->at(j);// * 10.0f;
	  float cz = z->at(j);// * 10.0f; 
	  hitsBoxSet->AddBox(cx, cy, cz, halfSize, halfSize, halfSize);
	  
	}
    }
  gEve->AddElement(hitsBoxSet);
}
//*********************************************//
//*********************************************//
void DrawEvent(Long64_t ientry)
{
  if (!gEve) {
    std::cerr << "TEveManager is not initialized!" << std::endl;
    return;
  }
  
  TEveEventManager* currEvent = gEve->GetCurrentEvent();
  if( currEvent ) currEvent->DestroyElements();
  
  DisplayEventHits(ientry, gEventTree);
  
  gEve->Redraw3D(kTRUE);
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void SetEventNumber()
{
      std::cout << "SetEventNumber() called!" << std::endl;
    if (!fNumberEntry) {
        std::cerr << "ERROR: fNumberEntry is null!" << std::endl;
        return;
    }
    EventNumber = fNumberEntry->GetNumber();
    std::cout << "Loading event " << EventNumber << std::endl;
    DrawEvent(EventNumber);
    
    //EventNumber = fNumberEntry->GetNumber();
    //std::cout << "Event Number is " << EventNumber << std::endl;
  //if(!ApplyIsolation) 
  //DrawEvent(EventNumber);
  //if(ApplyIsolation) DrawIsolated(EventNumber);
  
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void NextEvent()
{
  CurrentEventNumber = EventNumber;
  EventNumber = CurrentEventNumber + 1;
  std::cout<< "Next Event Number is " << EventNumber << std::endl;
  //if(!ApplyIsolation) 
  DrawEvent(EventNumber);
  //if(ApplyIsolation) DrawIsolated(EventNumber);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void PreviousEvent()
{
  CurrentEventNumber = EventNumber;
  EventNumber = CurrentEventNumber - 1;
  std::cout<< "Previous Event Number is " << EventNumber << std::endl;
  //if(!ApplyIsolation) 
  DrawEvent(EventNumber);
  //if(ApplyIsolation) DrawIsolated(EventNumber);
  
}


//*********************************************//
//*********************************************//
void DoExit()
{
  //cout << "Exit application..." << endl;
  gROOT->Reset();
  gApplication->Terminate(0);
}
//*********************************************//
//*********************************************//
void GetDisplay()
{

  //MyGuiHelper* gHelper = nullptr;

  // Event Display Manager
  TEveManager::Create(kTRUE,"V");
  TEveViewer *ev;
  TGLViewer *gv;
  ev = gEve->GetDefaultViewer();
  gEve->GetWindowManager()->HideAllEveDecorations();
  gv = ev->GetGLViewer();
  gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
  gv->SetClearColor(kBlack);

  TEveBrowser * browser = gEve->GetBrowser();
  browser->StartEmbedding(TRootBrowser::kLeft);
  
  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  frmMain->SetWindowName("XX GUI");
  frmMain->SetCleanup(kDeepCleanup);

  //gHelper = new MyGuiHelper();

  TGTextButton *fb = nullptr;
  int posy=0;
  TGVerticalFrame* hf = new TGVerticalFrame(frmMain);
  {
    TGGroupFrame *fGroupFrame2 = new TGGroupFrame(hf,"Event Display");
    fGroupFrame2->SetLayoutBroken(kTRUE);
    {
      posy = 20;
      fNumberEntry = new TGNumberEntry(fGroupFrame2, (Double_t) 0,6,-1,(TGNumberFormat::EStyle) 5,(TGNumberFormat::EAttribute) 1,(TGNumberFormat::ELimit) 2,0,10000);
      fGroupFrame2->AddFrame(fNumberEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,62,2));
      fNumberEntry->MoveResize(20,posy,90,18);
      fb = new TGTextButton(fGroupFrame2, "Show Event#");
      fb->SetCommand("SetEventNumber()");
      //fb->Connect("Clicked()", 0, 0, "SetEventNumber()");

      fGroupFrame2->AddFrame(fb, new TGLayoutHints(kLHintsExpandX));
      fb->MoveResize(120,posy,90,18);
      //
      posy+=28;
      fb = new TGTextButton(fGroupFrame2,"Next");
      //fb->SetCommand("NextEvent()");
      fb->Connect("Clicked()", 0, 0, "NextEvent()");

      fb->MoveResize(10,posy,95,18);
      fb->SetToolTipText("Show next event");
      //
      fb = new TGTextButton(fGroupFrame2,"Previous");
      fb->MoveResize(110,posy,95,18);
      //fb->SetCommand("PreviousEvent()");
      fb->Connect("Clicked()", 0, 0, "PreviousEvent()");

      fb->SetToolTipText("Show previous event");
    }
    hf->AddFrame(fGroupFrame2, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fGroupFrame2->MoveResize(0,200,240,105);
    // Group4
    TGGroupFrame *fGroupFrame4 = new TGGroupFrame(hf,"Display");
    fGroupFrame4->SetLayoutBroken(kTRUE);
    {
      posy+=21;
      fb = new TGTextButton(fGroupFrame4,"Exit");
      fb->MoveResize(55,posy,100,20);
      fb->Connect("Clicked()", "TApplication", gApplication, "Terminate()");

      //fb->SetCommand("DoExit()");
    }
    hf->AddFrame(fGroupFrame4, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fGroupFrame4->MoveResize(0,260,240,240);
  }
  frmMain->AddFrame(hf);
  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();

  browser->StopEmbedding();
  browser->SetTabTitle("Main", 0);
  //
  // Draw Display
  gEve->Redraw3D(kTRUE);	

  GetDetector("detector_geometry.gdml");
  AddAxes(); 

}

//*********************************************//
//*********************************************//
TTree* ReadDataFile(const std::string& rootFileName)
{
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
        file->Close();
        return nullptr;
    }
    
  return tree;
}
//*********************************************//
//*********************************************//
