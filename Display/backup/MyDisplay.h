#ifndef MYDISPLAY_H
#define MYDISPLAY_H

#include <TROOT.h>
#include <TSystem.h>

#include <TEveManager.h>
#include <TTree.h>
#include <vector>
#include <TApplication.h>
#include <TGeoManager.h>
#include <TFile.h>
#include <TBranch.h>

#include <TEveManager.h>
#include "TEveWindowManager.h"
#include "TEveEventManager.h"
#include <TEveGeoNode.h>
#include <TEveGeoShape.h>
#include <TEveBoxSet.h>
#include <TEveViewer.h>
#include "TEveBrowser.h"
#include "TEveViewer.h"
#include <TEvePointSet.h>
#include <TEveLine.h>


#include <TGNumberEntry.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGFrame.h>
#include <TGLViewer.h>
#include "TGLUtil.h"
#include "TGTab.h"

#include "TColor.h"

#include <iostream>
#include <vector>
#include <string>
#include <RQ_OBJECT.h>


class MyDisplay {
public:
    MyDisplay();
    virtual ~MyDisplay();
    
    void Initialize(const std::string& rootFile, const std::string& gdmlFile);
    void Run();
    
    // Button handlers
    void SetEventNumber();
    void NextEvent();
    void PreviousEvent();
    void DoExit();

private:
    void LoadGeometry(const std::string& gdmlFile);
    void GetDetector(const std::string& gdmlFile);
    void AddAxes();
    void DisplayEventHits(int targetEventID);
    void DrawEvent(Long64_t ientry);
    TTree* ReadDataFile(const std::string& rootFileName);
    void CreateGUI();
  void SetupViewer();
    // Data members
    TTree* fEventTree;
    TGNumberEntry* fNumberEntry;
    Int_t fEventNumber;
    TEveElementList* fHitElements;

  ClassDef(MyDisplay, 1);
};

#endif
