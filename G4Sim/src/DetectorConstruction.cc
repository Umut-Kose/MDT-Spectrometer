#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4GDMLParser.hh"


// Magnetic field includes
#include "MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"

//MDT detector (Monitored Drift Tubes)
#include "MDTSD.hh"
#include "G4SDManager.hh"

#include <cmath>
#include <cstdio>

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fTubeInnerLV(nullptr)
{}

DetectorConstruction::~DetectorConstruction() {}

void DetectorConstruction::ConstructSDandField()
{
  auto sdManager = G4SDManager::GetSDMpointer();

  auto mdtSD = new MDTSD("MDTSD");
  sdManager->AddNewDetector(mdtSD);

  MDTSD::SetVerbose(1); // 0=silent
                        // 1=stored hits
                        // 2=every ProcessHits
                        // 3=full dump

sdManager->AddNewDetector(mdtSD);

  if (fTubeInnerLV) {
    fTubeInnerLV->SetSensitiveDetector(mdtSD);
    G4cout << "MDTSD attached to TubeInnerLV" << G4endl;
  }

}

// Geometry construction
// ATLAS MDT Configuration: 4 stations × 3 planes each + 3 magnets (400mm thick iron)

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  auto nist = G4NistManager::Instance();

  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* steel = nist->FindOrBuildMaterial("G4_Fe");
  G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");
  G4Material* argon = nist->FindOrBuildMaterial("G4_Ar");

  G4bool checkOverlaps = true;

  // ============================================================
  // World
  // ============================================================

  G4double worldZ = 3000. * mm;  // Increased for 400mm magnets
  G4Box* worldBox = new G4Box("World", 1.3*m, 1.3*m, worldZ/2);
  G4LogicalVolume* worldLV = new G4LogicalVolume(worldBox, air, "WorldLV");
  G4VPhysicalVolume* worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, checkOverlaps);

// ============================================================
  // MDT configuration
  // ============================================================
  const int nStations = 4;   // 4 tracking stations
  const int nLayers   = 3;   // 3 tube layers per station (honeycomb)
  const int nMagnets  = 3;   // 3 magnets between 4 stations

  // MDT tube parameters (ATLAS-style)
  G4double tubeOuterDiameter = 30.0 * mm;
  G4double tubeOuterRadius   = tubeOuterDiameter / 2.0;
  G4double tubeWallThickness = 0.4  * mm;   // 400 µm aluminium walls
  G4double tubeInnerRadius   = tubeOuterRadius - tubeWallThickness;
  G4double tubeLength        = 860. * mm;

  // Honeycomb packing: tubes touching within a layer (Y-pitch = diameter)
  //   and close-packed between adjacent layers (Z-pitch = sqrt(3)/2 × d)
  G4double tubeYPitch  = tubeOuterDiameter;                      // 30 mm
  G4double layerZPitch = tubeOuterDiameter * std::sqrt(3.) / 2.; // ≈ 25.98 mm
  G4double layerYShift = tubeOuterRadius;                        // 15 mm  (odd-layer Y shift)

  const int nRows = 34;  // tubes per layer → Y coverage ≈ ±510 mm

  // Thin aluminium endplates on both Z-faces of each station (the gray panels)
  G4double endplateThickness = 10.0 * mm;

  // Magnet / gap parameters
  G4double magnetThickness = 400. * mm;
  G4double gapBeforeMagnet = 15.  * mm;
  G4double gapAfterMagnet  = 15.  * mm;

  // Per-station Z size: 2 endplates + tube envelope
  //   tube envelope = (nLayers-1)*layerZPitch + tubeOuterDiameter
  G4double tubeEnvelopeZ = (nLayers - 1) * layerZPitch + tubeOuterDiameter;
  G4double stationZSize  = 2. * endplateThickness + tubeEnvelopeZ;

  // Total detector length
  G4double totalLength = nStations * stationZSize
                       + nMagnets  * (gapBeforeMagnet + magnetThickness + gapAfterMagnet);
  
  G4cout << "========================================" << G4endl;
  G4cout << "MDT Spectrometer — 3-layer Honeycomb" << G4endl;
  G4cout << "  " << nStations << " stations × " << nLayers << " layers × " << nRows << " tubes" << G4endl;
  G4cout << "  Tube diameter = " << tubeOuterDiameter/mm << " mm  Tube wall = " << tubeWallThickness/mm << " mm" << G4endl;
  G4cout << "  Y-pitch (within layer) = " << tubeYPitch/mm << " mm" << G4endl;
  G4cout << "  Z-pitch (layer-to-layer) = " << layerZPitch/mm << " mm" << G4endl;
  G4cout << "  Y-shift (odd layers) = " << layerYShift/mm << " mm" << G4endl;
  G4cout << "  Station Z size = " << stationZSize/mm << " mm  (incl. endplates)" << G4endl;
  G4cout << "  Total tubes = " << nStations*nLayers*nRows << G4endl;
  G4cout << "  Total length = " << totalLength/mm << " mm" << G4endl;
  G4cout << "========================================" << G4endl;

  // ============================================================
  // Create MDT tube geometry 
  //
  // G4Tubs axis is local Z.
  // Later we rotate tubes by 90 deg around Y:
  //
  // local Z --> global X
  //
  // Therefore MDT wires are along global X.
  // ============================================================
  // Outer aluminium tube
  // tubeOuterRadius = 15 mm, tubeLength = 860 mm
  // Inner gas volume (sensitive)
  // tubeOuterRadius= 15 mm, tubeInnerRadius = 14.6 mm, tubeLength = 860 mm
  
  G4Tubs* tubeOuterSolid = new G4Tubs("TubeOuter", 0, tubeOuterRadius, tubeLength/2, 0, 360*deg);
  G4Tubs* tubeInnerSolid = new G4Tubs("TubeInner", 0, tubeInnerRadius, tubeLength/2, 0, 360*deg);
  
  G4LogicalVolume* tubeOuterLV = new G4LogicalVolume(tubeOuterSolid, aluminum, "TubeOuterLV");
  //G4LogicalVolume* tubeInnerLV = new G4LogicalVolume(tubeInnerSolid, argon, "TubeInnerLV");
  fTubeInnerLV = new G4LogicalVolume(tubeInnerSolid, argon, "TubeInnerLV");
  
  // Place gas volume inside aluminum tube
  new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), fTubeInnerLV, "TubeInner", tubeOuterLV, false, 0, checkOverlaps);
  
  // Visualization
  auto tubeOuterVis = new G4VisAttributes(G4Colour(0.7, 0.7, 0.8));  // Light blue-gray
  tubeOuterVis->SetForceWireframe(true);
  tubeOuterLV->SetVisAttributes(tubeOuterVis);

  auto tubeInnerVis = new G4VisAttributes(G4Colour(0.0, 0.8, 1.0, 0.3));  // Cyan, transparent
  tubeInnerVis->SetForceWireframe(true);
  fTubeInnerLV->SetVisAttributes(tubeInnerVis);

  // Rotation matrix for horizontal tube placement
  // All tubes along X-direction: rotate 90° around Y-axis
  G4RotationMatrix* rotX = new G4RotationMatrix();
  rotX->rotateY(90.*deg);

  // ============================================================
  // Endplate geometry
  // ============================================================
  // Aluminium endplate geometry (reused per station) 
  // Y half-extent: covers nRows tubes (even layer) + margin
  G4double       epHalfY = (nRows - 1) / 2.0 * tubeYPitch + tubeOuterRadius + 8.*mm;
  G4Box*           epBox = new G4Box("Endplate", tubeLength/2., epHalfY, endplateThickness/2.);
  G4LogicalVolume* epLV  = new G4LogicalVolume(epBox, aluminum, "EndplateLV");
  
  auto epVis = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
  epVis->SetForceSolid(true);
  epLV->SetVisAttributes(epVis);

  // ============================================================
  // Build full detector
  //
  // Build 4 Stations (3-layer honeycomb)
  // ============================================================
  G4double zStart = -totalLength / 2.;
  G4double zPos   = zStart;

  for (int station = 0; station < nStations; ++station) {
    G4cout << "Building Station " << (station + 1) << " at z = " << zPos/mm << " mm" << G4endl;

    // ---- front endplate ----
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zPos + endplateThickness/2.),
                      epLV, "Endplate", worldLV, false, station * 10 + 0, checkOverlaps);
    zPos += endplateThickness;

    // ---- 3 MDT tube layers in honeycomb arrangement ----
    for (int layer = 0; layer < nLayers; ++layer) {
      // Z centre of this layer
      G4double layerZ = zPos + tubeOuterRadius + layer * layerZPitch;

      // Odd layers are shifted by half a tube diameter in Y → honeycomb offset
      G4double yOffset = (layer % 2 == 0) ? 0.0 : layerYShift;

      // Centre the array in Y around 0
      G4double yStart = -(nRows - 1) / 2.0 * tubeYPitch + yOffset;

      for (int row = 0; row < nRows; ++row) {
        G4double tubeY    = yStart + row * tubeYPitch;
        // stationID = 1,2,3,4
        // planeID   = 0,1,2
        // tubeID    = 0,...,33
        G4int copyNumber  = (station + 1) * 10000 + layer * 1000 + row;
        new G4PVPlacement(rotX, G4ThreeVector(0., tubeY, layerZ),
                          tubeOuterLV, "MDTTube", worldLV, false, copyNumber, checkOverlaps);
      }

      G4cout << "  Layer " << layer << " (" << (layer%2?"odd":"even") << "): "
             << nRows << " tubes at Z=" << layerZ/mm << " mm, Ystart="
             << (-(nRows-1)/2.0*tubeYPitch + yOffset)/mm << " mm" << G4endl;
    }
    // Advance past the tube envelope
    zPos += tubeEnvelopeZ;

    // ---- back endplate ----
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zPos + endplateThickness/2.),
                      epLV, "Endplate", worldLV, false, station * 10 + 1, checkOverlaps);
    zPos += endplateThickness;
    
    // After last station, don't place a magnet
    if (station == nStations - 1) break;
    
    // ============================================================
    // Magnet
    // ============================================================
    zPos += gapBeforeMagnet;
    
    // Place magnet with slits (same as BabayMIND/FASERCAL)
    G4Box* solidMagnet = new G4Box("Magnet", 500.*mm, 500.*mm, magnetThickness/2);
    G4Box* slitBox = new G4Box("Slit", 250.*mm, 10.*mm, magnetThickness/2);

    G4SubtractionSolid* magnetMinusSlit1 = new G4SubtractionSolid(
      "MagnetMinusSlit1", solidMagnet, slitBox, nullptr, G4ThreeVector(0, 250.*mm, 0));
    G4SubtractionSolid* magnetWithSlits = new G4SubtractionSolid(
      "MagnetWithSlits", magnetMinusSlit1, slitBox, nullptr, G4ThreeVector(0, -250.*mm, 0));

    G4LogicalVolume* magnetLV = new G4LogicalVolume(magnetWithSlits, steel, "MagnetLV");
    G4double magnetZ = zPos + magnetThickness/2;
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, magnetZ), magnetLV, "Magnet", worldLV, false, station, checkOverlaps);

    // Magnetic field on the magnet volume
    auto field = new MagneticField();
    field->SetSlitPosition(250.*mm);
    auto fieldManager = new G4FieldManager(field);
    fieldManager->CreateChordFinder(field);
    magnetLV->SetFieldManager(fieldManager, true);

    // Magnet visualization (gray)
    auto magnetVis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    magnetVis->SetForceSolid(true);
    magnetLV->SetVisAttributes(magnetVis);

    // ============================================================
    // Aluminum conductor strips
    // ============================================================
    int nTopStrips = 11;
    G4double stripWidth = 25.*mm;
    G4double stripHalfLength = 125.*mm;
    G4double stripThickness = 2.*mm;

    for (int j = 0; j < nTopStrips; ++j) {
      G4bool placeFront = (j % 2 == 0);
      G4double xOffset = -250.*mm + j * 50.*mm;
      G4double yOffset = 375.*mm;
      G4Box* aluStripSolid = new G4Box("AluStrip", stripWidth, stripHalfLength, stripThickness);
      G4LogicalVolume* aluStripLV = new G4LogicalVolume(aluStripSolid, aluminum, "AluStripLV");
      
      auto aluVis = new G4VisAttributes(G4Colour(0.8, 0.6, 0.4));
      aluVis->SetForceSolid(true);
      aluStripLV->SetVisAttributes(aluVis);
      
      G4double zPlacement = magnetZ + (placeFront ? -magnetThickness/2 - 2.*mm : magnetThickness/2 + 2.*mm);
      G4ThreeVector pos(xOffset, yOffset, zPlacement);
      new G4PVPlacement(nullptr, pos, aluStripLV, "AluStrip", worldLV, false, station * 1000 + j, checkOverlaps);
    }
    for (int j = 0; j < nTopStrips; ++j) {
      G4bool placeFront = (j % 2 == 0);
      G4double xOffset = -250.*mm + j * 50.*mm;
      G4double yOffset = -375.*mm;
      G4Box* aluStripSolid = new G4Box("AluStrip", stripWidth, stripHalfLength, stripThickness);
      G4LogicalVolume* aluStripLV = new G4LogicalVolume(aluStripSolid, aluminum, "AluStripLV");
      
      auto aluVis = new G4VisAttributes(G4Colour(0.8, 0.6, 0.4));
      aluVis->SetForceSolid(true);
      aluStripLV->SetVisAttributes(aluVis);
      
      G4double zPlacement = magnetZ + (placeFront ? -magnetThickness/2 - 2.*mm : magnetThickness/2 + 2.*mm);
      G4ThreeVector pos(xOffset, yOffset, zPlacement);
      new G4PVPlacement(nullptr, pos, aluStripLV, "AluStrip", worldLV, false, station * 1000 + j + 100, checkOverlaps);
    }
    int nMiddleStrips = 12;
    G4double MiddlestripHalfLength = 250.*mm;
    for (int j = 1; j < nMiddleStrips; ++j) {
      G4bool placeFront = (j % 2 == 0);
      G4double xOffset = -300.*mm + j * 50.*mm;
      G4double yOffset = 0;
      G4Box* aluStripSolid = new G4Box("AluStrip", stripWidth, MiddlestripHalfLength, stripThickness);
      G4LogicalVolume* aluStripLV = new G4LogicalVolume(aluStripSolid, aluminum, "AluStripLV");
      
      auto aluVis = new G4VisAttributes(G4Colour(0.8, 0.6, 0.4));
      aluVis->SetForceSolid(true);
      aluStripLV->SetVisAttributes(aluVis);
      
      G4double zPlacement = magnetZ + (placeFront ? -magnetThickness/2 - 2.*mm : magnetThickness/2 + 2.*mm);
      G4ThreeVector pos(xOffset, yOffset, zPlacement);
      new G4PVPlacement(nullptr, pos, aluStripLV, "AluStrip", worldLV, false, station * 1000 + j + 200, checkOverlaps);
    }

    // Advance past magnet
    zPos += magnetThickness + gapAfterMagnet;
    
    G4cout << "Magnet " << (station + 1) << " placed at z = " << magnetZ/mm << " mm" << G4endl;
  }

  // World (invisible)
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  // ============================================================
  // Export GDML
  // ============================================================
  G4GDMLParser parser;
  G4String geometryFilename = "detector_geometry.gdml";

  std::remove(geometryFilename.c_str());
  parser.Write(geometryFilename, worldPV);

  G4cout << "Geometry construction complete!" << G4endl;
  
  return worldPV;
}

