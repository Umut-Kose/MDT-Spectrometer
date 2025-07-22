#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
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

//SciFi detector smearing
#include "SciFiSD.hh"
#include "G4SDManager.hh"

DetectorConstruction::DetectorConstruction() {}
DetectorConstruction::~DetectorConstruction() {}

// Geometry construction
// 4layers_SciFi + 10modules( 150mm_Magnet + 40mm_AirGap)


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  auto nist = G4NistManager::Instance();
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* steel = nist->FindOrBuildMaterial("G4_Fe");
  G4Material* plastic = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
  G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");

  G4double worldZ = 2000. * mm;
  G4Box* worldBox = new G4Box("World", 1.3*m, 1.3*m, worldZ/2);
  G4LogicalVolume* worldLV = new G4LogicalVolume(worldBox, air, "WorldLV");
  G4VPhysicalVolume* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0);

  // Parameters
  const int nMagnets = 10;
  const int nSciFiGroups = nMagnets + 1;
  const int nSciFiPerGroup = 4;
  const int nSciFiPlanes = nSciFiGroups * nSciFiPerGroup;
  G4double scifilayerThickness = 2.5 * mm;
  G4double magnetThickness = 150. * mm; // Thickness of each magnet changed from 100 to 150
  G4double gapBeforeMagnet = 10. * mm;
  G4double gapAfterMagnet = 10. * mm;

  // Calculate total length needed
  G4double totalSciFiThickness = nSciFiPlanes * scifilayerThickness;
  G4double totalMagnetThickness = nMagnets * magnetThickness;
  G4double totalGapBefore = nMagnets * gapBeforeMagnet;
  G4double totalGapAfter = nMagnets * gapAfterMagnet;
  G4double totalLength = totalSciFiThickness + totalMagnetThickness + totalGapBefore + totalGapAfter;

  // Center everything in the world volume
  G4double zStart = -totalLength/2;

  // Register SciFi as sensitive
  G4SDManager* sdManager = G4SDManager::GetSDMpointer();
  static SciFiSD* scifiSD = new SciFiSD("SciFiSD");
  sdManager->AddNewDetector(scifiSD);

  int scifiLayerID = 0;
  G4double zPos = zStart + scifilayerThickness/2;

  for (int group = 0; group < nSciFiGroups; ++group) {
    // Place 4 SciFi layers
    for (int l = 0; l < nSciFiPerGroup; ++l) {
      G4Box* scifiBox = new G4Box("SciFiLayer", 500.*mm, 500.*mm, scifilayerThickness / 2);
      G4LogicalVolume* scifiLV = new G4LogicalVolume(scifiBox, plastic, "SciFiLayerLV");
      G4String physName = "SciFiLayer_" + std::to_string(scifiLayerID);
      new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zPos), scifiLV, physName, worldLV, false, scifiLayerID);

      auto scifiVis = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // green
      scifiVis->SetForceSolid(true);
      scifiLV->SetVisAttributes(scifiVis);
      scifiLV->SetSensitiveDetector(scifiSD);

      zPos += scifilayerThickness;
      scifiLayerID++;
    }

    // After last group, don't place a magnet
    if (group == nMagnets) break;

    // Gap before magnet
    zPos += gapBeforeMagnet;

    // Place magnet
    G4Box* solidMagnet = new G4Box("Magnet", 500.*mm, 500.*mm, magnetThickness/2);
    G4Box* slitBox = new G4Box("Slit", 250.*mm, 10.*mm, magnetThickness/2);

    G4SubtractionSolid* magnetMinusSlit1 = new G4SubtractionSolid(
      "MagnetMinusSlit1", solidMagnet, slitBox, nullptr, G4ThreeVector(0, 250.*mm, 0));
    G4SubtractionSolid* magnetWithSlits = new G4SubtractionSolid(
      "MagnetWithSlits", magnetMinusSlit1, slitBox, nullptr, G4ThreeVector(0, -250.*mm, 0));

    G4LogicalVolume* magnetLV = new G4LogicalVolume(magnetWithSlits, steel, "MagnetLV");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zPos + magnetThickness/2 - scifilayerThickness/2), magnetLV, "Magnet", worldLV, false, group);

    // Magnetic field on the magnet volume
    auto field = new MagneticField();
    auto fieldManager = new G4FieldManager(field);
    fieldManager->CreateChordFinder(field);
    magnetLV->SetFieldManager(fieldManager, true);

    // Magnet (gray)
    auto magnetVis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)); // RGB: gray
    magnetVis->SetForceSolid(true);
    magnetLV->SetVisAttributes(magnetVis);

    // Place aluminum strips as before, using the magnet center z
    G4double magnetZ = zPos + magnetThickness/2 - scifilayerThickness/2;

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
      G4ThreeVector pos(xOffset, yOffset , zPlacement);
      new G4PVPlacement(nullptr, pos, aluStripLV, "AluStrip", worldLV, false, group * 1000 + j);
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
      G4ThreeVector pos(xOffset, yOffset , zPlacement);
      new G4PVPlacement(nullptr, pos, aluStripLV, "AluStrip", worldLV, false, group * 1000 + j + 100);
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
      G4ThreeVector pos(xOffset, yOffset , zPlacement);
      new G4PVPlacement(nullptr, pos, aluStripLV, "AluStrip", worldLV, false, group * 1000 + j + 200);
    }

    // Advance zPos past the magnet
    zPos += magnetThickness;

    // Gap after magnet
    zPos += gapAfterMagnet;
  }

  // World (invisible)
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  // Export GDML file
  G4GDMLParser parser;
  parser.Write("detector_geometry.gdml", worldPV);

  return worldPV;
}