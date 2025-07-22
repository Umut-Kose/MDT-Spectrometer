{
  TGeoManager::Import("detector_geometry.gdml");
  gGeoManager->GetTopVolume()->Draw("ogl");
}
