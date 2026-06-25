#include "PrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    fParticleGun(new G4ParticleGun(1)),
    fEnergyIndex(0)
{
  // Define energy steps from 10 GeV to 100 GeV
  for (G4double e = 10 * GeV; e <= 100 * GeV; e += 10 * GeV)
    fEnergies.push_back(e);
  // From 100 GeV to 1 TeV (200 GeV step)
  for (G4double e = 100 * GeV + 100 * GeV; e <= 1 * TeV; e += 200 * GeV)
    fEnergies.push_back(e);
//

  // Set particle type to mu-
  G4ParticleDefinition* muon = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
  fParticleGun->SetParticleDefinition(muon);
  
  // Start position at origin
  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -1500. * mm));
  //fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, -worldZ/2 - 10*mm));

  // Direction along +Z
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //if (fEnergyIndex >= fEnergies.size())
  //  return; // stop firing when out of energy steps
  /*
  the following code give duplicate eventID numbers!!
  if (fEnergyIndex >= fEnergies.size()) {
    fEnergyIndex = 0; // Reset to loop through energies again
  }
 
  G4double energy = fEnergies[fEnergyIndex];
  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->GeneratePrimaryVertex(anEvent);
  
  G4cout << "Fired muon with energy: " << energy / GeV << " GeV" << G4endl;
  */
// Use modulo to cycle through energies safely
  G4int currentIndex = fEnergyIndex % fEnergies.size();
  G4double energy = fEnergies[currentIndex];
  
  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->GeneratePrimaryVertex(anEvent);
  
  G4cout << "Event " << anEvent->GetEventID() 
         << ": Fired muon with energy: " << energy / GeV << " GeV" << G4endl;
  

  ++fEnergyIndex;
}
