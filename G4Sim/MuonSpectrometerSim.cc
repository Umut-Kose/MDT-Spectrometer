#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "QBBC.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "QGSP_BERT.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "ActionInitialization.hh"

int main(int argc, char** argv) {
    G4UIExecutive* ui = nullptr;
    if (argc == 1) {
        ui = new G4UIExecutive(argc, argv);
    }

    //auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    runManager->SetUserInitialization(new DetectorConstruction());
    //runManager->SetUserInitialization(new QBBC);
    // Geant4 natively loads its high-energy standard electromagnetic package when QGSP_BERT initialization is used.
    runManager->SetUserInitialization(new QGSP_BERT());

    //runManager->SetUserAction(new PrimaryGeneratorAction());
    runManager->SetUserInitialization(new ActionInitialization());

    runManager->Initialize();

    G4VisManager* visManager = new G4VisExecutive();
    visManager->Initialize();

    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (!ui) {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    } else {
        UImanager->ApplyCommand("/control/execute macro/run.mac");
        ui->SessionStart();
        delete ui;
    }

    delete visManager;
    delete runManager;
    return 0;
}
