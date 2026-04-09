#include <iostream>

// 1 = prototype scintillator + tessellated light guide (CBDsimDetectorConstructionProto)
// 0 = legacy calorimeter tower geometry (CBDsimDetectorConstruction)
#ifndef CBDsim_USE_PROTO_GEOMETRY
#define CBDsim_USE_PROTO_GEOMETRY 1
#endif

#if CBDsim_USE_PROTO_GEOMETRY
#include "CBDsimDetectorConstructionProto.hh"
#else
#include "CBDsimDetectorConstruction.hh"
#endif
#include "CBDsimActionInitialization.hh"

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh"	
#include "FTFP_BERT.hh"
#include "FTFP_BERT_HP.hh"
#include "QGSP_BERT_HP.hh"
#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main(int argc, char** argv) {
  // Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) ui = new G4UIExecutive(argc, argv);


  G4int seed = 0;
  G4String filename;
  if (argc > 2) seed = std::atoi(argv[2]);
  if (argc > 3) filename = argv[3];

  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(seed);

  // Interactive mode uses Serial for stable UI/vis command handling.03.20
  // Batch mode keeps Default (typically MT/Tasking) for throughput.
  const auto runManagerType = (argc == 1)
    ? G4RunManagerType::Serial
    : G4RunManagerType::Default;
  auto runManager = G4RunManagerFactory::CreateRunManager(runManagerType);


  // Mandatory user initialization classes
#if CBDsim_USE_PROTO_GEOMETRY
  runManager->SetUserInitialization(new CBDsimDetectorConstructionProto());
#else
  runManager->SetUserInitialization(new CBDsimDetectorConstruction());
#endif

  // physics module
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  physicsList->RegisterPhysics(opticalPhysics);
  

  auto opt = G4OpticalParameters::Instance();
  // opt->SetCerenkovStackPhotons(true); 변경 점이 많아서 일단은 주석 으로 남겨 놓기?
  // opt->SetScintillationStackPhotons(true);
  opt->SetCerenkovTrackSecondariesFirst(true);
  opt->SetScintTrackSecondariesFirst(true);
  runManager->SetUserInitialization(physicsList);

  

  // User action initialization
  runManager->SetUserInitialization(new CBDsimActionInitialization(seed,filename));


  // Visualization manager construction
  G4VisManager* visManager = new G4VisExecutive(argc, argv);
  visManager->Initialize();
  G4UImanager* UImanager = G4UImanager::GetUIpointer();


  if ( argc != 1 ) {
    // execute an argument macro file if exist
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  } else {
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    UImanager->ApplyCommand("/control/execute init.mac");
    if (ui->IsGUI()) { UImanager->ApplyCommand("/control/execute gui.mac"); }
    // start interactive session
    ui->SessionStart();
    delete ui;
  }


  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;

  return 0;
}
