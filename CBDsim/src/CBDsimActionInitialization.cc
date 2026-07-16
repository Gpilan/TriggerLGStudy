#include "CBDsimActionInitialization.hh"
#include "CBDsimRunAction.hh"
#include "CBDsimPrimaryGeneratorAction.hh"
#include "CBDsimEventAction.hh"
#include "CBDsimSteppingAction.hh"
#include "CBDsimTrackingAction.hh"

CBDsimActionInitialization::CBDsimActionInitialization(G4int seed, G4String filename)
: G4VUserActionInitialization()
{
  fSeed = seed;
  fFilename = filename;
}

CBDsimActionInitialization::~CBDsimActionInitialization() {}

void CBDsimActionInitialization::BuildForMaster() const {
  SetUserAction(new CBDsimRunAction(fSeed,fFilename));
}

void CBDsimActionInitialization::Build() const {
  SetUserAction(new CBDsimPrimaryGeneratorAction(fSeed));
  SetUserAction(new CBDsimRunAction(fSeed,fFilename));

  CBDsimEventAction* evtAct = new CBDsimEventAction();
  SetUserAction(evtAct);

  SetUserAction(new CBDsimSteppingAction(evtAct));
  SetUserAction(new CBDsimTrackingAction());
}
