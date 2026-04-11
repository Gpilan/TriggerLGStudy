#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"
#include "Randomize.hh"
#include "CBDsimPrimaryGeneratorAction.hh"

namespace { G4Mutex CBDsimPrimaryGeneratorActionMutex = G4MUTEX_INITIALIZER; }
int CBDsimPrimaryGeneratorAction::sNumEvt = 0;
G4ThreadLocal int CBDsimPrimaryGeneratorAction::sIdxEvt = 0;
G4ThreadLocal G4double CBDsimPrimaryGeneratorAction::sLastPrimaryEkin = 0.;
G4ThreadLocal G4double CBDsimPrimaryGeneratorAction::sLastPrimaryVx = 0.;
G4ThreadLocal G4double CBDsimPrimaryGeneratorAction::sLastPrimaryVy = 0.;
G4ThreadLocal G4double CBDsimPrimaryGeneratorAction::sLastPrimaryVz = 0.;
G4ThreadLocal G4double CBDsimPrimaryGeneratorAction::sLastPrimaryDirX = 0.;
G4ThreadLocal G4double CBDsimPrimaryGeneratorAction::sLastPrimaryDirY = 0.;
G4ThreadLocal G4double CBDsimPrimaryGeneratorAction::sLastPrimaryDirZ = 1.;

void CBDsimPrimaryGeneratorAction::GetLastPrimaryKinematics(G4double& ekinMeV, G4double& vx, G4double& vy,
                                                            G4double& vz, G4double& dx, G4double& dy,
                                                            G4double& dz) {
  ekinMeV = sLastPrimaryEkin;
  vx = sLastPrimaryVx;
  vy = sLastPrimaryVy;
  vz = sLastPrimaryVz;
  dx = sLastPrimaryDirX;
  dy = sLastPrimaryDirY;
  dz = sLastPrimaryDirZ;
}

CBDsimPrimaryGeneratorAction::CBDsimPrimaryGeneratorAction(G4int seed)
: G4VUserPrimaryGeneratorAction() {
  fSeed = seed;
  fNumPtc = 1;
  fTheta = 0.;
  fPhi = 0.;
  fRandX = 0.*mm;
  fRandY = 0.*mm;
  fY_0 = 0.*cm;
  fZ_0 = 0.*cm;

  fParticleGun = new G4ParticleGun(fNumPtc);

  G4ParticleTable* ptcTable = G4ParticleTable::GetParticleTable();
  G4String ptcName;
  fElectron = ptcTable->FindParticle(ptcName="e-");
  fMuon = ptcTable->FindParticle(ptcName="mu-");
  fPion = ptcTable->FindParticle(ptcName="pi+");
  fKaon0L = ptcTable->FindParticle(ptcName="kaon0L");
  fProton = ptcTable->FindParticle(ptcName="proton");
  fOptPhoton = ptcTable->FindParticle(ptcName="opticalphoton");
  fGamma = ptcTable->FindParticle(ptcName="gamma");

  DefineCommands();
}

CBDsimPrimaryGeneratorAction::~CBDsimPrimaryGeneratorAction() {
  if (fMessenger) delete fMessenger;
}

void CBDsimPrimaryGeneratorAction::GeneratePrimaries(G4Event* evt) {
  fY = (G4UniformRand()-0.5)*fRandX + fY_0;
  fZ = (G4UniformRand()-0.5)*fRandY + fZ_0;
  fOrigin.set(0.,fY,fZ);

  fParticleGun->SetParticlePosition(fOrigin);

  // Base +z (Proto trigger thin axis, typical /gun/direction 0 0 1). For legacy +x use
  // /CBDsim/generator/theta 1.570796327 rad (rotateY: +z -> +x).
  fDirection.set(0., 0., 1.);
  fDirection.rotateY(fTheta);
  fDirection.rotateZ(fPhi);

  fParticleGun->SetParticleMomentumDirection(fDirection);

  G4AutoLock lock(&CBDsimPrimaryGeneratorActionMutex);
  fParticleGun->GeneratePrimaryVertex(evt);
  sLastPrimaryEkin = fParticleGun->GetParticleEnergy();
  sLastPrimaryVx = fOrigin.x() / mm;
  sLastPrimaryVy = fOrigin.y() / mm;
  sLastPrimaryVz = fOrigin.z() / mm;
  {
    G4ThreeVector d = fDirection;
    if (d.mag2() > 0.) d = d.unit();
    sLastPrimaryDirX = d.x();
    sLastPrimaryDirY = d.y();
    sLastPrimaryDirZ = d.z();
  }
  sIdxEvt = sNumEvt;
  sNumEvt++;
}

void CBDsimPrimaryGeneratorAction::DefineCommands() {
  // Define /CBDsim/generator command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/CBDsim/generator/", "Primary generator control");

  G4GenericMessenger::Command& etaCmd = fMessenger->DeclareMethodWithUnit("theta","rad",&CBDsimPrimaryGeneratorAction::SetTheta,"theta of beam");
  etaCmd.SetParameterName("theta",true);
  etaCmd.SetDefaultValue("0.");

  G4GenericMessenger::Command& phiCmd = fMessenger->DeclareMethodWithUnit("phi","rad",&CBDsimPrimaryGeneratorAction::SetPhi,"phi of beam");
  phiCmd.SetParameterName("phi",true);
  phiCmd.SetDefaultValue("0.");

  G4GenericMessenger::Command& y0Cmd = fMessenger->DeclareMethodWithUnit("y0","cm",&CBDsimPrimaryGeneratorAction::SetY0,"y_0 of beam");
  y0Cmd.SetParameterName("y0",true);
  y0Cmd.SetDefaultValue("0.");

  G4GenericMessenger::Command& z0Cmd = fMessenger->DeclareMethodWithUnit("z0","cm",&CBDsimPrimaryGeneratorAction::SetZ0,"z_0 of beam");
  z0Cmd.SetParameterName("z0",true);
  z0Cmd.SetDefaultValue("0.");

  G4GenericMessenger::Command& randxCmd = fMessenger->DeclareMethodWithUnit("randx","mm",&CBDsimPrimaryGeneratorAction::SetRandX,"x width of beam");
  randxCmd.SetParameterName("randx",true);
  randxCmd.SetDefaultValue("10.");

  G4GenericMessenger::Command& randyCmd = fMessenger->DeclareMethodWithUnit("randy","mm",&CBDsimPrimaryGeneratorAction::SetRandY,"y width of beam");
  randyCmd.SetParameterName("randy",true);
  randyCmd.SetDefaultValue("10.");
}
