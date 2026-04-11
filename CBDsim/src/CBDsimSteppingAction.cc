#include "CBDsimSteppingAction.hh"

#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4OpticalPhoton.hh"
#include "G4UnitsTable.hh"
#include <iomanip>
#include <string>

namespace {
constexpr G4bool kVerboseStepping = false;

/** Proto: copy number on shared LVs (0=T1, 1=T2). Legacy: former GetCopyNumber(depth-1). */
G4int TriggerIndexFromTouchable(const G4TouchableHandle& th) {
  const G4int depth = th->GetHistoryDepth();
  for (G4int d = 0; d <= depth; ++d) {
    auto* vol = th->GetVolume(d);
    if (!vol) continue;
    const G4String& name = vol->GetName();
    if (name == "protoScintPhys" || name == "protoLightGuidePhys" || name == "protoSipmEnvPhys") {
      return vol->GetCopyNo();
    }
    if (name.find("protoFoil") != G4String::npos) return vol->GetCopyNo();
  }
  if (depth >= 1) return th->GetCopyNumber(depth - 1);
  return 0;
}
}  // namespace

CBDsimSteppingAction::CBDsimSteppingAction(CBDsimEventAction* evtAct)
: G4UserSteppingAction(), fEventAction(evtAct) {
}

CBDsimSteppingAction::~CBDsimSteppingAction() {}
void CBDsimSteppingAction::UserSteppingAction(const G4Step* step)
{
  if ( step->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() ) return;
  // G4cout<<"where are you?"<<G4endl;
  G4Track* track = step->GetTrack();
  (void)track;
  // std::cout<<"22222222"<<std::endl;
  
  // 변수 초기화
  num_test = 0;

  G4StepPoint* presteppoint = step->GetPreStepPoint();
  // std::cout<<"33333333"<<std::endl;
  G4LogicalVolume* preVol = presteppoint->GetPhysicalVolume()->GetLogicalVolume();
  // std::cout<<"444444444"<<std::endl;
  G4TouchableHandle theTouchable = presteppoint->GetTouchableHandle();
  // std::cout<<"55555555"<<std::endl;

  G4String matName = preVol->GetMaterial()->GetName();
  // std::cout<<"66666666"<<std::endl;

  if ( matName=="G4_Galactic" || matName=="G4_AIR" ) return;
  // std::cout<<"77777777"<<std::endl;

  fEdep.Edep = step->GetTotalEnergyDeposit();
  fEdep.triggerNum = TriggerIndexFromTouchable(theTouchable);
  //std::cout<<"Edep is "<<fEdep.Edep<<std::endl;

  // std::cout<<"99999999999"<<std::endl;
  if ( fEdep.Edep > 0. ) fEventAction->fillEdeps(fEdep);
  //////////////////
  // std::cout<<"10101010101010"<<std::endl;
  G4StepStatus stepStatus = fpSteppingManager->GetfStepStatus();
  // std::cout<<"12121212121212"<<std::endl;

  //G4cout << "The step status is " << stepStatus << G4endl;

  // std::cout<<"13131313131313"<<std::endl;
  G4int nSecAtRest = fpSteppingManager->GetfN2ndariesAtRestDoIt();
  // std::cout<<"14141414141414"<<std::endl;
  G4int nSecAlong  = fpSteppingManager->GetfN2ndariesAlongStepDoIt();
  G4int nSecPost   = fpSteppingManager->GetfN2ndariesPostStepDoIt();
  G4int nSecTotal  = nSecAtRest+nSecAlong+nSecPost;
  G4TrackVector* secVec = fpSteppingManager->GetfSecondary();

  // get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  (void)volume;

  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();
  auto deltaE = step->GetDeltaEnergy();

  const G4ParticleDefinition* particle;
  particle = step->GetTrack()->GetParticleDefinition();

  // Process name
  G4String proc_name = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  (void)proc_name;
  (void)edep;
  (void)deltaE;
  (void)stepStatus;

  G4double sum_all_ke = 0;
  G4double sum_all_e = 0;
  G4int oPnumber = 0;

  if (nSecTotal > 0) {
    if (kVerboseStepping) {
      G4cout << "******************************" << G4endl;
      G4cout << "Step is limited by " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
      G4cout << "Particle Name is " << particle->GetParticleName()
             << " KE = " << step->GetPreStepPoint()->GetKineticEnergy() / MeV << " MeV" << G4endl;
    }

    for (size_t lp1 = (*secVec).size() - nSecTotal; lp1 < (*secVec).size(); lp1++) {
      if ((*secVec)[lp1]->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
        oPnumber += 1;
        continue;
      }
      G4double x = (*secVec)[lp1]->GetPosition().getX();
      G4double y = (*secVec)[lp1]->GetPosition().getY();
      G4double z = (*secVec)[lp1]->GetPosition().getZ();
      G4double energy = (*secVec)[lp1]->GetTotalEnergy();
      G4String PartName = (*secVec)[lp1]->GetDefinition()->GetParticleName();
      G4String PhysicName = (*secVec)[lp1]->GetCreatorProcess()->GetProcessName();

      physical = CBDsimInterface::CBDsimPhysicalevent();
      fEventAction->fillPhysics(physical, x, y, z, energy, PartName, PhysicName);
      num_test += 1;

      if (kVerboseStepping) {
        G4cout << "    : " << G4BestUnit((*secVec)[lp1]->GetPosition(), "Length") << " "
               << G4BestUnit((*secVec)[lp1]->GetKineticEnergy(), "Energy") << " "
               << (*secVec)[lp1]->GetDefinition()->GetParticleName() << G4endl;
      }
      sum_all_ke += (*secVec)[lp1]->GetKineticEnergy();
      sum_all_e += (*secVec)[lp1]->GetTotalEnergy();

      if (kVerboseStepping && (*secVec)[lp1]->GetDefinition()->GetParticleType() == "nucleus") {
        G4cout << (*secVec)[lp1]->GetDefinition()->GetParticleName() << G4endl;
      }
    }
    photon = CBDsimInterface::CBDsimPhoton();
    fEventAction->fillOpticalPhoton(photon, oPnumber);
  }
  (void)sum_all_ke;
  (void)sum_all_e;
}
