#include "CBDsimSteppingAction.hh"

#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4OpticalPhoton.hh"
#include "G4UnitsTable.hh"
#include <iomanip>



CBDsimSteppingAction::CBDsimSteppingAction(CBDsimEventAction* evtAct)
: G4UserSteppingAction(), fEventAction(evtAct) {
}

CBDsimSteppingAction::~CBDsimSteppingAction() {}
void CBDsimSteppingAction::UserSteppingAction(const G4Step* step)
{
  if ( step->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() ) return;
  // G4cout<<"where are you?"<<G4endl;
  G4Track* track = step->GetTrack();
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
  // std::cout<<"8888888888"<<std::endl;
  fEdep.towerNum = theTouchable->GetCopyNumber(theTouchable->GetHistoryDepth()-1);
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
  // G4cout << volume->GetName() << G4endl;

  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();
  auto deltaE = step->GetDeltaEnergy();

  const G4ParticleDefinition* particle;
  particle = step->GetTrack()->GetParticleDefinition();

  // Process name
  G4String proc_name = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();


  G4cout << "******************************" << G4endl;
  G4cout << "Step is limited by "
  << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
  << G4endl;
  G4cout << "Processes involved to the step" << G4endl;
  G4cout << "******************************" << G4endl;


  G4cout << "Particle Name is " << particle->GetParticleName() <<
  " and its kinetic energy is " << step->GetPreStepPoint()->GetKineticEnergy()/MeV <<
  " and the total energy is " << step->GetPreStepPoint()->GetTotalEnergy()/MeV <<
  " Energy deposit is " << edep <<
  " Track ID is " << step->GetTrack()->GetTrackID() <<
  " parent ID is " << step->GetTrack()->GetParentID() <<
  " lepton # " << particle->GetLeptonNumber() <<
  " Baryon # " << particle->GetBaryonNumber() <<
  " Type " << particle->GetParticleType() << G4endl;

  // G4cout << "----------------------------------------" << G4endl;

  G4double sum_all_ke = 0;
  G4double sum_all_e = 0;
  //G4int num_test=0;
  G4int oPnumber = 0;  // 이벤트마다 리셋

  if(nSecTotal>0)
  {
    G4cout << "  :----- List of 2ndaries - " << std::setw(3) << nSecTotal
    << " (Rest=" << std::setw(2) << nSecAtRest
    << ",Along=" << std::setw(2) << nSecAlong
    << ",Post="  << std::setw(2) << nSecPost << ")" << G4endl;

      for(size_t lp1=(*secVec).size()-nSecTotal; lp1<(*secVec).size(); lp1++)
      {
          // Scintillation can create O(10^3–10^5) optical photons per step; per-photon G4cout/fillPhysics
          // freezes the session. Count only; transport still runs (optical steps return at line 19).
          if ( (*secVec)[lp1]->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() ) {
            oPnumber += 1;
            continue;
          }
          G4double x= (*secVec)[lp1]->GetPosition().getX();
          G4double y=(*secVec)[lp1]->GetPosition().getY();
          G4double z=(*secVec)[lp1]->GetPosition().getZ();
          G4double energy=(*secVec)[lp1]->GetTotalEnergy();
          G4String PartName=(*secVec)[lp1]->GetDefinition()->GetParticleName();
          G4String PhysicName=(*secVec)[lp1]->GetCreatorProcess()->GetProcessName();
          
          // physical 객체 초기화
          physical = CBDsimInterface::CBDsimPhysicalevent();
          fEventAction->fillPhysics(physical,x,y,z,energy,PartName,PhysicName);
          num_test+=1;
          // G4cout<<"##########################"<<G4endl;
          // G4cout<<"this is debug line "<<G4endl;
          // G4cout<<"test number is "<<num_test<<G4endl;
          // G4cout<<"##########################"<<G4endl;
          G4cout << "    : "
          << G4BestUnit((*secVec)[lp1]->GetPosition(), "Length") << " "
          << std::setw( 9) << G4BestUnit((*secVec)[lp1]->GetKineticEnergy() , "Energy") << " "
          << std::setw( 9) << G4BestUnit((*secVec)[lp1]->GetTotalEnergy() , "Energy") << " "
          << std::setw(18) << (*secVec)[lp1]->GetDefinition()->GetParticleName()
          << " generated by " << (*secVec)[lp1]->GetCreatorProcess()->GetProcessName() << "  "
          << " Track ID " << (*secVec)[lp1]->GetTrackID() << G4endl;
          sum_all_ke += (*secVec)[lp1]->GetKineticEnergy();
          sum_all_e += (*secVec)[lp1]->GetTotalEnergy();

          if ( (*secVec)[lp1]->GetDefinition()->GetParticleType() == "nucleus" ){

              G4cout << (*secVec)[lp1]->GetDefinition()->GetParticleName() << "   "
              << (*secVec)[lp1]->GetDefinition()->GetAtomicMass() << "   "
              << (*secVec)[lp1]->GetDefinition()->GetAtomicNumber() << "   "
              << (*secVec)[lp1]->GetDefinition()->GetParticleSubType() << G4endl;
          }
      }
      // photon 객체 초기화
      photon = CBDsimInterface::CBDsimPhoton();
      fEventAction->fillOpticalPhoton(photon,oPnumber);
      G4cout<<"infinity loop?"<<G4endl;
    }

}
