// Exercises the production solids, materials, surfaces and SD; no substitute geometry.
#include "CBDsimDetectorConstructionProto.hh"
#include "CBDsimOpticalDiagnostics.hh"
#include "G4RunManager.hh"
#include "G4Navigator.hh"
#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4TessellatedSolid.hh"
#include "G4PhysicalConstants.hh"
#include <iomanip>
#include "G4VUserPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpAbsorption.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <cmath>

namespace {
bool noLG = false;
double wavelength=450., angle=0.; bool inside=false; long waferSteps=0, relocation=0, relocationKilled=0;
long firstBad=0, firstGood=0, loops=0;
G4OpBoundaryProcess* boundary=nullptr;
class Physics : public G4VUserPhysicsList {
  void ConstructParticle() override { G4OpticalPhoton::Definition(); }
  void ConstructProcess() override {
    AddTransportation(); auto* pm=G4OpticalPhoton::Definition()->GetProcessManager();
    pm->AddDiscreteProcess(new G4OpAbsorption);
    boundary=new G4OpBoundaryProcess;pm->AddDiscreteProcess(boundary);
  }
  void SetCuts() override {}
};
class Gun : public G4VUserPrimaryGeneratorAction {
  G4ParticleGun gun{1};
  void GeneratePrimaries(G4Event* event) override {
    G4VPhysicalVolume* window=nullptr;
    for(auto* pv:*G4PhysicalVolumeStore::GetInstance())
      if(pv->GetName()=="protoSipmWindowPhys" && pv->GetCopyNo()==0) window=pv;
    const auto rot=window->GetObjectRotationValue();
    // Window local z points along package axis; derive front-to-back using rotation.
    // In production local +z is world +y, so incoming photons travel local -z.
    G4ThreeVector local(0,0, inside ? -0.145*mm-0.001*mm : 0.145*mm-0.001*mm);
    gun.SetParticleDefinition(G4OpticalPhoton::Definition());
    gun.SetParticleEnergy(h_Planck*c_light/(wavelength*nm));
    gun.SetParticlePosition(rot*local+window->GetObjectTranslation());
    gun.SetParticleMomentumDirection(rot*G4ThreeVector(std::sin(angle),0,-std::cos(angle)));
    gun.SetParticlePolarization(rot*G4ThreeVector(0,1,0));gun.GeneratePrimaryVertex(event);
  }
};
class Events : public G4UserEventAction {
  void BeginOfEventAction(const G4Event* e) override {CBDsimOpticalDiagnostics::BeginEvent(e->GetEventID());}
  void EndOfEventAction(const G4Event* e) override {CBDsimOpticalDiagnostics::EndEvent(e->IsAborted());}
};
class Tracks : public G4UserTrackingAction {
  void PreUserTrackingAction(const G4Track* t) override {CBDsimOpticalDiagnostics::PreUserTrackingAction(t);}
  void PostUserTrackingAction(const G4Track* t) override {CBDsimOpticalDiagnostics::PostUserTrackingAction(t);}
};
class Steps : public G4UserSteppingAction {
  void UserSteppingAction(const G4Step* s) override {
    CBDsimOpticalDiagnostics::UserSteppingAction(s);
    auto* pv=s->GetPreStepPoint()->GetPhysicalVolume();
    if(pv && pv->GetName()=="protoSipmWaferPhys") {
      if(s->GetPostStepPoint()->GetStepStatus()==fGeomBoundary && boundary->GetStatus()==StepTooSmall) {
        ++relocation;if(s->GetTrack()->GetTrackStatus()==fStopAndKill)++relocationKilled;
      } else ++waferSteps;
    }
    if(G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()<10) {
      auto* post=s->GetPostStepPoint()->GetPhysicalVolume();
      std::cout<<"TRACE event="<<G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()
       <<" pre="<<(pv?pv->GetName():"none")<<" post="<<(post?post->GetName():"none")
       <<" length="<<s->GetStepLength()/mm<<" status="<<boundary->GetStatus()<<" killed="<<(s->GetTrack()->GetTrackStatus()==fStopAndKill)<<std::endl;
    }
    if(s->GetTrack()->GetCurrentStepNumber()>10000) {++loops;s->GetTrack()->SetTrackStatus(fStopAndKill);}
  }
};
}
int main(int argc,char** argv) {
  if(argc!=8) return 2;
  const int n=std::stoi(argv[1]); const long seed=std::stol(argv[2]);
  wavelength=std::stod(argv[4]);angle=std::stod(argv[5])*deg;inside=std::stoi(argv[6]);noLG=std::stoi(argv[7]);
  setenv("CBDsim_PROTO_NO_LG",noLG?"1":"0",1);
  setenv("CBDsim_OPTICAL_DIAG","1",1);setenv("CBDsim_OPTICAL_DIAG_HIST","0",1);
  CLHEP::HepRandom::setTheSeed(seed);
  auto* rm=new G4RunManager;rm->SetUserInitialization(new CBDsimDetectorConstructionProto);
  rm->SetUserInitialization(new Physics);rm->SetUserAction(new Gun);rm->SetUserAction(new Events);
  rm->SetUserAction(new Tracks);rm->SetUserAction(new Steps);rm->Initialize();
  CBDsimOpticalDiagnostics::ResetForRun();rm->BeamOn(n);CBDsimOpticalDiagnostics::MergeWorkerIntoMaster();
  CBDsimOpticalDiagnostics::WriteSummaryFile(argv[3],n);
  std::ifstream in(argv[3]);std::map<std::string,double> c;std::string k;double v;while(in>>k>>v)c[k]=v;
  const auto attempts=c["budget_fate_detected"]+c["budget_fate_qe_reject"];
  bool ok=relocationKilled==0 && loops==0 && c["budget_closed"]==1 && c["budget_fate_unknown"]==0 &&
    c["budget_fate_no_rindex"]==0 && waferSteps==attempts && attempts<=n;
  std::cout<<"SENSOR attempts="<<attempts<<" wafer_steps="<<waferSteps<<" detected="<<c["budget_fate_detected"]<<" relocation="<<relocation<<" relocation_killed="<<relocationKilled<<" pass="<<ok<<std::endl;
  delete rm;return ok?0:1;
}
