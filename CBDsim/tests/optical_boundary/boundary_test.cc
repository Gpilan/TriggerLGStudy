// Real Geant4 transport through one planar interface, exercising production diagnostics.
#include "CBDsimOpticalDiagnostics.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpAbsorption.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleGun.hh"
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

namespace {
std::string mode;
G4OpBoundaryProcess* boundary = nullptr;
class Detector : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override {
    auto makeMaterial = [](const char* name, bool absorb) {
      auto* m = new G4Material(name, 1., 1.01*g/mole, 1.*g/cm3);
      auto* p = new G4MaterialPropertiesTable;
      G4double e[] = {2.*eV, 4.*eV}, n[] = {1.5, 1.5}, a[] = {1.e-3*mm, 1.e-3*mm};
      p->AddProperty("RINDEX", e, n, 2);
      if (absorb) p->AddProperty("ABSLENGTH", e, a, 2);
      m->SetMaterialPropertiesTable(p);
      return m;
    };
    auto* w = new G4LogicalVolume(new G4Box("world",100*mm,100*mm,100*mm),
                                 makeMaterial("outside",false),"world");
    auto* wp = new G4PVPlacement(nullptr,{},w,"protoWorldPhys",nullptr,false,0);
    auto* slab = new G4LogicalVolume(new G4Box("slab",10*mm,10*mm,5*mm),
                                     makeMaterial("inside",mode=="bulk"),"slab");
    auto* sp = new G4PVPlacement(nullptr,{0,0,-5*mm},slab,
        mode=="other_absorb" ? "protoScintPhys" : "protoLightGuidePhys",w,false,0,true);
    auto* surface = new G4OpticalSurface("test",glisur,polished,
        mode=="transmit" || mode=="bulk" ? dielectric_dielectric : dielectric_metal);
    auto* p = new G4MaterialPropertiesTable;
    G4double e[] = {2.*eV,4.*eV};
    double r = mode=="reflect" ? 1. : mode=="mixed" ? 0.5 : 0.;
    G4double refl[] = {r,r}, eff[] = {mode=="detect" ? 1. : 0.,mode=="detect" ? 1. : 0.};
    if (mode!="transmit" && mode!="bulk") {
      p->AddProperty("REFLECTIVITY",e,refl,2);
      p->AddProperty("EFFICIENCY",e,eff,2);
    }
    surface->SetMaterialPropertiesTable(p);
    new G4LogicalBorderSurface("interface",sp,wp,surface);
    return wp;
  }
};
class Physics : public G4VUserPhysicsList {
  void ConstructParticle() override { G4OpticalPhoton::OpticalPhotonDefinition(); }
  void ConstructProcess() override {
    AddTransportation();
    auto* pm=G4OpticalPhoton::OpticalPhotonDefinition()->GetProcessManager();
    pm->AddDiscreteProcess(new G4OpAbsorption);
    boundary=new G4OpBoundaryProcess;
    pm->AddDiscreteProcess(boundary);
  }
  void SetCuts() override {}
};
class Gun : public G4VUserPrimaryGeneratorAction {
  G4ParticleGun gun{1};
  void GeneratePrimaries(G4Event* event) override {
    gun.SetParticleDefinition(G4OpticalPhoton::OpticalPhotonDefinition());
    gun.SetParticleEnergy(3.*eV);
    gun.SetParticlePosition({0,0,-5*mm});
    gun.SetParticleMomentumDirection({0,0,1});
    gun.SetParticlePolarization({1,0,0});
    gun.GeneratePrimaryVertex(event);
  }
};
class Steps : public G4UserSteppingAction {
 public:
  std::map<int,long> observed;
  long nonBoundary=0;
  void UserSteppingAction(const G4Step* s) override {
    CBDsimOpticalDiagnostics::UserSteppingAction(s);
    if (s->GetPostStepPoint()->GetStepStatus()==fGeomBoundary) {
      ++observed[static_cast<int>(boundary->GetStatus())];
      // This test ends after ONE encounter; these artificial terminations are not a fate-budget test.
      s->GetTrack()->SetTrackStatus(fStopAndKill);
    } else ++nonBoundary;
  }
};
class Tracks : public G4UserTrackingAction {
  void PostUserTrackingAction(const G4Track* t) override {
    CBDsimOpticalDiagnostics::PostUserTrackingAction(t);
  }
};
}
int main(int argc,char** argv) {
  if(argc!=5) return 2;
  mode=argv[1]; const int n=std::stoi(argv[2]); const long seed=std::stol(argv[3]);
  if(mode!="reflect" && mode!="absorb" && mode!="other_absorb" && mode!="transmit" && mode!="mixed" && mode!="bulk" && mode!="detect") return 2;
  setenv("CBDsim_OPTICAL_DIAG","1",1); setenv("CBDsim_OPTICAL_DIAG_HIST","1",1);
  CLHEP::HepRandom::setTheSeed(seed);
  auto* rm=new G4RunManager;
  rm->SetUserInitialization(new Detector); rm->SetUserInitialization(new Physics);
  rm->SetUserAction(new Gun); auto* steps=new Steps; rm->SetUserAction(steps);
  rm->SetUserAction(new Tracks); rm->Initialize();
  CBDsimOpticalDiagnostics::ResetForRun(); rm->BeamOn(n);
  CBDsimOpticalDiagnostics::MergeWorkerIntoMaster();
  CBDsimOpticalDiagnostics::WriteSummaryFile(argv[4],n);
  CBDsimOpticalDiagnostics::WritePathHistogramFile(std::string(argv[4])+".hist");
  std::ifstream file(argv[4]); std::map<std::string,double> c; std::string k; double v;
  while(file>>k>>v) c[k]=v;
  bool ok=true;
  for(const auto& p: steps->observed)
    ok &= c["boundary_status_"+std::to_string(p.first)]==p.second;
  const double absorption=c["boundary_absorb_lg"]+c["boundary_absorb_other"];
  const double bulk=c["bulk_absorb_lg"];
  const long reflect=steps->observed[SpikeReflection];
  if(mode=="absorb" || mode=="other_absorb") ok &= absorption==n && bulk==0;
  if(mode=="other_absorb") ok &= c["boundary_absorb_other"]==n;
  if(mode=="reflect") ok &= reflect==n && absorption==0 && bulk==0;
  if(mode=="transmit") ok &= steps->observed[FresnelRefraction]==n && absorption==0 && bulk==0;
  if(mode=="detect") ok &= steps->observed[Detection]==n && absorption==0 && bulk==0;
  if(mode=="bulk") ok &= bulk==n && absorption==0 && c["boundary_total"]==0;
  if(mode=="mixed") ok &= absorption+reflect==n && std::abs(absorption-0.5*n)<=5.*std::sqrt(0.25*n) && bulk==0;
  std::cout<<"VALIDATION "<<mode<<" seed="<<seed<<" n="<<n<<" absorb="<<absorption
           <<" reflect="<<reflect<<" bulk="<<bulk<<" pass="<<ok<<std::endl;
  delete rm; return ok ? 0 : 1;
}
