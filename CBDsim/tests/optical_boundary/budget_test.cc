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
#include "G4OpWLS.hh"
#include "G4Event.hh"
#include "G4UserEventAction.hh"
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
      G4double e[] = {0.5*eV, 4.*eV}, n[] = {1.5, 1.5}, a[] = {1.e-3*mm, 1.e-3*mm};
      if (!(mode=="no_rindex" && std::string(name)=="outside")) p->AddProperty("RINDEX", e, n, 2);
      if (absorb) p->AddProperty("ABSLENGTH", e, a, 2);
      if (mode=="wls" && std::string(name)=="inside") {
        G4double we[]={1.*eV,2.*eV,2.5*eV,4.*eV}, wa[]={1.e9*mm,1.e9*mm,1.e-3*mm,1.e-3*mm};
        G4double se[]={1.*eV,1.5*eV,2.*eV}, intensity[]={0.,1.,0.};
        p->AddProperty("WLSABSLENGTH",we,wa,4);
        p->AddProperty("WLSCOMPONENT",se,intensity,3);
        p->AddConstProperty("WLSTIMECONSTANT",1.*ns);
      }
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
        mode=="no_rindex" || mode=="transmit" || mode=="bulk" || mode=="wls" || mode=="qe" || mode=="sd_detect" || mode=="unknown" ? dielectric_dielectric : dielectric_metal);
    auto* p = new G4MaterialPropertiesTable;
    G4double e[] = {2.*eV,4.*eV};
    double r = mode=="reflect" ? 1. : mode=="mixed" ? 0.5 : 0.;
    G4double refl[] = {r,r}, eff[] = {mode=="detect" ? 1. : 0.,mode=="detect" ? 1. : 0.};
    if (mode!="no_rindex" && mode!="transmit" && mode!="bulk" && mode!="wls" && mode!="qe" && mode!="sd_detect" && mode!="unknown") {
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
    if (mode=="wls") pm->AddDiscreteProcess(new G4OpWLS);
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
    if (mode=="qe" || mode=="sd_detect") {
      if (mode=="qe") CBDsimOpticalDiagnostics::RecordSipmQeReject(s);
      else CBDsimOpticalDiagnostics::RecordSipmDetect(s);
      s->GetTrack()->SetTrackStatus(fStopAndKill);
    }
    CBDsimOpticalDiagnostics::UserSteppingAction(s);
    if (mode=="unknown") s->GetTrack()->SetTrackStatus(fStopAndKill);

  }
};
class Events : public G4UserEventAction {
  void BeginOfEventAction(const G4Event* e) override { CBDsimOpticalDiagnostics::BeginEvent(e->GetEventID()); }
  void EndOfEventAction(const G4Event* e) override { CBDsimOpticalDiagnostics::EndEvent(e->IsAborted()); }
};
class Tracks : public G4UserTrackingAction {
  void PreUserTrackingAction(const G4Track* t) override { CBDsimOpticalDiagnostics::PreUserTrackingAction(t); }
  void PostUserTrackingAction(const G4Track* t) override {
    CBDsimOpticalDiagnostics::PostUserTrackingAction(t);
  }
};
}
int main(int argc,char** argv) {
  if(argc!=6) return 2;
  mode=argv[1]; const int n=std::stoi(argv[2]); const long seed=std::stol(argv[3]);
  if (mode!="absorb" && mode!="other_absorb" && mode!="transmit" && mode!="bulk" &&
      mode!="detect" && mode!="qe" && mode!="sd_detect" && mode!="unknown" && mode!="wls" && mode!="no_rindex") return 2;
  setenv("CBDsim_OPTICAL_DIAG","1",1); setenv("CBDsim_OPTICAL_DIAG_HIST",argv[5],1);
  CLHEP::HepRandom::setTheSeed(seed);
  auto* rm=new G4RunManager;
  rm->SetUserInitialization(new Detector); rm->SetUserInitialization(new Physics);
  rm->SetUserAction(new Gun); auto* steps=new Steps; rm->SetUserAction(steps);
  rm->SetUserAction(new Tracks); rm->SetUserAction(new Events); rm->Initialize();
  CBDsimOpticalDiagnostics::ResetForRun(); rm->BeamOn(n);
  CBDsimOpticalDiagnostics::MergeWorkerIntoMaster();
  CBDsimOpticalDiagnostics::WriteSummaryFile(argv[4],n);
  CBDsimOpticalDiagnostics::WritePathHistogramFile(std::string(argv[4])+".hist");
  std::ifstream file(argv[4]); std::map<std::string,double> c; std::string k; double v;
  while(file>>k>>v) c[k]=v;
  std::string fate = mode=="no_rindex" ? "no_rindex" : mode=="transmit" || mode=="reflect" ? "world_exit" :
    mode=="bulk" ? "absorb_lg" : mode=="detect" || mode=="sd_detect" ? "detected" :
    mode=="qe" ? "qe_reject" : mode=="unknown" ? "unknown" :
    mode=="other_absorb" ? "absorb_other_surface" : mode=="wls" ? "wls_conversion" : "absorb_lg_surface";
  bool ok = c["budget_fate_"+fate]==n && c["budget_started"]==(mode=="wls" ? 2*n : n)
    && c["budget_terminal"]==c["budget_started"] && c["budget_fate_unfinished"]==0
    && c["budget_balance_delta"]==0 && c["budget_events"]==n
    && c["budget_source_primary_optical"]==n && c["budget_terminal_without_start"]==0;
  if (mode=="wls") ok &= c["budget_source_OpWLS"]==n && c["budget_fate_world_exit"]==n;
  std::cout << "BUDGET " << mode << " hist=" << argv[5] << " started=" << c["budget_started"]
    << " terminal=" << c["budget_terminal"] << " fate=" << fate << " count=" << c["budget_fate_"+fate] << " pass=" << ok << std::endl;
  delete rm; return ok ? 0 : 1;
}
