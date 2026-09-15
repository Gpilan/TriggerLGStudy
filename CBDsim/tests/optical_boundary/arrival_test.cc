// D02 integration test: known, nonoverlapping layered geometry and repeated window entries.
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
#include "G4UserEventAction.hh"
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
G4OpBoundaryProcess* boundary=nullptr;
class Detector : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override {
    auto material=[](const char* name) {
      auto* m=new G4Material(name,1.,1.01*g/mole,1.*g/cm3);
      auto* p=new G4MaterialPropertiesTable;
      G4double e[]={2*eV,4*eV}, n[]={1.5,1.5};
      p->AddProperty("RINDEX",e,n,2); m->SetMaterialPropertiesTable(p); return m;
    };
    auto* w=new G4LogicalVolume(new G4Box("world",100*mm,100*mm,100*mm),material("outside"),"world");
    auto* wp=new G4PVPlacement(nullptr,{},w,"protoWorldPhys",nullptr,false,0);
    auto layer=[&](const char* name,double low,double high) {
      auto* lv=new G4LogicalVolume(new G4Box(name,10*mm,10*mm,(high-low)/2),material(name),name);
      return new G4PVPlacement(nullptr,{0,0,(high+low)/2},lv,name,w,false,0,true);
    };
    auto* scint=layer("protoScintPhys",-20*mm,mode=="nolg" ? 0 : -12*mm);
    if(mode!="nolg") {
      layer("protoScintLGGelPhys",-12*mm,-10*mm);
      layer("protoLightGuidePhys",-10*mm,0);
    }
    auto* gel=layer("protoSipmGelPhys",0,2*mm);
    auto* glass=layer("protoSipmWindowPhys",2*mm,4*mm);
    auto mirror=[&](const char* name,G4VPhysicalVolume* a,G4VPhysicalVolume* b) {
      auto* s=new G4OpticalSurface(name,glisur,polished,dielectric_metal);
      auto* p=new G4MaterialPropertiesTable;
      G4double e[]={2*eV,4*eV}, r[]={1.,1.};p->AddProperty("REFLECTIVITY",e,r,2);
      s->SetMaterialPropertiesTable(p);new G4LogicalBorderSurface(name,a,b,s);
    };
    if(mode=="reflect") mirror("front",gel,glass);
    if(mode=="reenter") { mirror("back",glass,wp);mirror("left",scint,wp); }
    return wp;
  }
};
class Physics : public G4VUserPhysicsList {
  void ConstructParticle() override {G4OpticalPhoton::OpticalPhotonDefinition();}
  void ConstructProcess() override {
    AddTransportation();auto* pm=G4OpticalPhoton::OpticalPhotonDefinition()->GetProcessManager();
    pm->AddDiscreteProcess(new G4OpAbsorption);boundary=new G4OpBoundaryProcess;pm->AddDiscreteProcess(boundary);
  }
  void SetCuts() override {}
};
class Gun : public G4VUserPrimaryGeneratorAction {
  G4ParticleGun gun{2}; // Two track IDs per event; IDs repeat across events.
  void GeneratePrimaries(G4Event* e) override {
    gun.SetParticleDefinition(G4OpticalPhoton::OpticalPhotonDefinition());gun.SetParticleEnergy(3*eV);
    gun.SetParticlePosition({0,0,-16*mm});gun.SetParticleMomentumDirection({0,0,1});
    gun.SetParticlePolarization({1,0,0});gun.GeneratePrimaryVertex(e);
  }
};
class Events : public G4UserEventAction {
  void BeginOfEventAction(const G4Event*) override {CBDsimOpticalDiagnostics::BeginEvent();}
};
class Steps : public G4UserSteppingAction {
 public:
  long entries=0, reflections=0;
  std::map<int,int> trackEntries;
  void UserSteppingAction(const G4Step* s) override {
    CBDsimOpticalDiagnostics::UserSteppingAction(s);
    auto* t=s->GetTrack();bool stop=false;
    if(s->GetPostStepPoint()->GetStepStatus()==fGeomBoundary) {
      const auto status=boundary->GetStatus();
      const double z=s->GetPostStepPoint()->GetPosition().z()/mm;
      // Independent oracle uses known plane coordinate/direction, not production volume classification.
      if(std::abs(z-2.)<1.e-6 && s->GetPreStepPoint()->GetMomentumDirection().z()>0) {
        if(status==FresnelRefraction) {
          ++entries;int count=++trackEntries[t->GetTrackID()];
          stop=mode!="reenter" || count==2;
        } else if(status==SpikeReflection) {++reflections;stop=mode=="reflect";}
      }
      if(mode=="inlet" && std::abs(z+10.)<1.e-6) stop=true;
    }
    if(t->GetCurrentStepNumber()>100) {std::cerr<<"Unexpected transport loop\n";std::exit(3);}
    if(stop) {trackEntries.erase(t->GetTrackID());t->SetTrackStatus(fStopAndKill);}
  }
};
class Tracks : public G4UserTrackingAction {
  void PostUserTrackingAction(const G4Track* t) override {CBDsimOpticalDiagnostics::PostUserTrackingAction(t);}
};
bool mappings() {
  const std::map<std::string,std::string> cases={
    {"protoScintLGGelPhys","inlet_gel"},{"protoScintLGGelLog","inlet_gel"},
    {"protoSipmGelPhys","outlet_gel"},{"protoSipmGelLog","outlet_gel"},
    {"protoScintPhys","scint"},{"protoScintLog","scint"},
    {"protoLightGuidePhys","lg"},{"protoLightGuideLog","lg"},
    {"protoSipmWindowPhys","sipm_glass"},{"protoSipmWindowLog","sipm_glass"},
    {"protoSipmWaferPhys","sipm_wafer"},{"protoSipmWaferLog","sipm_wafer"},
    {"protoSipmEnvPhys","sipm_env"},{"protoSipmEnvLog","sipm_env"},
    {"protoFoilWrapPhys","foil"},{"protoFoilCornerPhys","foil"},
    {"protoFoilTipRingPhys","foil"},{"protoFoilYmPhys","foil"},
    {"protoWorldPhys","world"},{"protoWorldLog","world"},
    {"worldPhysical","world"},{"worldLogical","world"},{"unmapped","other"},{"none","other"}};
  for(const auto& p:cases) if(CBDsimOpticalDiagnostics::RegionForVolume(p.first)!=p.second) return false;
  return true;
}
}
int main(int argc,char** argv) {
  if(argc!=6) return 2;
  mode=argv[1];const int n=std::stoi(argv[2]);const long seed=std::stol(argv[3]);
  if(mode!="forward" && mode!="inlet" && mode!="reflect" && mode!="reenter" && mode!="nolg") return 2;
  setenv("CBDsim_OPTICAL_DIAG","1",1);setenv("CBDsim_OPTICAL_DIAG_HIST",argv[5],1);
  CLHEP::HepRandom::setTheSeed(seed);
  auto* rm=new G4RunManager;rm->SetUserInitialization(new Detector);rm->SetUserInitialization(new Physics);
  rm->SetUserAction(new Gun);rm->SetUserAction(new Events);auto* s=new Steps;rm->SetUserAction(s);
  rm->SetUserAction(new Tracks);rm->Initialize();CBDsimOpticalDiagnostics::ResetForRun();rm->BeamOn(n);
  CBDsimOpticalDiagnostics::MergeWorkerIntoMaster();CBDsimOpticalDiagnostics::WriteSummaryFile(argv[4],n);
  std::ifstream f(argv[4]);std::map<std::string,double> c;std::string k;double v;while(f>>k>>v)c[k]=v;
  const long photons=2*n;const bool reach=mode!="inlet" && mode!="reflect";
  const long expectedEntries=reach ? photons*(mode=="reenter" ? 2 : 1) : 0;
  bool ok=mappings() && c["unique_window_arrivals"]==(reach ? photons : 0)
    && c["window_entry_steps"]==expectedEntries && s->entries==expectedEntries;
  ok &= c["boundary_transmit_outlet_gel_to_sipm_glass"]==expectedEntries;
  if(mode=="reflect") ok &= s->reflections==photons && c["boundary_encounter_outlet_gel_to_sipm_glass"]==photons;
  if(mode=="inlet" || mode=="forward" || mode=="reflect")
    ok &= c["boundary_transmit_scint_to_inlet_gel"]==photons
       && c["boundary_transmit_inlet_gel_to_lg"]==photons
       && c["scint_to_lg"]==0 && std::abs(c["path_inlet_gel_mm_per_event"]-4.)<1.e-6;
  if(mode=="inlet") ok &= c["lg_to_sipm"]==0 && c["path_outlet_gel_mm_per_event"]==0;
  if(mode=="forward") ok &= c["lg_to_sipm"]==2*photons; // explicitly legacy sum, not arrival
  if(mode=="nolg") ok &= c["path_inlet_gel_mm_per_event"]==0 && c["path_lg_mm_per_event"]==0;
  std::cout<<"D02 "<<mode<<" hist="<<argv[5]<<" photons="<<photons<<" entries="<<c["window_entry_steps"]
    <<" unique="<<c["unique_window_arrivals"]<<" pass="<<ok<<std::endl;
  delete rm;return ok ? 0 : 1;
}
