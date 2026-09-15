// Exercises the production solids, materials, surfaces and SD; no substitute geometry.
#include "CBDsimDetectorConstructionProto.hh"
#include "CBDsimOpticalDiagnostics.hh"
#include "G4RunManager.hh"
#include "G4Navigator.hh"
#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4TessellatedSolid.hh"
#include "G4Polyhedron.hh"
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
    const int i=event->GetEventID();
    // Grid contains near-edge and near-corner points, both assemblies.
    const int trigger=i%2, cell=(i/2)%100;
    const double hx=noLG ? 7.5 : 20.;
    const double x=(-1.+2.*(cell%10)/9.)*(hx-0.001)*mm;
    const double z=(-1.+2.*(cell/10)/9.)*(2.5-0.001)*mm;
    G4VPhysicalVolume* scint=nullptr;
    for(auto* pv:*G4PhysicalVolumeStore::GetInstance())
      if(pv->GetName()=="protoScintPhys" && pv->GetCopyNo()==trigger) scint=pv;
    const auto rot=scint->GetObjectRotationValue();
    gun.SetParticleDefinition(G4OpticalPhoton::Definition());gun.SetParticleEnergy(3.*eV);
    gun.SetParticlePosition(rot*G4ThreeVector(x,-29.999*mm,z)+scint->GetObjectTranslation());
    gun.SetParticleMomentumDirection(rot*G4ThreeVector(0,-1,0));
    gun.SetParticlePolarization(rot*G4ThreeVector(1,0,0));gun.GeneratePrimaryVertex(event);
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
    if(s->GetTrack()->GetCurrentStepNumber()==1) {
      auto* a=s->GetPreStepPoint()->GetPhysicalVolume();auto* b=s->GetPostStepPoint()->GetPhysicalVolume();
      if(a && b && a->GetName()=="protoScintPhys" && b->GetName()==(noLG?"protoSipmGelPhys":"protoScintLGGelPhys")) ++firstGood;
      else {++firstBad;std::cout<<"BAD_FIRST "<<(a?a->GetName():"none")<<" -> "<<(b?b->GetName():"none")<<std::endl;}
    }
    if(s->GetTrack()->GetCurrentStepNumber()>10000) {++loops;s->GetTrack()->SetTrackStatus(fStopAndKill);}
  }
};
}
int main(int argc,char** argv) {
  if(argc!=6 && argc!=7) return 2;
  const std::string mode=argv[1]; const int n=std::stoi(argv[2]);const long seed=std::stol(argv[3]);
  noLG=std::stoi(argv[5])!=0;
  setenv("CBDsim_PROTO_NO_LG",noLG?"1":"0",1);
  setenv("CBDsim_OPTICAL_DIAG","1",1);setenv("CBDsim_OPTICAL_DIAG_HIST","1",1);
  CLHEP::HepRandom::setTheSeed(seed);
  // Validation-only tolerance sensitivity; never changes production defaults.
  if(argc==7) G4GeometryManager::GetInstance()->SetWorldMaximumExtent(std::stod(argv[6])*mm);
  std::cout<<"SURFACE_TOLERANCE_MM "<<G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm<<std::endl;
  auto* rm=new G4RunManager;rm->SetUserInitialization(new CBDsimDetectorConstructionProto);
  rm->SetUserInitialization(new Physics);rm->SetUserAction(new Gun);rm->SetUserAction(new Events);
  rm->SetUserAction(new Tracks);rm->SetUserAction(new Steps);rm->Initialize();
  bool ok=true;
  if(mode=="overlap") {
    int checked=0, failures=0;
    for(auto* pv:*G4PhysicalVolumeStore::GetInstance()) {
      G4ThreeVector lo,hi;pv->GetLogicalVolume()->GetSolid()->BoundingLimits(lo,hi);
      std::cout<<"PLACEMENT "<<pv->GetName()<<" copy="<<pv->GetCopyNo()<<" local_center="<<pv->GetObjectTranslation()
        <<" local_bbox="<<lo<<":"<<hi<<std::endl;
      if(!pv->GetMotherLogical()) continue;
      ++checked;if(pv->CheckOverlaps(n,0.,true,20)) ++failures;
    }
    ok=failures==0;std::cout<<"OVERLAP checked="<<checked<<" failing_placements="<<failures<<" samples="<<n<<" tolerance_mm=0 seed="<<seed<<" pass="<<ok<<std::endl;
  } else if(mode=="geometry") {
    G4VPhysicalVolume* world=nullptr;
    for(auto* pv:*G4PhysicalVolumeStore::GetInstance()) if(!pv->GetMotherLogical()) world=pv;
    G4Navigator nav;nav.SetWorldVolume(world);
    long checked=0, bad=0;
    for(auto* scint:*G4PhysicalVolumeStore::GetInstance()) {
      if(scint->GetName()!="protoScintPhys") continue;
      const auto rot=scint->GetObjectRotationValue();
      auto check=[&](G4ThreeVector local,const char* expected) {
        auto* pv=nav.LocateGlobalPointAndSetup(rot*local+scint->GetObjectTranslation(),nullptr,false);
        ++checked;
        if(!pv || pv->GetName()!=expected) {++bad;std::cout<<"BAD_POINT local="<<local<<" expected="<<expected<<" got="<<(pv?pv->GetName():"none")<<std::endl;}
      };
      for(double eps : {1.e-6*mm,1.e-5*mm,1.e-4*mm}) {
        const char* ext = std::getenv("CBDsim_PROTO_TAPE_EXTENSIONS");
        const char* tapeOpt = std::getenv("CBDsim_PROTO_GEL_TAPE");
        const bool extendedTape = (!ext || std::string(ext)=="1") && (!tapeOpt || std::string(tapeOpt)=="1");
        if (extendedTape) {
          const double hx=noLG?7.5:20.;
          const double cy=noLG?-30.25:-30.01;
          for (int sign : {-1,1}) {
            check({sign*(hx+eps)*mm,cy*mm,0.},"protoGelTapePhys");
            check({0.,cy*mm,sign*(2.5*mm+eps)},"protoGelTapePhys");
            check({sign*(hx+eps)*mm,cy*mm,sign*(2.5*mm+eps)},"protoGelTapePhys");
            check({sign*(hx-eps)*mm,cy*mm,0.},noLG?"protoSipmWindowPhys":"protoScintLGGelPhys");
            check({sign*(hx+.05)*mm+sign*eps,cy*mm,0.},"protoWorldPhys");
          }
        }
        const char* seal = std::getenv("CBDsim_PROTO_CORNER_SEAL");
        if (!seal || std::string(seal)=="1") {
          for (int sign : {-1,1}) {
            check({sign*(20.*mm+eps),-29.99*mm,0.},"protoFoilCornerSealPhys");
            check({0.,-29.99*mm,sign*(2.5*mm+eps)},"protoFoilCornerSealPhys");
            if (!noLG) check({0.,-30.01*mm,sign*(2.5*mm+eps)},extendedTape?"protoGelTapePhys":"protoWorldPhys");
          }
        }
        const char* tape = std::getenv("CBDsim_PROTO_GEL_TAPE");
        if (noLG && (!tape || std::string(tape)=="1")) {
          for (int sign : {-1,1}) {
            check({sign*(7.5*mm-eps),-30.05*mm,0.},"protoSipmGelPhys");
            check({sign*(7.5*mm+eps),-30.05*mm,0.},"protoGelTapePhys");
            check({0.,-30.05*mm,sign*(2.5*mm+eps)},"protoGelTapePhys");
            check({0.,-30.05*mm,sign*(2.55*mm+eps)},"protoWorldPhys");
            // Former foil-strip overlap and the 1 um wrap lip are occupied by tape.
            check({sign*(7.5*mm+eps),-30.02*mm,0.},"protoGelTapePhys");
            check({sign*(7.55*mm+eps),-30.02*mm,0.},"protoFoilYmPhys");
            check({0.,-30.0005*mm,sign*(2.52*mm)},"protoGelTapePhys");
          }
        }

        for(int ix=0;ix<10;++ix) for(int iz=0;iz<10;++iz) {
          const double x=(-1.+2.*ix/9.)*((noLG?7.5:20.)-0.001)*mm;
          const double z=(-1.+2.*iz/9.)*2.499*mm;
          check({x,-30.*mm+eps,z},"protoScintPhys");
          check({x,-30.*mm-eps,z},noLG?"protoSipmGelPhys":"protoScintLGGelPhys");
          if(!noLG) {
            check({x,-30.02*mm+eps,z},"protoScintLGGelPhys");
            check({x,-30.02*mm-eps,z},"protoLightGuidePhys");
          }
        }
        const double tip=(noLG?-30.:-60.02)*mm;
        if (!noLG) {
          const char* contact = std::getenv("CBDsim_PROTO_TIP_CONTACT");
          const bool touching = !contact || std::string(contact)=="1";
          for (int j=0;j<32;++j) {
            const double phi=twopi*j/32.;
            for (double y : {-60.07,-60.25,-60.415}) {
              const char* inside = y>-60.12 ? "protoSipmGelPhys" : (y>-60.41 ? "protoSipmWindowPhys" : "protoSipmWaferPhys");
              check({(7.5*mm-eps)*std::cos(phi),y*mm,(7.5*mm-eps)*std::sin(phi)},inside);
              check({(7.5*mm+eps)*std::cos(phi),y*mm,(7.5*mm+eps)*std::sin(phi)},touching?"protoFoilTipRingPhys":"protoWorldPhys");
            }
          }
        }
        for(int j=0;j<100;++j) {
          const double phi=twopi*j/100.;
          const double x=noLG?7.499*std::cos(phi)*mm:7.49*std::cos(phi)*mm;
          const double z=noLG?2.499*std::sin(phi)*mm:7.49*std::sin(phi)*mm;
          check({x,tip+eps,z},noLG?"protoScintPhys":"protoLightGuidePhys");
          check({x,tip-eps,z},"protoSipmGelPhys");
          check({x,tip-0.1*mm+eps,z},"protoSipmGelPhys");
          check({x,tip-0.1*mm-eps,z},"protoSipmWindowPhys");
          check({x,tip-0.39*mm+eps,z},"protoSipmWindowPhys");
          check({x,tip-0.39*mm-eps,z},"protoSipmWaferPhys");
        }
      }
    }
    ok=bad==0;
    if(!noLG) {
      for(auto* pv:*G4PhysicalVolumeStore::GetInstance()) {
        if(pv->GetName()!="protoLightGuidePhys" || pv->GetCopyNo()!=0) continue;
        auto* ts=dynamic_cast<G4TessellatedSolid*>(pv->GetLogicalVolume()->GetSolid());
        if (!ts) {
          const char* round=std::getenv("CBDsim_PROTO_ROUND_TIP");
          ok &= !round || std::string(round)=="1";
          auto* solid=pv->GetLogicalVolume()->GetSolid();
          auto* poly=solid->GetPolyhedron();
          ok &= poly && poly->GetNoFacets()>0;
          std::cout<<"ROUND_TIP visual_facets="<<(poly?poly->GetNoFacets():0)<<std::endl;
          for (int j=0;j<256;++j) {
            double phi=twopi*(j+.5)/256.;
            G4ThreeVector q(7.4998*mm*std::cos(phi),-60.019*mm,7.4998*mm*std::sin(phi));
            ok &= solid->Inside(q)==kInside;
            ok &= std::abs(solid->DistanceToOut(q,G4ThreeVector(0,-1,0))-.001*mm)<1.e-7*mm;
          }
          ok &= solid->Inside(G4ThreeVector(0,-60.021*mm,0))==kOutside;
          ok &= solid->Inside(G4ThreeVector(0,-29.999*mm,0))==kOutside;
          std::cout<<"ROUND_TIP outlet rays=256 total_length=30 pass="<<ok<<std::endl;
          continue;
        }
        double inletArea=0.,outletArea=0.,volume=0.,sag=0.;
        for(int j=0;j<ts->GetNumberOfFacets();++j) {
          const auto* f=ts->GetFacet(j);auto a=f->GetVertex(0),b=f->GetVertex(1),c=f->GetVertex(2);
          volume+=a.dot(b.cross(c))/6.;
          if(std::abs(a.y()+30.02*mm)<1.e-9 && std::abs(b.y()-a.y())<1.e-9 && std::abs(c.y()-a.y())<1.e-9) inletArea+=f->GetArea();
          if(std::abs(a.y()+60.02*mm)<1.e-9 && std::abs(b.y()-a.y())<1.e-9 && std::abs(c.y()-a.y())<1.e-9) {
            outletArea+=f->GetArea();const auto mid=(b+c)*0.5;sag=std::max(sag,7.5*mm-std::hypot(mid.x(),mid.z()));
          }
        }
        // Before looking at results: exact rectangle, outlet error <=0.05%, sag <=3 um.
        ok &= std::abs(inletArea/(mm*mm)-200.)<1.e-8 &&
          std::abs(outletArea/(pi*7.5*7.5*mm*mm)-1.)<0.0005 && sag<0.003*mm && volume>0.;
        std::cout<<std::setprecision(12)<<"MESH inlet_mm2="<<inletArea/(mm*mm)<<" outlet_mm2="<<outletArea/(mm*mm)
          <<" volume_mm3="<<volume/(mm*mm*mm)<<" max_sag_mm="<<sag/mm<<" facets="<<ts->GetNumberOfFacets()<<std::endl;
      }
    }
    std::cout<<"GEOMETRY checked="<<checked<<" bad="<<bad<<" eps_mm=1e-6,1e-5,1e-4 pass="<<ok<<std::endl;
  } else if(mode=="photons") {
    CBDsimOpticalDiagnostics::ResetForRun();rm->BeamOn(n);CBDsimOpticalDiagnostics::MergeWorkerIntoMaster();
    CBDsimOpticalDiagnostics::WriteSummaryFile(argv[4],n);CBDsimOpticalDiagnostics::WritePathHistogramFile(std::string(argv[4])+".hist");
    std::ifstream in(argv[4]);std::map<std::string,double> c;std::string k;double v;while(in>>k>>v)c[k]=v;
    ok=firstBad==0 && firstGood==n && loops==0 && c["budget_closed"]==1 && c["budget_started"]==n &&
      c["budget_fate_unknown"]==0 && c["budget_fate_no_rindex"]==0 && c["budget_fate_unfinished"]==0;
    std::cout<<"PHOTONS first_good="<<firstGood<<" first_bad="<<firstBad<<" no_rindex="<<c["budget_fate_no_rindex"]
      <<" loops="<<loops<<" pass="<<ok<<std::endl;
  } else return 2;
  delete rm;return ok?0:1;
}
