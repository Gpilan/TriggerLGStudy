#include "CBDsimPrimaryGeneratorAction.hh"
#include "CBDsimDetectorConstructionProto.hh"
#include "G4RunManager.hh"
#include "G4VUserPhysicsList.hh"
#include "G4Electron.hh"
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4UImanager.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <cmath>
class Physics: public G4VUserPhysicsList {
 void ConstructParticle() override {G4Electron::Definition();}
 void ConstructProcess() override {AddTransportation();}
 void SetCuts() override {}
};
int main() {
 auto* rm=new G4RunManager;auto* detector=new CBDsimDetectorConstructionProto;
 rm->SetUserInitialization(detector);rm->SetUserInitialization(new Physics);
 auto* gun=new CBDsimPrimaryGeneratorAction(42);rm->SetUserAction(gun);rm->Initialize();
 auto* ui=G4UImanager::GetUIpointer();
 auto command=[&](const char* s){if(ui->ApplyCommand(s)!=0) throw std::runtime_error(s);};
 command("/gun/particle e-");command("/gun/energy 60 GeV");
 const char* positions[]={"/gun/position 0 0 0 mm","/gun/position 18 0 0 mm","/gun/position 0 -18 0 mm","/gun/position 18 18 -2 mm"};
 G4ThreeVector expected[]={{0,0,0},{18,0,0},{0,-18,0},{18,18,-2}};
 G4Navigator nav;nav.SetWorldVolume(G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());
 for(int i=0;i<4;i++) for(int tilt=0;tilt<2;tilt++) {
  command(positions[i]);command(tilt?"/gun/direction 0.1 0.05 1":"/gun/direction 0 0 1");
  G4Event event(i);gun->GeneratePrimaries(&event);auto* v=event.GetPrimaryVertex();
  const auto p=v->GetPosition()/mm;const auto d=v->GetPrimary()->GetMomentumDirection();
  const auto want=tilt?G4ThreeVector(.1,.05,1).unit():G4ThreeVector(0,0,1);
  if((p-expected[i]).mag()>1e-12 || (d-want).mag()>1e-12) return 1;
  double k,x,y,z,dx,dy,dz;gun->GetLastPrimaryKinematics(k,x,y,z,dx,dy,dz);
  if(std::abs(k-60000)>1e-9 || (G4ThreeVector(x,y,z)-p).mag()>1e-12 || (G4ThreeVector(dx,dy,dz)-d).mag()>1e-12)return 2;
  for(int trig=0;trig<2;trig++) {
   const auto hit=p+d*((8.6+6*trig-p.z())/d.z());auto* pv=nav.LocateGlobalPointAndSetup(hit*mm);
   if(!pv || pv->GetName()!="protoScintPhys" || pv->GetCopyNo()!=trig)return 3;
   std::cout<<"PRIMARY case="<<i<<" tilt="<<tilt<<" trigger="<<trig<<" origin="<<p<<" dir="<<d<<" tile_center_cross="<<hit<<std::endl;
  }
 }
 command("/gun/position 2 3 0 mm");command("/CBDsim/generator/randx 2 mm");command("/CBDsim/generator/randy 4 mm");
 for(int i=0;i<100;i++) {G4Event e(i);gun->GeneratePrimaries(&e);auto p=e.GetPrimaryVertex()->GetPosition()/mm;
  if(p.x()<1||p.x()>3||p.y()<1||p.y()>5||p.z()!=0)return 4;}
 std::cout<<"PASS 8 kinematics cases, 16 tile intersections, 100 spread events"<<std::endl;delete rm;
}
