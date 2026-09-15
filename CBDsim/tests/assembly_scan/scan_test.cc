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
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
class Physics: public G4VUserPhysicsList {
 void ConstructParticle() override {G4Electron::Definition();}
 void ConstructProcess() override {AddTransportation();}
 void SetCuts() override {}
};
int main() {
 const double s=std::atof(std::getenv("CBDsim_PROTO_SCAN_S_MM"));
 auto* rm=new G4RunManager;
 rm->SetUserInitialization(new CBDsimDetectorConstructionProto);
 rm->SetUserInitialization(new Physics);
 auto* gun=new CBDsimPrimaryGeneratorAction(42);rm->SetUserAction(gun);rm->Initialize();
 auto* ui=G4UImanager::GetUIpointer();
 for(auto command:{"/gun/particle e-","/gun/energy 60 GeV","/gun/position 0 0 0 mm","/gun/direction 0 0 1"})
  if(ui->ApplyCommand(command)!=0)return 1;
 { // Destroy the event and local navigator before the run manager owns teardown.
 G4Event event(0);gun->GeneratePrimaries(&event);auto* v=event.GetPrimaryVertex();
 if(v->GetPosition().mag()>1e-12 || (v->GetPrimary()->GetMomentumDirection()-G4ThreeVector(0,0,1)).mag()>1e-12)return 2;
 auto* world=G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();
 std::cout<<std::setprecision(14);
 for(int i=0;i<world->GetLogicalVolume()->GetNoDaughters();++i){
  auto* p=world->GetLogicalVolume()->GetDaughter(i);
  const auto r=p->GetObjectRotationValue();
  std::cout<<"PLACEMENT "<<p->GetName()<<" "<<p->GetCopyNo()<<" "<<p->GetObjectTranslation()/mm
           <<" "<<r.xx()<<","<<r.xy()<<","<<r.xz()<<","<<r.yx()<<","<<r.yy()<<","<<r.yz()<<","<<r.zx()<<","<<r.zy()<<","<<r.zz()<<std::endl;
  if(p->GetName()=="protoScintPhys"){
   G4ThreeVector lo,hi;p->GetLogicalVolume()->GetSolid()->BoundingLimits(lo,hi);
   std::cout<<"TILE_BOUNDS "<<p->GetCopyNo()<<" "<<lo/mm<<" "<<hi/mm<<std::endl;
  }
 }
 G4Navigator nav;nav.SetWorldVolume(world);
 for(int trig=0;trig<2;trig++)for(double dz:{-2.5+1e-5,0.,2.5-1e-5}){
  const G4ThreeVector hit(0,0,8.6+6*trig+dz);auto* pv=nav.LocateGlobalPointAndSetup(hit*mm);
  if(!pv || pv->GetName()!="protoScintPhys" || pv->GetCopyNo()!=trig)return 3;
  const auto local=pv->GetObjectRotationValue().inverse()*(hit*mm-pv->GetObjectTranslation())/mm;
  if((local-G4ThreeVector(0,s,dz)).mag()>1e-9)return 4;
  std::cout<<"INTERSECTION trigger="<<trig<<" world_mm="<<hit<<" local_mm="<<local<<" far_end_distance_mm="<<30-local.y()<<std::endl;
 }
 std::cout<<"PASS fixed production primary; six tile intersections"<<std::endl;
 }
 delete rm;
}
