#ifndef CBDsimBoundaryProcess_hh
#define CBDsimBoundaryProcess_hh
// Recovery of genuine opaque entries skipped as StepTooSmall in Geant4 11.2.0/11.4.0.
// See docs/BOUNDARY-RECOVERY.md for scope, dispatch-marker semantics and validation.
#include "G4Version.hh"
#include "G4Box.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4GeometryTolerance.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4TouchableHistory.hh"
#include "G4ParticleChange.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4SystemOfUnits.hh"
#include <memory>
#include <vector>
#include <cmath>
#include <iomanip>
static_assert(G4VERSION_NUMBER == 1120 || G4VERSION_NUMBER == 1140,
              "Review boundary adapter before changing Geant4 version");

class CBDsimBoundaryProcess final : public G4OpBoundaryProcess {
 public:
  CBDsimBoundaryProcess():G4OpBoundaryProcess("OpBoundary"), tolerance_(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) {}
  void StartTracking(G4Track* track) override {
    expected_.reset(); proxyStep_.reset(); proxyTrack_.reset();
    G4OpBoundaryProcess::StartTracking(track);
  }
  void EndTracking() override {
    // Release touchable handles while the worker's geometry/allocators are alive.
    G4OpBoundaryProcess::EndTracking();
    expected_.reset(); proxyStep_.reset(); proxyTrack_.reset();
  }
  G4VParticleChange* PostStepDoIt(const G4Track& track,const G4Step& step) override {
    const auto* pre=step.GetPreStepPoint(); const auto* post=step.GetPostStepPoint();
    const G4StepPoint* incident=pre;
    bool remembered=false;
    if (!HasRindex(pre->GetMaterial()) && expected_ &&
        (post->GetPosition()-lastBoundary_).mag()<=8.*tolerance_) {
      incident=expected_.get(); remembered=true;
    }
    G4OpticalSurface* surface=nullptr; G4ThreeVector normal;
    bool correct=step.GetPostStepPoint()->GetStepStatus()==fGeomBoundary &&
      track.GetStepLength()<=tolerance_ && G4ParallelWorldProcess::GetHyperStep()==nullptr &&
      HasRindex(incident->GetMaterial()) && !HasRindex(post->GetMaterial()) &&
      SupportedSurface(incident,post,track.GetKineticEnergy(),surface) &&
      ForwardInside(post,track.GetMomentumDirection(),8.*tolerance_) &&
      ForwardInside(post,track.GetMomentumDirection(),4.*tolerance_);
    if(correct) {
      G4bool valid=false;
      normal=G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()
        ->GetGlobalExitNormal(post->GetPosition(),&valid);
      correct=valid && normal.mag2()>0. && std::abs(normal.dot(track.GetMomentumDirection()))>1.e-12;
    }
    if (correct && surface->GetMaterialPropertiesTable()->GetProperty(kREFLECTIVITY)
                       ->Value(track.GetKineticEnergy()) > 0.) {
      const auto n = normal.unit();
      const auto reflected = track.GetMomentumDirection() - 2.*track.GetMomentumDirection().dot(n)*n;
      if (ForwardInside(post, reflected, 8.*tolerance_))
        G4Exception("CBDsimBoundaryProcess", "AmbiguousOpaqueNormal", FatalException,
                    "Opaque entry has a normal that would reflect into its bulk; repair the geometry.");
    }
    // Bounded edge ambiguity: the scint exit face can differ from the unique
    // receiving rim-box face within tolerance. Do not reflect into opaque bulk.
    G4ThreeVector receivingExit;
    const bool overrideNormal = track.GetStepLength()>tolerance_ &&
      post->GetStepStatus()==fGeomBoundary && !G4ParallelWorldProcess::GetHyperStep() &&
      ReceivingRimNormal(pre,post,track.GetMomentumDirection(),track.GetKineticEnergy(),receivingExit);
    ScopedNormal normalScope(overrideNormal,receivingExit);
    G4VParticleChange* result=nullptr;
    if(correct) {
      // G4VParticleChange retains a track pointer until its proposals are applied.
      // Both proxies therefore live as process members through the next invocation.
      // The proxy length is a dispatch marker above both short-step guards.
      // Only the track's short-step flag is changed; actual path length/time/points
      // remain in the unchanged copied G4Step used by ParticleChange initialization.
      proxyTrack_=std::make_unique<G4Track>(track);
      proxyStep_=std::make_unique<G4Step>(step);
      proxyTrack_->SetTrackID(track.GetTrackID());
      proxyTrack_->SetParentID(track.GetParentID());
      proxyTrack_->SetStep(proxyStep_.get());
      proxyStep_->SetTrack(proxyTrack_.get());
      proxyTrack_->SetStepLength(16.*tolerance_);
      auto* proxyPre=proxyStep_->GetPreStepPoint();
      proxyPre->SetMaterial(incident->GetMaterial());
      proxyPre->SetMaterialCutsCouple(incident->GetMaterialCutsCouple());
      proxyPre->SetTouchableHandle(incident->GetTouchableHandle());
      proxyTrack_->SetTouchableHandle(incident->GetTouchableHandle());
      result=G4OpBoundaryProcess::PostStepDoIt(*proxyTrack_,*proxyStep_);
      if (result->GetTrueStepLength() != step.GetStepLength())
        G4Exception("CBDsimBoundaryProcess", "BoundaryTrueLength", FatalException,
                    "Boundary dispatch must preserve the actual geometrical path length.");
      auto* change=dynamic_cast<G4ParticleChange*>(result);
      const G4ThreeVector outgoing=change ? *change->GetMomentumDirection():track.GetMomentumDirection();
      G4cout<<std::setprecision(17)<<"[Boundary tiny opaque] track="<<track.GetTrackID()
        <<" step="<<track.GetCurrentStepNumber()<<" actual_pre="<<Volume(pre)
        <<" effective_pre="<<Volume(incident)<<" post="<<Volume(post)
        <<" remembered="<<remembered<<" surface="<<surface->GetName()
        <<" original_length_mm="<<step.GetStepLength()/mm
        <<" proxy_flag_length_mm="<<proxyTrack_->GetStepLength()/mm
        <<" preserved_true_length_mm="<<result->GetTrueStepLength()/mm
        <<" position_mm="<<post->GetPosition()/mm<<" incoming="<<track.GetMomentumDirection()
        <<" navigator_normal="<<normal<<" outgoing="<<outgoing
        <<" reflected_forward_inside="<<ForwardInside(post,outgoing,4.*tolerance_)
        <<" status="<<static_cast<int>(GetStatus())<<G4endl;
    } else {
      result=G4OpBoundaryProcess::PostStepDoIt(track,step);
    }
    const auto status=GetStatus();
    if(status==StepTooSmall) {
      // Keep expected medium only within the local relocation chain.
      if(expected_ && (post->GetPosition()-lastBoundary_).mag()>8.*tolerance_) expected_.reset();
    } else if(IsReflection(status)) {
      if(HasRindex(incident->GetMaterial())) {
        expected_=std::make_unique<G4StepPoint>(*incident);lastBoundary_=post->GetPosition();
      } else expected_.reset();
    } else if(status==Transmission || status==FresnelRefraction || status==SameMaterial ||
              status==CoatedDielectricRefraction || status==CoatedDielectricFrustratedTransmission) {
      if(HasRindex(post->GetMaterial())) {
        expected_=std::make_unique<G4StepPoint>(*post);lastBoundary_=post->GetPosition();
      } else expected_.reset();
    } else expected_.reset();
    return result;
  }
 private:
  class ExitNormalNavigator final : public G4Navigator {
   public:
    G4ThreeVector normal;
    G4ThreeVector GetGlobalExitNormal(const G4ThreeVector&,G4bool* valid) override {
      *valid=true;return normal;
    }
  };
  // Only G4OpBoundaryProcess's synchronous normal query sees this facade.
  // Tracking navigator and transport geometry are never replaced.
  class ScopedNormal {
    G4Navigator** slot_=nullptr;G4Navigator* old_=nullptr;
    std::unique_ptr<ExitNormalNavigator> facade_;
   public:
    ScopedNormal(bool active,const G4ThreeVector& normal) {
      if(!active)return;
      auto* tm=G4TransportationManager::GetTransportationManager();
      if(tm->GetNoActiveNavigators()!=1 || G4ParallelWorldProcess::GetHypNavigatorID()!=0)
        G4Exception("CBDsimBoundaryProcess","RimNormalNavigator",FatalException,"Expected one mass-world navigator");
      slot_=&*tm->GetActiveNavigatorsIterator();old_=*slot_;facade_=std::make_unique<ExitNormalNavigator>();facade_->normal=normal;*slot_=facade_.get();
    }
    ~ScopedNormal(){if(slot_)*slot_=old_;}
  };
  bool ReceivingRimNormal(const G4StepPoint* pre,const G4StepPoint* post,
                          const G4ThreeVector& direction,G4double energy,G4ThreeVector& exit) {
    if(Volume(pre)!="protoScintPhys" || Volume(post)!="protoFoilCornerSealPhys" ||
       !HasRindex(pre->GetMaterial()) || HasRindex(post->GetMaterial()))return false;
    auto* box=dynamic_cast<G4Box*>(post->GetPhysicalVolume()->GetLogicalVolume()->GetSolid());
    auto* touch=dynamic_cast<const G4TouchableHistory*>(post->GetTouchable());
    G4OpticalSurface* surface=nullptr;
    if(!box || box->GetName()!="protoNoLGCornerSealBox" || !touch ||
       !SupportedSurface(pre,post,energy,surface) ||
       !ForwardInside(post,direction,4.*tolerance_) || !ForwardInside(post,direction,8.*tolerance_) ||
       !IncidentBefore(pre,post,direction,8.*tolerance_))return false;
    G4bool valid=false;
    const auto nav=G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()
      ->GetGlobalExitNormal(post->GetPosition(),&valid);
    if(!valid || nav.mag2()==0.)return false;
    const auto n=nav.unit();const auto wrong=direction-2.*direction.dot(n)*n;
    if(!ForwardInside(post,wrong,4.*tolerance_) || !ForwardInside(post,wrong,8.*tolerance_))return false;
    const auto tr=touch->GetHistory()->GetTopTransform();
    const auto q=tr.TransformPoint(post->GetPosition());
    const G4double half[3]={box->GetXHalfLength(),box->GetYHalfLength(),box->GetZHalfLength()};
    G4ThreeVector local;int faces=0;
    for(int a=0;a<3;++a) {
      const auto distance=std::abs(std::abs(q[a])-half[a]);
      if(distance<=tolerance_) {local[a]=q[a]>=0.?1.:-1.;++faces;}
      else if(distance<=8.*tolerance_)return false;
    }
    if(faces!=1)return false;
    const auto outward=tr.Inverse().TransformAxis(local).unit();
    if(direction.dot(outward)>=-1.e-12)return false;
    const auto reflected=direction-2.*direction.dot(outward)*outward;
    if(ForwardInside(post,reflected,4.*tolerance_) || ForwardInside(post,reflected,8.*tolerance_))return false;
    exit=-outward;
    G4cout<<"[Boundary rim normal] old_exit="<<nav<<" receiving_exit="<<exit
      <<" position_mm="<<post->GetPosition()/mm<<G4endl;
    return true;
  }
  static bool IncidentBefore(const G4StepPoint* pre,const G4StepPoint* post,
                             const G4ThreeVector& dir,G4double distance) {
    auto* touch=dynamic_cast<const G4TouchableHistory*>(pre->GetTouchable());
    auto* pv=pre->GetPhysicalVolume();if(!touch || !pv)return false;
    const auto local=touch->GetHistory()->GetTopTransform().TransformPoint(post->GetPosition()-distance*dir);
    return pv->GetLogicalVolume()->GetSolid()->Inside(local)==kInside;
  }
  static bool HasRindex(const G4Material* m) {
    return m && m->GetMaterialPropertiesTable() && m->GetMaterialPropertiesTable()->GetProperty(kRINDEX);
  }
  static G4String Volume(const G4StepPoint* p) {return p->GetPhysicalVolume()?p->GetPhysicalVolume()->GetName():"none";}
  static bool IsReflection(G4OpBoundaryProcessStatus s) {
    return s==FresnelReflection || s==TotalInternalReflection || s==LambertianReflection ||
      s==LobeReflection || s==SpikeReflection || s==BackScattering || s==CoatedDielectricReflection;
  }
  static bool SupportedSurface(const G4StepPoint* a,const G4StepPoint* b,G4double energy,G4OpticalSurface*& out) {
    auto* av=a->GetPhysicalVolume();auto* bv=b->GetPhysicalVolume();if(!av || !bv)return false;
    G4LogicalSurface* s=G4LogicalBorderSurface::GetSurface(av,bv);
    if(!s) {
      if(bv->GetMotherLogical()==av->GetLogicalVolume()) {
        s=G4LogicalSkinSurface::GetSurface(bv->GetLogicalVolume());
        if(!s)s=G4LogicalSkinSurface::GetSurface(av->GetLogicalVolume());
      } else {
        s=G4LogicalSkinSurface::GetSurface(av->GetLogicalVolume());
        if(!s)s=G4LogicalSkinSurface::GetSurface(bv->GetLogicalVolume());
      }
    }
    if(!s)return false;
    out=dynamic_cast<G4OpticalSurface*>(s->GetSurfaceProperty());
    if(!out || (out->GetName()!="AluminumSurf" && out->GetName()!="protoGelTapeAbsorber") ||
       out->GetType()!=dielectric_metal || out->GetFinish()!=polished || out->GetModel()!=unified)return false;
    const auto* m=out->GetMaterialPropertiesTable();if(!m)return false;
    auto* r=m->GetProperty(kREFLECTIVITY);auto* t=m->GetProperty(kTRANSMITTANCE);auto* e=m->GetProperty(kEFFICIENCY);
    return r && r->Value(energy)>=0. && r->Value(energy)<=1. && (!t || t->Value(energy)==0.) && (!e || e->Value(energy)==0.);
  }
  static bool ForwardInside(const G4StepPoint* p,const G4ThreeVector& dir,G4double distance) {
    auto* pv=p->GetPhysicalVolume();
    auto* touch=dynamic_cast<const G4TouchableHistory*>(p->GetTouchable());
    if(!pv || !touch)return false;
    const auto local=touch->GetHistory()->GetTopTransform().TransformPoint(p->GetPosition()+distance*dir);
    return pv->GetLogicalVolume()->GetSolid()->Inside(local)==kInside;
  }
  G4double tolerance_;
  std::unique_ptr<G4StepPoint> expected_;
  G4ThreeVector lastBoundary_;
  std::unique_ptr<G4Track> proxyTrack_;
  std::unique_ptr<G4Step> proxyStep_;
};

class CBDsimBoundaryPhysics final : public G4VPhysicsConstructor {
 public:
  CBDsimBoundaryPhysics():G4VPhysicsConstructor("CBDsimBoundaryPhysics",0) {}
  void ConstructParticle() override {}
  void ConstructProcess() override {
    auto* pm=G4OpticalPhoton::Definition()->GetProcessManager();
    auto capture=[&](G4ProcessVectorTypeIndex kind) {
      std::vector<G4VProcess*> v;auto* p=pm->GetPostStepProcessVector(kind);
      for(std::size_t i=0;i<p->size();++i)v.push_back((*p)[i]);return v;
    };
    auto oldGPIL=capture(typeGPIL),oldDoIt=capture(typeDoIt);
    G4OpBoundaryProcess* old=nullptr;
    for(auto* p:oldGPIL)if(auto* b=dynamic_cast<G4OpBoundaryProcess*>(p)) {old=b;break;}
    if(!old) {G4Exception("CBDsimBoundaryPhysics","MissingOpBoundary",FatalException,"Expected standard optical boundary process");return;}
    auto* replacement=new CBDsimBoundaryProcess;
    pm->RemoveProcess(old);pm->AddDiscreteProcess(replacement);
    // Preserve the exact effective GPIL/DoIt ordering, including previous default-order ties.
    for(std::size_t i=0;i<oldDoIt.size();++i)
      pm->SetProcessOrdering(oldDoIt[i]==old ? replacement:oldDoIt[i],idxPostStep,static_cast<int>(i)+1);
    auto newGPIL=capture(typeGPIL),newDoIt=capture(typeDoIt);
    auto same=[&](const auto& a,const auto& b){
      if(a.size()!=b.size())return false;
      for(std::size_t i=0;i<a.size();++i)if((a[i]==old?replacement:a[i])!=b[i])return false;return true;
    };
    if(!same(oldGPIL,newGPIL) || !same(oldDoIt,newDoIt))
      G4Exception("CBDsimBoundaryPhysics","BoundaryProcessOrdering",FatalException,"Optical GPIL/DoIt ordering changed");
    G4cout<<"[Boundary installer] GPIL/DoIt ordering preserved; tiny opaque entry recovery Geant4="<<G4VERSION_NUMBER<<G4endl;
    delete old;
  }
};
#endif
