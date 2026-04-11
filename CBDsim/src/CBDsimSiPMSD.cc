#include "CBDsimSiPMSD.hh"
#include "CBDsimSiPMHit.hh"

#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VPhysicalVolume.hh"

#include <atomic>

using namespace std;

namespace {
constexpr int kSiPMSD_DbgMaxPrint = 25;
}

namespace {
/** Legacy: depth-2 = SiPM cell copy, depth-3 = tower parent. Proto: depth-1 = protoSipmEnvPhys (copy 0=T1,1=T2), single SiPM pixel. */
void SiPMAndTowerFromTouchable(G4Step* step, G4int& SiPMnum, G4int& towernum) {
  auto* touch = step->GetPostStepPoint()->GetTouchable();
  G4VPhysicalVolume* vol1 = touch->GetVolume(1);
  G4VPhysicalVolume* vol2 = touch->GetVolume(2);
  G4VPhysicalVolume* vol3 = touch->GetVolume(3);
  if (vol1 && vol1->GetName() == "protoSipmEnvPhys") {
    towernum = vol1->GetCopyNo();
    SiPMnum = 0;
    return;
  }
  SiPMnum = vol2 ? vol2->GetCopyNo() : 0;
  towernum = vol3 ? vol3->GetCopyNo() : 0;
}
}  // namespace

CBDsimSiPMSD::CBDsimSiPMSD(const G4String& name, const G4String& hitsCollectionName, std::pair<int,int> xy)
: G4VSensitiveDetector(name), fHitCollection(0), fHCID(-1), fTowerXY(xy),fWavBin(600),fTimeBin(24000),fWavlenStart(900.),fWavlenEnd(300.),fTimeStart(0.),fTimeEnd(240.)
{
  collectionName.insert(hitsCollectionName);
  fWavlenStep = (fWavlenStart-fWavlenEnd)/(float)fWavBin;
  fTimeStep = (fTimeEnd-fTimeStart)/(float)fTimeBin;
}

CBDsimSiPMSD::~CBDsimSiPMSD() {}

void CBDsimSiPMSD::Initialize(G4HCofThisEvent* hce) {
  fHitCollection = new CBDsimSiPMHitsCollection(SensitiveDetectorName,collectionName[0]);
  if (fHCID<0) { fHCID = GetCollectionID(0); }
  hce->AddHitsCollection(fHCID,fHitCollection);
}

G4bool CBDsimSiPMSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
  if (step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return false;
  G4int SiPMnum = 0;
  G4int towernum = 0;
  SiPMAndTowerFromTouchable(step, SiPMnum, towernum);

  // 진단용: 광학 광자가 웨이퍼 SD로 들어올 때만 호출됨. 확인 후 제거하거나 kSiPMSD_DbgMaxPrint 를 0 으로.
  {
    static std::atomic<int> sDbgProcessHitsCount{0};
    const int n = ++sDbgProcessHitsCount;
    if (n <= kSiPMSD_DbgMaxPrint) {
      G4String postName = "(no post PV)";
      if (auto* pv = step->GetPostStepPoint()->GetPhysicalVolume()) postName = pv->GetName();
      G4int evId = -1;
      if (auto* em = G4EventManager::GetEventManager()) {
        if (auto* ev = em->GetConstCurrentEvent()) evId = ev->GetEventID();
      }
      G4cout << "[SiPM SD] ProcessHits #" << n << " ev=" << evId << " tower=" << towernum
             << " SiPMnum=" << SiPMnum << " postPV=" << postName << G4endl;
    }
  }

  G4double hitTime  = step->GetPostStepPoint()->GetGlobalTime();
  G4double energy = step->GetTrack()->GetTotalEnergy();
  G4int towerX = fTowerXY.first;
  G4int towerY = fTowerXY.second;
  G4int sipmX = SiPMnum/towerY;
  G4int sipmY = SiPMnum%towerY;
  G4int nofHits = fHitCollection->entries();
  CBDsimSiPMHit* hit = NULL;


  for (G4int i = 0; i < nofHits; i++) {
    if ( ((*fHitCollection)[i]->GetSiPMnum()==SiPMnum) && ((*fHitCollection)[i]->GetTowernum()==towernum) ) {
      hit = (*fHitCollection)[i];
      break;
    }
  }
  if (hit==NULL) {

    hit = new CBDsimSiPMHit(fWavBin,fTimeBin);
    hit->SetSiPMnum(SiPMnum);
    hit->SetTowernum(towernum);
    hit->SetX(sipmX);
    hit->SetY(sipmY);
    hit->SetTowerX(towerX);
    hit->SetTowerY(towerY);

    fHitCollection->insert(hit);
  }

  hit->photonCount();

  CBDsimInterface::hitRange wavRange = findWavRange(energy);
  hit->CountWavlenSpectrum(wavRange);


  CBDsimInterface::hitRange timeRange = findTimeRange(hitTime);
  hit->CountTimeStruct(timeRange);


  return true;
}

void CBDsimSiPMSD::EndOfEvent(G4HCofThisEvent*) {
  if ( verboseLevel>1 ) {
    G4int nofHits = fHitCollection->entries();
    G4cout
    << G4endl
    << "-------->Hits Collection: in this event they are " << nofHits
    << " hits in the tracker chambers: " << G4endl;
    for ( G4int i=0; i<nofHits; i++ ) (*fHitCollection)[i]->Print();
  }
}

CBDsimInterface::hitRange CBDsimSiPMSD::findWavRange(G4double en) {
  int i = 0;
  for ( ; i < fWavBin+1; i++) {
    if ( en < wavToE( (fWavlenStart - (float)i*fWavlenStep)*nm ) ) break;
  }

  if (i==0) return std::make_pair(fWavlenStart,99999.);
  else if (i==fWavBin+1) return std::make_pair(0.,fWavlenEnd);

  return std::make_pair( fWavlenStart-(float)i*fWavlenStep, fWavlenStart-(float)(i-1)*fWavlenStep );
}

CBDsimInterface::hitRange CBDsimSiPMSD::findTimeRange(G4double stepTime) {
  int i = 0;
  for ( ; i < fTimeBin+1; i++) {
    if ( stepTime < ( (fTimeStart + (float)i*fTimeStep)*ns ) ) break;
  }

  if (i==0) return std::make_pair(0.,fTimeStart);
  else if (i==fTimeBin+1) return std::make_pair(fTimeEnd,99999.);

  return std::make_pair( fTimeStart+(float)(i-1)*fTimeStep, fTimeStart+(float)i*fTimeStep );
}
