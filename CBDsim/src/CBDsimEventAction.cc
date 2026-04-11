#include "CBDsimRunAction.hh"
#include "CBDsimEventAction.hh"
#include "CBDsimPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4VHitsCollection.hh"
#include "G4AutoLock.hh"
#include "G4Threading.hh"
#include <algorithm>

namespace {
  G4Mutex CBDsimEventActionMutex = G4MUTEX_INITIALIZER;
  G4Condition CBDsimEventActionCV = G4CONDITION_INITIALIZER;

  void FlattenTimeStruct(CBDsimInterface::CBDsimSiPMData& d) {
    d.timeBinEdgeLow.clear();
    d.timeBinEdgeHigh.clear();
    d.timeBinCounts.clear();
    for (const auto& kv : d.timeStruct) {
      d.timeBinEdgeLow.push_back(kv.first.first);
      d.timeBinEdgeHigh.push_back(kv.first.second);
      d.timeBinCounts.push_back(kv.second);
    }
  }

  void MergeSiPMTimeInto(std::vector<int>& accCounts, std::vector<float>& accLow,
                         std::vector<float>& accHigh, const CBDsimInterface::CBDsimSiPMData& s) {
    if (s.timeBinCounts.empty()) return;
    if (accCounts.empty()) {
      accCounts = s.timeBinCounts;
      accLow = s.timeBinEdgeLow;
      accHigh = s.timeBinEdgeHigh;
      return;
    }
    const size_t n = std::min(accCounts.size(), s.timeBinCounts.size());
    for (size_t i = 0; i < n; ++i) accCounts[i] += s.timeBinCounts[i];
  }

  void FillTriggerSummaries(CBDsimInterface::CBDsimEventData& ev) {
    ev.siPMPhotonSumTrig0 = ev.siPMPhotonSumTrig1 = 0;
    ev.nSiPMChannelsTrig0 = ev.nSiPMChannelsTrig1 = 0;
    ev.hasSiPMTrig0 = ev.hasSiPMTrig1 = 0;
    ev.timeMergedCountsTrig0.clear();
    ev.timeMergedEdgeLowTrig0.clear();
    ev.timeMergedEdgeHighTrig0.clear();
    ev.timeMergedCountsTrig1.clear();
    ev.timeMergedEdgeLowTrig1.clear();
    ev.timeMergedEdgeHighTrig1.clear();

    auto accumulateTower = [&](const CBDsimInterface::CBDsimTowerData& tw, int trig) {
      int* sumPtr = (trig == 0) ? &ev.siPMPhotonSumTrig0 : &ev.siPMPhotonSumTrig1;
      int* nChPtr = (trig == 0) ? &ev.nSiPMChannelsTrig0 : &ev.nSiPMChannelsTrig1;
      int* hasPtr = (trig == 0) ? &ev.hasSiPMTrig0 : &ev.hasSiPMTrig1;
      std::vector<int>* mCnt = (trig == 0) ? &ev.timeMergedCountsTrig0 : &ev.timeMergedCountsTrig1;
      std::vector<float>* mLo = (trig == 0) ? &ev.timeMergedEdgeLowTrig0 : &ev.timeMergedEdgeLowTrig1;
      std::vector<float>* mHi = (trig == 0) ? &ev.timeMergedEdgeHighTrig0 : &ev.timeMergedEdgeHighTrig1;

      auto accumulateSide = [&](const std::vector<CBDsimInterface::CBDsimSiPMData>& side) {
        for (const auto& s : side) {
          *sumPtr += s.count;
          (*nChPtr)++;
          if (s.count > 0) *hasPtr = 1;
          MergeSiPMTimeInto(*mCnt, *mLo, *mHi, s);
        }
      };
      accumulateSide(tw.SiPMs);
    };
    accumulateTower(ev.towerT1, 0);
    accumulateTower(ev.towerT2, 1);
  }
}

CBDsimEventAction::CBDsimEventAction() {
  fEventData = nullptr;
}

CBDsimEventAction::~CBDsimEventAction() {
  if (fEventData) {
    delete fEventData;
    fEventData = nullptr;
  }
}

void CBDsimEventAction::BeginOfEventAction(const G4Event* evt) {
  clear();
  int evtNo = evt->GetEventID();

  
  fEventData = new CBDsimInterface::CBDsimEventData();


}

void CBDsimEventAction::clear() {
  fTowerMap.clear();
  fEdepMap.clear();
//memset(fPhysical,0,sizeof(CBDsimInterface::CBDsimTotalPhysical));
  fPhysicalMap.clear();
  fPhotonVector.clear();
  
  // fEventData 메모리 누수 방지
  if (fEventData) {
    delete fEventData;
    fEventData = nullptr;
  }
}

void CBDsimEventAction::EndOfEventAction(const G4Event* evt) {
    G4HCofThisEvent* hce = evt->GetHCofThisEvent();
  if (!hce) {
    std::cout << "No hit collection!" << std::endl;
    return;
  }
  // GetHC(i): 이 이벤트에 쌓인 HC의 순서 인덱스(0 … n-1).
  // SiPM은 G4THitsCollection<CBDsimSiPMHit> 이지만, dynamic_cast 가 RTTI/DSO 이슈로
  // nullptr 이 되는 환경이 있어, 컬렉션 이름 확인 후 static_cast 사용.
  const int nhc = hce->GetNumberOfCollections();
  const G4int evtId = evt->GetEventID();

  // 이벤트 0 한 번: HC 이름/ID·이후 fTowerMap 키 확인용 (디버그)
  if (evtId == 0) {
    G4cout << "\n=== [CBDsim HC debug] EventID=0, GetNumberOfCollections=" << nhc << " ===" << G4endl;
    for (int i = 0; i < nhc; ++i) {
      G4VHitsCollection* vhc = hce->GetHC(i);
      if (!vhc) {
        G4cout << "  HC[" << i << "] (null)" << G4endl;
        continue;
      }
      G4cout << "  HC[" << i << "] GetName=\"" << vhc->GetName() << "\"" << G4endl;
    }
  }

  for (int i = 0; i < nhc; ++i) {
    G4VHitsCollection* vhc = hce->GetHC(i);
    if (!vhc) continue;
    const G4String& nm = vhc->GetName();
    if (nm != "SiPMSDBC") continue;
    auto* sipmHC = static_cast<CBDsimSiPMHitsCollection*>(vhc);
    const int sipms = sipmHC->entries();
    if (evtId == 0 && sipms > 0) {
      G4cout << "  SiPM HC matched, entries=" << sipms << " first hit towernum=" << (*sipmHC)[0]->GetTowernum()
             << " SiPMnum=" << (*sipmHC)[0]->GetSiPMnum() << G4endl;
    }
    for (int iHC = 0; iHC < sipms; ++iHC) {
      fillHits((*sipmHC)[iHC]);
    }
  }

  if (evtId == 0) {
    G4cout << "=== [CBDsim HC debug] fTowerMap.size=" << fTowerMap.size() << " ===" << G4endl;
    for (const auto& kv : fTowerMap) {
      G4cout << "  tower key=" << kv.first << "  SiPMs=" << kv.second.SiPMs.size() << G4endl;
    }
    G4cout << "=== [CBDsim HC debug] end ===\n" << G4endl;
  }

  fEventData->event_number = CBDsimPrimaryGeneratorAction::sIdxEvt;
  {
    G4double ekin = 0., vx = 0., vy = 0., vz = 0., dx = 0., dy = 0., dz = 0.;
    CBDsimPrimaryGeneratorAction::GetLastPrimaryKinematics(ekin, vx, vy, vz, dx, dy, dz);
    fEventData->primaryEkin = static_cast<float>(ekin);
    fEventData->primaryVx = static_cast<float>(vx);
    fEventData->primaryVy = static_cast<float>(vy);
    fEventData->primaryVz = static_cast<float>(vz);
    fEventData->primaryDirX = static_cast<float>(dx);
    fEventData->primaryDirY = static_cast<float>(dy);
    fEventData->primaryDirZ = static_cast<float>(dz);
  }

  {
    fEventData->towerT1 = CBDsimInterface::CBDsimTowerData();
    fEventData->towerT1.triggerNum = 0;
    fEventData->towerT2 = CBDsimInterface::CBDsimTowerData();
    fEventData->towerT2.triggerNum = 1;
    if (fTowerMap.count(0)) fEventData->towerT1 = fTowerMap[0];
    if (fTowerMap.count(1)) fEventData->towerT2 = fTowerMap[1];
    fEventData->towerT1.triggerNum = 0;
    fEventData->towerT2.triggerNum = 1;
  }

  for (const auto& edepMap : fEdepMap) {
    fEventData->Edeps.push_back(edepMap.second);
  }
  for (const auto& physicalMap : fPhysicalMap){
    fEventData->totPhysicals.push_back(physicalMap);
  }
  for (const auto& photonVector : fPhotonVector){
    fEventData->opticalPhotons.push_back(photonVector);
  }

  FillTriggerSummaries(*fEventData);

  queue();
}
void CBDsimEventAction::fillOpticalPhoton(CBDsimInterface::CBDsimPhoton& Photondata,G4int oPnumber) {
  Photondata.opticalPhotonNumber = oPnumber;
  fPhotonVector.push_back(Photondata);
}

void CBDsimEventAction::fillHits(CBDsimSiPMHit* hit) {
  CBDsimInterface::CBDsimSiPMData sipmData;
  sipmData.count = hit->GetPhotonCount();
  sipmData.SiPMnum = hit->GetSiPMnum();
  sipmData.x = hit->GetX();
  sipmData.y = hit->GetY();
  sipmData.timeStruct = hit->GetTimeStruct();
  sipmData.wavlenSpectrum = hit->GetWavlenSpectrum();
  FlattenTimeStruct(sipmData);
  int towernum = hit->GetTowernum();
  auto towerIter = fTowerMap.find(towernum);

  if (towerIter == fTowerMap.end()) {
    CBDsimInterface::CBDsimTowerData towerData;
    towerData.triggerNum = hit->GetTowernum();
    towerData.numx = hit->GetTowerX();
    towerData.numy = hit->GetTowerY();
    towerData.SiPMs.push_back(sipmData);
    fTowerMap.insert(std::make_pair(towernum, towerData));
  } else {
    towerIter->second.SiPMs.push_back(sipmData);
  }
}

void CBDsimEventAction::fillEdeps(CBDsimInterface::CBDsimEdepData& edepData) {
  auto towerIter = fEdepMap.find(edepData.triggerNum);

  if ( towerIter==fEdepMap.end() ) {
    fEdepMap.insert(std::make_pair(edepData.triggerNum,edepData));
  } else {
    towerIter->second.Edep += edepData.Edep;
  }
}

void CBDsimEventAction::fillPhysics(CBDsimInterface::CBDsimPhysicalevent& physical, G4double x,G4double y, G4double z, G4double energy,G4String PartName,G4String PhysicName){
  physical.x=x;
  physical.y=y;
  physical.z=z;
  if (PartName=="gamma"){
    physical.particleID=0;
    if (PhysicName=="msc")physical.physicalID=4;
    else if (PhysicName=="eIoni")physical.physicalID=5;
    else if (PhysicName=="eBrem")physical.physicalID=6;
    else if (PhysicName=="annihil")physical.physicalID=7;

    //else G4cout<<"kkkkkkkkkjkjkjkjkjk"<<G4endl;


  }
  else if (PartName=="e-"){
    physical.particleID=1;
    if (PhysicName=="phot")physical.physicalID=0;
    else if (PhysicName=="compt")physical.physicalID=1;
    else if (PhysicName=="conv")physical.physicalID=2;
    else if (PhysicName=="Rayl")physical.physicalID=3;
    else if (PhysicName=="eIoni")physical.physicalID=5;
    else G4cout<<"kkkkkihugbffbdjlfkkkkjkjkjkjkjk"<<G4endl;

  }
  else if (PartName=="e+"){
    physical.particleID=2;
    if (PhysicName=="phot")physical.physicalID=0;
    else if (PhysicName=="compt")physical.physicalID=1;
    else if (PhysicName=="conv")physical.physicalID=2;
    else if (PhysicName=="Rayl")physical.physicalID=3;
    else if (PhysicName=="eIoni")physical.physicalID=5;
    else G4cout<<"kkkkqqqqqqqqqqqqkkkkkjkjkjkjkjk"<<G4endl;

  }


  physical.energy=energy;


  fPhysicalMap.push_back(physical);


}
//void CBDsimEventAction::fillVexs(CBDsimInterface::CBDsimVertax& vexData){



void CBDsimEventAction::queue() {
  while ( CBDsimRunAction::sNumEvt != CBDsimPrimaryGeneratorAction::sIdxEvt ) {
    G4AutoLock lock(&CBDsimEventActionMutex);
    std::cout << "thread = " << G4Threading::G4GetThreadId() << " | sNumEvt = " << CBDsimRunAction::sNumEvt << " | sIdxEvt = " << CBDsimPrimaryGeneratorAction::sIdxEvt << std::endl;
    if ( CBDsimRunAction::sNumEvt == CBDsimPrimaryGeneratorAction::sIdxEvt ) break;
    G4CONDITIONWAIT(&CBDsimEventActionCV, &lock);
  }
  G4AutoLock lock(&CBDsimEventActionMutex);
  if (CBDsimRunAction::sRootIO)
    CBDsimRunAction::sRootIO->fill(fEventData);
  CBDsimRunAction::sNumEvt++;
  G4CONDITIONBROADCAST(&CBDsimEventActionCV);
}
