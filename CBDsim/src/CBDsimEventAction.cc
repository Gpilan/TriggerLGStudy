#include "CBDsimRunAction.hh"
#include "CBDsimEventAction.hh"
#include "CBDsimPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4AutoLock.hh"
#include "G4Threading.hh"


namespace {
  G4Mutex CBDsimEventActionMutex = G4MUTEX_INITIALIZER;
  G4Condition CBDsimEventActionCV = G4CONDITION_INITIALIZER;
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

  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  fSiPMCID.push_back(SDman->GetCollectionID("SiPMSDBC"));
  fSiPMCID.push_back(SDman->GetCollectionID("SiPMSDFC"));


  fEventData = new CBDsimInterface::CBDsimEventData();


}

void CBDsimEventAction::clear() {
  fSiPMCID.clear();
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
  int totSDnum = hce->GetNumberOfCollections(); // which is 1, hence no need to run lopp on SD in this case
  CBDsimSiPMHitsCollection* sipmHC = 0;
  for (int iSD=0;iSD<totSDnum;iSD++){
    if (hce) {
      if (fSiPMCID[iSD]>=0) sipmHC = (CBDsimSiPMHitsCollection*)(hce->GetHC(fSiPMCID[iSD]));
    }
    if (sipmHC) {
      int sipms = sipmHC->entries();
      for (int iHC = 0; iHC < sipms; iHC++) {
        fillHits((*sipmHC)[iHC]);
      }
    }
  }

  fEventData->event_number = CBDsimPrimaryGeneratorAction::sIdxEvt;

  for (const auto& towerMap : fTowerMap) {
    fEventData->towers.push_back(towerMap.second);
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
  // fEventData->opticalPhotons.push_back();
  //fEventData->totPhysicals.push_back(fPhysical);


  queue();
}
void CBDsimEventAction::fillOpticalPhoton(CBDsimInterface::CBDsimPhoton& Photondata,G4int oPnumber) {
  Photondata.opticalPhotonNumber = oPnumber;
  fPhotonVector.push_back(Photondata);
  G4cout<<"the photon number is "<<Photondata.opticalPhotonNumber<<G4endl;
}

void CBDsimEventAction::fillHits(CBDsimSiPMHit* hit) {
  std::cout << "nPhotons = " << hit->GetPhotonCount() << ", XY = " << hit->GetX() << hit->GetY() << ", tower XY = " << hit->GetTowerX() << hit->GetTowerY() <<"Front="<<hit->GetIsFront()<< std::endl;
  CBDsimInterface::CBDsimSiPMData sipmData;
  CBDsimInterface::CBDsimSiPMData sipmDatafront;
  if(hit->GetIsFront()==0){
    sipmData.count = hit->GetPhotonCount();
    sipmData.SiPMnum = hit->GetSiPMnum();
    sipmData.x = hit->GetX();
    sipmData.y = hit->GetY();
    sipmData.timeStruct = hit->GetTimeStruct();
    sipmData.wavlenSpectrum = hit->GetWavlenSpectrum();
  }
  else if(hit->GetIsFront()==1){
    sipmDatafront.count = hit->GetPhotonCount();
    sipmDatafront.SiPMnum = hit->GetSiPMnum();
    sipmDatafront.x = hit->GetX();
    sipmDatafront.y = hit->GetY();
    sipmDatafront.timeStruct = hit->GetTimeStruct();
    sipmDatafront.wavlenSpectrum = hit->GetWavlenSpectrum();
  }
  int towernum = hit->GetTowernum();
  auto towerIter = fTowerMap.find(towernum);


  if ( towerIter==fTowerMap.end() ) {
    CBDsimInterface::CBDsimTowerData towerData;
    towerData.towerNum = hit->GetTowernum();
    towerData.numx = hit->GetTowerX();
    towerData.numy = hit->GetTowerY();


    if (hit->GetIsFront()==1) {

      towerData.SiPMFronts.push_back(sipmDatafront);
    }
    else if (hit->GetIsFront()==0){


      towerData.SiPMs.push_back(sipmData);
    }

    fTowerMap.insert(std::make_pair(towernum,towerData));
  } else {
    if (hit->GetIsFront()==1) {

      towerIter->second.SiPMFronts.push_back(sipmDatafront);
    }
    else if (hit->GetIsFront()==0){


      towerIter->second.SiPMs.push_back(sipmData);
    }

  }
}

void CBDsimEventAction::fillEdeps(CBDsimInterface::CBDsimEdepData& edepData) {
  auto towerIter = fEdepMap.find(edepData.towerNum);

  if ( towerIter==fEdepMap.end() ) {
    fEdepMap.insert(std::make_pair(edepData.towerNum,edepData));
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
