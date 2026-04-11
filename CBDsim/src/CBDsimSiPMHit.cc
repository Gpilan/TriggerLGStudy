#include "CBDsimSiPMHit.hh"

G4ThreadLocal G4Allocator<CBDsimSiPMHit>* CBDsimSiPMHitAllocator = 0;

CBDsimSiPMHit::CBDsimSiPMHit(G4int wavBin, G4int timeBin)
: G4VHit(),
fPhotons(0),
fSiPMnum(-1),
fTowernum(-1),
fX(-1),
fY(-1),
fTowerX(-1),
fTowerY(-1),
fWavBin(wavBin),
fTimeBin(timeBin)
{}

CBDsimSiPMHit::~CBDsimSiPMHit() {}

CBDsimSiPMHit::CBDsimSiPMHit(const CBDsimSiPMHit &right)
: G4VHit() {
  fSiPMnum = right.fSiPMnum;
  fPhotons = right.fPhotons;
  fTowernum = right.fTowernum;
  fX = right.fX;
  fY = right.fY;
  fTowerX = right.fTowerX;
  fTowerY = right.fTowerY;
  fTimeStruct = right.fTimeStruct;
  fWavlenSpectrum=right.fWavlenSpectrum;

}

const CBDsimSiPMHit& CBDsimSiPMHit::operator=(const CBDsimSiPMHit &right) {
  fSiPMnum = right.fSiPMnum;
  fPhotons = right.fPhotons;
  fTowernum = right.fTowernum;
  fX = right.fX;
  fY = right.fY;
  fTowerX = right.fTowerX;
  fTowerY = right.fTowerY;
  fTimeStruct = right.fTimeStruct;
  fWavlenSpectrum=right.fWavlenSpectrum;

  return *this;
}

G4bool CBDsimSiPMHit::operator==(const CBDsimSiPMHit &right) const {
  return (fSiPMnum==right.fSiPMnum && fTowernum==right.fTowernum);
}

void CBDsimSiPMHit::CountWavlenSpectrum(CBDsimInterface::hitRange range) {
  auto it = fWavlenSpectrum.find(range);
  if (it==fWavlenSpectrum.end()) fWavlenSpectrum.insert(std::make_pair(range,1));
  else it->second++;
}

void  CBDsimSiPMHit::CountTimeStruct(CBDsimInterface::hitRange range){
  auto it=fTimeStruct.find(range);
  if (it==fTimeStruct.end())  fTimeStruct.insert(std::make_pair(range,1));
  else it->second++;

}
