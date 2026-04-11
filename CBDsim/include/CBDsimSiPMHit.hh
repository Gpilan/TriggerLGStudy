#ifndef CBDsimSiPMHit_h
#define CBDsimSiPMHit_h 1

#include "CBDsimInterface.h"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class CBDsimSiPMHit : public G4VHit {
public:
  CBDsimSiPMHit(G4int wavBin, G4int timeBin);
  CBDsimSiPMHit(const CBDsimSiPMHit& right);
  virtual ~CBDsimSiPMHit();

  const CBDsimSiPMHit& operator=(const CBDsimSiPMHit &right);
  G4bool operator==(const CBDsimSiPMHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void* aHit);

  G4int GetSiPMnum() const { return fSiPMnum; }
  void SetSiPMnum(G4int n) { fSiPMnum = n; }
  G4int GetTowernum() const { return fTowernum; }
  void SetTowernum(G4int n) { fTowernum = n; }

  G4int GetX() const { return fX; }
  void SetX(G4int n) { fX = n; }
  G4int GetY() const { return fY; }
  void SetY(G4int n) { fY = n; }

  G4int GetTowerX() const { return fTowerX; }
  void SetTowerX(G4int n) { fTowerX = n; }
  G4int GetTowerY() const { return fTowerY; }
  void SetTowerY(G4int n) { fTowerY = n; }

  void photonCount() { fPhotons++; }
  G4int GetPhotonCount() const { return fPhotons; }

  void CountWavlenSpectrum(CBDsimInterface::hitRange range);
  CBDsimInterface::CBDsimWavlenSpectrum GetWavlenSpectrum() const { return fWavlenSpectrum; }

  void CountTimeStruct(CBDsimInterface::hitRange range);
  CBDsimInterface::CBDsimTimeStruct GetTimeStruct() const { return fTimeStruct; }
private:
  G4int fPhotons;
  G4int fSiPMnum;
  G4int fTowernum;
  G4int fX;
  G4int fY;
  G4int fTowerX;
  G4int fTowerY;
  CBDsimInterface::CBDsimWavlenSpectrum fWavlenSpectrum;
  CBDsimInterface::CBDsimTimeStruct fTimeStruct;

  G4int fWavBin;
  G4int fTimeBin;
};

typedef G4THitsCollection<CBDsimSiPMHit> CBDsimSiPMHitsCollection;
extern G4ThreadLocal G4Allocator<CBDsimSiPMHit>* CBDsimSiPMHitAllocator;

inline void* CBDsimSiPMHit::operator new(size_t) {
  if (!CBDsimSiPMHitAllocator) CBDsimSiPMHitAllocator = new G4Allocator<CBDsimSiPMHit>;
  return (void*)CBDsimSiPMHitAllocator->MallocSingle();
}

inline void CBDsimSiPMHit::operator delete(void*aHit) {
  CBDsimSiPMHitAllocator->FreeSingle((CBDsimSiPMHit*) aHit);
}

#endif
