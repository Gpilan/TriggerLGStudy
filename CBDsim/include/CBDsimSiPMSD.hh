#ifndef CBDsimSiPMSD_h
#define CBDsimSiPMSD_h 1

#include "CBDsimSiPMHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include <utility>

class CBDsimSiPMSD : public G4VSensitiveDetector {
public:
  CBDsimSiPMSD(const G4String& name, const G4String& hitsCollectionName, std::pair<int,int> xy);
  virtual ~CBDsimSiPMSD();

  virtual void Initialize(G4HCofThisEvent* HCE);
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);
  virtual void EndOfEvent(G4HCofThisEvent* HCE);

private:
  CBDsimSiPMHitsCollection* fHitCollection;
  G4int fHCID;
  std::pair<int,int> fTowerXY;
  G4int fWavBin;
  G4int fTimeBin;
  G4float fWavlenStart;
  G4float fWavlenEnd;
  G4float fTimeStart;
  G4float fTimeEnd;
  G4float fWavlenStep;
  G4float fTimeStep;

    G4double wavToE(G4double wav) { return h_Planck*c_light/wav; }

  CBDsimInterface::hitRange findWavRange(G4double en);
  CBDsimInterface::hitRange findTimeRange(G4double stepTime);
};

#endif
