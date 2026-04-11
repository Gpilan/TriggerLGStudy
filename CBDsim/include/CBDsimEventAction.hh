#ifndef CBDsimEventAction_h
#define CBDsimEventAction_h 1

#include "CBDsimInterface.h"
#include "CBDsimSiPMHit.hh"
#include "CBDsimSiPMSD.hh"

#include "G4UserEventAction.hh"

class CBDsimRunAction;

class CBDsimEventAction : public G4UserEventAction {
public:
  CBDsimEventAction();
  virtual ~CBDsimEventAction();

  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

  void fillEdeps(CBDsimInterface::CBDsimEdepData& edepData);
  void fillOpticalPhoton(CBDsimInterface::CBDsimPhoton& photon, G4int oPnumber);
  void fillPhysics(CBDsimInterface::CBDsimPhysicalevent& physical,G4double x,G4double y,G4double z,G4double energy,G4String PartName,G4String PhysicName);
  //void fillVexs(CBDsimInterface::CBDsimVertax& vexData);

private:
  void clear();
  void fillHits(CBDsimSiPMHit* hit);
  void queue();

  CBDsimInterface::CBDsimEventData* fEventData;
  std::map<int,CBDsimInterface::CBDsimTowerData> fTowerMap;
  std::map<int,CBDsimInterface::CBDsimEdepData> fEdepMap;
  std::vector<CBDsimInterface::CBDsimPhysicalevent> fPhysicalMap;
  std::vector<CBDsimInterface::CBDsimPhoton> fPhotonVector;
  // std::map<int,CBDsimInterface::CBDsimTotalPhysical> fPhysicalMap;
};

#endif
