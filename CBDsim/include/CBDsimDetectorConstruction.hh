#ifndef CBDsimDetectorConstruction_h
#define CBDsimDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "CBDsimMaterials.hh"
#include "CBDsimSiPMSD.hh"

class fiberSpec;

class CBDsimDetectorConstruction : public G4VUserDetectorConstruction {
public:
  CBDsimDetectorConstruction();
  virtual ~CBDsimDetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();

private:
  void DefineMaterials();
  G4Material* FindMaterial(G4String matName) { return fMaterials->GetMaterial(matName); }
  G4OpticalSurface* FindSurface(G4String surfName) { return fMaterials->GetOpticalSurface(surfName); }

  CBDsimMaterials* fMaterials;

  G4double fTotalX;
  G4double fAirX;
  G4double fFoilX;
  G4double fTowerH;
  G4double fSiPMH;
  G4double fTowerX;
  G4double fSiPMX;
  G4double fFilterT;
  G4double fEcalH;
  G4double fEcalX;
  G4double fEcalPMTH;
  G4double fEcalPMTX;
  G4int fTowerN;
  G4int fEcalN;
  G4int fTotalN;
  const fiberSpec* fFiberSpec;

  G4VisAttributes* fVisAttrOrange;
  G4VisAttributes* fVisAttrBlue;
  G4VisAttributes* fVisAttrTransparent;
  G4VisAttributes* fVisAttrGreen;

  std::vector<CBDsimSiPMSD*> fSiPMSD;
  G4LogicalVolume* fWaferlogical;
};

class fiberSpec {
public:
  fiberSpec(G4double rMin_cladC, G4double rMin_cladS, G4double rMax, G4double towerH,
    G4Material* matPS, G4Material* matPMMA, G4Material* matFluoPoly, G4VisAttributes* visAttrOrange, G4VisAttributes* visAttrBlue);
  virtual ~fiberSpec();

  G4double rMin_cladC;
  G4double rMin_cladS;
  G4double rMax;
  G4double towerH;

  G4Material* matPS;
  G4Material* matPMMA;
  G4Material* matFluoPoly;

  G4VisAttributes* visAttrOrange;
  G4VisAttributes* visAttrBlue;
};

#endif
