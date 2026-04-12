#ifndef CBDsimMaterials_h
#define CBDsimMaterials_h 1

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"

class CBDsimMaterials {
public:
  virtual ~CBDsimMaterials();
  static CBDsimMaterials* GetInstance();
  G4Material* GetMaterial(const G4String);
  G4OpticalSurface* GetOpticalSurface(const G4String);

private:
  CBDsimMaterials();
  void CreateMaterials();

  static CBDsimMaterials* fInstance;
  G4NistManager* fNistMan;
  G4Material* fVacuum;
  G4Material* fAir;
  G4Material* fFluoPoly;
  G4Material* fPMMA;
  /** PMMA-like bulk, RINDEX matched to scint (legacy name `ProtoLG_MatchScint`; current proto has no LG volume). */
  G4Material* fProtoLGMatchScint;
  G4Material* fPS;
  G4Material* fCu;
  G4Material* fSi;
  G4Material* fGlass;
  G4Material* fGelatin;
  G4Material* fLYSO;
  G4Material* fLSO;
  G4Material* fPWO;
  G4Material* fAl;
  G4OpticalSurface* fSiPMSurf;
  G4OpticalSurface* fFilterSurf;
  G4OpticalSurface* fAluminumSurf;
  G4OpticalSurface* fAirSurf;
};

#endif
