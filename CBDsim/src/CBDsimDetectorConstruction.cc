#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4NistManager.hh"
#include "G4PVParameterised.hh"
//#include "CBDsimCellParam.hh"

#include "CBDsimDetectorConstruction.hh"
#include "CBDsimFiberParam.hh"

CBDsimDetectorConstruction::CBDsimDetectorConstruction()
: fMaterials(NULL) {
  DefineMaterials();

  fVisAttrOrange = new G4VisAttributes(G4Colour(1.0,0.5,0.,1.0));
  fVisAttrOrange->SetVisibility(true);
  fVisAttrBlue = new G4VisAttributes(G4Colour(0.,0.,1.0,1.0));
  fVisAttrBlue->SetVisibility(true);

  fFiberSpec = new fiberSpec(
    0.49*mm,  // rMin_cladC
    0.485*mm, // rMin_cladS
    0.5*mm,   // rMax
    2.*m,     // towerH
    FindMaterial("Polystyrene"),
    FindMaterial("PMMA"),
    FindMaterial("FluorinatedPolymer"),
    fVisAttrOrange,
    fVisAttrBlue
  );


  fTowerH = 5.5*cm;
  fSiPMH = 0.3*mm;
  fTowerX = 5.*mm;
  fFoilX = 0.016*mm;
  fAirX= 0.01*mm;
  fSiPMX = 5.*mm;
  fFilterT = 0.01*mm;
  fTotalN=1;
  fTowerN =1;
  fTotalX= fTowerN*5.052*mm;



  fVisAttrTransparent = new G4VisAttributes();
  fVisAttrTransparent->SetVisibility(false);
  fVisAttrGreen = new G4VisAttributes(G4Colour(0.3,0.7,0.3));
  fVisAttrGreen->SetVisibility(true);
}

CBDsimDetectorConstruction::~CBDsimDetectorConstruction() {
  if (fMaterials) delete fMaterials;
  if (fFiberSpec) delete fFiberSpec;
  if (fWaferlogical) delete fWaferlogical;
}

void CBDsimDetectorConstruction::DefineMaterials() {
  fMaterials = CBDsimMaterials::GetInstance();
}

G4VPhysicalVolume* CBDsimDetectorConstruction::Construct() {
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();
  G4LogicalBorderSurface::CleanSurfaceTable();

  G4VSolid* worldSolid = new G4Box("worldBox",10./2.*m,10./2.*m,10./2.*m);
  G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid,FindMaterial("G4_Galactic"),"worldLogical");
  G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,false,0);
  worldLogical->SetVisAttributes(fVisAttrTransparent);

  G4VSolid* towertotalEnvSolid = new G4Box("towertotalEnvSolid",fTotalX/2.,fTotalX/2.,(fTowerH+fSiPMH)/2.);
  G4RotationMatrix* towerRot = new G4RotationMatrix(90.*deg,0.*deg,0.);
  G4LogicalVolume* towertotalEnvLogical = new G4LogicalVolume(towertotalEnvSolid,FindMaterial("G4_AIR"),"towertotalEnvLogical");
  //여기다 하나씩 배치할 것

  new G4PVPlacement(towerRot,G4ThreeVector(0.5*m,0,0),towertotalEnvLogical,"towertotalEnvPhysical",worldLogical,false,0);

  G4VSolid* towerEnvSolid = new G4Box("towerEnvSolid",fFoilX+fAirX+fTowerX/2.,fFoilX+fAirX+fTowerX/2.,(fTowerH+fSiPMH)/2.);
  G4LogicalVolume* towerEnvLogical = new G4LogicalVolume(towerEnvSolid,FindMaterial("G4_AIR"),"towerEnvLogical");

  G4VSolid* towerSolid = new G4Box("towerSolid",fTowerX/2.,fTowerX/2.,fTowerH/2.);
  G4LogicalVolume* towerLogical = new G4LogicalVolume(towerSolid,FindMaterial("LYSO"),"towerLogical");
  G4VPhysicalVolume* towerPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),towerLogical,"towerPhysical",towerEnvLogical,false,0);

  towerLogical->SetVisAttributes(fVisAttrBlue);

  G4VSolid* foilSolid1 = new G4Box("foilSolid1",fFoilX/2.,fTowerX/2.+fAirX+fFoilX/2.,(fTowerH+fSiPMH)/2.);
  G4LogicalVolume* foilLogical1 = new G4LogicalVolume(foilSolid1,FindMaterial("Aluminum"),"foilLogical1");

  G4VSolid* foilSolid2 = new G4Box("foilSolid2",fTowerX/2.+fAirX+fFoilX/2.,fFoilX/2.,(fTowerH+fSiPMH)/2.);
  G4LogicalVolume* foilLogical2 = new G4LogicalVolume(foilSolid2,FindMaterial("Aluminum"),"foilLogical2");

  G4VPhysicalVolume* foilPhysical1 = new G4PVPlacement(0,G4ThreeVector(fTowerX/2.+fFoilX/2.+fAirX,fFoilX/2.,0.),foilLogical1,"foilPhysical1",towerEnvLogical,false,0);
  G4VPhysicalVolume* foilPhysical3 = new G4PVPlacement(0,G4ThreeVector(-(fTowerX/2.+fFoilX/2.+fAirX),-(fFoilX/2.),0.),foilLogical1,"foilPhysical3",towerEnvLogical,false,0);
  G4VPhysicalVolume* foilPhysical2 = new G4PVPlacement(0,G4ThreeVector(-(fFoilX/2.),fTowerX/2.+fFoilX/2.+fAirX,0.),foilLogical2,"foilPhysical2",towerEnvLogical,false,0);
  G4VPhysicalVolume* foilPhysical4 = new G4PVPlacement(0,G4ThreeVector(fFoilX/2.,-(fTowerX/2.+fFoilX/2.+fAirX),0.),foilLogical2,"foilPhysical4",towerEnvLogical,false,0);

  new G4LogicalSkinSurface("AlSurf1",foilLogical1,FindSurface("AluminumSurf"));
  new G4LogicalSkinSurface("AlSurf2",foilLogical2,FindSurface("AluminumSurf"));

  foilLogical1->SetVisAttributes(fVisAttrOrange);
  foilLogical2->SetVisAttributes(fVisAttrOrange);

  G4VSolid* sipmEnvSolid = new G4Box("sipmEnvSolid",fSiPMX/2.,fSiPMX/2.,fSiPMH/2.);
  G4LogicalVolume* sipmEnvLogical = new G4LogicalVolume(sipmEnvSolid,FindMaterial("G4_Galactic"),"sipmEnvLogical");

  G4VSolid* waferSolid = new G4Box("waferSolid",fSiPMX/2.,fSiPMX/2.,fFilterT/2.);
  fWaferlogical = new G4LogicalVolume(waferSolid,FindMaterial("Silicon"),"waferLogical");
  G4VPhysicalVolume* waferPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(fSiPMH-fFilterT)/2.),fWaferlogical,"waferPhysical",sipmEnvLogical,false,0);
  G4LogicalSkinSurface* waferSurface = new G4LogicalSkinSurface("SiPMSurf",fWaferlogical,FindSurface("SiPMSurf"));

  fWaferlogical->SetVisAttributes(fVisAttrGreen);

  G4VSolid* windowSolid = new G4Box("windowSolid",fSiPMX/2.,fSiPMX/2.,(fSiPMH-fFilterT)/2.);
  G4LogicalVolume* windowLogical = new G4LogicalVolume(windowSolid,FindMaterial("Glass"),"windowLogical");
  G4VPhysicalVolume* windowPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,-fFilterT/2.),windowLogical,"windowPhysical",sipmEnvLogical,false,0);

  G4VPhysicalVolume* sipmEnvPhysical = new G4PVPlacement(0,G4ThreeVector(0.,0.,(fTowerH+fSiPMH)/2.),sipmEnvLogical,"sipmEnvPhysical",towerEnvLogical,false,0);

  CBDsimFiberParam* fiberParam = new CBDsimFiberParam(fTowerN,fTowerN,fFiberSpec);
  G4VPVParameterisation* cellParam = new CBDsimCellParam(fiberParam);

  G4PVParameterised* towerEnvPhysical = new G4PVParameterised("cell",towerEnvLogical,towertotalEnvLogical,kZAxis,fiberParam->numx*fiberParam->numy,cellParam);

  return worldPhysical;
}

void CBDsimDetectorConstruction::ConstructSDandField() {
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SiPMname = "SiPMSD";

  CBDsimSiPMSD* sipmSD = new CBDsimSiPMSD(SiPMname+"B",SiPMname+"BC",std::make_pair(fTowerN,fTowerN));
  SDman->AddNewDetector(sipmSD);
  fWaferlogical->SetSensitiveDetector(sipmSD);
}

fiberSpec::fiberSpec(G4double rMin_cladC_, G4double rMin_cladS_, G4double rMax_, G4double towerH_,
  G4Material* matPS_, G4Material* matPMMA_, G4Material* matFluoPoly_, G4VisAttributes* visAttrOrange_, G4VisAttributes* visAttrBlue_) {

  rMin_cladC = rMin_cladC_;
  rMin_cladS = rMin_cladS_;
  rMax = rMax_;
  towerH = towerH_;
  matPS = matPS_;
  matPMMA = matPMMA_;
  matFluoPoly = matFluoPoly_;
  visAttrOrange = visAttrOrange_;
  visAttrBlue = visAttrBlue_;
}

fiberSpec::~fiberSpec() {}
