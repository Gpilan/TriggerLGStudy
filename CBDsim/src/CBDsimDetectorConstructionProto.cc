#include "CBDsimDetectorConstructionProto.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "CBDsimSiPMSD.hh"

namespace {
// Beam +z through thin z (5 mm). Scint 10×60 mm (half 5×30 mm in x,y), thin z half = 2.5 mm (5 mm full).
constexpr G4double kHxWide = 5.0 * mm;
constexpr G4double kHyLong = 30.0 * mm;
constexpr G4double kHzThin = 2.5 * mm;
/** Rectangular SiPM + window on -y face: 9.9 mm (x) × 5 mm (z), centered on tile. */
constexpr G4double kSipmHalfX = 4.95 * mm;
constexpr G4double kSipmHalfZ = kHzThin;
/** Vacuum gap between scint -y face and SiPM window (0 = flush). */
constexpr G4double kSiPMAirGap = 0.0 * mm;

/** Z-offset margin: keeps trigger assemblies off z=0 in world coordinates. */
constexpr G4double kOuterAirSafetyMargin = 1.0 * mm;

constexpr G4double kFoilT = 0.016 * mm;
constexpr G4double kFoilCornerInset = 0.02 * mm;
constexpr G4double kAirGap = 0.01 * mm;
constexpr G4double kEnvMarginZ = 0.1 * mm;
constexpr G4double kTrig12TileFaceGapZ = 1.0 * mm;
constexpr G4double kTrig12TileCenterSeparationZ = 2.0 * kHzThin + kTrig12TileFaceGapZ;
constexpr G4double kSiPMH = 0.3 * mm;
constexpr G4double kFilterT = 0.01 * mm;

/** Half-width in x of each Al strip on -y beside the SiPM (one gap = kHxWide − kSipmHalfX → half = half gap). */
constexpr G4double kFoilYmStripHalfX = 0.5 * (kHxWide - kSipmHalfX);

}  // namespace

CBDsimDetectorConstructionProto::CBDsimDetectorConstructionProto() {
  DefineMaterials();
  fVisWorld = new G4VisAttributes(false);
  fVisScint = new G4VisAttributes(G4Colour(0.2, 0.6, 1.0, 1.0));
  fVisScint->SetVisibility(true);
  fVisFoil = new G4VisAttributes(G4Colour(1.0, 0.5, 0.0, 1.0));
  fVisFoil->SetVisibility(true);
  fVisSiPM = new G4VisAttributes(G4Colour(0.3, 0.7, 0.3));
  fVisSiPM->SetVisibility(true);
}

CBDsimDetectorConstructionProto::~CBDsimDetectorConstructionProto() {
  delete fVisWorld;
  delete fVisScint;
  delete fVisFoil;
  delete fVisSiPM;
}

void CBDsimDetectorConstructionProto::DefineMaterials() {
  fMaterials = CBDsimMaterials::GetInstance();
  G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
}

G4VPhysicalVolume* CBDsimDetectorConstructionProto::Construct() {
  fProtoWaferLog = nullptr;
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();
  G4LogicalBorderSurface::CleanSurfaceTable();

  auto* worldSolid = new G4Box("protoWorld", 0.5 * m, 0.5 * m, 0.5 * m);
  auto* worldLog = new G4LogicalVolume(worldSolid, FindMaterial("G4_Galactic"), "protoWorldLog");
  auto* worldPhys = new G4PVPlacement(nullptr, {}, worldLog, "protoWorldPhys", nullptr, false, 0);
  worldLog->SetVisAttributes(fVisWorld);

  const G4double hzEnvThin = kHzThin + kAirGap + kFoilT + kEnvMarginZ;
  const G4double hzEnv = hzEnvThin;

  // Half the z distance between T1 and T2 tile (scint) centers: |kZTrig2 - kZTrig1| = 2 * kEnvZHalfSep.
  const G4double kEnvZHalfSep = 0.5 * kTrig12TileCenterSeparationZ;

  /** Same absolute z positions as before (outer at kProtoAssemblyZ0, env at ±kEnvZHalfSep). World is vacuum (G4_Galactic). */
  const G4double kProtoAssemblyZ0 = kEnvZHalfSep + hzEnv + kOuterAirSafetyMargin;
  const G4double kZTrig1 = kProtoAssemblyZ0 - kEnvZHalfSep;
  const G4double kZTrig2 = kProtoAssemblyZ0 + kEnvZHalfSep;

  G4RotationMatrix rotTrig2;
  rotTrig2.rotateZ(halfpi);

  auto trWorld1 = [&](const G4Transform3D& localInAssembly) {
    return G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., kZTrig1)) * localInAssembly;
  };
  auto trWorld2 = [&](const G4Transform3D& localInAssembly) {
    return G4Transform3D(rotTrig2, G4ThreeVector(0., 0., kZTrig2)) * localInAssembly;
  };

  auto* scintSolid = new G4Box("protoScint", kHxWide, kHyLong, kHzThin);
  auto* scintLog =
      new G4LogicalVolume(scintSolid, FindMaterial("Polystyrene"), "protoScintLog");
  new G4PVPlacement(trWorld1(G4Transform3D()), scintLog, "protoScintPhys", worldLog, false, 0);
  new G4PVPlacement(trWorld2(G4Transform3D()), scintLog, "protoScintPhys", worldLog, false, 1);
  scintLog->SetVisAttributes(fVisScint);

  // SiPM: rectangular box on -y face; rotateX(-90°) maps local +z (stack axis) to world -y (toward scint +y).
  G4RotationMatrix sipmRot;
  sipmRot.rotateX(-halfpi);
  const G4double ySipmCenter = -kHyLong - kSiPMAirGap - kSiPMH * 0.5;

  auto* sipmEnvS =
      new G4Box("protoSipmEnv", kSipmHalfX, kSipmHalfZ, kSiPMH * 0.5);
  auto* sipmEnvLog =
      new G4LogicalVolume(sipmEnvS, FindMaterial("G4_Galactic"), "protoSipmEnvLog");
  sipmEnvLog->SetVisAttributes(fVisWorld);

  auto* sipmWindowS = new G4Box(
      "protoSipmWindow", kSipmHalfX, kSipmHalfZ, (kSiPMH - kFilterT) * 0.5);
  auto* sipmWindowLog = new G4LogicalVolume(sipmWindowS, FindMaterial("Glass"), "protoSipmWindowLog");
  new G4PVPlacement(nullptr, {0., 0., kFilterT * 0.5}, sipmWindowLog, "protoSipmWindowPhys", sipmEnvLog,
                    false, 0);

  auto* sipmWaferS = new G4Box("protoSipmWafer", kSipmHalfX, kSipmHalfZ, kFilterT * 0.5);
  auto* sipmWaferLog =
      new G4LogicalVolume(sipmWaferS, FindMaterial("SiPM_WaferSilicon"), "protoSipmWaferLog");
  new G4PVPlacement(nullptr, {0., 0., -(kSiPMH - kFilterT) * 0.5}, sipmWaferLog, "protoSipmWaferPhys",
                    sipmEnvLog, false, 0);
  sipmWaferLog->SetVisAttributes(fVisSiPM);
  fProtoWaferLog = sipmWaferLog;

  const G4double xh = kHxWide + kAirGap + 0.5 * kFoilT;
  const G4double yh = kHyLong + kAirGap + 0.5 * kFoilT;
  const G4double zh = kHzThin + kAirGap + 0.5 * kFoilT;
  const G4double ft2 = 0.5 * kFoilT;
  const G4double xhF = xh - kFoilCornerInset;
  const G4double yhF = yh - kFoilCornerInset;
  const G4double zhF = zh - kFoilCornerInset;

  auto* foilXmS = new G4Box("protoFoilXm", ft2, yhF, zhF);
  auto* foilXmLog = new G4LogicalVolume(foilXmS, FindMaterial("Aluminum"), "protoFoilXmLog");
  foilXmLog->SetVisAttributes(fVisFoil);
  new G4LogicalSkinSurface("protoAlSurfXm", foilXmLog, FindSurface("AluminumSurf"));

  auto* foilXpS = new G4Box("protoFoilXp", ft2, yhF, zhF);
  auto* foilXpLog = new G4LogicalVolume(foilXpS, FindMaterial("Aluminum"), "protoFoilXpLog");
  foilXpLog->SetVisAttributes(fVisFoil);
  new G4LogicalSkinSurface("protoAlSurfXp", foilXpLog, FindSurface("AluminumSurf"));

  auto* foilZmS = new G4Box("protoFoilZm", xhF, yhF, ft2);
  auto* foilZmLog = new G4LogicalVolume(foilZmS, FindMaterial("Aluminum"), "protoFoilZmLog");
  foilZmLog->SetVisAttributes(fVisFoil);
  new G4LogicalSkinSurface("protoAlSurfZm", foilZmLog, FindSurface("AluminumSurf"));

  auto* foilZpS = new G4Box("protoFoilZp", xhF, yhF, ft2);
  auto* foilZpLog = new G4LogicalVolume(foilZpS, FindMaterial("Aluminum"), "protoFoilZpLog");
  foilZpLog->SetVisAttributes(fVisFoil);
  new G4LogicalSkinSurface("protoAlSurfZp", foilZpLog, FindSurface("AluminumSurf"));

  const G4double yFoilPlus = kHyLong + kAirGap + ft2;
  auto* foilYpS = new G4Box("protoFoilYp", xhF, ft2, zhF);
  auto* foilYpLog = new G4LogicalVolume(foilYpS, FindMaterial("Aluminum"), "protoFoilYpLog");
  foilYpLog->SetVisAttributes(fVisFoil);
  new G4LogicalSkinSurface("protoAlSurfYp", foilYpLog, FindSurface("AluminumSurf"));

  const G4double yFoilMinus = -kHyLong - kAirGap - ft2;
  const G4double xLeftYm = -0.5 * (kHxWide + kSipmHalfX);
  const G4double xRightYm = 0.5 * (kHxWide + kSipmHalfX);
  auto* foilYmStripS = new G4Box("protoFoilYmStrip", kFoilYmStripHalfX, ft2, zhF);
  auto* foilYmStripLog = new G4LogicalVolume(foilYmStripS, FindMaterial("Aluminum"), "protoFoilYmStripLog");
  foilYmStripLog->SetVisAttributes(fVisFoil);
  new G4LogicalSkinSurface("protoAlSurfYmStrip", foilYmStripLog, FindSurface("AluminumSurf"));

  {
    const G4Transform3D localSipm(G4Transform3D(sipmRot, G4ThreeVector(0., ySipmCenter, 0.)));
    new G4PVPlacement(trWorld1(localSipm), sipmEnvLog, "protoSipmEnvPhys", worldLog, false, 0);
    new G4PVPlacement(trWorld2(localSipm), sipmEnvLog, "protoSipmEnvPhys", worldLog, false, 1);
  }
  {
    const G4Transform3D tx(G4RotationMatrix(), G4ThreeVector(-xh, 0., 0.));
    const G4Transform3D tpx(G4RotationMatrix(), G4ThreeVector(xh, 0., 0.));
    const G4Transform3D tzm(G4RotationMatrix(), G4ThreeVector(0., 0., -zh));
    const G4Transform3D tzp(G4RotationMatrix(), G4ThreeVector(0., 0., zh));
    const G4Transform3D typ(G4RotationMatrix(), G4ThreeVector(0., yFoilPlus, 0.));
    new G4PVPlacement(trWorld1(tx), foilXmLog, "protoFoilXmPhys", worldLog, false, 0);
    new G4PVPlacement(trWorld2(tx), foilXmLog, "protoFoilXmPhys", worldLog, false, 1);
    new G4PVPlacement(trWorld1(tpx), foilXpLog, "protoFoilXpPhys", worldLog, false, 0);
    new G4PVPlacement(trWorld2(tpx), foilXpLog, "protoFoilXpPhys", worldLog, false, 1);
    new G4PVPlacement(trWorld1(tzm), foilZmLog, "protoFoilZmPhys", worldLog, false, 0);
    new G4PVPlacement(trWorld2(tzm), foilZmLog, "protoFoilZmPhys", worldLog, false, 1);
    new G4PVPlacement(trWorld1(tzp), foilZpLog, "protoFoilZpPhys", worldLog, false, 0);
    new G4PVPlacement(trWorld2(tzp), foilZpLog, "protoFoilZpPhys", worldLog, false, 1);
    new G4PVPlacement(trWorld1(typ), foilYpLog, "protoFoilYpPhys", worldLog, false, 0);
    new G4PVPlacement(trWorld2(typ), foilYpLog, "protoFoilYpPhys", worldLog, false, 1);
    const G4Transform3D tymL(G4RotationMatrix(), G4ThreeVector(xLeftYm, yFoilMinus, 0.));
    const G4Transform3D tymR(G4RotationMatrix(), G4ThreeVector(xRightYm, yFoilMinus, 0.));
    new G4PVPlacement(trWorld1(tymL), foilYmStripLog, "protoFoilYmPhys", worldLog, false, 0);
    new G4PVPlacement(trWorld2(tymL), foilYmStripLog, "protoFoilYmPhys", worldLog, false, 1);
    new G4PVPlacement(trWorld1(tymR), foilYmStripLog, "protoFoilYmPhys", worldLog, false, 2);
    new G4PVPlacement(trWorld2(tymR), foilYmStripLog, "protoFoilYmPhys", worldLog, false, 3);
  }

  return worldPhys;
}

void CBDsimDetectorConstructionProto::ConstructSDandField() {
  // After CloseGeometry(); reports volume overlaps that often precede GeomNav1002 at runtime.
  if (auto* w = G4PhysicalVolumeStore::GetInstance()->GetVolume("protoWorldPhys")) {
    w->CheckOverlaps();
  }
  if (!fProtoWaferLog) return;
  auto* SDman = G4SDManager::GetSDMpointer();
  auto* sipmSD = new CBDsimSiPMSD("SiPMSDB", "SiPMSDBC", std::make_pair(1, 1));
  SDman->AddNewDetector(sipmSD);
  fProtoWaferLog->SetSensitiveDetector(sipmSD);
}
