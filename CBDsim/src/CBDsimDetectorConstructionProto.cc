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
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4Tubs.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "CBDsimSiPMSD.hh"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {
// Beam +z through thin z (5 mm). LG couples to 40x5 mm face (x-z plane, -y face), extends to -y.
constexpr G4double kHxWide = 20.0 * mm;
constexpr G4double kHyLong = 30.0 * mm;
constexpr G4double kHzThin = 2.5 * mm;
constexpr G4double kLguide = 30.0 * mm;
/** LG outlet circle radius (fits inside 40x60 mm coupling face). */
constexpr G4double kRguide = 7.5 * mm;
constexpr G4int kNPhi = 64;
constexpr G4int kNSlice = 20;
/**
 * y-offset between PS bottom and LG top. Set to 0 (flush) with vacuum world: no n≠1 layer between
 * scint and LG for optics. If GeomNav1002 increases at coplanar boundaries, restore a μm gap.
 */
constexpr G4double kLGZGap = 0.0 * mm;
/** y gap between LG tip and SiPM package (LG tessellated cap vs G4Tubs). 0 = flush in vacuum. */
constexpr G4double kSiPMLGAirGap = 0.0 * mm;

/** Z-offset margin (same role as before outer-air removal): keeps trigger assemblies off z=0 in world coordinates. */
constexpr G4double kOuterAirSafetyMargin = 1.0 * mm;

// Match legacy tower wrapping (CBDsimDetectorConstruction)
constexpr G4double kFoilT = 0.016 * mm;
/** Inset from nominal foil half-extents so adjacent foil boxes do not share corner volume (avoids ~um overlaps). */
constexpr G4double kFoilCornerInset = 0.02 * mm;
constexpr G4double kAirGap = 0.01 * mm;
constexpr G4double kEnvMarginXY = 2.0 * mm;
/** Used with hzEnvThin to set hzEnv (LG tip radius in local z). */
constexpr G4double kEnvMarginZ = 0.1 * mm;
/**
 * Clear gap along world z between the two scintillator tiles (beam +z through thin z).
 * Tile z full thickness = 2*kHzThin. Center separation must be at least that; we use
 * 2*kHzThin + kTrig12TileFaceGapZ so the slabs do not overlap.
 */
constexpr G4double kTrig12TileFaceGapZ = 1.0 * mm;
/** |kZTrig2 - kZTrig1| = distance between scintillator centers along z. */
constexpr G4double kTrig12TileCenterSeparationZ = 2.0 * kHzThin + kTrig12TileFaceGapZ;
// SiPM package (match CBDsimDetectorConstruction front stack thicknesses)
constexpr G4double kSiPMH = 0.3 * mm;
constexpr G4double kFilterT = 0.01 * mm;

void rectBoundaryXZ(G4double phi, G4double hx, G4double hz, G4double& px, G4double& pz) {
  const G4double c = std::cos(phi);
  const G4double s = std::sin(phi);
  G4double scale;
  if (std::fabs(c) < 1e-12) {
    scale = hz / std::fabs(s);
  } else if (std::fabs(s) < 1e-12) {
    scale = hx / std::fabs(c);
  } else {
    scale = std::min(hx / std::fabs(c), hz / std::fabs(s));
  }
  px = scale * c;
  pz = scale * s;
}

G4TessellatedSolid* BuildLightGuideTessellated() {
  auto* ts = new G4TessellatedSolid("ProtoLightGuide");

  const G4int nRings = kNSlice + 1;
  std::vector<G4ThreeVector> v;
  v.reserve(static_cast<size_t>(nRings * kNPhi));

  for (G4int i = 0; i < nRings; ++i) {
    const G4double t = static_cast<G4double>(i) / static_cast<G4double>(kNSlice);
    const G4double y = -kHyLong - kLGZGap - t * kLguide;
    for (G4int j = 0; j < kNPhi; ++j) {
      const G4double phi = twopi * static_cast<G4double>(j) / static_cast<G4double>(kNPhi);
      G4double rx, rz;
      rectBoundaryXZ(phi, kHxWide, kHzThin, rx, rz);
      const G4double cx = kRguide * std::cos(phi);
      const G4double cz = kRguide * std::sin(phi);
      const G4double px = (1.0 - t) * rx + t * cx;
      const G4double pz = (1.0 - t) * rz + t * cz;
      v.emplace_back(px, y, pz);
    }
  }

  auto idx = [](G4int ring, G4int j) { return static_cast<size_t>(ring * kNPhi + j); };

  for (G4int i = 0; i < kNSlice; ++i) {
    for (G4int j = 0; j < kNPhi; ++j) {
      const G4int jp = (j + 1) % kNPhi;
      const G4ThreeVector& a = v[idx(i, j)];
      const G4ThreeVector& b = v[idx(i, jp)];
      const G4ThreeVector& c = v[idx(i + 1, jp)];
      const G4ThreeVector& d = v[idx(i + 1, j)];
      // Quad a–b / d–c: outward normals for positive enclosed volume (GeomSolids1001 if wrong).
      ts->AddFacet(new G4TriangularFacet(a, b, c, ABSOLUTE));
      ts->AddFacet(new G4TriangularFacet(a, c, d, ABSOLUTE));
    }
  }

  const G4ThreeVector cbot(0.0, -kHyLong - kLGZGap, 0.0);
  for (G4int j = 0; j < kNPhi; ++j) {
    const G4int jp = (j + 1) % kNPhi;
    ts->AddFacet(new G4TriangularFacet(cbot, v[idx(0, jp)], v[idx(0, j)], ABSOLUTE));
  }

  const G4ThreeVector ctop(0.0, -kHyLong - kLGZGap - kLguide, 0.0);
  for (G4int j = 0; j < kNPhi; ++j) {
    const G4int jp = (j + 1) % kNPhi;
    ts->AddFacet(new G4TriangularFacet(ctop, v[idx(kNSlice, j)], v[idx(kNSlice, jp)], ABSOLUTE));
  }

  ts->SetSolidClosed(true);
  return ts;
}

}  // namespace

CBDsimDetectorConstructionProto::CBDsimDetectorConstructionProto() {
  DefineMaterials();
  fVisWorld = new G4VisAttributes(false);
  fVisScint = new G4VisAttributes(G4Colour(0.2, 0.6, 1.0, 1.0));
  fVisScint->SetVisibility(true);
  fVisLG = new G4VisAttributes(G4Colour(1.0, 0.65, 0.2, 1.0));
  fVisLG->SetVisibility(true);
  fVisFoil = new G4VisAttributes(G4Colour(1.0, 0.5, 0.0, 1.0));
  fVisFoil->SetVisibility(true);
  fVisSiPM = new G4VisAttributes(G4Colour(0.3, 0.7, 0.3));
  fVisSiPM->SetVisibility(true);
}

CBDsimDetectorConstructionProto::~CBDsimDetectorConstructionProto() {
  delete fVisWorld;
  delete fVisScint;
  delete fVisLG;
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
  const G4double hzEnv = std::max(hzEnvThin, kRguide + kEnvMarginZ);

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

  G4TessellatedSolid* lgSolid = BuildLightGuideTessellated();
  auto* lgLog =
      new G4LogicalVolume(lgSolid, FindMaterial("ProtoLG_MatchScint"), "protoLightGuideLog");
  auto* lgPV1 = new G4PVPlacement(trWorld1(G4Transform3D()), lgLog, "protoLightGuidePhys", worldLog, false, 0);
  auto* lgPV2 = new G4PVPlacement(trWorld2(G4Transform3D()), lgLog, "protoLightGuidePhys", worldLog, false, 1);
  lgLog->SetVisAttributes(fVisLG);

  // SiPM: G4Tubs on z; rotateX(-90°) maps local +z to world +y (window toward LG at +y).
  G4RotationMatrix sipmRot;
  sipmRot.rotateX(-halfpi);
  const G4double yLgTip = -kHyLong - kLGZGap - kLguide;
  const G4double ySipmCenter = yLgTip - kSiPMH * 0.5 - kSiPMLGAirGap;

  auto* sipmEnvS = new G4Tubs("protoSipmEnv", 0., kRguide, kSiPMH * 0.5, 0., twopi);
  auto* sipmEnvLog =
      new G4LogicalVolume(sipmEnvS, FindMaterial("G4_Galactic"), "protoSipmEnvLog");
  sipmEnvLog->SetVisAttributes(fVisWorld);

  auto* sipmWindowS =
      new G4Tubs("protoSipmWindow", 0., kRguide, (kSiPMH - kFilterT) * 0.5, 0., twopi);
  auto* sipmWindowLog = new G4LogicalVolume(sipmWindowS, FindMaterial("Glass"), "protoSipmWindowLog");
  new G4PVPlacement(nullptr, {0., 0., kFilterT * 0.5}, sipmWindowLog, "protoSipmWindowPhys", sipmEnvLog,
                    false, 0);

  auto* sipmWaferS = new G4Tubs("protoSipmWafer", 0., kRguide, kFilterT * 0.5, 0., twopi);
  auto* sipmWaferLog =
      new G4LogicalVolume(sipmWaferS, FindMaterial("Silicon"), "protoSipmWaferLog");
  new G4PVPlacement(nullptr, {0., 0., -(kSiPMH - kFilterT) * 0.5}, sipmWaferLog, "protoSipmWaferPhys",
                    sipmEnvLog, false, 0);
  new G4LogicalSkinSurface("protoSiPMSurf", sipmWaferLog, FindSurface("SiPMSurf"));
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
  }

  // LG | world vacuum (protoWorldLog = G4_Galactic).
  new G4LogicalBorderSurface("protoLGFoilBorder_protoWorldPhys_T1", lgPV1, worldPhys, FindSurface("AluminumSurf"));
  new G4LogicalBorderSurface("protoLGFoilBorder_protoWorldPhys_T2", lgPV2, worldPhys, FindSurface("AluminumSurf"));

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
