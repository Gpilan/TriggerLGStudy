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
 * Thin gap between PS bottom and LG top (y). Coplanar scint box + tessellated LG can trigger
 * GeomNav1002 "stuck track" when many optical secondaries hug boundaries; μm-scale air is negligible
 * for bulk physics but stabilizes navigation. Tune down only if you verify overlaps are clean.
 */
constexpr G4double kLGZGap = 0.001 * mm;
/** Same idea at LG tip vs SiPM window (coplanar tessellated cap vs G4Tubs). */
constexpr G4double kSiPMLGAirGap = 0.001 * mm;

/** Extra half-thickness on protoOuterAirEnv so both triggers (incl. T2 Rz) sit comfortably inside one air mother. */
constexpr G4double kOuterAirSafetyMargin = 1.0 * mm;

// Match legacy tower wrapping (CBDsimDetectorConstruction)
constexpr G4double kFoilT = 0.016 * mm;
/** Inset from nominal foil half-extents so adjacent foil boxes do not share corner volume (avoids ~um overlaps). */
constexpr G4double kFoilCornerInset = 0.02 * mm;
constexpr G4double kAirGap = 0.01 * mm;
/** Extra half-width on protoAirEnv in x,y (outer air already clears world; keep modest). */
constexpr G4double kEnvMarginXY = 2.0 * mm;
/** Extra half-length on protoAirEnv in z (LG tip spans ~kRguide in local z — must fit inside mother). */
constexpr G4double kEnvMarginZ = 0.1 * mm;
/** Gap between outer faces of the two protoAirEnv boxes along world z (env centers at ±(hzEnv + gap/2)). */
constexpr G4double kGapBetweenEnvBoxes = 1.0 * mm;
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
      ts->AddFacet(new G4TriangularFacet(a, d, c, ABSOLUTE));
      ts->AddFacet(new G4TriangularFacet(a, c, b, ABSOLUTE));
    }
  }

  const G4ThreeVector cbot(0.0, -kHyLong - kLGZGap, 0.0);
  for (G4int j = 0; j < kNPhi; ++j) {
    const G4int jp = (j + 1) % kNPhi;
    ts->AddFacet(new G4TriangularFacet(cbot, v[idx(0, j)], v[idx(0, jp)], ABSOLUTE));
  }

  const G4ThreeVector ctop(0.0, -kHyLong - kLGZGap - kLguide, 0.0);
  for (G4int j = 0; j < kNPhi; ++j) {
    const G4int jp = (j + 1) % kNPhi;
    ts->AddFacet(new G4TriangularFacet(ctop, v[idx(kNSlice, jp)], v[idx(kNSlice, j)], ABSOLUTE));
  }

  ts->SetSolidClosed(true);
  return ts;
}

/** Descendant search without G4VPhysicalVolume::GetMother / G4PVPlacement::GetMotherPhysical (API varies by G4 version; see repo _G4_version_release). */
G4VPhysicalVolume* FindDaughterPVByName(G4VPhysicalVolume* parent, const G4String& name) {
  auto* lv = parent->GetLogicalVolume();
  for (G4int i = 0; i < lv->GetNoDaughters(); ++i) {
    auto* d = lv->GetDaughter(i);
    if (d->GetName() == name) return d;
    if (auto* found = FindDaughterPVByName(d, name)) return found;
  }
  return nullptr;
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
  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
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
  auto* worldLog = new G4LogicalVolume(worldSolid, FindMaterial("G4_AIR"), "protoWorldLog");
  auto* worldPhys = new G4PVPlacement(nullptr, {}, worldLog, "protoWorldPhys", nullptr, false, 0);
  worldLog->SetVisAttributes(fVisWorld);

  const G4double hxEnv = kHxWide + kAirGap + kFoilT + kEnvMarginXY;
  const G4double hzEnvThin = kHzThin + kAirGap + kFoilT + kEnvMarginZ;
  // Tessellated LG tip is a circle of radius kRguide in local x-z; mother must cover ±kRguide in z.
  const G4double hzEnv = std::max(hzEnvThin, kRguide + kEnvMarginZ);
  const G4double hyExtentUp = kHyLong + kAirGap + kFoilT;
  const G4double hyExtentDown = kHyLong + kLGZGap + kLguide + kSiPMH;
  const G4double hyEnv = std::max(hyExtentUp, hyExtentDown) + kEnvMarginXY;

  // Env centers: hzEnv must cover LG (±kRguide in z); spacing hzEnv + gap/2 avoids T1/T2 env overlap.
  // Scintillator centers are placed at ±kEnvZHalfSep (same as env centers) — no local z shift (shift would push LG past mother).
  // Effective gap between scint inner faces = 2*kEnvZHalfSep - 2*kHzThin (>> 1 mm if kRguide ~ LG tip radius).
  const G4double kEnvZHalfSep = hzEnv + 0.5 * kGapBetweenEnvBoxes;

  // One outer air envelope (protoOuterAirEnv) fully contains both trigger env placements.
  const G4double outerHx = std::max(hxEnv, hyEnv) + kOuterAirSafetyMargin;
  const G4double outerHy = std::max(hxEnv, hyEnv) + kOuterAirSafetyMargin;
  const G4double outerHz = kEnvZHalfSep + hzEnv + kOuterAirSafetyMargin;

  /** World +z shift so T1/T2 (centers at ±kEnvZHalfSep inside outer) both lie in z>0; avoids straddling z=0. */
  const G4double kProtoAssemblyZ0 = kEnvZHalfSep + hzEnv + kOuterAirSafetyMargin;

  auto* envSolid = new G4Box("protoAirEnv", hxEnv, hyEnv, hzEnv);
  auto* envLog = new G4LogicalVolume(envSolid, FindMaterial("G4_AIR"), "protoAirEnvLog");

  auto* outerEnvSolid = new G4Box("protoOuterAirEnv", outerHx, outerHy, outerHz);
  auto* outerEnvLog = new G4LogicalVolume(outerEnvSolid, FindMaterial("G4_AIR"), "protoOuterAirEnvLog");
  new G4PVPlacement(nullptr, G4ThreeVector(0., 0., kProtoAssemblyZ0), outerEnvLog, "protoOuterAirEnvPhys", worldLog,
                    false, 0);
  outerEnvLog->SetVisAttributes(fVisWorld);

  auto* scintSolid = new G4Box("protoScint", kHxWide, kHyLong, kHzThin);
  auto* scintLog =
      new G4LogicalVolume(scintSolid, FindMaterial("Polystyrene"), "protoScintLog");
  new G4PVPlacement(nullptr, {}, scintLog, "protoScintPhys", envLog, false, 0);
  scintLog->SetVisAttributes(fVisScint);

  G4TessellatedSolid* lgSolid = BuildLightGuideTessellated();
  auto* lgLog =
      new G4LogicalVolume(lgSolid, FindMaterial("ProtoLG_MatchScint"), "protoLightGuideLog");
  new G4PVPlacement(nullptr, {}, lgLog, "protoLightGuidePhys", envLog, false, 0);
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

  new G4PVPlacement(G4Transform3D(sipmRot, G4ThreeVector(0., ySipmCenter, 0.)), sipmEnvLog, "protoSipmEnvPhys",
                    envLog, false, 0);
  new G4PVPlacement(nullptr, {-xh, 0., 0.}, foilXmLog, "protoFoilXmPhys", envLog, false, 0);
  new G4PVPlacement(nullptr, {xh, 0., 0.}, foilXpLog, "protoFoilXpPhys", envLog, false, 0);
  new G4PVPlacement(nullptr, {0., 0., -zh}, foilZmLog, "protoFoilZmPhys", envLog, false, 0);
  new G4PVPlacement(nullptr, {0., 0., zh}, foilZpLog, "protoFoilZpPhys", envLog, false, 0);
  new G4PVPlacement(nullptr, {0., yFoilPlus, 0.}, foilYpLog, "protoFoilYpPhys", envLog, false, 0);

  // Two trigger assemblies share envLog; scint centers at ±kEnvZHalfSep (see kEnvZHalfSep comment above).
  // Trigger 2: +90° about +z vs trigger 1. Env copy numbers 0 / 1 → SiPM SD uses volume depth 2 (env) as in legacy proto.
  auto* envPhysT1 =
      new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -kEnvZHalfSep), envLog, "protoAirEnvPhys_T1", outerEnvLog,
                        false, 0);
  G4RotationMatrix rotTrig2;
  rotTrig2.rotateZ(halfpi);
  auto* envPhysT2 = new G4PVPlacement(G4Transform3D(rotTrig2, G4ThreeVector(0., 0., kEnvZHalfSep)), envLog,
                                        "protoAirEnvPhys_T2", outerEnvLog, false, 1);

  // One border surface per env instance: navigate from each env physical volume (no parent getters; G4 API differs by version).
  if (auto* lgPV1 = FindDaughterPVByName(envPhysT1, "protoLightGuidePhys")) {
    new G4LogicalBorderSurface("protoLGFoilBorder_protoAirEnvPhys_T1", lgPV1, envPhysT1,
                               FindSurface("AluminumSurf"));
  }
  if (auto* lgPV2 = FindDaughterPVByName(envPhysT2, "protoLightGuidePhys")) {
    new G4LogicalBorderSurface("protoLGFoilBorder_protoAirEnvPhys_T2", lgPV2, envPhysT2,
                               FindSurface("AluminumSurf"));
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
