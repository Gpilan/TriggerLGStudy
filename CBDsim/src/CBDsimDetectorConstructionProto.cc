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
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "CBDsimSiPMSD.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <vector>

namespace {
// Beam +z through thin z (5 mm). LG couples to 10x5 mm face (x-z plane, -y face), extends to -y.
constexpr G4double kHxWide = 20 * mm;
constexpr G4double kHyLong = 30.0 * mm;
constexpr G4double kHzThin = 2.5 * mm;
constexpr G4double kLguide = 30.0 * mm;
/** LG outlet circle radius (kept from LG-v3 baseline). */
constexpr G4double kRguide = 7.5 * mm;
constexpr G4int kNPhi = 256;
constexpr G4int kNSlice = 20;
/**
 * y-offset between PS bottom and LG top. Set to 0 (flush) with vacuum world: no n!=1 layer between
 * scint and LG for optics. If GeomNav1002 increases at coplanar boundaries, restore a um gap.
 */
constexpr G4double kLGZGap = 0.0 * mm;
/** y gap between LG tip and SiPM package (LG tessellated cap vs G4Tubs). 0 = flush in vacuum. */
constexpr G4double kSiPMLGAirGap = 0.0 * mm;
/** Index-matching grease (Gelatin n≈1.52) between LG tip and SiPM window. */
constexpr G4double kSiPMGelT = 0.10 * mm;
/** Gel overlaps into LG tip opening for flush navigation (no vacuum sliver). */
constexpr G4double kLGTipGelOverlap = 0.02 * mm;
/** Axial half-length of tip Al ring (covers LG|SiPM junction, not just 8 um foil). */
constexpr G4double kTipRingHalfY = 0.30 * mm;

/** Z-offset margin (same role as before outer-air removal): keeps trigger assemblies off z=0 in world coordinates. */
constexpr G4double kOuterAirSafetyMargin = 1.0 * mm;

// Match legacy tower wrapping (CBDsimDetectorConstruction)
constexpr G4double kFoilT = 0.016 * mm;
/** Corner foil pad on scint -y face (covers scint-only wedge outside LG polygon). */
constexpr G4double kFoilCornerPadMax = 4.0 * mm;
/** LG tip ring: outer radius extension beyond kRguide (covers polygon–circle sliver at outlet). */
constexpr G4double kTipRingOuterExtra = 0.05 * mm;
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
/** no-LG rectangular SiPM/gel aperture half-width in x (= LG outlet radius). */
constexpr G4double kSipmRectHalfX = kRguide;
/** no-LG aperture half-height in z (full scint thin dimension). */
constexpr G4double kSipmRectHalfZ = kHzThin;

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
      // Quad a-b / d-c: outward normals for positive enclosed volume (GeomSolids1001 if wrong).
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

/** When set (CBDsim_PROTO_NO_LG=1), skip LG; rectangular SiPM on scint -y with foil strips. */
bool ProtoBuildLightGuide() {
  const char* v = std::getenv("CBDsim_PROTO_NO_LG");
  if (!v || v[0] == '\0') return true;
  switch (v[0]) {
    case '0':
    case 'f':
    case 'F':
    case 'n':
    case 'N':
      return true;
    default:
      return false;
  }
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
  fVisGel = new G4VisAttributes(G4Colour(0.95, 0.85, 0.2, 0.45));
  fVisGel->SetVisibility(true);
}

CBDsimDetectorConstructionProto::~CBDsimDetectorConstructionProto() {
  delete fVisWorld;
  delete fVisScint;
  delete fVisLG;
  delete fVisFoil;
  delete fVisSiPM;
  delete fVisGel;
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

  /** Same absolute z positions as before (outer at kProtoAssemblyZ0, env at +/-kEnvZHalfSep). World is vacuum (G4_Galactic). */
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
  auto* scintPV1 =
      new G4PVPlacement(trWorld1(G4Transform3D()), scintLog, "protoScintPhys", worldLog, false, 0);
  auto* scintPV2 =
      new G4PVPlacement(trWorld2(G4Transform3D()), scintLog, "protoScintPhys", worldLog, false, 1);
  scintLog->SetVisAttributes(fVisScint);

  const bool withLG = ProtoBuildLightGuide();
  {
    const char* envNoLg = std::getenv("CBDsim_PROTO_NO_LG");
    G4cout << "\n=== [Proto geometry] CBDsim_PROTO_NO_LG="
           << (envNoLg ? envNoLg : "(unset)") << " => " << (withLG ? "LG" : "no-LG")
           << " ===\n"
           << G4endl;
  }
  G4PVPlacement* lgPV1 = nullptr;
  G4PVPlacement* lgPV2 = nullptr;
  if (withLG) {
    G4TessellatedSolid* lgSolid = BuildLightGuideTessellated();
    auto* lgLog =
        new G4LogicalVolume(lgSolid, FindMaterial("ProtoLG_MatchScint"), "protoLightGuideLog");
    lgPV1 = new G4PVPlacement(trWorld1(G4Transform3D()), lgLog, "protoLightGuidePhys", worldLog, false, 0);
    lgPV2 = new G4PVPlacement(trWorld2(G4Transform3D()), lgLog, "protoLightGuidePhys", worldLog, false, 1);
    lgLog->SetVisAttributes(fVisLG);
  } else {
    G4cout << "[Proto geometry] no-LG mode: rect SiPM " << (2.0 * kSipmRectHalfX) / mm << " x "
           << (2.0 * kSipmRectHalfZ) / mm << " mm on scint -y (kHxWide=" << kHxWide / mm
           << " mm), v5 gel+foil" << G4endl;
  }

  // SiPM stack: LG mode = circular tubs at LG tip; no-LG = rectangular boxes on scint -y.
  G4RotationMatrix sipmRot;
  sipmRot.rotateX(-halfpi);
  const G4double yScintFace = -kHyLong - kLGZGap;
  const G4double yCouplingFace = withLG ? (yScintFace - kLguide) : yScintFace;
  const G4double gelHalfY = kSiPMGelT * 0.5;
  const G4double windowHalfY = (kSiPMH - kFilterT) * 0.5;
  const G4double waferHalfY = kFilterT * 0.5;
  const G4double yGelCenter = yCouplingFace - gelHalfY + kLGTipGelOverlap;
  const G4double yWindowCenter = yGelCenter - gelHalfY - windowHalfY;
  const G4double yWaferCenter = yWindowCenter - windowHalfY - waferHalfY;

  G4cout << "[Proto geometry] coupling face y=" << yCouplingFace / mm << " mm, gel center y="
         << yGelCenter / mm << " mm (gel top y=" << (yGelCenter + gelHalfY) / mm << " mm, overlap "
         << kLGTipGelOverlap / mm << " mm)" << G4endl;

  G4LogicalVolume* sipmGelLog = nullptr;
  G4LogicalVolume* sipmWindowLog = nullptr;
  G4LogicalVolume* sipmEnvLog = nullptr;
  G4LogicalVolume* sipmWaferLog = nullptr;
  G4PVPlacement* sipmWaferPV = nullptr;

  if (withLG) {
    auto* sipmGelS = new G4Tubs("protoSipmGel", 0., kRguide, gelHalfY, 0., twopi);
    sipmGelLog = new G4LogicalVolume(sipmGelS, FindMaterial("Gelatin"), "protoSipmGelLog");
    sipmGelLog->SetVisAttributes(fVisGel);

    auto* sipmWindowS = new G4Tubs("protoSipmWindow", 0., kRguide, windowHalfY, 0., twopi);
    sipmWindowLog =
        new G4LogicalVolume(sipmWindowS, FindMaterial("Glass"), "protoSipmWindowLog");
    sipmWindowLog->SetVisAttributes(fVisSiPM);

    auto* sipmEnvS = new G4Tubs("protoSipmEnv", 0., kRguide, waferHalfY, 0., twopi);
    sipmEnvLog = new G4LogicalVolume(sipmEnvS, FindMaterial("G4_Galactic"), "protoSipmEnvLog");
    sipmEnvLog->SetVisAttributes(fVisWorld);

    auto* sipmWaferS = new G4Tubs("protoSipmWafer", 0., kRguide, waferHalfY, 0., twopi);
    sipmWaferLog =
        new G4LogicalVolume(sipmWaferS, FindMaterial("SiPM_WaferSilicon"), "protoSipmWaferLog");
    sipmWaferPV =
        new G4PVPlacement(nullptr, G4ThreeVector(), sipmWaferLog, "protoSipmWaferPhys", sipmEnvLog,
                          false, 0);
  } else {
    auto* sipmGelS =
        new G4Box("protoSipmGel", kSipmRectHalfX, kSipmRectHalfZ, gelHalfY);
    sipmGelLog = new G4LogicalVolume(sipmGelS, FindMaterial("Gelatin"), "protoSipmGelLog");
    sipmGelLog->SetVisAttributes(fVisGel);

    auto* sipmWindowS =
        new G4Box("protoSipmWindow", kSipmRectHalfX, kSipmRectHalfZ, windowHalfY);
    sipmWindowLog =
        new G4LogicalVolume(sipmWindowS, FindMaterial("Glass"), "protoSipmWindowLog");
    sipmWindowLog->SetVisAttributes(fVisSiPM);

    auto* sipmEnvS =
        new G4Box("protoSipmEnv", kSipmRectHalfX, kSipmRectHalfZ, waferHalfY);
    sipmEnvLog = new G4LogicalVolume(sipmEnvS, FindMaterial("G4_Galactic"), "protoSipmEnvLog");
    sipmEnvLog->SetVisAttributes(fVisWorld);

    auto* sipmWaferS =
        new G4Box("protoSipmWafer", kSipmRectHalfX, kSipmRectHalfZ, waferHalfY);
    sipmWaferLog =
        new G4LogicalVolume(sipmWaferS, FindMaterial("SiPM_WaferSilicon"), "protoSipmWaferLog");
    sipmWaferPV =
        new G4PVPlacement(nullptr, G4ThreeVector(), sipmWaferLog, "protoSipmWaferPhys", sipmEnvLog,
                          false, 0);
  }
  sipmWaferLog->SetVisAttributes(fVisSiPM);
  fProtoWaferLog = sipmWaferLog;

  const G4double xi = kHxWide + kAirGap;
  const G4double yi = kHyLong + kAirGap;
  const G4double zi = kHzThin + kAirGap;
  const G4double xo = xi + kFoilT;
  const G4double yo = yi + kFoilT;
  const G4double zo = zi + kFoilT;
  const G4double ft2 = 0.5 * kFoilT;

  // Continuous 5-face Al wrap (sides + beam +y top), open on -y for LG coupling.
  auto* foilOuterS = new G4Box("protoFoilOuter", xo, yo, zo);
  auto* foilInnerS = new G4Box("protoFoilInner", xi, yi, zi);
  auto* foilShellS = new G4SubtractionSolid("protoFoilShell", foilOuterS, foilInnerS);
  auto* foilBottomCutS =
      new G4Box("protoFoilBottomCut", xo + 1.e-3 * mm, kFoilT + 1.e-3 * mm, zo + 1.e-3 * mm);
  auto* foilWrapS = new G4SubtractionSolid("protoFoilWrap", foilShellS, foilBottomCutS, nullptr,
                                           G4ThreeVector(0., -yo + ft2, 0.));
  auto* foilWrapLog = new G4LogicalVolume(foilWrapS, FindMaterial("Aluminum"), "protoFoilWrapLog");
  foilWrapLog->SetVisAttributes(fVisFoil);
  new G4LogicalSkinSurface("protoAlSurfWrap", foilWrapLog, FindSurface("AluminumSurf"));

  // Wedge leak patch: thin Al pads at the 4 corners of the scint -y face (LG polygon misses these).
  const G4double cornerPadHalf =
      std::min(kFoilCornerPadMax, 0.45 * std::min(kHxWide, kHzThin));
  auto* foilCornerS = new G4Box("protoFoilCorner", cornerPadHalf, ft2, cornerPadHalf);
  auto* foilCornerLog =
      new G4LogicalVolume(foilCornerS, FindMaterial("Aluminum"), "protoFoilCornerLog");
  foilCornerLog->SetVisAttributes(fVisFoil);
  new G4LogicalSkinSurface("protoAlSurfCorner", foilCornerLog, FindSurface("AluminumSurf"));
  const G4double yCorner = -kHyLong - kAirGap - ft2;
  const G4double xCorner = kHxWide - cornerPadHalf;
  const G4double zCorner = kHzThin - cornerPadHalf;

  // LG tip annulus (LG mode only): blocks radial leak at outlet.
  const G4double tipRingInnerR =
      kRguide * std::cos(pi / static_cast<G4double>(kNPhi)) - 0.01 * mm;
  const G4double tipRingOuterR = kRguide + kAirGap + kFoilT + kTipRingOuterExtra;
  const G4double yTipRing = yCouplingFace - kTipRingHalfY + kLGTipGelOverlap;
  auto* tipRingS =
      new G4Tubs("protoFoilTipRing", tipRingInnerR, tipRingOuterR, kTipRingHalfY, 0., twopi);
  auto* tipRingLog =
      new G4LogicalVolume(tipRingS, FindMaterial("Aluminum"), "protoFoilTipRingLog");
  tipRingLog->SetVisAttributes(fVisFoil);
  new G4LogicalSkinSurface("protoAlSurfTipRing", tipRingLog, FindSurface("AluminumSurf"));
  G4RotationMatrix tipRingRot;
  tipRingRot.rotateX(-halfpi);

  // no-LG: Al strips on scint -y face beside rectangular SiPM (v3 idea, v5 Al surface).
  G4LogicalVolume* foilYmStripLog = nullptr;
  G4double xLeftYm = 0.;
  G4double xRightYm = 0.;
  if (!withLG) {
    const G4double foilYmStripHalfX = 0.5 * (kHxWide - kSipmRectHalfX);
    auto* foilYmStripS =
        new G4Box("protoFoilYmStrip", foilYmStripHalfX, ft2, zi);
    foilYmStripLog =
        new G4LogicalVolume(foilYmStripS, FindMaterial("Aluminum"), "protoFoilYmStripLog");
    foilYmStripLog->SetVisAttributes(fVisFoil);
    new G4LogicalSkinSurface("protoAlSurfYmStrip", foilYmStripLog, FindSurface("AluminumSurf"));
    xLeftYm = -0.5 * (kHxWide + kSipmRectHalfX);
    xRightYm = 0.5 * (kHxWide + kSipmRectHalfX);
  }

  {
    const G4Transform3D localGel(G4Transform3D(sipmRot, G4ThreeVector(0., yGelCenter, 0.)));
    const G4Transform3D localWindow(G4Transform3D(sipmRot, G4ThreeVector(0., yWindowCenter, 0.)));
    const G4Transform3D localEnv(G4Transform3D(sipmRot, G4ThreeVector(0., yWaferCenter, 0.)));
    auto* sipmGelPV1 = new G4PVPlacement(trWorld1(localGel), sipmGelLog, "protoSipmGelPhys", worldLog,
                                         false, 0);
    auto* sipmGelPV2 = new G4PVPlacement(trWorld2(localGel), sipmGelLog, "protoSipmGelPhys", worldLog,
                                         false, 1);
    auto* sipmWindowPV1 = new G4PVPlacement(trWorld1(localWindow), sipmWindowLog, "protoSipmWindowPhys",
                                            worldLog, false, 0);
    auto* sipmWindowPV2 = new G4PVPlacement(trWorld2(localWindow), sipmWindowLog, "protoSipmWindowPhys",
                                            worldLog, false, 1);
    auto* sipmPV1 =
        new G4PVPlacement(trWorld1(localEnv), sipmEnvLog, "protoSipmEnvPhys", worldLog, false, 0);
    auto* sipmPV2 =
        new G4PVPlacement(trWorld2(localEnv), sipmEnvLog, "protoSipmEnvPhys", worldLog, false, 1);

    if (withLG) {
      new G4LogicalBorderSurface("protoLG1ToWorldReflect", lgPV1, worldPhys, FindSurface("AluminumSurf"));
      new G4LogicalBorderSurface("protoWorldToLG1Reflect", worldPhys, lgPV1, FindSurface("AluminumSurf"));
      new G4LogicalBorderSurface("protoLG2ToWorldReflect", lgPV2, worldPhys, FindSurface("AluminumSurf"));
      new G4LogicalBorderSurface("protoWorldToLG2Reflect", worldPhys, lgPV2, FindSurface("AluminumSurf"));

      new G4LogicalBorderSurface("protoLG1ToGel1Trans", lgPV1, sipmGelPV1, FindSurface("AirSurf"));
      new G4LogicalBorderSurface("protoGel1ToLG1Trans", sipmGelPV1, lgPV1, FindSurface("AirSurf"));
      new G4LogicalBorderSurface("protoLG2ToGel2Trans", lgPV2, sipmGelPV2, FindSurface("AirSurf"));
      new G4LogicalBorderSurface("protoGel2ToLG2Trans", sipmGelPV2, lgPV2, FindSurface("AirSurf"));
    } else {
      new G4LogicalBorderSurface("protoScint1ToGel1Trans", scintPV1, sipmGelPV1, FindSurface("AirSurf"));
      new G4LogicalBorderSurface("protoGel1ToScint1Trans", sipmGelPV1, scintPV1, FindSurface("AirSurf"));
      new G4LogicalBorderSurface("protoScint2ToGel2Trans", scintPV2, sipmGelPV2, FindSurface("AirSurf"));
      new G4LogicalBorderSurface("protoGel2ToScint2Trans", sipmGelPV2, scintPV2, FindSurface("AirSurf"));
    }

    new G4LogicalBorderSurface("protoGel1ToWindow1Trans", sipmGelPV1, sipmWindowPV1, FindSurface("AirSurf"));
    new G4LogicalBorderSurface("protoWindow1ToGel1Trans", sipmWindowPV1, sipmGelPV1, FindSurface("AirSurf"));
    new G4LogicalBorderSurface("protoGel2ToWindow2Trans", sipmGelPV2, sipmWindowPV2, FindSurface("AirSurf"));
    new G4LogicalBorderSurface("protoWindow2ToGel2Trans", sipmWindowPV2, sipmGelPV2, FindSurface("AirSurf"));

    new G4LogicalBorderSurface("protoWaferToSiPM1Absorb", sipmWaferPV, sipmPV1, FindSurface("SiPMSurf"));
    new G4LogicalBorderSurface("protoWaferToSiPM2Absorb", sipmWaferPV, sipmPV2, FindSurface("SiPMSurf"));
  }
  {
    const G4Transform3D tWrap(G4RotationMatrix(), G4ThreeVector(0., 0., 0.));
    new G4PVPlacement(trWorld1(tWrap), foilWrapLog, "protoFoilWrapPhys", worldLog, false, 0);
    new G4PVPlacement(trWorld2(tWrap), foilWrapLog, "protoFoilWrapPhys", worldLog, false, 1);

    const std::array<G4ThreeVector, 4> cornerPos = {
        G4ThreeVector(xCorner, yCorner, zCorner),
        G4ThreeVector(-xCorner, yCorner, zCorner),
        G4ThreeVector(xCorner, yCorner, -zCorner),
        G4ThreeVector(-xCorner, yCorner, -zCorner),
    };
    for (const auto& pos : cornerPos) {
      const G4Transform3D tc(G4RotationMatrix(), pos);
      new G4PVPlacement(trWorld1(tc), foilCornerLog, "protoFoilCornerPhys", worldLog, false, 0);
      new G4PVPlacement(trWorld2(tc), foilCornerLog, "protoFoilCornerPhys", worldLog, false, 1);
    }

    const G4Transform3D tTipRing(G4Transform3D(tipRingRot, G4ThreeVector(0., yTipRing, 0.)));
    if (withLG) {
      new G4PVPlacement(trWorld1(tTipRing), tipRingLog, "protoFoilTipRingPhys", worldLog, false, 0);
      new G4PVPlacement(trWorld2(tTipRing), tipRingLog, "protoFoilTipRingPhys", worldLog, false, 1);
    }

    if (!withLG && foilYmStripLog) {
      const G4double yFoilMinus = -kHyLong - kAirGap - ft2;
      const G4Transform3D tymL(G4RotationMatrix(), G4ThreeVector(xLeftYm, yFoilMinus, 0.));
      const G4Transform3D tymR(G4RotationMatrix(), G4ThreeVector(xRightYm, yFoilMinus, 0.));
      new G4PVPlacement(trWorld1(tymL), foilYmStripLog, "protoFoilYmPhys", worldLog, false, 0);
      new G4PVPlacement(trWorld2(tymL), foilYmStripLog, "protoFoilYmPhys", worldLog, false, 1);
      new G4PVPlacement(trWorld1(tymR), foilYmStripLog, "protoFoilYmPhys", worldLog, false, 2);
      new G4PVPlacement(trWorld2(tymR), foilYmStripLog, "protoFoilYmPhys", worldLog, false, 3);
    }
  }

  // LG|world: AluminumSurf on taper sides (LG mode). Tip ring + window flush block outlet bypass.

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
