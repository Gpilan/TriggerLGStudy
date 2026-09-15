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
#include "G4UnionSolid.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SDManager.hh"
#include "G4Exception.hh"

#include "CBDsimSiPMSD.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

namespace {
// Beam +z through thin z (5 mm). LG couples to the 40x5 mm x-z face and extends toward -y.
constexpr G4double kHxWide = 20 * mm;
constexpr G4double kHyLong = 30.0 * mm;
constexpr G4double kHzThin = 2.5 * mm;
constexpr G4double kLguide = 30.0 * mm;
/** LG outlet circle radius (kept from LG-v3 baseline). */
constexpr G4double kRguide = 7.5 * mm;
constexpr G4int kNPhi = 256;
constexpr G4int kNSlice = 20;
/** Distinct inlet grease layer: scint ends at -kHyLong; LG starts one layer below. */
constexpr G4double kScintLGGelT = 0.02 * mm;
constexpr G4double kLGZGap = kScintLGGelT;
/** Distinct outlet grease layer, flush with LG (or scint in no-LG mode). */
constexpr G4double kSiPMGelT = 0.10 * mm;
/** Axial half-length of tip Al ring (covers LG|SiPM junction, not just 8 um foil). */
constexpr G4double kTipRingHalfY = 0.30 * mm;

/** Z-offset margin (same role as before outer-air removal): keeps trigger assemblies off z=0 in world coordinates. */
constexpr G4double kOuterAirSafetyMargin = 1.0 * mm;

// Match legacy tower wrapping (CBDsimDetectorConstruction)
constexpr G4double kFoilT = 0.016 * mm;
/** Radial foil sleeve outside the circular sensor stack; does not enter its aperture. */
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

G4TessellatedSolid* BuildLightGuideTessellated(G4double endInset = 0., G4double endRadius = kRguide) {
  auto* ts = new G4TessellatedSolid("ProtoLightGuide");

  // Uniform azimuths alone cut off the rectangular inlet corners. Include the four
  // exact corner directions in EVERY ring so adjacent rings keep matching topology.
  std::vector<G4double> angles;
  for (G4int j = 0; j < kNPhi; ++j) angles.push_back(twopi * j / kNPhi);
  const G4double corner = std::atan2(kHzThin, kHxWide);
  for (const auto angle : {corner, pi-corner, pi+corner, twopi-corner}) angles.push_back(angle);
  std::sort(angles.begin(), angles.end());
  angles.erase(std::unique(angles.begin(), angles.end(),
      [](G4double a, G4double b) { return std::abs(a-b) < 1.e-12; }), angles.end());
  const G4int nPhi = static_cast<G4int>(angles.size());
  const G4int nRings = kNSlice + 1;
  std::vector<G4ThreeVector> v;
  v.reserve(static_cast<size_t>(nRings * nPhi));

  for (G4int i = 0; i < nRings; ++i) {
    const G4double t = static_cast<G4double>(i) / static_cast<G4double>(kNSlice);
    const G4double y = -kHyLong - kLGZGap - t * (kLguide-endInset);
    for (G4int j = 0; j < nPhi; ++j) {
      const G4double phi = angles[j];
      G4double rx, rz;
      rectBoundaryXZ(phi, kHxWide, kHzThin, rx, rz);
      const G4double cx = endRadius * std::cos(phi);
      const G4double cz = endRadius * std::sin(phi);
      const G4double px = (1.0 - t) * rx + t * cx;
      const G4double pz = (1.0 - t) * rz + t * cz;
      v.emplace_back(px, y, pz);
    }
  }

  auto idx = [nPhi](G4int ring, G4int j) { return static_cast<size_t>(ring * nPhi + j); };

  for (G4int i = 0; i < kNSlice; ++i) {
    for (G4int j = 0; j < nPhi; ++j) {
      const G4int jp = (j + 1) % nPhi;
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
  for (G4int j = 0; j < nPhi; ++j) {
    const G4int jp = (j + 1) % nPhi;
    ts->AddFacet(new G4TriangularFacet(cbot, v[idx(0, jp)], v[idx(0, j)], ABSOLUTE));
  }

  const G4ThreeVector ctop(0.0, -kHyLong - kLGZGap - kLguide + endInset, 0.0);
  for (G4int j = 0; j < nPhi; ++j) {
    const G4int jp = (j + 1) % nPhi;
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

  // Fixed world beam, rigid assembly motion toward each local window (-y).
  // Thus the beam crosses local +y=s, approaching the opposite end (+30 mm).
  G4double scanS = 0.;
  if (const char* value = std::getenv("CBDsim_PROTO_SCAN_S_MM")) {
    char* end = nullptr;
    scanS = std::strtod(value, &end) * mm;
    if (end == value || *end != '\0' || !std::isfinite(scanS) || scanS < 0. || scanS > 29.5*mm) {
      G4Exception("CBDsimDetectorConstructionProto::Construct", "InvalidScanPosition", FatalException,
                  "CBDsim_PROTO_SCAN_S_MM must be finite and in [0,29.5] mm.");
    }
  }
  G4cout << "[Proto scan] s_mm=" << scanS/mm
         << " T1_translation=" << G4ThreeVector(0.,-scanS,kZTrig1)/mm
         << " T2_translation=" << G4ThreeVector(scanS,0.,kZTrig2)/mm << G4endl;

  auto trWorld1 = [&](const G4Transform3D& localInAssembly) {
    return G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., -scanS, kZTrig1)) * localInAssembly;
  };
  auto trWorld2 = [&](const G4Transform3D& localInAssembly) {
    return G4Transform3D(rotTrig2, G4ThreeVector(scanS, 0., kZTrig2)) * localInAssembly;
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
    const char* roundSetting = std::getenv("CBDsim_PROTO_ROUND_TIP");
    if (roundSetting && G4String(roundSetting)!="0" && G4String(roundSetting)!="1")
      G4Exception("CBDsimDetectorConstructionProto::Construct", "InvalidRoundTip", FatalException,
                  "CBDsim_PROTO_ROUND_TIP must be 0 or 1.");
    const bool roundTip = !roundSetting || G4String(roundSetting)=="1";
    G4VSolid* lgSolid = BuildLightGuideTessellated(roundTip ? 0.1*mm : 0., roundTip ? kRguide-0.001*mm : kRguide);
    if (roundTip) {
      // Trial: 0.2 mm cylinder, 0.1 mm positive Boolean overlap; total LG remains 30 mm.
      auto* cylinder = new G4Tubs("ProtoRoundTip",0.,kRguide,0.1*mm,0.,twopi);
      G4RotationMatrix rotation; rotation.rotateX(-halfpi);
      lgSolid = new G4UnionSolid("ProtoLightGuideRoundTip",lgSolid,cylinder,
          G4Transform3D(rotation,G4ThreeVector(0.,-kHyLong-kLGZGap-kLguide+0.1*mm,0.)));
    }
    G4cout << "[Proto round tip] active=" << roundTip << " total_length_mm=" << kLguide/mm
           << " cylinder_length_mm=" << (roundTip ? 0.2 : 0.) << " boolean_overlap_mm=" << (roundTip ? 0.1 : 0.) << G4endl;
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
  // Scint -y face is always y=-kHyLong; LG inlet is shifted by kLGZGap (= inlet gel thickness).
  G4RotationMatrix sipmRot;
  sipmRot.rotateX(-halfpi);
  const G4double yScintFace = -kHyLong;
  const G4double yLgInlet = yScintFace - kLGZGap;
  const G4double yCouplingFace = withLG ? (yLgInlet - kLguide) : yScintFace;
  const G4double gelHalfY = kSiPMGelT * 0.5;
  const G4double windowHalfY = (kSiPMH - kFilterT) * 0.5;
  const G4double waferHalfY = kFilterT * 0.5;
  const G4double yGelCenter = yCouplingFace - gelHalfY;
  const G4double yWindowCenter = yGelCenter - gelHalfY - windowHalfY;
  const G4double yWaferCenter = yWindowCenter - windowHalfY - waferHalfY;

  // Absorbing sleeves: no-LG gel + window sides, LG inlet gel sides.
  // Extensions can be disabled to reproduce the previous gel-only model.
  // 50 um thickness, direct contact, R=T=0 are explicit idealizations, not measured tape data.
  const char* tapeSetting = std::getenv("CBDsim_PROTO_GEL_TAPE");
  if (tapeSetting && G4String(tapeSetting) != "0" && G4String(tapeSetting) != "1")
    G4Exception("CBDsimDetectorConstructionProto::Construct", "InvalidGelTape", FatalException,
                "CBDsim_PROTO_GEL_TAPE must be 0 (legacy) or 1 (absorbing tape sleeve).");
  const char* extensionSetting = std::getenv("CBDsim_PROTO_TAPE_EXTENSIONS");
  if (extensionSetting && G4String(extensionSetting)!="0" && G4String(extensionSetting)!="1")
    G4Exception("CBDsimDetectorConstructionProto::Construct", "InvalidTapeExtensions", FatalException,
                "CBDsim_PROTO_TAPE_EXTENSIONS must be 0 or 1.");
  const bool tapeExtensions = !extensionSetting || G4String(extensionSetting)=="1";
  const bool gelTape = (!withLG || tapeExtensions) && (!tapeSetting || G4String(tapeSetting) == "1");
  const G4double tapeHalfY = withLG ? kScintLGGelT*0.5 : gelHalfY + (tapeExtensions ? windowHalfY : 0.);
  const G4double tapeCenterY = yScintFace - tapeHalfY;
  const G4double tapeHalfX = withLG ? kHxWide : kSipmRectHalfX;
  const G4double tapeHalfZ = withLG ? kHzThin : kSipmRectHalfZ;
  const G4double tapeT = 0.05*mm;
  G4SubtractionSolid* gelTapeSolid = nullptr;
  if (gelTape) {
    auto* outer = new G4Box("protoGelTapeOuter", tapeHalfX+tapeT, tapeHalfY, tapeHalfZ+tapeT);
    auto* aperture = new G4Box("protoGelTapeAperture", tapeHalfX, tapeHalfY+1.*mm, tapeHalfZ);
    gelTapeSolid = new G4SubtractionSolid("protoGelTape", outer, aperture);
    auto* tapeLog = new G4LogicalVolume(gelTapeSolid,
        G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE"), "protoGelTapeLog");
    auto* surface = new G4OpticalSurface("protoGelTapeAbsorber", unified, polished, dielectric_metal);
    auto* properties = new G4MaterialPropertiesTable;
    G4double energies[] = {1.*eV, 10.*eV};
    G4double zero[] = {0., 0.};
    properties->AddProperty("REFLECTIVITY", energies, zero, 2);
    properties->AddProperty("TRANSMITTANCE", energies, zero, 2);
    properties->AddProperty("EFFICIENCY", energies, zero, 2);
    surface->SetMaterialPropertiesTable(properties);
    new G4LogicalSkinSurface("protoGelTapeAbsorberSkin", tapeLog, surface);
    const G4Transform3D localTape(G4RotationMatrix(), G4ThreeVector(0., tapeCenterY, 0.));
    new G4PVPlacement(trWorld1(localTape), tapeLog, "protoGelTapePhys", worldLog, false, 0);
    new G4PVPlacement(trWorld2(localTape), tapeLog, "protoGelTapePhys", worldLog, false, 1);
    tapeLog->SetVisAttributes(fVisFoil);
  }
  G4cout << "[Proto gel tape] active=" << gelTape << " extensions=" << tapeExtensions
         << " scope=" << (withLG ? "LG_inlet_gel" : (tapeExtensions ? "noLG_gel_and_window" : "noLG_gel_only"))
         << " thickness_mm=" << tapeT/mm << " direct_contact=1 R=0 T=0 efficiency=0" << G4endl;

  G4cout << "[Proto geometry] coupling face y=" << yCouplingFace / mm << " mm, gel center y="
         << yGelCenter / mm << " mm (gel top y=" << (yGelCenter + gelHalfY) / mm << " mm, overlap "
         << 0.0 << " mm)" << G4endl;
  if (withLG) {
    G4cout << "[Proto geometry] scint-LG inlet gel: T=" << kScintLGGelT / mm
           << " mm, overlap=" << 0.0 << " mm, LG inlet y=" << yLgInlet / mm
           << " mm" << G4endl;
  }

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

  // Inlet optical grease between scint -y and LG (LG mode only).
  G4LogicalVolume* scintLGGelLog = nullptr;
  if (withLG) {
    const G4double inletGelHalfY = 0.5 * kScintLGGelT;
    auto* scintLGGelS = new G4Box("protoScintLGGel", kHxWide, inletGelHalfY, kHzThin);
    scintLGGelLog =
        new G4LogicalVolume(scintLGGelS, FindMaterial("Gelatin"), "protoScintLGGelLog");
    scintLGGelLog->SetVisAttributes(fVisGel);
  }

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
  // Tape occupies this small part of the former foil lip: subtract it rather than overlap.
  if (gelTapeSolid)
    foilWrapS = new G4SubtractionSolid("protoFoilWrapTapeCut", foilWrapS, gelTapeSolid, nullptr,
                                     G4ThreeVector(0., tapeCenterY, 0.));
  auto* foilWrapLog = new G4LogicalVolume(foilWrapS, FindMaterial("Aluminum"), "protoFoilWrapLog");
  foilWrapLog->SetVisAttributes(fVisFoil);
  new G4LogicalSkinSurface("protoAlSurfWrap", foilWrapLog, FindSurface("AluminumSurf"));

  // No corner pads: the LG inlet now covers the exact rectangle. Pads occupied
  // the inlet grease/LG volume (and duplicated the no-LG foil strips).

  // LG tip annulus (LG mode only): blocks radial leak at outlet.
  // Independent sensor radial clearance: contact is the default optical model after validation.
  // Do not change the scintillator foil air gap.
  const char* contactSetting = std::getenv("CBDsim_PROTO_TIP_CONTACT");
  if (contactSetting && G4String(contactSetting)!="0" && G4String(contactSetting)!="1")
    G4Exception("CBDsimDetectorConstructionProto::Construct", "InvalidTipContact", FatalException,
                "CBDsim_PROTO_TIP_CONTACT must be 0 or 1.");
  const bool tipContact = !contactSetting || G4String(contactSetting)=="1";
  const G4double tipRingInnerR = kRguide + (tipContact ? 0. : kAirGap);
  G4cout << "[Proto tip contact] active=" << (withLG && tipContact)
         << " radial_clearance_mm=" << (tipRingInnerR-kRguide)/mm << G4endl;
  const G4double tipRingOuterR = kRguide + kAirGap + kFoilT + kTipRingOuterExtra;
  const G4double yTipRing = yCouplingFace - kTipRingHalfY;
  auto* tipRingS =
      new G4Tubs("protoFoilTipRing", tipRingInnerR, tipRingOuterR, kTipRingHalfY, 0., twopi);
  auto* tipRingLog =
      new G4LogicalVolume(tipRingS, FindMaterial("Aluminum"), "protoFoilTipRingLog");
  tipRingLog->SetVisAttributes(fVisFoil);
  new G4LogicalSkinSurface("protoAlSurfTipRing", tipRingLog, FindSurface("AluminumSurf"));
  G4RotationMatrix tipRingRot;
  tipRingRot.rotateX(-halfpi);

  // no-LG: Al strips on scint -y face beside rectangular SiPM (v3 idea, v5 Al surface).
  G4VSolid* foilYmStripSolid = nullptr;
  G4LogicalVolume* foilYmStripLog = nullptr;
  G4double xLeftYm = 0.;
  G4double xRightYm = 0.;
  if (!withLG) {
    const G4double stripInnerX = kSipmRectHalfX + (gelTape ? tapeT : 0.);
    const G4double foilYmStripHalfX = 0.5 * (kHxWide - stripInnerX);
    auto* foilYmStripS =
        new G4Box("protoFoilYmStrip", foilYmStripHalfX, ft2, zi);
    foilYmStripSolid = foilYmStripS;
    foilYmStripLog =
        new G4LogicalVolume(foilYmStripS, FindMaterial("Aluminum"), "protoFoilYmStripLog");
    foilYmStripLog->SetVisAttributes(fVisFoil);
    new G4LogicalSkinSurface("protoAlSurfYmStrip", foilYmStripLog, FindSurface("AluminumSurf"));
    xLeftYm = -0.5 * (kHxWide + stripInnerX);
    xRightYm = 0.5 * (kHxWide + stripInnerX);
  }

  // Folded Al rim closes only the open tile/foil air-gap perimeter.
  // In LG mode it stops exactly at the scint face: no reflector is added to inlet gel sides.
  const char* sealSetting = std::getenv("CBDsim_PROTO_CORNER_SEAL");
  if (sealSetting && G4String(sealSetting)!="0" && G4String(sealSetting)!="1")
    G4Exception("CBDsimDetectorConstructionProto::Construct", "InvalidCornerSeal", FatalException,
                "CBDsim_PROTO_CORNER_SEAL must be 0 or 1.");
  const bool cornerSeal = !sealSetting || G4String(sealSetting)=="1";
  if (cornerSeal) {
    const G4double sealTop = -kHyLong + 0.02*mm;
    const G4double sealBottom = withLG ? -kHyLong : -yo;
    const G4double sealY = 0.5*(sealTop+sealBottom);
    const G4double sealHalfY = 0.5*(sealTop-sealBottom);
    auto* rimOuter = new G4Box("protoCornerSealOuter", xo, sealHalfY, zo);
    auto* rimAperture = new G4Box("protoCornerSealAperture", kHxWide, sealHalfY+1.*mm, kHzThin);
    G4VSolid* rim = new G4SubtractionSolid("protoCornerSealRing", rimOuter, rimAperture);
    // Existing volumes keep their ownership; the new rim fills only the missing material.
    rim = new G4SubtractionSolid("protoCornerSealWrapCut", rim, foilWrapS, nullptr,
                                G4ThreeVector(0.,-sealY,0.));
    if (gelTapeSolid)
      rim = new G4SubtractionSolid("protoCornerSealTapeCut", rim, gelTapeSolid, nullptr,
                                  G4ThreeVector(0.,tapeCenterY-sealY,0.));
    if (foilYmStripSolid) {
      const G4double foilY = -kHyLong-kAirGap-ft2;
      rim = new G4SubtractionSolid("protoCornerSealLeftCut", rim, foilYmStripSolid, nullptr,
                                  G4ThreeVector(xLeftYm,foilY-sealY,0.));
      rim = new G4SubtractionSolid("protoCornerSealRightCut", rim, foilYmStripSolid, nullptr,
                                  G4ThreeVector(xRightYm,foilY-sealY,0.));
    }
    auto* rimLog = new G4LogicalVolume(rim, FindMaterial("Aluminum"), "protoFoilCornerSealLog");
    new G4LogicalSkinSurface("protoAlSurfCornerSeal", rimLog, FindSurface("AluminumSurf"));
    const G4Transform3D localRim(G4RotationMatrix(),G4ThreeVector(0.,sealY,0.));
    new G4PVPlacement(trWorld1(localRim),rimLog,"protoFoilCornerSealPhys",worldLog,false,0);
    new G4PVPlacement(trWorld2(localRim),rimLog,"protoFoilCornerSealPhys",worldLog,false,1);
    rimLog->SetVisAttributes(fVisFoil);
  }
  G4cout << "[Proto corner seal] active=" << cornerSeal
         << " surface=AluminumSurf scope=tile_foil_perimeter" << G4endl;

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

    G4PVPlacement* scintLGGelPV1 = nullptr;
    G4PVPlacement* scintLGGelPV2 = nullptr;
    if (withLG && scintLGGelLog) {
      const G4double yInletGel = -kHyLong - 0.5 * kScintLGGelT;
      const G4Transform3D localInletGel(G4RotationMatrix(), G4ThreeVector(0., yInletGel, 0.));
      scintLGGelPV1 = new G4PVPlacement(trWorld1(localInletGel), scintLGGelLog, "protoScintLGGelPhys",
                                        worldLog, false, 0);
      scintLGGelPV2 = new G4PVPlacement(trWorld2(localInletGel), scintLGGelLog, "protoScintLGGelPhys",
                                        worldLog, false, 1);
    }

    if (withLG) {
      new G4LogicalBorderSurface("protoLG1ToWorldReflect", lgPV1, worldPhys, FindSurface("AluminumSurf"));
      new G4LogicalBorderSurface("protoWorldToLG1Reflect", worldPhys, lgPV1, FindSurface("AluminumSurf"));
      new G4LogicalBorderSurface("protoLG2ToWorldReflect", lgPV2, worldPhys, FindSurface("AluminumSurf"));
      new G4LogicalBorderSurface("protoWorldToLG2Reflect", worldPhys, lgPV2, FindSurface("AluminumSurf"));

      new G4LogicalBorderSurface("protoLG1ToGel1Trans", lgPV1, sipmGelPV1, FindSurface("AirSurf"));
      new G4LogicalBorderSurface("protoGel1ToLG1Trans", sipmGelPV1, lgPV1, FindSurface("AirSurf"));
      new G4LogicalBorderSurface("protoLG2ToGel2Trans", lgPV2, sipmGelPV2, FindSurface("AirSurf"));
      new G4LogicalBorderSurface("protoGel2ToLG2Trans", sipmGelPV2, lgPV2, FindSurface("AirSurf"));

      if (scintLGGelPV1 && scintLGGelPV2) {
        new G4LogicalBorderSurface("protoScint1ToInletGel1Trans", scintPV1, scintLGGelPV1,
                                   FindSurface("AirSurf"));
        new G4LogicalBorderSurface("protoInletGel1ToScint1Trans", scintLGGelPV1, scintPV1,
                                   FindSurface("AirSurf"));
        new G4LogicalBorderSurface("protoScint2ToInletGel2Trans", scintPV2, scintLGGelPV2,
                                   FindSurface("AirSurf"));
        new G4LogicalBorderSurface("protoInletGel2ToScint2Trans", scintLGGelPV2, scintPV2,
                                   FindSurface("AirSurf"));
        new G4LogicalBorderSurface("protoLG1ToInletGel1Trans", lgPV1, scintLGGelPV1,
                                   FindSurface("AirSurf"));
        new G4LogicalBorderSurface("protoInletGel1ToLG1Trans", scintLGGelPV1, lgPV1,
                                   FindSurface("AirSurf"));
        new G4LogicalBorderSurface("protoLG2ToInletGel2Trans", lgPV2, scintLGGelPV2,
                                   FindSurface("AirSurf"));
        new G4LogicalBorderSurface("protoInletGel2ToLG2Trans", scintLGGelPV2, lgPV2,
                                   FindSurface("AirSurf"));
      }
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

  // LG|world retains the existing reflective side surface; the sleeve is outside the sensor stack.

  return worldPhys;
}

void CBDsimDetectorConstructionProto::ConstructSDandField() {
  // Checking only worldPhys does not check its daughters. Full placement checks
  // and coupling probes are provided by tests/proto_geometry (tolerance = 0).
  if (!fProtoWaferLog) return;
  auto* SDman = G4SDManager::GetSDMpointer();
  auto* sipmSD = new CBDsimSiPMSD("SiPMSDB", "SiPMSDBC", std::make_pair(1, 1));
  SDman->AddNewDetector(sipmSD);
  fProtoWaferLog->SetSensitiveDetector(sipmSD);
}
