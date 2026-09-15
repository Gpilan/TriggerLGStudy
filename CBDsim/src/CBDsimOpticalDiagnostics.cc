#include "CBDsimOpticalDiagnostics.hh"

#include "G4AutoLock.hh"
#include "G4Threading.hh"
#include "G4Material.hh"
#include "G4Exception.hh"
#include <sstream>
#include <vector>
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <map>
#include <unordered_map>
#include <unordered_set>

namespace {
enum class TrackFate : int {
  kUnknown = 0,
  kDetected,
  kQeReject,
  kAbsorbScint,
  kAbsorbLg,
  kAbsorbLgSurface,  // OpBoundary Absorption on LG (was lumped into kill_other)
  kAbsorbGlass,
  kAbsorbSi,
  kAbsorbOther,
  kKillWorld,
  kKillNav,  // navigation/stuck heuristic (split from kill_other / some kill_world)
  kKillOther,
  kAbsorbOtherSurface,
  kWorldExit,
  kWls,
  kNoRindex,
};

constexpr int kNFateHists = 15;
constexpr int kNBins = 400;
constexpr G4double kXmaxMm = 2000.0;
/** Step length below this (mm) counts toward a tiny-step streak. */
constexpr G4double kTinyStepMm = 1.0e-6;
/** Consecutive tiny steps (or max streak) at/above this → kill_nav candidate. */
constexpr G4int kTinyStreakThresh = 3;

constexpr std::array<const char*, kNFateHists> kFateNames = {
    "detected",          "qe_reject",    "absorb_scint", "absorb_lg",
    "absorb_lg_surface", "absorb_glass", "absorb_si",    "absorb_other",
    "legacy_kill_world", "legacy_kill_nav", "unknown",
    "absorb_other_surface", "world_exit", "wls_conversion", "no_rindex",
};

struct PathHist {
  G4long totalBins[kNBins] = {};
  G4long totalUnder = 0;
  G4long totalOver = 0;
  G4long fateBins[kNFateHists][kNBins] = {};
  G4long fateUnder[kNFateHists] = {};
  G4long fateOver[kNFateHists] = {};
  G4long fateCounts[kNFateHists] = {};
};

/** Last-volume region counts for one kill cause. */
struct KillLoc {
  G4long scint = 0;
  G4long lg = 0;
  G4long gel = 0;
  G4long inletGel = 0;
  G4long outletGel = 0;
  G4long world = 0;
  G4long sipm_env = 0;
  G4long foil = 0;
  G4long sipm = 0;  // wafer + glass
  G4long other = 0;
};

struct TrackStepState {
  G4int tinyStreak = 0;
  G4int maxTinyStreak = 0;
  G4double lastStepMm = 0.;
  G4String lastProc;
};

struct BudgetTrack {
  G4String source;
  G4String terminal;
  G4String region = "none";
  G4String process = "none";
  G4int parent = 0;
  G4int status = 0;
  G4bool tiny = false;
  G4String preVolume = "none", postVolume = "none", postMaterial = "none";
  G4int boundaryStatus = -1;
  G4ThreeVector position;

};
struct EventBudget {
  G4int id = -1;
  std::unordered_map<G4int, BudgetTrack> tracks;
};
G4ThreadLocal EventBudget* eventBudget = nullptr;
EventBudget& Budget() {
  if (!eventBudget) eventBudget = new EventBudget;
  return *eventBudget;
}
struct Counters {
  std::map<G4String, G4long> budgetCounts;
  std::vector<G4String> eventRows;
  std::vector<G4String> exceptionRows;
  std::map<G4int, G4long> boundaryStatuses;
  std::map<G4String, G4long> encounters;
  std::map<G4String, G4long> transmissions;
  G4long windowEntrySteps = 0;
  G4long uniqueWindowArrivals = 0;
  G4double pathInletGelMm = 0.;
  G4double pathOutletGelMm = 0.;
  G4long nBoundarySteps = 0;
  G4long nScintToLg = 0;
  G4long nLgToScint = 0;
  G4long nScintToWorld = 0;
  G4long nWorldToScint = 0;
  G4long nLgToWorld = 0;
  G4long nWorldToLg = 0;
  G4long nLgToSipm = 0;
  G4long nSipmToLg = 0;
  G4long nOtherBoundary = 0;
  G4long nBulkAbsorbScint = 0;
  G4long nBulkAbsorbLg = 0;
  G4long nBulkAbsorbGlass = 0;
  G4long nBulkAbsorbSi = 0;
  G4long nBulkAbsorbOther = 0;
  G4long nBoundaryAbsorbLg = 0;     // OpBoundary Absorption, pre-vol = LG
  G4long nBoundaryAbsorbOther = 0;  // OpBoundary Absorption, other pre-vol
  G4long nKillOtherResidual = 0;    // still kill_other after nav/surface reclass
  G4long nKillNav = 0;
  G4long nTracksKilled = 0;
  G4long nKillLastScint = 0;
  G4long nKillLastLg = 0;
  G4long nKillLastWorld = 0;
  G4long nKillLastSipm = 0;
  G4long nKillLastOther = 0;
  KillLoc killNavLoc{};
  KillLoc killWorldLoc{};
  KillLoc killOtherLoc{};
  G4long nSipmDetect = 0;
  G4long nSipmQeReject = 0;
  G4double pathScintMm = 0.;
  G4double pathLgMm = 0.;
  G4double pathGlassMm = 0.;
  G4double pathSiMm = 0.;
  G4double pathWorldMm = 0.;
  PathHist pathHist{};
};

G4Mutex diagMutex = G4MUTEX_INITIALIZER;
Counters master;
G4ThreadLocal Counters* worker = nullptr;
G4ThreadLocal std::unordered_map<G4int, TrackFate>* trackFates = nullptr;
G4ThreadLocal std::unordered_map<G4int, TrackStepState>* trackSteps = nullptr;
G4ThreadLocal std::unordered_set<G4int>* windowArrivals = nullptr;

std::unordered_set<G4int>& WindowArrivals() {
  if (!windowArrivals) windowArrivals = new std::unordered_set<G4int>();
  return *windowArrivals;
}

G4bool EnvEnabled() {
  const char* v = std::getenv("CBDsim_OPTICAL_DIAG");
  return v && (v[0] == '1' || v[0] == 'y' || v[0] == 'Y' || v[0] == 't' || v[0] == 'T');
}

G4bool DiagHistEnabled() {
  if (!EnvEnabled()) return false;
  const char* v = std::getenv("CBDsim_OPTICAL_DIAG_HIST");
  return v && (v[0] == '1' || v[0] == 'y' || v[0] == 'Y' || v[0] == 't' || v[0] == 'T');
}

Counters& Worker() {
  if (!worker) worker = new Counters();
  return *worker;
}

std::unordered_map<G4int, TrackFate>& Fates() {
  if (!trackFates) trackFates = new std::unordered_map<G4int, TrackFate>();
  return *trackFates;
}

std::unordered_map<G4int, TrackStepState>& Steps() {
  if (!trackSteps) trackSteps = new std::unordered_map<G4int, TrackStepState>();
  return *trackSteps;
}

void AddKillLoc(KillLoc& a, const KillLoc& b) {
  a.scint += b.scint;
  a.lg += b.lg;
  a.gel += b.gel;
  a.inletGel += b.inletGel;
  a.outletGel += b.outletGel;
  a.world += b.world;
  a.sipm_env += b.sipm_env;
  a.foil += b.foil;
  a.sipm += b.sipm;
  a.other += b.other;
}

void AddPathHist(PathHist& a, const PathHist& b) {
  a.totalUnder += b.totalUnder;
  a.totalOver += b.totalOver;
  for (int i = 0; i < kNBins; ++i) a.totalBins[i] += b.totalBins[i];
  for (int f = 0; f < kNFateHists; ++f) {
    a.fateUnder[f] += b.fateUnder[f];
    a.fateOver[f] += b.fateOver[f];
    a.fateCounts[f] += b.fateCounts[f];
    for (int i = 0; i < kNBins; ++i) a.fateBins[f][i] += b.fateBins[f][i];
  }
}

void Add(Counters& a, const Counters& b) {
  for (const auto& e : b.budgetCounts) a.budgetCounts[e.first] += e.second;
  a.eventRows.insert(a.eventRows.end(), b.eventRows.begin(), b.eventRows.end());
  a.exceptionRows.insert(a.exceptionRows.end(), b.exceptionRows.begin(), b.exceptionRows.end());
  for (const auto& entry : b.boundaryStatuses) a.boundaryStatuses[entry.first] += entry.second;
  for (const auto& entry : b.encounters) a.encounters[entry.first] += entry.second;
  for (const auto& entry : b.transmissions) a.transmissions[entry.first] += entry.second;
  a.windowEntrySteps += b.windowEntrySteps;
  a.uniqueWindowArrivals += b.uniqueWindowArrivals;
  a.pathInletGelMm += b.pathInletGelMm;
  a.pathOutletGelMm += b.pathOutletGelMm;
  a.nBoundarySteps += b.nBoundarySteps;
  a.nScintToLg += b.nScintToLg;
  a.nLgToScint += b.nLgToScint;
  a.nScintToWorld += b.nScintToWorld;
  a.nWorldToScint += b.nWorldToScint;
  a.nLgToWorld += b.nLgToWorld;
  a.nWorldToLg += b.nWorldToLg;
  a.nLgToSipm += b.nLgToSipm;
  a.nSipmToLg += b.nSipmToLg;
  a.nOtherBoundary += b.nOtherBoundary;
  a.nBulkAbsorbScint += b.nBulkAbsorbScint;
  a.nBulkAbsorbLg += b.nBulkAbsorbLg;
  a.nBulkAbsorbGlass += b.nBulkAbsorbGlass;
  a.nBulkAbsorbSi += b.nBulkAbsorbSi;
  a.nBulkAbsorbOther += b.nBulkAbsorbOther;
  a.nBoundaryAbsorbLg += b.nBoundaryAbsorbLg;
  a.nBoundaryAbsorbOther += b.nBoundaryAbsorbOther;
  a.nKillOtherResidual += b.nKillOtherResidual;
  a.nKillNav += b.nKillNav;
  a.nTracksKilled += b.nTracksKilled;
  a.nKillLastScint += b.nKillLastScint;
  a.nKillLastLg += b.nKillLastLg;
  a.nKillLastWorld += b.nKillLastWorld;
  a.nKillLastSipm += b.nKillLastSipm;
  a.nKillLastOther += b.nKillLastOther;
  AddKillLoc(a.killNavLoc, b.killNavLoc);
  AddKillLoc(a.killWorldLoc, b.killWorldLoc);
  AddKillLoc(a.killOtherLoc, b.killOtherLoc);
  a.nSipmDetect += b.nSipmDetect;
  a.nSipmQeReject += b.nSipmQeReject;
  a.pathScintMm += b.pathScintMm;
  a.pathLgMm += b.pathLgMm;
  a.pathGlassMm += b.pathGlassMm;
  a.pathSiMm += b.pathSiMm;
  a.pathWorldMm += b.pathWorldMm;
  AddPathHist(a.pathHist, b.pathHist);
}

G4String VolRegion(const G4String& volName) {
  // Specific names precede the scint prefix: inlet grease is not scintillator.
  if (volName.find("protoScintLGGel") == 0) return "inlet_gel";
  if (volName.find("protoSipmGel") == 0) return "outlet_gel";
  if (volName.find("protoScint") != G4String::npos) return "scint";
  if (volName.find("protoLightGuide") != G4String::npos) return "lg";
  if (volName.find("protoSipmWafer") != G4String::npos) return "sipm_wafer";
  if (volName.find("protoSipmWindow") != G4String::npos) return "sipm_glass";
  if (volName.find("protoSipmEnv") != G4String::npos) return "sipm_env";
  if (volName.find("protoFoil") != G4String::npos) return "foil";
  if (volName.find("protoWorld") != G4String::npos || volName == "worldPhysical" ||
      volName == "worldLogical")
    return "world";
  return "other";
}

void BumpKillLoc(KillLoc& loc, const G4String& region) {
  if (region == "scint") loc.scint++;
  else if (region == "lg") loc.lg++;
  else if (region == "inlet_gel") { loc.gel++; loc.inletGel++; }
  else if (region == "outlet_gel") { loc.gel++; loc.outletGel++; }
  else if (region == "world") loc.world++;
  else if (region == "sipm_env") loc.sipm_env++;
  else if (region == "foil") loc.foil++;
  else if (region == "sipm_wafer" || region == "sipm_glass") loc.sipm++;
  else loc.other++;
}

void RecordPath(Counters& c, const G4String& region, G4double mm) {
  if (region == "scint") c.pathScintMm += mm;
  else if (region == "lg") c.pathLgMm += mm;
  else if (region == "sipm_glass") c.pathGlassMm += mm;
  else if (region == "sipm_wafer") c.pathSiMm += mm;
  else if (region == "inlet_gel") c.pathInletGelMm += mm;
  else if (region == "outlet_gel") c.pathOutletGelMm += mm;
  else if (region == "world" || region == "sipm_env") c.pathWorldMm += mm;
}

void RecordBoundary(Counters& c, const G4String& pre, const G4String& post) {
  c.nBoundarySteps++;
  const G4String a = VolRegion(pre);
  const G4String b = VolRegion(post);
  ++c.encounters[a + "_to_" + b];
  if (a == "scint" && b == "lg") c.nScintToLg++;
  else if (a == "lg" && b == "scint") c.nLgToScint++;
  else if (a == "scint" && b == "world") c.nScintToWorld++;
  else if (a == "world" && b == "scint") c.nWorldToScint++;
  else if (a == "lg" && b == "world") c.nLgToWorld++;
  else if (a == "world" && b == "lg") c.nWorldToLg++;
  // Legacy aggregate retained for old readers; NOT a unique photon arrival count.
  else if (a == "lg" && (b == "sipm_env" || b == "sipm_glass" || b == "outlet_gel")) c.nLgToSipm++;
  else if ((a == "sipm_env" || a == "sipm_glass" || a == "outlet_gel") && b == "lg") c.nSipmToLg++;
  else if (a == "outlet_gel" && b == "sipm_glass") c.nLgToSipm++;
  else if (a == "sipm_glass" && b == "outlet_gel") c.nSipmToLg++;
  else c.nOtherBoundary++;
}

void RecordKillRegion(Counters& c, const G4String& volName) {
  c.nTracksKilled++;
  const G4String r = VolRegion(volName);
  if (r == "scint") c.nKillLastScint++;
  else if (r == "lg") c.nKillLastLg++;
  else if (r == "world" || r == "sipm_env") c.nKillLastWorld++;
  else if (r == "sipm_wafer" || r == "sipm_glass" || r == "outlet_gel") c.nKillLastSipm++;
  else c.nKillLastOther++;
}

void UpdateTinyStep(TrackStepState& st, G4double stepMm, const G4VProcess* proc) {
  st.lastStepMm = stepMm;
  if (proc) st.lastProc = proc->GetProcessName();
  else st.lastProc = "";
  if (stepMm < kTinyStepMm) {
    st.tinyStreak++;
    st.maxTinyStreak = std::max(st.maxTinyStreak, st.tinyStreak);
  } else {
    st.tinyStreak = 0;
  }
}

G4bool IsNavKillCandidate(const TrackStepState& st) {
  if (st.maxTinyStreak >= kTinyStreakThresh) return true;
  if (st.lastStepMm < kTinyStepMm && st.lastProc == "Transportation") return true;
  return false;
}

int FateHistIndex(TrackFate fate) {
  switch (fate) {
    case TrackFate::kDetected: return 0;
    case TrackFate::kQeReject: return 1;
    case TrackFate::kAbsorbScint: return 2;
    case TrackFate::kAbsorbLg: return 3;
    case TrackFate::kAbsorbLgSurface: return 4;
    case TrackFate::kAbsorbGlass: return 5;
    case TrackFate::kAbsorbSi: return 6;
    case TrackFate::kAbsorbOther: return 7;
    case TrackFate::kKillWorld: return 8;
    case TrackFate::kKillNav: return 9;
    case TrackFate::kKillOther: return 10;
    case TrackFate::kAbsorbOtherSurface: return 11;
    case TrackFate::kWorldExit: return 12;
    case TrackFate::kWls: return 13;
    case TrackFate::kNoRindex: return 14;
    default: return 10;
  }
}

G4OpBoundaryProcess* FindOpBoundaryProcess(const G4Track* track) {
  G4ProcessManager* pm = track->GetDefinition()->GetProcessManager();
  if (!pm) return nullptr;
  G4ProcessVector* pv = pm->GetPostStepProcessVector();
  if (!pv) return nullptr;
  for (G4int i = 0; i < pv->entries(); ++i) {
    G4VProcess* p = (*pv)[i];
    if (p && p->GetProcessName() == "OpBoundary") {
      return dynamic_cast<G4OpBoundaryProcess*>(p);
    }
  }
  return nullptr;
}

void WriteKillLoc(std::ostream& out, const char* prefix, const KillLoc& loc) {
  out << prefix << "_scint " << loc.scint << "\n";
  out << prefix << "_lg " << loc.lg << "\n";
  out << prefix << "_gel " << loc.gel << "\n";
  out << prefix << "_inlet_gel " << loc.inletGel << "\n";
  out << prefix << "_outlet_gel " << loc.outletGel << "\n";
  out << prefix << "_world " << loc.world << "\n";
  out << prefix << "_sipm_env " << loc.sipm_env << "\n";
  out << prefix << "_foil " << loc.foil << "\n";
  out << prefix << "_sipm " << loc.sipm << "\n";
  out << prefix << "_other " << loc.other << "\n";
}

void PrintKillLoc(const char* label, const KillLoc& loc) {
  G4cout << "    " << label << ": scint=" << loc.scint << " lg=" << loc.lg << " gel=" << loc.gel
         << " world=" << loc.world << " sipm_env=" << loc.sipm_env << " foil=" << loc.foil
         << " sipm=" << loc.sipm << " other=" << loc.other << G4endl;
}

void FillPathHist(PathHist& h, TrackFate fate, G4double mm) {
  const G4double x = std::max(0.0, mm);
  int bin = static_cast<int>(x / kXmaxMm * static_cast<G4double>(kNBins));
  if (bin < 0) {
  } else if (bin >= kNBins) {
    h.totalOver++;
  } else {
    h.totalBins[bin]++;
  }
  if (x < 0.) h.totalUnder++;

  const int fi = FateHistIndex(fate);
  h.fateCounts[fi]++;
  if (x < 0.) {
    h.fateUnder[fi]++;
    return;
  }
  if (bin >= kNBins) {
    h.fateOver[fi]++;
    return;
  }
  h.fateBins[fi][bin]++;
}

G4double MeanFromHist(const G4long* bins, G4long under, G4long over, G4long count) {
  if (count <= 0) return 0.;
  const G4double w = kXmaxMm / static_cast<G4double>(kNBins);
  G4double sum = 0.;
  for (int i = 0; i < kNBins; ++i) sum += static_cast<G4double>(bins[i]) * (i + 0.5) * w;
  sum += static_cast<G4double>(over) * (kXmaxMm + 0.5 * w);
  sum += static_cast<G4double>(under) * (0.5 * w);
  return sum / static_cast<G4double>(count);
}

void PrintCounters(const Counters& c, G4int nEvents) {
  const G4double inv = nEvents > 0 ? 1.0 / static_cast<G4double>(nEvents) : 0.0;
  G4cout << "\n=== [CBDsim optical diagnostics] events=" << nEvents << " ===" << G4endl;
  for (const auto& e : c.budgetCounts) G4cout << "  budget_" << e.first << ": " << e.second << G4endl;
  G4cout << "  boundary steps (total): " << c.nBoundarySteps << G4endl;
  G4cout << "    scint->LG: " << c.nScintToLg << "  LG->scint: " << c.nLgToScint << G4endl;
  G4cout << "    scint->world: " << c.nScintToWorld << "  world->scint: " << c.nWorldToScint << G4endl;
  G4cout << "    LG->world: " << c.nLgToWorld << "  world->LG: " << c.nWorldToLg << G4endl;
  G4cout << "    legacy coupling encounters (not arrivals): " << c.nLgToSipm
         << "  reverse: " << c.nSipmToLg << G4endl;
  G4cout << "    window coupling entries: " << c.windowEntrySteps
         << "  unique photons: " << c.uniqueWindowArrivals << G4endl;
  G4cout << "    other boundary: " << c.nOtherBoundary << G4endl;
  G4cout << "  bulk OpAbsorption: scint=" << c.nBulkAbsorbScint << " lg=" << c.nBulkAbsorbLg
         << " glass=" << c.nBulkAbsorbGlass << " si=" << c.nBulkAbsorbSi
         << " other=" << c.nBulkAbsorbOther << G4endl;
  G4cout << "  OpBoundary Absorption: lg=" << c.nBoundaryAbsorbLg
         << " other=" << c.nBoundaryAbsorbOther << G4endl;
  G4cout << "  tracks killed: " << c.nTracksKilled
         << " (last vol: scint=" << c.nKillLastScint << " lg=" << c.nKillLastLg
         << " world/env=" << c.nKillLastWorld << " sipm=" << c.nKillLastSipm
         << " other=" << c.nKillLastOther << ")" << G4endl;
  if (DiagHistEnabled()) {
    G4cout << "  kill_nav: " << c.nKillNav << "  kill_other residual: " << c.nKillOtherResidual
           << G4endl;
    G4cout << "  kill cause x last region:" << G4endl;
    PrintKillLoc("kill_nav", c.killNavLoc);
    PrintKillLoc("kill_world", c.killWorldLoc);
    PrintKillLoc("kill_other", c.killOtherLoc);
  }
  G4cout << "  SiPM detect steps: " << c.nSipmDetect << "  QE reject: " << c.nSipmQeReject << G4endl;
  G4cout << "  path length sum (mm/event): scint=" << c.pathScintMm * inv << " lg=" << c.pathLgMm * inv
         << " glass=" << c.pathGlassMm * inv << " si=" << c.pathSiMm * inv
         << " world/env=" << c.pathWorldMm * inv << G4endl;
  if (DiagHistEnabled()) {
    G4long totalPhotons = 0;
    for (int f = 0; f < kNFateHists; ++f) totalPhotons += c.pathHist.fateCounts[f];
    G4cout << "  optical photon fate counts (path hist): total_tracks=" << totalPhotons << G4endl;
    for (int f = 0; f < kNFateHists; ++f) {
      const G4long n = c.pathHist.fateCounts[f];
      if (n == 0) continue;
      const G4double frac = totalPhotons > 0 ? 100.0 * static_cast<G4double>(n) / static_cast<G4double>(totalPhotons)
                                             : 0.0;
      const G4double mean =
          MeanFromHist(c.pathHist.fateBins[f], c.pathHist.fateUnder[f], c.pathHist.fateOver[f], n);
      G4cout << "    " << kFateNames[f] << ": n=" << n << " (" << frac << "%) mean_path_mm=" << mean
             << G4endl;
    }
    const G4long nt = totalPhotons;
    const G4double meanAll =
        MeanFromHist(c.pathHist.totalBins, c.pathHist.totalUnder, c.pathHist.totalOver, nt);
    G4cout << "    all: mean_path_mm=" << meanAll << G4endl;
  }
  G4cout << "=== [end optical diagnostics] ===\n" << G4endl;
}

void WriteCountersFile(const Counters& c, const G4String& path, G4int nEvents) {
  const G4double inv = nEvents > 0 ? 1.0 / static_cast<G4double>(nEvents) : 0.0;
  std::ofstream out(path.c_str());
  if (!out) return;
  out << "events " << nEvents << "\n";
  out << "diagnostics_schema_version 3\n";
  auto count = [&](const char* key) { auto it = c.budgetCounts.find(key); return it == c.budgetCounts.end() ? 0L : it->second; };
  out << "budget_missing_event_records " << nEvents - count("events") << "\n";
  out << "budget_closed " << (count("events") == nEvents && count("aborted_events") == 0 &&
      count("fate_unfinished") == 0 && count("balance_delta") == 0 &&
      count("terminal_without_start") == 0 && count("duplicate_terminal_callbacks") == 0) << "\n";
  for (const auto& e : c.budgetCounts) out << "budget_" << e.first << " " << e.second << "\n";
  out << "window_entry_steps " << c.windowEntrySteps << "\n";
  out << "unique_window_arrivals " << c.uniqueWindowArrivals << "\n";
  for (const auto& entry : c.encounters)
    out << "boundary_encounter_" << entry.first << " " << entry.second << "\n";
  for (const auto& entry : c.transmissions)
    out << "boundary_transmit_" << entry.first << " " << entry.second << "\n";
  out << "boundary_total " << c.nBoundarySteps << "\n";
  // Numeric status values correspond to G4OpBoundaryProcessStatus in this Geant4 build.
  for (const auto& entry : c.boundaryStatuses)
    out << "boundary_status_" << entry.first << " " << entry.second << "\n";
  out << "scint_to_lg " << c.nScintToLg << "\n";
  out << "lg_to_scint " << c.nLgToScint << "\n";
  out << "scint_to_world " << c.nScintToWorld << "\n";
  out << "world_to_scint " << c.nWorldToScint << "\n";
  out << "lg_to_world " << c.nLgToWorld << "\n";
  out << "world_to_lg " << c.nWorldToLg << "\n";
  out << "lg_to_sipm " << c.nLgToSipm << "\n";
  out << "sipm_to_lg " << c.nSipmToLg << "\n";
  out << "other_boundary " << c.nOtherBoundary << "\n";
  out << "bulk_absorb_scint " << c.nBulkAbsorbScint << "\n";
  out << "bulk_absorb_lg " << c.nBulkAbsorbLg << "\n";
  out << "bulk_absorb_glass " << c.nBulkAbsorbGlass << "\n";
  out << "bulk_absorb_si " << c.nBulkAbsorbSi << "\n";
  out << "boundary_absorb_lg " << c.nBoundaryAbsorbLg << "\n";
  out << "boundary_absorb_other " << c.nBoundaryAbsorbOther << "\n";
  out << "kill_nav " << c.nKillNav << "\n";
  out << "kill_other_residual " << c.nKillOtherResidual << "\n";
  WriteKillLoc(out, "kill_nav", c.killNavLoc);
  WriteKillLoc(out, "kill_world", c.killWorldLoc);
  WriteKillLoc(out, "kill_other", c.killOtherLoc);
  out << "sipm_detect " << c.nSipmDetect << "\n";
  out << "sipm_qe_reject " << c.nSipmQeReject << "\n";
  out << "path_scint_mm_per_event " << c.pathScintMm * inv << "\n";
  out << "path_lg_mm_per_event " << c.pathLgMm * inv << "\n";
  out << "path_inlet_gel_mm_per_event " << c.pathInletGelMm * inv << "\n";
  out << "path_outlet_gel_mm_per_event " << c.pathOutletGelMm * inv << "\n";
  out << "path_glass_mm_per_event " << c.pathGlassMm * inv << "\n";
  out << "path_si_mm_per_event " << c.pathSiMm * inv << "\n";
  out << "path_world_mm_per_event " << c.pathWorldMm * inv << "\n";
}

void WritePathHistToFile(const Counters& c, const G4String& path) {
  std::ofstream out(path.c_str());
  if (!out) return;
  out << "# CBDsim optical photon total track-length histogram\n";
  out << "nbins " << kNBins << "\n";
  out << "xmax_mm " << kXmaxMm << "\n";
  out << "diagnostics_schema_version 3\n";
  out << "nfates " << kNFateHists << "\n";
  G4long totalPhotons = 0;
  for (int f = 0; f < kNFateHists; ++f) totalPhotons += c.pathHist.fateCounts[f];
  out << "total_tracks " << totalPhotons << "\n";
  const G4double meanAll =
      MeanFromHist(c.pathHist.totalBins, c.pathHist.totalUnder, c.pathHist.totalOver, totalPhotons);
  out << "mean_path_mm_all " << meanAll << "\n";
  out << "hist total under " << c.pathHist.totalUnder << " over " << c.pathHist.totalOver << "\n";
  out << "counts";
  for (int i = 0; i < kNBins; ++i) out << " " << c.pathHist.totalBins[i];
  out << "\n";
  for (int f = 0; f < kNFateHists; ++f) {
    const G4long n = c.pathHist.fateCounts[f];
    const G4double mean =
        MeanFromHist(c.pathHist.fateBins[f], c.pathHist.fateUnder[f], c.pathHist.fateOver[f], n);
    const G4double frac =
        totalPhotons > 0 ? static_cast<G4double>(n) / static_cast<G4double>(totalPhotons) : 0.0;
    out << "fate " << kFateNames[f] << " count " << n << " frac " << frac << " mean_path_mm " << mean
        << " under " << c.pathHist.fateUnder[f] << " over " << c.pathHist.fateOver[f] << "\n";
    out << "counts";
    for (int i = 0; i < kNBins; ++i) out << " " << c.pathHist.fateBins[f][i];
    out << "\n";
  }
}
}  // namespace

bool CBDsimOpticalDiagnostics::Enabled() { return EnvEnabled(); }

bool CBDsimOpticalDiagnostics::HistEnabled() { return DiagHistEnabled(); }

G4String CBDsimOpticalDiagnostics::RegionForVolume(const G4String& name) { return VolRegion(name); }

void CBDsimOpticalDiagnostics::BeginEvent(G4int eventId) {
  if (!EnvEnabled()) return;
  if (Budget().id >= 0) EndEvent(true);
  Budget() = EventBudget{};
  Budget().id = eventId;
  WindowArrivals().clear();
  Fates().clear();
  Steps().clear();
}

void CBDsimOpticalDiagnostics::PreUserTrackingAction(const G4Track* track) {
  if (!EnvEnabled() || track->GetDefinition() != G4OpticalPhoton::Definition()) return;
  // A suspended track can enter tracking again; count it only on first entry.
  auto& tracks = Budget().tracks;
  if (tracks.count(track->GetTrackID())) return;
  BudgetTrack t;
  t.parent = track->GetParentID();
  const auto* creator = track->GetCreatorProcess();
  t.source = t.parent == 0 ? "primary_optical" : creator ? creator->GetProcessName() : "unknown_creator";
  tracks.emplace(track->GetTrackID(), t);
}

void CBDsimOpticalDiagnostics::EndEvent(G4bool aborted) {
  if (!EnvEnabled() || Budget().id < 0) return;
  auto& c = Worker();
  std::map<G4String, G4long> counts;
  counts["started"] = 0;
  counts["fate_unfinished"] = 0;
  for (int i = 0; i < kNFateHists; ++i)
    if (i != 8 && i != 9) counts[G4String("fate_") + kFateNames[i]] = 0;
  for (const auto& e : Budget().tracks) {
    const auto& t = e.second;
    ++counts["started"];
    ++counts["source_" + t.source];
    const G4String terminal = t.terminal.empty() ? "unfinished" : t.terminal;
    ++counts["fate_" + terminal];
    if (terminal == "unknown" || terminal == "unfinished" || terminal == "no_rindex") {
      std::ostringstream row;
      row << G4Threading::G4GetThreadId() << " " << Budget().id << " " << e.first
          << " " << t.parent << " " << t.source << " " << terminal << " " << t.region
          << " " << t.process << " " << t.status << " " << t.tiny
          << " " << t.preVolume << " " << t.postVolume << " " << t.postMaterial
          << " " << t.boundaryStatus << " " << t.position.x()/CLHEP::mm
          << " " << t.position.y()/CLHEP::mm << " " << t.position.z()/CLHEP::mm;
      c.exceptionRows.push_back(row.str());
    }
  }
  G4long terminalSum = 0;
  for (const auto& e : counts) if (e.first.find("fate_") == 0 && e.first != "fate_unfinished") terminalSum += e.second;
  counts["terminal"] = terminalSum;
  counts["balance_delta"] = counts["started"] - terminalSum - counts["fate_unfinished"];
  counts["events"] = 1;
  counts["aborted_events"] = aborted ? 1 : 0;
  std::ostringstream row;
  row << G4Threading::G4GetThreadId() << " " << Budget().id;
  for (const auto& e : counts) { c.budgetCounts[e.first] += e.second; row << " " << e.first << "=" << e.second; }
  const double n = counts["started"];
  for (const auto& e : counts) if (e.first.find("fate_") == 0)
    row << " fraction_" << e.first.substr(5) << "=" << (n > 0 ? e.second/n : 0.);
  c.eventRows.push_back(row.str());
  Budget() = EventBudget{};
  Fates().clear(); Steps().clear(); WindowArrivals().clear();
}

void CBDsimOpticalDiagnostics::ResetForRun() {
  G4AutoLock lock(&diagMutex);
  if (!G4Threading::IsWorkerThread()) master = Counters{};
  Worker() = Counters{};
  Budget() = EventBudget{};
  Fates().clear();
  Steps().clear();
  WindowArrivals().clear();
}

void CBDsimOpticalDiagnostics::MergeWorkerIntoMaster() {
  if (!EnvEnabled()) return;
  G4AutoLock lock(&diagMutex);
  Add(master, Worker());
  Worker() = Counters{};
  Budget() = EventBudget{};
  Fates().clear();
  Steps().clear();
  WindowArrivals().clear();
}

void CBDsimOpticalDiagnostics::UserSteppingAction(const G4Step* step) {
  if (!EnvEnabled()) return;
  if (step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return;

  auto& c = Worker();
  const G4double mm = step->GetStepLength() / CLHEP::mm;
  G4String preVol = "none";
  if (auto* pv = step->GetPreStepPoint()->GetPhysicalVolume()) preVol = pv->GetName();
  RecordPath(c, VolRegion(preVol), mm);

  const G4StepPoint* post = step->GetPostStepPoint();
  if (post->GetStepStatus() == fGeomBoundary) {
    G4String postVol = "none";
    if (auto* pv = post->GetPhysicalVolume()) postVol = pv->GetName();
    RecordBoundary(c, preVol, postVol);
  }

  const G4VProcess* proc = post->GetProcessDefinedStep();
  if (EnvEnabled()) {
    UpdateTinyStep(Steps()[step->GetTrack()->GetTrackID()], mm, proc);
  }

  if (proc && proc->GetProcessName() == "OpAbsorption") {
    const G4String r = VolRegion(preVol);
    if (r == "scint") c.nBulkAbsorbScint++;
    else if (r == "lg") c.nBulkAbsorbLg++;
    else if (r == "sipm_glass") c.nBulkAbsorbGlass++;
    else if (r == "sipm_wafer") c.nBulkAbsorbSi++;
    else c.nBulkAbsorbOther++;

    if (EnvEnabled()) {
      const G4int tid = step->GetTrack()->GetTrackID();
      auto& fmap = Fates();
      TrackFate& fate = fmap[tid];
      if (fate == TrackFate::kUnknown) {
        if (r == "scint") fate = TrackFate::kAbsorbScint;
        else if (r == "lg") fate = TrackFate::kAbsorbLg;
        else if (r == "sipm_glass") fate = TrackFate::kAbsorbGlass;
        else if (r == "sipm_wafer") fate = TrackFate::kAbsorbSi;
        else fate = TrackFate::kAbsorbOther;
      }
    }
  }
  auto& fate = Fates()[step->GetTrack()->GetTrackID()];
  if (fate == TrackFate::kUnknown) {
    if (post->GetStepStatus() == fWorldBoundary && !post->GetPhysicalVolume()) fate = TrackFate::kWorldExit;
    else if (proc && (proc->GetProcessName() == "OpWLS" || proc->GetProcessName() == "OpWLS2")) fate = TrackFate::kWls;
  }
  // OpBoundary is forced, and does not limit the step (normally Transportation does).
  // Its status is meaningful here only for the current geometrical boundary step.
  if (post->GetStepStatus() == fGeomBoundary) {
    G4OpBoundaryProcess* opb = FindOpBoundaryProcess(step->GetTrack());
    if (opb) ++c.boundaryStatuses[static_cast<G4int>(opb->GetStatus())];
    if (opb && (opb->GetStatus() == Transmission || opb->GetStatus() == FresnelRefraction ||
                opb->GetStatus() == SameMaterial)) {
      const auto* pv = post->GetPhysicalVolume();
      const G4String a = VolRegion(preVol);
      const G4String b = VolRegion(pv ? pv->GetName() : "none");
      if (a != b) ++c.transmissions[a + "_to_" + b];
      // Reference arrival surface: coupling-side entry into the optical window.
      // Excludes inlet grease, world-side leaks, wafer returns and reflection encounters.
      if (b == "sipm_glass" && (a == "outlet_gel" || a == "lg" || a == "scint")) {
        ++c.windowEntrySteps;
        if (WindowArrivals().insert(step->GetTrack()->GetTrackID()).second)
          ++c.uniqueWindowArrivals;
      }
    }
    if (opb && fate == TrackFate::kUnknown) {
      if (opb->GetStatus() == Detection) fate = TrackFate::kDetected;
      else if (opb->GetStatus() == NoRINDEX) fate = TrackFate::kNoRindex;
    }
    if (opb && opb->GetStatus() == Absorption) {
      const G4String r = VolRegion(preVol);
      if (r == "lg") {
        c.nBoundaryAbsorbLg++;
        if (EnvEnabled()) {
          const G4int tid = step->GetTrack()->GetTrackID();
          auto& fmap = Fates();
          TrackFate& fate = fmap[tid];
          if (fate == TrackFate::kUnknown) fate = TrackFate::kAbsorbLgSurface;
        }
      } else {
        c.nBoundaryAbsorbOther++;
        if (EnvEnabled()) {
          auto& fate = Fates()[step->GetTrack()->GetTrackID()];
          if (fate == TrackFate::kUnknown) fate = TrackFate::kAbsorbOtherSurface;
        }
      }
    }
  }
}

void CBDsimOpticalDiagnostics::PostUserTrackingAction(const G4Track* track) {
  if (!EnvEnabled()) return;
  if (track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return;
  if (auto it = Budget().tracks.find(track->GetTrackID()); it != Budget().tracks.end()) {
    it->second.status = static_cast<G4int>(track->GetTrackStatus());
  }
  if (track->GetTrackStatus() != fStopAndKill && track->GetTrackStatus() != fKillTrackAndSecondaries) return;

  G4String vol = "none";
  if (auto* pv = track->GetVolume()) vol = pv->GetName();

  // Retain on suspension; release on actual termination, even with histogram mode off.
  WindowArrivals().erase(track->GetTrackID());

  const G4int tid = track->GetTrackID();
  auto& fmap = Fates();
  TrackFate fate = TrackFate::kUnknown;
  if (auto it = fmap.find(tid); it != fmap.end()) fate = it->second;

  const G4String region = VolRegion(vol);
  TrackStepState stepSt{};
  if (auto it = Steps().find(tid); it != Steps().end()) stepSt = it->second;

  if (fate == TrackFate::kUnknown) fate = TrackFate::kKillOther;
  auto it = Budget().tracks.find(tid);
  if (it != Budget().tracks.end()) {
    auto& t = it->second;
    if (!t.terminal.empty()) {
      ++Worker().budgetCounts["duplicate_terminal_callbacks"];
      return;
    }
    t.terminal = kFateNames[FateHistIndex(fate)];
    t.region = region; t.process = stepSt.lastProc.empty() ? "none" : stepSt.lastProc;
    t.status = static_cast<G4int>(track->GetTrackStatus());
    t.tiny = IsNavKillCandidate(stepSt);
    t.position = track->GetPosition();
    if (const auto* step = track->GetStep()) {
      if (auto* pv = step->GetPreStepPoint()->GetPhysicalVolume()) t.preVolume = pv->GetName();
      const auto* post = step->GetPostStepPoint();
      if (auto* pv = post->GetPhysicalVolume()) t.postVolume = pv->GetName();
      if (auto* material = post->GetMaterial()) t.postMaterial = material->GetName();
      if (post->GetStepStatus() == fGeomBoundary)
        if (auto* boundary = FindOpBoundaryProcess(track)) t.boundaryStatus = static_cast<int>(boundary->GetStatus());
    }
    if (t.tiny) ++Worker().budgetCounts["tiny_step_candidates"];
  } else {
    ++Worker().budgetCounts["terminal_without_start"];
    std::ostringstream row;
    row << G4Threading::G4GetThreadId() << " " << Budget().id << " " << tid
        << " " << track->GetParentID() << " missing_start unknown " << region
        << " " << (stepSt.lastProc.empty() ? "none" : stepSt.lastProc)
        << " " << static_cast<int>(track->GetTrackStatus()) << " " << IsNavKillCandidate(stepSt) << " none none none -1 0 0 0";
    Worker().exceptionRows.push_back(row.str());
  }

  RecordKillRegion(Worker(), vol);
  if (fate == TrackFate::kKillNav) {
    Worker().nKillNav++;
    BumpKillLoc(Worker().killNavLoc, region);
  } else if (fate == TrackFate::kKillWorld) {
    BumpKillLoc(Worker().killWorldLoc, region);
  } else if (fate == TrackFate::kKillOther) {
    Worker().nKillOtherResidual++;
    BumpKillLoc(Worker().killOtherLoc, region);
  }

  const G4double pathMm = track->GetTrackLength() / CLHEP::mm;
  if (DiagHistEnabled()) FillPathHist(Worker().pathHist, fate, pathMm);
  fmap.erase(tid);
  Steps().erase(tid);
}

void CBDsimOpticalDiagnostics::RecordSipmDetect(const G4Step* step) {
  if (!EnvEnabled()) return;
  Worker().nSipmDetect++;
  if (!step) return;
  Fates()[step->GetTrack()->GetTrackID()] = TrackFate::kDetected;
}

void CBDsimOpticalDiagnostics::RecordSipmQeReject(const G4Step* step) {
  if (!EnvEnabled()) return;
  Worker().nSipmQeReject++;
  if (!step) return;
  const G4int tid = step->GetTrack()->GetTrackID();
  auto& fmap = Fates();
  if (fmap[tid] != TrackFate::kDetected) fmap[tid] = TrackFate::kQeReject;
}

void CBDsimOpticalDiagnostics::PrintSummary(G4int nEvents) {
  if (!EnvEnabled()) return;
  PrintCounters(master, nEvents);
}

void CBDsimOpticalDiagnostics::WriteSummaryFile(const G4String& path, G4int nEvents) {
  if (!EnvEnabled() || path.empty()) return;
  WriteCountersFile(master, path, nEvents);
  std::ofstream events(path + ".events.txt"), exceptions(path + ".exceptions.txt");
  if (!events || !exceptions) {
    G4Exception("CBDsimOpticalDiagnostics", "BudgetOutput", FatalException, "Cannot write event budget sidecars");
    return;
  }
  events << "# thread event key=value; started means distinct optical tracks entering tracking, not all generated/stack-killed photons\n";
  exceptions << "# thread event track parent source fate last_region last_process track_status tiny_candidate pre_volume post_volume post_material boundary_status x_mm y_mm z_mm\n";
  for (const auto& row : master.eventRows) events << row << "\n";
  for (const auto& row : master.exceptionRows) exceptions << row << "\n";
}

void CBDsimOpticalDiagnostics::WritePathHistogramFile(const G4String& path) {
  if (!DiagHistEnabled() || path.empty()) return;
  WritePathHistToFile(master, path);
}
