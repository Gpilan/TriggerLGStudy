#include "CBDsimOpticalDiagnostics.hh"

#include "G4AutoLock.hh"
#include "G4OpticalPhoton.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

#include <array>
#include <cstdlib>
#include <fstream>
#include <unordered_map>

namespace {
enum class TrackFate : int {
  kUnknown = 0,
  kDetected,
  kQeReject,
  kAbsorbScint,
  kAbsorbLg,
  kAbsorbGlass,
  kAbsorbSi,
  kAbsorbOther,
  kKillWorld,
  kKillOther,
};

constexpr int kNFateHists = 9;
constexpr int kNBins = 400;
constexpr G4double kXmaxMm = 2000.0;

constexpr std::array<const char*, kNFateHists> kFateNames = {
    "detected",    "qe_reject",    "absorb_scint", "absorb_lg",  "absorb_glass",
    "absorb_si",   "absorb_other", "kill_world",   "kill_other",
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

struct Counters {
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
  G4long nTracksKilled = 0;
  G4long nKillLastScint = 0;
  G4long nKillLastLg = 0;
  G4long nKillLastWorld = 0;
  G4long nKillLastSipm = 0;
  G4long nKillLastOther = 0;
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
  a.nTracksKilled += b.nTracksKilled;
  a.nKillLastScint += b.nKillLastScint;
  a.nKillLastLg += b.nKillLastLg;
  a.nKillLastWorld += b.nKillLastWorld;
  a.nKillLastSipm += b.nKillLastSipm;
  a.nKillLastOther += b.nKillLastOther;
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
  if (volName.find("protoScint") != G4String::npos) return "scint";
  if (volName.find("protoLightGuide") != G4String::npos) return "lg";
  if (volName.find("protoSipmWafer") != G4String::npos) return "sipm_wafer";
  if (volName.find("protoSipmWindow") != G4String::npos) return "sipm_glass";
  if (volName.find("protoSipmEnv") != G4String::npos) return "sipm_env";
  if (volName.find("protoSipmGel") != G4String::npos) return "gel";
  if (volName.find("protoFoil") != G4String::npos) return "foil";
  if (volName.find("protoWorld") != G4String::npos || volName == "worldPhysical" ||
      volName == "worldLogical")
    return "world";
  return "other";
}

void RecordPath(Counters& c, const G4String& region, G4double mm) {
  if (region == "scint") c.pathScintMm += mm;
  else if (region == "lg") c.pathLgMm += mm;
  else if (region == "sipm_glass") c.pathGlassMm += mm;
  else if (region == "sipm_wafer") c.pathSiMm += mm;
  else if (region == "world" || region == "sipm_env") c.pathWorldMm += mm;
}

void RecordBoundary(Counters& c, const G4String& pre, const G4String& post) {
  c.nBoundarySteps++;
  const G4String a = VolRegion(pre);
  const G4String b = VolRegion(post);
  if (a == "scint" && b == "lg") c.nScintToLg++;
  else if (a == "lg" && b == "scint") c.nLgToScint++;
  else if (a == "scint" && b == "world") c.nScintToWorld++;
  else if (a == "world" && b == "scint") c.nWorldToScint++;
  else if (a == "lg" && b == "world") c.nLgToWorld++;
  else if (a == "world" && b == "lg") c.nWorldToLg++;
  else if (a == "lg" && (b == "sipm_env" || b == "sipm_glass" || b == "gel")) c.nLgToSipm++;
  else if ((a == "sipm_env" || a == "sipm_glass" || a == "gel") && b == "lg") c.nSipmToLg++;
  else if (a == "gel" && b == "sipm_glass") c.nLgToSipm++;
  else if (a == "sipm_glass" && b == "gel") c.nSipmToLg++;
  else c.nOtherBoundary++;
}

void RecordKillRegion(Counters& c, const G4String& volName) {
  c.nTracksKilled++;
  const G4String r = VolRegion(volName);
  if (r == "scint") c.nKillLastScint++;
  else if (r == "lg") c.nKillLastLg++;
  else if (r == "world" || r == "sipm_env") c.nKillLastWorld++;
  else if (r == "sipm_wafer" || r == "sipm_glass" || r == "gel") c.nKillLastSipm++;
  else c.nKillLastOther++;
}

int FateHistIndex(TrackFate fate) {
  switch (fate) {
    case TrackFate::kDetected: return 0;
    case TrackFate::kQeReject: return 1;
    case TrackFate::kAbsorbScint: return 2;
    case TrackFate::kAbsorbLg: return 3;
    case TrackFate::kAbsorbGlass: return 4;
    case TrackFate::kAbsorbSi: return 5;
    case TrackFate::kAbsorbOther: return 6;
    case TrackFate::kKillWorld: return 7;
    case TrackFate::kKillOther: return 8;
    default: return 8;
  }
}

TrackFate FateFromRegion(const G4String& region) {
  if (region == "world" || region == "sipm_env" || region == "foil") return TrackFate::kKillWorld;
  if (region == "scint") return TrackFate::kKillOther;
  if (region == "lg") return TrackFate::kKillOther;
  if (region == "gel" || region == "sipm_glass" || region == "sipm_wafer") return TrackFate::kKillOther;
  return TrackFate::kKillOther;
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
  G4cout << "  boundary steps (total): " << c.nBoundarySteps << G4endl;
  G4cout << "    scint->LG: " << c.nScintToLg << "  LG->scint: " << c.nLgToScint << G4endl;
  G4cout << "    scint->world: " << c.nScintToWorld << "  world->scint: " << c.nWorldToScint << G4endl;
  G4cout << "    LG->world: " << c.nLgToWorld << "  world->LG: " << c.nWorldToLg << G4endl;
  G4cout << "    LG->SiPM: " << c.nLgToSipm << "  SiPM->LG: " << c.nSipmToLg << G4endl;
  G4cout << "    other boundary: " << c.nOtherBoundary << G4endl;
  G4cout << "  bulk OpAbsorption: scint=" << c.nBulkAbsorbScint << " lg=" << c.nBulkAbsorbLg
         << " glass=" << c.nBulkAbsorbGlass << " si=" << c.nBulkAbsorbSi
         << " other=" << c.nBulkAbsorbOther << G4endl;
  G4cout << "  tracks killed: " << c.nTracksKilled
         << " (last vol: scint=" << c.nKillLastScint << " lg=" << c.nKillLastLg
         << " world/env=" << c.nKillLastWorld << " sipm=" << c.nKillLastSipm
         << " other=" << c.nKillLastOther << ")" << G4endl;
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
  out << "boundary_total " << c.nBoundarySteps << "\n";
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
  out << "sipm_detect " << c.nSipmDetect << "\n";
  out << "sipm_qe_reject " << c.nSipmQeReject << "\n";
  out << "path_scint_mm_per_event " << c.pathScintMm * inv << "\n";
  out << "path_lg_mm_per_event " << c.pathLgMm * inv << "\n";
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

void CBDsimOpticalDiagnostics::ResetForRun() {
  G4AutoLock lock(&diagMutex);
  master = Counters{};
  Worker() = Counters{};
  Fates().clear();
}

void CBDsimOpticalDiagnostics::MergeWorkerIntoMaster() {
  if (!EnvEnabled()) return;
  G4AutoLock lock(&diagMutex);
  Add(master, Worker());
  Worker() = Counters{};
  Fates().clear();
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
  if (proc && proc->GetProcessName() == "OpAbsorption") {
    const G4String r = VolRegion(preVol);
    if (r == "scint") c.nBulkAbsorbScint++;
    else if (r == "lg") c.nBulkAbsorbLg++;
    else if (r == "sipm_glass") c.nBulkAbsorbGlass++;
    else if (r == "sipm_wafer") c.nBulkAbsorbSi++;
    else c.nBulkAbsorbOther++;

    if (DiagHistEnabled()) {
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
}

void CBDsimOpticalDiagnostics::PostUserTrackingAction(const G4Track* track) {
  if (!EnvEnabled()) return;
  if (track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return;
  if (track->GetTrackStatus() != fStopAndKill) return;

  G4String vol = "none";
  if (auto* pv = track->GetVolume()) vol = pv->GetName();
  RecordKillRegion(Worker(), vol);

  if (!DiagHistEnabled()) return;

  const G4int tid = track->GetTrackID();
  auto& fmap = Fates();
  TrackFate fate = TrackFate::kUnknown;
  if (auto it = fmap.find(tid); it != fmap.end()) fate = it->second;
  if (fate == TrackFate::kUnknown) fate = FateFromRegion(VolRegion(vol));

  const G4double pathMm = track->GetTrackLength() / CLHEP::mm;
  FillPathHist(Worker().pathHist, fate, pathMm);
  fmap.erase(tid);
}

void CBDsimOpticalDiagnostics::RecordSipmDetect(const G4Step* step) {
  if (!EnvEnabled()) return;
  Worker().nSipmDetect++;
  if (!DiagHistEnabled() || !step) return;
  Fates()[step->GetTrack()->GetTrackID()] = TrackFate::kDetected;
}

void CBDsimOpticalDiagnostics::RecordSipmQeReject(const G4Step* step) {
  if (!EnvEnabled()) return;
  Worker().nSipmQeReject++;
  if (!DiagHistEnabled() || !step) return;
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
}

void CBDsimOpticalDiagnostics::WritePathHistogramFile(const G4String& path) {
  if (!DiagHistEnabled() || path.empty()) return;
  WritePathHistToFile(master, path);
}
