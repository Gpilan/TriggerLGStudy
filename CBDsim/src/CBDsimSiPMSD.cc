#include "CBDsimSiPMSD.hh"
#include "CBDsimSiPMHit.hh"
#include "CBDsimOpticalDiagnostics.hh"

#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4StepPoint.hh"
#include "G4Exception.hh"
#include "Randomize.hh"

#include <atomic>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace {
constexpr int kSiPMSD_DbgMaxPrint = 25;
constexpr const char* kQeCsvDefaultPath =
    "analysis/reference/pmt_r2076/R2076_quantum_efficiency_percent.csv";
}

namespace {
struct QeTable {
  std::vector<G4double> wavelengthNm;
  std::vector<G4double> qeFraction;
};

G4double Clamp01(G4double x) {
  if (x < 0.) return 0.;
  if (x > 1.) return 1.;
  return x;
}

QeTable LoadQeCsvOrThrow(const std::string& path) {
  QeTable table;
  std::ifstream fin(path);
  if (!fin.is_open()) {
    G4ExceptionDescription msg;
    msg << "Failed to open QE CSV: " << path;
    G4Exception("CBDsimSiPMSD::LoadQeCsvOrThrow", "QE_CSV_OPEN", FatalException, msg);
  }

  std::string line;
  // Skip header: wavelength_nm,quantum_efficiency_percent
  std::getline(fin, line);
  while (std::getline(fin, line)) {
    if (line.empty()) continue;

    std::stringstream ss(line);
    std::string wStr;
    std::string qStr;
    if (!std::getline(ss, wStr, ',')) {
      G4ExceptionDescription msg;
      msg << "Malformed QE CSV row (wavelength missing): " << line;
      G4Exception("CBDsimSiPMSD::LoadQeCsvOrThrow", "QE_CSV_PARSE", FatalException, msg);
    }
    if (!std::getline(ss, qStr, ',')) {
      G4ExceptionDescription msg;
      msg << "Malformed QE CSV row (QE missing): " << line;
      G4Exception("CBDsimSiPMSD::LoadQeCsvOrThrow", "QE_CSV_PARSE", FatalException, msg);
    }

    const G4double w = std::stod(wStr);
    const G4double qPercent = std::stod(qStr);
    table.wavelengthNm.push_back(w);
    table.qeFraction.push_back(Clamp01(qPercent * 0.01));
  }

  if (table.wavelengthNm.size() < 2) {
    G4ExceptionDescription msg;
    msg << "QE CSV has too few points: " << table.wavelengthNm.size();
    G4Exception("CBDsimSiPMSD::LoadQeCsvOrThrow", "QE_CSV_SIZE", FatalException, msg);
  }

  for (size_t i = 1; i < table.wavelengthNm.size(); ++i) {
    if (!(table.wavelengthNm[i - 1] < table.wavelengthNm[i])) {
      G4ExceptionDescription msg;
      msg << "QE CSV wavelength must be strictly increasing.";
      G4Exception("CBDsimSiPMSD::LoadQeCsvOrThrow", "QE_CSV_ORDER", FatalException, msg);
    }
  }

  G4cout << "[SiPM SD] Loaded QE CSV: " << path
         << " (points=" << table.wavelengthNm.size() << ")" << G4endl;
  return table;
}

const QeTable& GetQeTable() {
  static const QeTable table = LoadQeCsvOrThrow(kQeCsvDefaultPath);
  return table;
}

G4double InterpolateQeFraction(const QeTable& table, G4double wavelengthNm) {
  const auto& xs = table.wavelengthNm;
  const auto& ys = table.qeFraction;
  if (wavelengthNm <= xs.front()) return ys.front();
  if (wavelengthNm >= xs.back()) return ys.back();

  auto it = std::lower_bound(xs.begin(), xs.end(), wavelengthNm);
  const size_t i1 = static_cast<size_t>(it - xs.begin());
  const size_t i0 = i1 - 1;
  const G4double x0 = xs[i0];
  const G4double x1 = xs[i1];
  const G4double y0 = ys[i0];
  const G4double y1 = ys[i1];
  if (x1 <= x0) return y0;
  const G4double t = (wavelengthNm - x0) / (x1 - x0);
  return Clamp01(y0 + t * (y1 - y0));
}
}  // namespace

namespace {

G4bool PvIsSipmWafer(const G4VPhysicalVolume* pv) {
  if (!pv) return false;
  const G4String& nm = pv->GetName();
  return nm == "protoSipmWaferPhys" || nm == "waferPhysical";
}

G4bool TowerFromProtoSipmTouchable(const G4VTouchable* touch, G4int& towernum) {
  if (!touch) return false;
  const G4int depth = touch->GetHistoryDepth();
  for (G4int d = 0; d <= depth; ++d) {
    auto* vol = touch->GetVolume(d);
    if (vol && vol->GetName() == "protoSipmEnvPhys") {
      towernum = vol->GetCopyNo();
      return true;
    }
  }
  return false;
}

/** Legacy: depth-2 = SiPM cell copy, depth-3 = tower parent. Proto: walk touchable for protoSipmEnvPhys.
 *  Optical photons: use the step point whose PV is the wafer (pre or post); Post alone can mis-identify at boundaries. */
void SiPMAndTowerFromTouchable(G4Step* step, G4int& SiPMnum, G4int& towernum) {
  SiPMnum = 0;
  towernum = 0;
  const G4StepPoint* inWafer = nullptr;
  if (PvIsSipmWafer(step->GetPostStepPoint()->GetPhysicalVolume())) {
    inWafer = step->GetPostStepPoint();
  } else if (PvIsSipmWafer(step->GetPreStepPoint()->GetPhysicalVolume())) {
    inWafer = step->GetPreStepPoint();
  }
  if (inWafer) {
    if (TowerFromProtoSipmTouchable(inWafer->GetTouchable(), towernum)) return;
  }
  auto* touch = step->GetPostStepPoint()->GetTouchable();
  G4VPhysicalVolume* vol1 = touch->GetVolume(1);
  G4VPhysicalVolume* vol2 = touch->GetVolume(2);
  G4VPhysicalVolume* vol3 = touch->GetVolume(3);
  if (vol1 && vol1->GetName() == "protoSipmEnvPhys") {
    towernum = vol1->GetCopyNo();
    SiPMnum = 0;
    return;
  }
  SiPMnum = vol2 ? vol2->GetCopyNo() : 0;
  towernum = vol3 ? vol3->GetCopyNo() : 0;
}
}  // namespace

CBDsimSiPMSD::CBDsimSiPMSD(const G4String& name, const G4String& hitsCollectionName, std::pair<int,int> xy)
: G4VSensitiveDetector(name), fHitCollection(0), fHCID(-1), fTowerXY(xy),fWavBin(600),fTimeBin(24000),fWavlenStart(900.),fWavlenEnd(300.),fTimeStart(0.),fTimeEnd(240.)
{
  collectionName.insert(hitsCollectionName);
  fWavlenStep = (fWavlenStart-fWavlenEnd)/(float)fWavBin;
  fTimeStep = (fTimeEnd-fTimeStart)/(float)fTimeBin;
}

CBDsimSiPMSD::~CBDsimSiPMSD() {}

void CBDsimSiPMSD::Initialize(G4HCofThisEvent* hce) {
  fHitCollection = new CBDsimSiPMHitsCollection(SensitiveDetectorName,collectionName[0]);
  if (fHCID<0) { fHCID = GetCollectionID(0); }
  hce->AddHitsCollection(fHCID,fHitCollection);
}

G4bool CBDsimSiPMSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
  if (step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return false;
  const G4double energy = step->GetTrack()->GetTotalEnergy();
  const G4double wavelengthNm = (h_Planck * c_light / energy) / nm;
  const QeTable& qeTable = GetQeTable();
  const G4double qe = InterpolateQeFraction(qeTable, wavelengthNm);
  const bool detected = (G4UniformRand() < qe);
  if (!detected) {
  CBDsimOpticalDiagnostics::RecordSipmQeReject(step);
    // Detection model: photon is absorbed in Si wafer; only accepted fraction is counted.
    step->GetTrack()->SetTrackStatus(fStopAndKill);
    return false;
  }

  G4int SiPMnum = 0;
  G4int towernum = 0;
  SiPMAndTowerFromTouchable(step, SiPMnum, towernum);

  // 진단용: 광학 광자가 웨이퍼 SD로 들어올 때만 호출됨. 확인 후 제거하거나 kSiPMSD_DbgMaxPrint 를 0 으로.
  {
    static std::atomic<int> sDbgProcessHitsCount{0};
    const int n = ++sDbgProcessHitsCount;
    if (n <= kSiPMSD_DbgMaxPrint) {
      G4String postName = "(no post PV)";
      if (auto* pv = step->GetPostStepPoint()->GetPhysicalVolume()) postName = pv->GetName();
      G4int evId = -1;
      if (auto* em = G4EventManager::GetEventManager()) {
        if (auto* ev = em->GetConstCurrentEvent()) evId = ev->GetEventID();
      }
      G4cout << "[SiPM SD] ProcessHits #" << n << " ev=" << evId << " tower=" << towernum
             << " SiPMnum=" << SiPMnum << " postPV=" << postName << G4endl;
    }
  }

  G4double hitTime  = step->GetPostStepPoint()->GetGlobalTime();
  G4int towerX = fTowerXY.first;
  G4int towerY = fTowerXY.second;
  G4int sipmX = SiPMnum/towerY;
  G4int sipmY = SiPMnum%towerY;
  G4int nofHits = fHitCollection->entries();
  CBDsimSiPMHit* hit = NULL;


  for (G4int i = 0; i < nofHits; i++) {
    if ( ((*fHitCollection)[i]->GetSiPMnum()==SiPMnum) && ((*fHitCollection)[i]->GetTowernum()==towernum) ) {
      hit = (*fHitCollection)[i];
      break;
    }
  }
  if (hit==NULL) {

    hit = new CBDsimSiPMHit(fWavBin,fTimeBin);
    hit->SetSiPMnum(SiPMnum);
    hit->SetTowernum(towernum);
    hit->SetX(sipmX);
    hit->SetY(sipmY);
    hit->SetTowerX(towerX);
    hit->SetTowerY(towerY);

    fHitCollection->insert(hit);
  }

  hit->photonCount();
  CBDsimOpticalDiagnostics::RecordSipmDetect(step);

  CBDsimInterface::hitRange wavRange = findWavRange(energy);
  hit->CountWavlenSpectrum(wavRange);


  CBDsimInterface::hitRange timeRange = findTimeRange(hitTime);
  hit->CountTimeStruct(timeRange);

  // One optical photon should contribute once when it enters/steps in wafer.
  // Kill track here to prevent multiple SD hits from repeated internal stepping/reflections.
  step->GetTrack()->SetTrackStatus(fStopAndKill);


  return true;
}

void CBDsimSiPMSD::EndOfEvent(G4HCofThisEvent*) {
  if ( verboseLevel>1 ) {
    G4int nofHits = fHitCollection->entries();
    G4cout
    << G4endl
    << "-------->Hits Collection: in this event they are " << nofHits
    << " hits in the tracker chambers: " << G4endl;
    for ( G4int i=0; i<nofHits; i++ ) (*fHitCollection)[i]->Print();
  }
}

CBDsimInterface::hitRange CBDsimSiPMSD::findWavRange(G4double en) {
  int i = 0;
  for ( ; i < fWavBin+1; i++) {
    if ( en < wavToE( (fWavlenStart - (float)i*fWavlenStep)*nm ) ) break;
  }

  if (i==0) return std::make_pair(fWavlenStart,99999.);
  else if (i==fWavBin+1) return std::make_pair(0.,fWavlenEnd);

  return std::make_pair( fWavlenStart-(float)i*fWavlenStep, fWavlenStart-(float)(i-1)*fWavlenStep );
}

CBDsimInterface::hitRange CBDsimSiPMSD::findTimeRange(G4double stepTime) {
  int i = 0;
  for ( ; i < fTimeBin+1; i++) {
    if ( stepTime < ( (fTimeStart + (float)i*fTimeStep)*ns ) ) break;
  }

  if (i==0) return std::make_pair(0.,fTimeStart);
  else if (i==fTimeBin+1) return std::make_pair(fTimeEnd,99999.);

  return std::make_pair( fTimeStart+(float)(i-1)*fTimeStep, fTimeStart+(float)i*fTimeStep );
}
