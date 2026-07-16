#include "CBDsimRunAction.hh"
#include "CBDsimOpticalDiagnostics.hh"

#include "G4AutoLock.hh"
#include "G4Run.hh"
#include "G4Threading.hh"

#include <cstdlib>
#include <iostream>

namespace { G4Mutex CBDsimRunActionMutex = G4MUTEX_INITIALIZER; }
CBDsimRootInterface* CBDsimRunAction::sRootIO = 0;
int CBDsimRunAction::sNumEvt = 0;

CBDsimRunAction::CBDsimRunAction(G4int seed, G4String filename)
: G4UserRunAction() {
  fSeed = seed;
  fFilename = filename;

  G4AutoLock lock(&CBDsimRunActionMutex);

  // ROOT 출력: 파일명이 있을 때만 생성 (빈 문자열이면 GUI/지오메트리 확인용으로 .root 미생성)
  if (!fFilename.empty() && !sRootIO) {
    sRootIO = new CBDsimRootInterface(fFilename+"_"+std::to_string(fSeed)+".root");
    sRootIO->create();
  }
}

CBDsimRunAction::~CBDsimRunAction() {
  if (IsMaster()) {
    G4AutoLock lock(&CBDsimRunActionMutex);

    if (sRootIO) {
      sRootIO->write();
      sRootIO->close();
      delete sRootIO;
      sRootIO = 0;
    }
  }
}

void CBDsimRunAction::BeginOfRunAction(const G4Run*) {
  CBDsimOpticalDiagnostics::ResetForRun();
}

void CBDsimRunAction::EndOfRunAction(const G4Run* run) {
  CBDsimOpticalDiagnostics::MergeWorkerIntoMaster();
  if (IsMaster() && run) {
    const G4int nEvents = run->GetNumberOfEvent();
    CBDsimOpticalDiagnostics::PrintSummary(nEvents);
    if (const char* out = std::getenv("CBDsim_OPTICAL_DIAG_OUT")) {
      CBDsimOpticalDiagnostics::WriteSummaryFile(out, nEvents);
    }
    if (const char* histOut = std::getenv("CBDsim_OPTICAL_DIAG_HIST_OUT")) {
      CBDsimOpticalDiagnostics::WritePathHistogramFile(histOut);
    }
  }
}
