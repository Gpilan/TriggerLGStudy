// Callback lifecycle contract test; real transport is covered by budget_test.cc.
#include "CBDsimOpticalDiagnostics.hh"
#include "G4DynamicParticle.hh"
#include "G4OpticalPhoton.hh"
#include "G4Track.hh"
#include <cstdlib>
#include <fstream>
#include <map>
#include <iostream>
int main(int argc, char** argv) {
  if (argc != 3) return 2;
  setenv("CBDsim_OPTICAL_DIAG", "1", 1);
  setenv("CBDsim_OPTICAL_DIAG_HIST", argv[2], 1);
  CBDsimOpticalDiagnostics::ResetForRun();
  for (int event = 0; event < 3; ++event) {
    CBDsimOpticalDiagnostics::BeginEvent(event);
    G4Track track(new G4DynamicParticle(G4OpticalPhoton::Definition(), {0,0,1}, 3.e-6), 0., {});
    track.SetTrackID(1); track.SetParentID(0);
    CBDsimOpticalDiagnostics::PreUserTrackingAction(&track);
    track.SetTrackStatus(fSuspend);
    CBDsimOpticalDiagnostics::PostUserTrackingAction(&track);
    CBDsimOpticalDiagnostics::PreUserTrackingAction(&track);
    if (event < 2) {
      track.SetTrackStatus(event == 0 ? fStopAndKill : fKillTrackAndSecondaries);
      CBDsimOpticalDiagnostics::PostUserTrackingAction(&track);
      if (event == 0) CBDsimOpticalDiagnostics::PostUserTrackingAction(&track);
    }
    CBDsimOpticalDiagnostics::EndEvent(event == 2);
  }
  CBDsimOpticalDiagnostics::MergeWorkerIntoMaster();
  CBDsimOpticalDiagnostics::WriteSummaryFile(argv[1], 3);
  std::ifstream in(argv[1]); std::string k; double v; std::map<std::string,double> c;
  while (in >> k >> v) c[k] = v;
  const bool ok = c["budget_started"] == 3 && c["budget_terminal"] == 2 &&
    c["budget_fate_unknown"] == 2 && c["budget_fate_unfinished"] == 1 &&
    c["budget_aborted_events"] == 1 && c["budget_duplicate_terminal_callbacks"] == 1 &&
    c["budget_balance_delta"] == 0;
  std::cout << "LIFECYCLE hist=" << argv[2] << " pass=" << ok << std::endl;
  return ok ? 0 : 1;
}
