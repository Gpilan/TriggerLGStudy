#ifndef CBDsimOpticalDiagnostics_h
#define CBDsimOpticalDiagnostics_h 1

#include "globals.hh"

class G4Step;
class G4Track;

/**
 * Run-level optical photon boundary / fate counters (enable: CBDsim_OPTICAL_DIAG=1).
 * Per-photon total track-length histograms by fate (enable: CBDsim_OPTICAL_DIAG_HIST=1,
 * output: CBDsim_OPTICAL_DIAG_HIST_OUT).
 */
class CBDsimOpticalDiagnostics {
public:
  static bool Enabled();
  static bool HistEnabled();
  static void ResetForRun();
  static void MergeWorkerIntoMaster();
  static void UserSteppingAction(const G4Step* step);
  static void PostUserTrackingAction(const G4Track* track);
  static void RecordSipmDetect(const G4Step* step);
  static void RecordSipmQeReject(const G4Step* step);
  static void PrintSummary(G4int nEvents);
  static void WriteSummaryFile(const G4String& path, G4int nEvents);
  static void WritePathHistogramFile(const G4String& path);
};

#endif
