#ifndef CBDsimOpticalDiagnostics_h
#define CBDsimOpticalDiagnostics_h 1

#include "globals.hh"

class G4Step;
class G4Track;

/**
 * Run-level optical photon boundary / fate counters (enable: CBDsim_OPTICAL_DIAG=1).
 * Per-photon total track-length histograms by fate (enable: CBDsim_OPTICAL_DIAG_HIST=1,
 * output: CBDsim_OPTICAL_DIAG_HIST_OUT).
 * absorb_lg_surface: OpBoundary Absorption with pre-volume = LG (split out of kill_other).
 * absorb_other_surface: OpBoundary Absorption in other regions.
 * Summary boundary_status_<enum value> records current fGeomBoundary status encounters;
 * these are step counts, not unique photon arrivals or a complete photon budget.
 * Schema 2 separates inlet_gel/outlet_gel and adds per-region encounter/transmission counts.
 * unique_window_arrivals counts each event/track once across all prototype windows, on
 * transmission from outlet_gel, LG or scint into glass (not photocathode detection).
 * lg_to_sipm/sipm_to_lg are deprecated legacy encounter sums, NOT unique arrivals.
 * Schema 3: event budgets are independent of histograms. Only fWorldBoundary means
 * world_exit; tiny steps are auxiliary observations, never a navigation-loss cause.
 * Budget denominator: distinct optical tracks entering tracking (includes WLS children).
 * Untracked/stack-killed photons are outside this denominator; incomplete events stay explicit.
 */
class CBDsimOpticalDiagnostics {
public:
  static bool Enabled();
  static bool HistEnabled();
  static void ResetForRun();
  // Called on each worker at event start; track IDs are only unique within an event.
  static void BeginEvent(G4int eventId = -1);
  static void EndEvent(G4bool aborted = false);
  static void PreUserTrackingAction(const G4Track* track);
  static G4String RegionForVolume(const G4String& name);
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
