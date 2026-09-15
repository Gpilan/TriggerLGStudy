# D01 optical boundary validation

This standalone target compiles the production diagnostics with Geant4 11.2 transport,
without modifying the main CBDsim executable. It uses a single planar interface, one
3 eV photon per event at normal incidence, index 1.5 on both sides, and fixed polarization.
No production geometry, PMT model, ROOT output, or Condor submission is involved.

From the repository root:

```bash
source envset.sh
cmake -S CBDsim/tests/optical_boundary -B /tmp/trigger-d01-build
cmake --build /tmp/trigger-d01-build -j2
mkdir -p /tmp/trigger-d01-results
/tmp/trigger-d01-build/optical_boundary_test mixed 5000 42 /tmp/trigger-d01-results/mixed.txt
```

Modes: `reflect`, `absorb`, `other_absorb`, `transmit`, `mixed`, `bulk`, `detect`.
Use a fresh output path for each invocation. Each run exits nonzero on failed assertions.
The test compares production summary status counts to observed boundary statuses, and
checks analytical expectations: exact 5000 counts for deterministic modes, zero bulk/surface
cross-contamination, and binomial 5-sigma tolerance for R=0.5. Run `mixed` with seeds 42 and 43.
`detect` tests Geant4 surface Detection status, not the application's SiPM/QE model.

The stepping action deliberately terminates surviving photons after their first encounter.
Consequently the test validates boundary diagnostics only; it cannot certify D03 final fates,
world escape, the production geometry, or T1/T2 timing. `boundary_status_<integer>` uses the
G4OpBoundaryProcessStatus enum of the recorded Geant4 version. Legacy transition counters
remain boundary encounters, not unique arrivals (D02).

For before/after regression, configure a second build with
`-DDIAGNOSTICS_SOURCE=/absolute/path/to/before-diagnostics.cc`.

## D02 layered arrival validation

The same build also produces `optical_arrival_test`. Example:

```bash
/tmp/trigger-d01-build/optical_arrival_test reenter 100 42 /tmp/trigger-d01-results/reenter.txt 1
```

Arguments: mode, events, seed, fresh summary filename, histogram enabled (`0` or `1`).
Modes: `forward`, `inlet`, `reflect`, `reenter`, `nolg`.
There are two photons per event, with track IDs reused across events. Run each mode
with histogram off and on. The test exercises production BeginEvent/stepping/tracking/merge.
The oracle uses the known z=2 mm plane, travel direction and optical status rather than
the production region mapper. Twenty-four representative volume names are also checked.

- `forward`: scint → inlet gel → LG → outlet gel → window. One window entry per photon.
- `inlet`: stop after inlet gel → LG. No window entries or legacy LG→SiPM encounters.
- `reflect`: perfect mirror at the window front. Encounter count increases but arrivals remain zero.
- `reenter`: mirrors return each photon through the window front a second time. Two entries,
  one unique arrival per photon. This tests deduplication before track termination.
- `nolg`: scint → outlet gel → window. No inlet gel/LG contribution.

All expected counts are exact; inlet path sum is 4 mm/event for two photons crossing a
2 mm gel layer once. Tests use planar, nonoverlapping layers of distinct equal-index
materials. They do not validate production geometry or a full photon fate budget.

## Summary schema 2

- `boundary_encounter_<region>_to_<region>`: boundary pre/post-volume encounters,
  including reflection attempts; these names alone do not imply physical transmission.
- `boundary_transmit_<region>_to_<region>`: distinct-region boundaries with status
  Transmission, FresnelRefraction or SameMaterial. Reflection/absorption excluded.
- `window_entry_steps`: transmitted coupling-side entries into prototype glass from
  outlet_gel, LG or scint. World-side leaks and wafer returns are excluded by definition.
- `unique_window_arrivals`: each event/track at most once across all prototype windows.
  Not a per-trigger sum, not a QE acceptance count and not a photocathode arrival count.
- `lg_to_sipm`/`sipm_to_lg`: retained, deprecated legacy coupling encounter aggregates.
  They can include multiple interfaces per photon. Use the new entry/unique fields for arrivals.
  Fixing inlet classification changes historical region aggregates; compare schema versions explicitly.
- `path_inlet_gel_mm_per_event`, `path_outlet_gel_mm_per_event`: distinct path sums.
  Kill-location `_gel` remains the sum of `_inlet_gel` and `_outlet_gel` for compatibility.

New integrations must call `CBDsimOpticalDiagnostics::BeginEvent()` on each worker at
event start. The production EventAction does so. Dedup state is also cleared on track
termination, run reset and merge, and is retained when a track is suspended.

## D03: event photon budgets (schema 3)

`optical_budget_test <mode> <events> <seed> <summary.txt> <hist:0|1>`
uses real Geant4 transport for `absorb`, `other_absorb`, `transmit`, `bulk`,
`detect` (boundary EFFICIENCY), `wls`, and `no_rindex`. `qe` and `sd_detect`
exercise the production diagnostic SD callbacks with forced decisions; they do
not validate the PMT QE probability/model. `unknown` explicitly kills a live
track. No tiny-step heuristic is allowed to label a navigation cause.

For example, after building as above:

```sh
./build/optical_budget_test wls 100 42 /tmp/wls.txt 0
./build/optical_budget_lifecycle_test /tmp/lifecycle.txt 0
```

Each WLS primary converts once; its independently tracked child leaves the world:
100 primaries + 100 WLS children = 100 conversions + 100 world exits.
The lifecycle test checks suspension/resumption, reused event-local IDs,
kill-with-secondaries status, duplicate terminal callbacks, and an aborted event
with an unfinished track. This is a callback contract test, not transport physics.

Production hooks: `BeginEvent(eventId)`, `PreUserTrackingAction(track)`, stepping/SD,
`PostUserTrackingAction(track)`, `EndEvent(aborted)`. The denominator is **distinct
optical tracks that enter tracking**; stack-killed/untracked photons are excluded.
The current application registers no user stacking action. Creator process names
(including Scintillation, Cerenkov and OpWLS/OpWLS2) remain separate source columns.
Suspension never finalizes a photon. At event end, missing final callbacks remain
`unfinished`; a stopped photon without a known cause is `unknown`.

With `CBDsim_OPTICAL_DIAG=1` and `CBDsim_OPTICAL_DIAG_OUT=<path>`, the summary now
also writes `<path>.events.txt` (one row/event, counts and fate fractions) and
`<path>.exceptions.txt` (thread/event/track/parent/source/cause/last region/process,
track status and tiny-step candidate). Histograms are optional and do not affect
these budgets. `budget_closed=1` means all event records and final callbacks close;
it does **not** mean unknown=0 or that geometry/PMT physics is validated. Check
`budget_fate_unknown`, `budget_fate_no_rindex`, and the exception file separately.
Run-level rows are merged from workers; order is not event order. Output storage
scales with the number of events and exceptions, plus tracks in the current event.

Schema 3 histogram fates use `world_exit`, `unknown`, `wls_conversion`, `no_rindex`.
Old kill_world/kill_nav slots are retained as zero-valued legacy placeholders;
their historical heuristic meaning must not be compared as the same cause.
D01/D02 harnesses retain their narrower boundary/arrival assertions and do not
supply all budget lifecycle hooks; their budget sidecars are not closure tests.
Actual prototype budget validation must be repeated after G01 geometry fixes.
