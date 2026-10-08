# Standalone Trigger Geant4 11.4 migration — updated 2026-10-08

This baseline uses Geant4 11.4.0, LCG_107 / GCC 11.3.0,
ROOT 6.34.02 and C++17. Start a fresh shell and use a separate build directory;
do not reuse libraries or CMake caches from the 11.2 / GCC13 build.

The project goal is validated standalone Trigger → event phase-space handoff →
standalone DRC → event-by-event timing correlation. Full implant, common geometry,
SD/ROOT merging and forced QE conversion are out of scope, including fallbacks.
The October 8 development checkpoint below records the user's decision to
continue with output-plane design while retaining the measured version differences.

## 2026-10-08 development checkpoint

The user adopted the cleaned Geant4 11.4 source as the new `main` development
baseline and requested renaming `integration/drc-implant` to
`integration/drc-timing` for subsequent output-plane work. Output-plane code is
not part of this checkpoint. Research outputs and historical production bundles
are archived separately from the Git source tree.

The expanded pointer-initialization comparison used 1,000 events per condition
for 11.2/11.4, LG/noLG and initialization on/off: 8,000 stored events and 4,000
matched pairs. Stored payloads matched in every pair; valid paired CFD time
differences were zero. Uninitialized variants still crashed during shutdown.
The pointer fix is retained. These results do not identify a specific Geant4
class as the cause of the observed version-dependent timing distributions.
The historical limitations below remain recorded; adopting a development
baseline does not establish physics equivalence between versions.

The optional Precision executable and observer hooks were removed on 2026-10-06.
Initial precision-build results below are historical; current commands select
the normal simulation only. This architecture change does not restore Precision.

The removal passed a fresh default build and two-event LG/noLG regressions.
Normal-target preprocessed diagnostics and canonical events/diagnostics/path
histograms match the pre-removal candidate exactly. Boundary recovery remains
enabled; its removal requires a separate statistical and failure-case decision.

## Build requirements

`envset.sh` selects the installed Geant4 CMake config, Xerces-C 3.3.0 and all
15 dataset paths specified by that installation. The LCG view alone selects
Xerces-C 3.2.4 and leaves several Geant4 dataset variables old or unset.
The installation uses CLHEP 2.4.7.1; this records the actual deployed version,
not a claim that it equals every upstream validation configuration.

Geant4's config also requires EXPAT development headers, version >= 2.5.0.
This server has the matching 2.5.0 runtime but no expat-devel package installed.
For the checkpoint, the official AlmaLinux EL9 package
`expat-devel-2.5.0-6.el9_8.5.x86_64.rpm` was extracted inside the local study.
Its headers were paired with `/usr/lib64/libexpat.so.1`; no system installation
or CVMFS file was changed. The package's source, hash and MIT notice are retained
in that local dependency copy. Supply your own matching header path when building:

```bash
source envset.sh
cmake -S . -B build-g4-11.4-dev \
  -DCMAKE_BUILD_TYPE=Release -DWITH_GEANT4_UIVIS=ON \
  -DEXPAT_INCLUDE_DIR=/path/to/matching/expat/include \
  -DEXPAT_LIBRARY_RELEASE=/usr/lib64/libexpat.so.1
cmake --build build-g4-11.4-dev -j2
```

The Trigger boundary adapter remains enabled by default. For comparison builds
only, `-DCBDsim_USE_BOUNDARY_ADAPTER=OFF` selects standard Geant4 handling.
The adapter's version guard permits 11.2.0 and 11.4.0; its recovery algorithm,
true-path-length check and process-order checks are preserved. CMake requires
11.4.0 exactly for this development configuration. See
[boundary recovery](BOUNDARY-RECOVERY.md) for its limited validation scope.

## Checks completed

- The initial main and Precision targets compiled; after Precision removal,
  a fresh normal/default build including Reco passed.
  Runtime dependencies resolve to Geant4 11.4 and the selected ROOT/Xerces stack.
- Existing standard optical-boundary tests: seven modes, 5,000 photons per mode.
- Current shared LG/noLG geometry: overlap checks with 20,000 samples, navigation
  probes and 400 coupling photons per mode. No geometry/material/QE changes.
- Existing wafer-response tests: 3,000 photons per condition, both models,
  with photons launched at window/wafer, 450 nm and normal incidence.
- The 64-photon LG edge control reproduced 16 NoRINDEX losses with the standard
  boundary and zero with the adapter. Normal 200-photon controls in each model
  have identical diagnostics, path histograms and final RNG states.
- Serial 60 GeV pencil positron runs, two events per model, pass normal shutdown,
  primary audits, photon budgets, ROOT event/count/time-bin checks and re-reading.
  Standard/adapter canonical event digests match within each model and seed.

The initial full runs exposed a pre-existing uninitialized `fChain` during
ROOT shutdown. GDB localized the fault to `CBDsimRootInterface::close()`.
Initializing the four pointer members to null fixes the observed shutdown;
the event schema, output physics payload and ownership policy were not changed.

## 2026-10-07 distribution comparison and remaining gates

The full 12,000-incident-event diagnostic comparison is complete: 3,000 per
version and LG/noLG mode, W40 × L60 × T5 mm, central 60 GeV pencil positrons,
Serial, unchanged geometry/material/QE/boundary model. M1 remains incomplete.
The whole 11.2/LCG105/GCC13/ROOT6.30 versus 11.4/LCG107/GCC11/ROOT6.34 stack is
compared; this does not isolate Geant4 as the sole cause.

| Pair timing width | 11.2 | 11.4 | Relative change, bootstrap 95% CI |
|---|---:|---:|---|
| LG σ(T1−T2) | 105.75 ps | 96.04 ps | −9.18%, [−12.48%, −5.92%] |
| noLG σ(T1−T2) | 119.30 ps | 96.82 ps | −18.84%, [−21.68%, −15.61%] |

Mean detected photons decrease by roughly 2–4%. Sample standard deviations and
16–84% widths show the same direction. Physics preservation has not been
established. Explain the response shifts before accepting a new baseline; do
not tune QE/geometry simply to recover old numbers. Proposed meaningful-change
margins have not yet been agreed with the user.

The 11.2 LG sample at seed 107710056, event 21 produced GeomNav1002 in the light
guide. Its 50-event ROOT/budget checks passed, but navigation validation failed.
The diagnostic report includes this sample without promoting its failed gate.
All 12,000 ROOT events remain available. Ten events contain one photon in the
[240,99999] ns overflow bin; the unchanged timing algorithm rejects that channel.
There are 11,990 valid timing pairs, and all incident events remain in efficiency
denominators. Warning impact and this timing-selection rule require review.

Additional tests passed: both version builds/datasets, primary/Gaussian 20,000,
14 existing analysis tests, repeated-run ROOT output, geometry reinitialization,
two-thread MT ROOT/budget, and UI/Vis OFF build with identical two-event noLG
canonical payloads compared to ON. Repeated beamOn diagnostic summaries contain
only the last run. Interactive GUI and other dimension-specific models are
outside this validation. Precision remains removed.

The 2026-10-06 checkpoint did not submit or stop jobs. The 2026-10-07 comparison
used new DAG 9731883; it finished with exit 1 because of the navigation warning
(239 distribution nodes passed, one rejected, four smokes passed). Thirty-six
jobs delayed on one worker were rescheduled with unchanged config/seeds and all
attempts preserved. The 850 events available before rescheduling have identical
canonical payloads in the completed replays. Original protected Condor inputs
and directories remain unchanged. No commit or push was performed.

Detailed results and uncertainty/limitations are in the private local study
`.local-context/studies/2026-10-07-chained-simulation/comparison-diagnostic/README.md`.
The current workspace plan is `../docs/integration-preparation/CHAIN_PLAN.md`.
The earlier build/control evidence remains in
`.local-context/studies/2026-10-06-geant4-migration/`.
