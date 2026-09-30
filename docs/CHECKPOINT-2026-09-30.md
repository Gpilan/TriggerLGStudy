# September30 scan checkpoint

This checkpoint records the shared source and selected completed scans before further research. It supersedes the September16 repository snapshot as a code checkpoint, not as a uniform new physics baseline.

## Shared source changes

- noLG coupling aperture is min(15mm,W−0.2mm):9.8mm for W10 and14.8mm for W15, preserving tape/foil edge allowance.
- Signed assembly positions accept−29.5..29.5mm; negative positions approach the sensor. Frozen-run input validation uses the same range.
- LG triangular-facet intersections are checked against the finite facet to reject accepted hits outside the triangle without changing geometry tolerance.
- Corner packaging uses equivalent touching box components, avoiding nested-Boolean corner ambiguity.
- The Geant4 11.2-specific boundary adapter preserves process ordering and applies narrowly guarded opaque-entry/rim-normal recovery. See [boundary recovery](BOUNDARY-RECOVERY.md).
- A reusable local results-portal generator organizes explicitly selected validated figures by topic. Private portal records remain excluded from Git.

## Completed research

See [selected scan results](../analysis/results/scan_checkpoint_20260930/README.md). Position scan now has3000 valid timing pairs for every condition. Tables, representative ROOT-rendered Δt plots, timing/count trends and waveform diagnostics are included.

The September23 noLG−25 tape-entry extension remains a [study patch](../patches/README.md); this GitHub update does not promote it into shared simulation or rebuild the user's GUI binary. Production conditions used multiple frozen versions, documented in the results.

## Validation

Historical path-specific regression evidence and its limits are documented in BOUNDARY-RECOVERY.md. The release check compiles the exact shared CBDsim source in an isolated build and runs the existing Python tests. The isolated shared-source CBDsim build passed with LCG105/Geant4 11.2.0. All10 existing Python tests passed (1 scripts +9 analysis). No production jobs are submitted by this checkpoint.

## Position-range limitation found during checkpoint checks

The shared default W40 geometry at s=−25mm reports real inter-assembly scintillator/foil–LG overlaps (up to about1.5mm). This is outside the validated W10 crossing-position campaign. The accepted signed input range is a parser bound, not a promise of valid geometry for every width/position combination. Keep shared W40 at its validated nominal configuration unless a new placement is explicitly checked. This publication does not change that geometry or claim this test passed.

Nominal W40 at s=0: both LG and noLG passed initialization and the10000-point-per-volume overlap check. These are initialization/geometry checks; no new production sample was generated.
