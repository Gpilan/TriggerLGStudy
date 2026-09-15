# Prototype geometry validation (G01)

This standalone executable compiles the **production** prototype geometry,
materials, optical diagnostics, and sensor SD. It does not create ROOT output or
rebuild the usual `build/CBDsim/CBDsim` executable.

```sh
source envset.sh
cmake -S CBDsim/tests/proto_geometry -B /tmp/trigger-geometry-test
cmake --build /tmp/trigger-geometry-test -j2
# Arguments: mode, samples/events, seed, summary path, noLG (0 or 1)
/tmp/trigger-geometry-test overlap 20000 42 /tmp/unused.txt 0
/tmp/trigger-geometry-test geometry 0 42 /tmp/unused.txt 0
/tmp/trigger-geometry-test photons 2000 42 /tmp/lg-photons.txt 0
# Repeat all three with last argument 1 for no-LG; repeat seed 43.
```

- `overlap`: checks every physical placement with a mother, including nested wafer,
  both assembly placements and all foil pieces. Calls `CheckOverlaps(samples, 0,
  true, 20)`; any reported overlapping placement fails. Logs transforms/bounds.
  Enumerating the physical-volume store checks each placed instance, unlike only
  checking the world volume. Shared daughter placements are checked relative to
  their common mother solid.
- `geometry`: uses a separate Geant4 navigator to sample both sides of the actual
  coupling interfaces at 1/10/100 nm offsets, including edges/corners in both
  assemblies. Scint, inlet gel, LG, outlet gel, window, and wafer must be correctly
  identified. The LG cap must have the full 200 mm² rectangular inlet; the outlet
  polygon area must be within 0.05% of the 15 mm circle and sag below 3 μm. Signed
  mesh volume is reported for independent resolution comparisons.
- `photons`: shoots 3 eV optical photons along the local coupling axis from a
  10×10 near-edge grid in both scintillators (LG: full inlet; noLG: sensor aperture).
  Uses real absorption/boundary processes and the production sensor SD/QE.
  First boundary must be scint→gel. Every event budget must close, with zero
  unknown, NoRINDEX, unfinished tracks or >10,000-step guard terminations.
  This is a transport correctness test, not a sensor calibration or timing scan.

An optional seventh argument supplies the world maximum extent in mm **only to
this test**, before any solid is constructed. It exercises Geant4 surface-tolerance
sensitivity (e.g. 1000/10000 mm give 1e-8/1e-7 mm instead of default 1e-9 mm).
The overlap reporting tolerance remains 0 mm, and production defaults are unchanged.

For source-before regression or mesh resolution studies, configure a separate
build with `-DPROTO_SOURCE=/absolute/path/to/snapshot.cc`. Record the snapshot,
seed, log, binary and source hashes. Do not modify the working geometry while a
comparison is running. A finite overlap sample cannot prove correctness for all
future dimensions: rerun these checks whenever dimensions/topology are changed.
