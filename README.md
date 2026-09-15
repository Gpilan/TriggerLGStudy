# TriggerLGStudy

Geant4 simulation and ROOT/Python analysis of optical photon collection and timing in a two-scintillator trigger. **LG and noLG share one implementation** and are selected through run configuration.

## Start here

- [Simulation](CBDsim/) and [geometry validation](CBDsim/tests/proto_geometry/README.md)
- [Frozen build and run workflow](scripts/README-frozen-baseline.md): explicit LG/noLG mode, seeds, beam settings, source and binary hashes
- [Optical diagnostics](CBDsim/tests/optical_boundary/README.md) and [sensor response assumptions](CBDsim/tests/sensor_response/README.md)
- [SPE waveform model](analysis/README-spe-response.md) and [analysis code](analysis/)
- [Version management](docs/VERSIONING.md) and [historical snapshot index](docs/SNAPSHOTS.md)

On the KNU environment, load the existing CERN LCG view with `source envset.sh`. This setup depends on CVMFS; it is not a portable dependency installer. Use the frozen workflow to build and run in new output directories without replacing a shared executable.

[2026-09-15 technical baseline validation](docs/BASELINE-2026-09-15.md).

## Versions and run conditions

`main` is the integrated code line. Short-lived `fix/`, `study/`, and `docs/` branches hold ongoing work. Completed historical versions live under `snapshot/*` tags; tested reference points use `baseline/*` tags. Tags preserve code, while run manifests identify geometry, LG/noLG mode, beam, seeds, sample size, binary and analysis settings.

A validation baseline records which technical checks passed. It does not certify the complete physical detector model. The current sensor response is an idealized conditional wafer-QE model. Absolute PMT timing, collection efficiency, electronics and other real-device effects require further validation.

For comparisons, distinguish Gaussian fit sigma from quantile widths, retain timing failures in efficiency denominators, and verify identical conditions beyond the LG/noLG selection. Existing historical figures and logs retain their original assumptions. Plot output defaults to PNG; add PDF only when explicitly needed.
