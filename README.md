# TriggerLGStudy

Geant4 simulation of a two-scintillator trigger, comparing light collection and timing with and without light guides (LG/noLG).

Environment: **Geant4 11.4.0, ROOT 6.34.02, C++17**. `envset.sh` selects the CERN LCG/CVMFS environment used on the KNU server.

| Path | Purpose |
|---|---|
| `CBDsim/` | Simulation, geometry and optical transport; validation in `CBDsim/tests/` |
| `rootIO/` | ROOT event data and file I/O |
| `Reco/` | ROOT event reader executable |
| `analysis/` | Photon-time analysis, SPE reconstruction and CFD timing |
| `scripts/` | Reproducible runs, ROOT validation and Condor template |

- [Build requirements](docs/GEANT4-MIGRATION.md#build-requirements) · [Frozen run workflow](scripts/README-frozen-baseline.md) · [Condor template](scripts/condor_template.sub)
- [SPE/CFD model](analysis/README-spe-response.md) · [LED selection](analysis/README-LED.md)
- [Sensor model](CBDsim/tests/sensor_response/README.md) · [QE data source](analysis/reference/pmt_r2076/README.md)

The sensor response and SPE waveform are idealized models; their assumptions are documented above. Research results and run histories are maintained separately from the code documentation.
