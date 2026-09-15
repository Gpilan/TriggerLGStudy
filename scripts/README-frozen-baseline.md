# Frozen technical baselines (R03)

Run from the repository after `source envset.sh`. This builds from a copied source
snapshot, including dirty source, and never rebuilds the shared Condor executable.
Choose new paths/IDs: existing bundles or runs are rejected rather than overwritten.

```sh
python3 scripts/frozen_baseline.py prepare /absolute/new-bundle
python3 /absolute/new-bundle/source/scripts/frozen_baseline.py run /absolute/new-bundle lg42 --mode lg --seed 42 --events 2
python3 /absolute/new-bundle/source/scripts/frozen_baseline.py run /absolute/new-bundle lg42-replay --mode lg --seed 42 --events 2
python3 /absolute/new-bundle/source/scripts/frozen_baseline.py run /absolute/new-bundle nolg42 --mode nolg --seed 42 --events 2
```

The supported request is deliberately narrow: 60 GeV electrons, one primary per
event, nominal G01 geometry, central +z pencil beam, one MT worker. General beam
control and final physics/timing analysis remain separate tasks. This wrapper does
not submit Condor jobs or select a physics winner.

Each bundle records HEAD, binary dirty/index patches, source-file hashes, build
commands/logs, environment/data-directory paths and checksums of linked libraries.
Each run has its own copied executable/rootIO dictionary, QE CSV at the expected
relative path, macro, explicit LG mode, diagnostic outputs, ROOT output and status
manifest. The runtime must resolve the private rootIO library; other runtime library
hashes must belong to the frozen build environment, and Geant4 data paths must match.
CVMFS libraries/data are recorded external dependencies, not vendored: the same LCG
installation and data paths are still required. Source and executable inputs are
checksum-verified; this is not an OS-enforced immutable storage system.

Use the **copied runner** under `bundle/source/scripts` for replay. Input changes
fail verification. Each event must be present in ROOT and the diagnostics, with
matching detected counts, exact stored central-beam kinematics and a closed optical
budget. NoRINDEX/unknown/unfinished or diagnostic exceptions fail the run. ROOT and
nonzero time-bin payloads are serialized canonically to validate reproducibility;
ROOT file checksums themselves need not agree because ROOT headers contain metadata.
Compare `runs/<id>/validation.json:event_digest_sha256` for same-seed replay.

`manifest.json` progresses through preparing/running to complete or failed. A
process killed externally may leave a running manifest; that is not a completed run.
Keep failure logs and choose a new run ID for retries. Analysis CLI is explicitly
null while only this technical validation is performed; the validation CLI is saved.
Final baseline statistics and timing claims require the remaining checklist work.

B01 이후 `run`은 `--position-mm x y z`와 `--direction dx dy dz`를 지원합니다.
방향을 정규화하고 요청값을 macro/manifest에 남깁니다. ROOT schema는 float32이므로
최신 validator는 그 저장 정밀도를 반영합니다. 이전 frozen bundle의 validator를
자동 교체하지 않으므로 일반 기울기 실행에는 최신 소스로 새 bundle을 prepare하세요.
