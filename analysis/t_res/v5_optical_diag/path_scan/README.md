# Optical photon path-length scan (v5 geometry)

Trigger size **1×1 … 4×4** (`kHxWide` = 5 … 20 mm) 에서 optical photon **total track length** 분포와 **fate** 비율을 비교합니다.

## 활성화

```bash
export CBDsim_OPTICAL_DIAG=1
export CBDsim_OPTICAL_DIAG_HIST=1
export CBDsim_OPTICAL_DIAG_HIST_OUT=.../path_4x4_100ev.path_hist.txt
```

## 파이프라인

```bash
analysis/t_res/v5_optical_diag/path_scan/run_path_scan.sh
```

- 100 ev × 7 sizes, 각 size마다 rebuild (`kHxWide` patch)
- 산출: `path_<tag>_100ev.path_hist.txt`, `path_<tag>_100ev.txt`
- 플롯: `analysis/plot_optical_path_scan.py` → `figures/`, `path_vs_size_summary.csv`

## Fate 정의

| fate | 의미 |
|------|------|
| `detected` | SiPM QE 통과 |
| `qe_reject` | wafer 도달, QE 탈락 |
| `absorb_scint` / `absorb_lg` | bulk OpAbsorption |
| `kill_world` | world/foil/env 에서 소멸 |
| `kill_other` | scint/lg/gel 등 기타 |

Histogram: **400 bins, 0–2000 mm** (total track length per photon).

## 결과 요약

스캔 완료 후 `PATH_SCAN_RESULTS.md` 참고.
