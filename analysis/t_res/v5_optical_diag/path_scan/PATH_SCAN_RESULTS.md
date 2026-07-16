# Optical photon path-length scan (v5 geometry)

**일시**: 2026-07-06  
**조건**: 60 GeV e⁻, 100 ev/size, v5 geometry (foil shell + tip ring + gel + SiPM coupling)  
**스캔**: trigger size 1×1 … 4×4 (`kHxWide` = 5 … 20 mm)

## 산출물

| 파일 | 설명 |
|------|------|
| `path_<tag>_100ev.path_hist.txt` | fate별 total track-length histogram (400 bins, 0–2000 mm) |
| `path_<tag>_100ev.txt` | 기존 boundary/path 합계 diag |
| `path_vs_size_summary.csv` | size별 fate 비율·mean path |
| `figures/fate_fraction_vs_size.png` | fate 비율 stacked bar |
| `figures/mean_path_vs_size.png` | mean path vs size |
| `figures/total_path_overlay_1x1_2x2_4x4.png` | total path 분포 overlay |

## 활성화 방법

```bash
export CBDsim_OPTICAL_DIAG=1
export CBDsim_OPTICAL_DIAG_HIST=1
export CBDsim_OPTICAL_DIAG_HIST_OUT=.../path_4x4_100ev.path_hist.txt
```

## 핵심 결과 (100 ev 평균)

### 1. 전체 mean path ↑ (trigger 커질수록 transport 부담 증가)

| size | mean path all (mm) | total tracks /100 ev |
|------|-------------------|----------------------|
| 1×1 | **130** | 1.90M |
| 2×2 | **158** | 1.98M |
| 3×3 | **209** | 1.95M |
| 4×4 | **253** | 1.91M |

1×1 → 4×4: mean path **약 1.9×**. photon 수는 비슷하지만 **한 광자당 이동 거리**가 크게 늘어남.

### 2. Fate 비율 — “어디로 가는가”

| fate | 1×1 | 4×4 | 해석 |
|------|-----|-----|------|
| **detected** | 10.4% | 5.5% | SiPM QE 통과 — **비율 절반** |
| **qe_reject** | 38.4% | 20.3% | wafer 도달 후 QE 탈락 |
| **absorb_scint** | 4.2% | 7.2% | PS bulk 흡수 ↑ |
| **kill_world** | 34.7% | 35.3% | world leak — **비율 거의 일정** |
| **kill_other** | 12.2% | 31.6% | LG/scint 등 — **4×4에서 급증** |

**요약**: trigger가 커지면 **detected + qe_reject 합(≈wafer 도달)** 이 49% → 25%로 줄고, **kill_other(LG/scint boundary 등)** 가 12% → 32%로 늘어남. world leak 비율은 ~35%로 유지.

### 3. Fate별 mean path (mm)

| fate | 1×1 | 4×4 |
|------|-----|-----|
| detected | 148 | **218** |
| kill_world | **79** | 208 |
| kill_other | 181 | **317** |
| absorb_scint | 179 | **316** |

- **detected photon**도 4×4에서 path **~50% 증가** (timing tail에 영향 가능).
- **kill_world**는 1×1에서 mean 79 mm (짧은 leak) → 4×4에서 208 mm (긴 ping-pong).
- **kill_other** tail이 4×4에서 매우 김 (facet/wedge grazing).

### 4. v4 LG scan과의 연결

- v4에서 4×4는 SiPM 수집·timing이 1×1 대비 불리했음.
- 이번 scan: **detected fraction 10% → 5%**, **mean detected path +47%** → LG 크기 증가 시 **수집 효율 + transport 길이** 동시 악화.
- `kill_world` ~35% 유지 → v5 foil patch 후에도 **world fate 비율**은 size에 거의 무관 (절대 개수는 증가).

## Condor 전 시사점

- v5 geometry로 **world path 합은 대폭 감소**(이전 diag)했으나, **fate 분해**에서는 여전히 ~35%가 `kill_world`.
- 4×4 Condor는 **accept** 가능하나, MPV/σ trend는 **detected path lengthening**과 **detected fraction 감소**를 함께 봐야 함.
- 잔여 vis leak (PS 2 + facet 3–4)는 `kill_world` / `kill_other` tail의 일부.

## 재실행

```bash
analysis/t_res/v5_optical_diag/path_scan/run_path_scan.sh
```

기존 `path_*_100ev.path_hist.txt`가 있으면 해당 size는 skip.
