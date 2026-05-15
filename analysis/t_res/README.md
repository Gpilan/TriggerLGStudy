# t_res (타이밍 분석 베이스라인)

- **`data/`** — 시뮬레이션 산출 `.root` (깃에는 안 올림; 루트 `.gitignore`의 `*.root`)
- **`figures/`** — `plot_trigger_timing.py` 기본 PNG 출력
- **`archive/version_1_2026-05-15/`** — 기존 v1 분석 결과 백업(`data/`, `figures/`)

## 새 ROOT 파일 재분석 준비 (v2)

1. 새 시뮬 결과 `.root` 파일을 `analysis/t_res/data/` 에 복사
2. 저장소 루트에서 환경 설정
3. 기존 분석 스크립트를 동일하게 실행

```bash
source envset.sh

# 기본 타이밍 플롯
python3 analysis/plot_trigger_timing.py

# LG vs no-LG 비교 (CFD 포함)
python3 analysis/compare_timing_resolution.py

# 필요 시 추가 분석
python3 analysis/plot_event_photon_time.py analysis/t_res/data/<new_file>.root --event 0
python3 analysis/plot_overlay_photon_time.py analysis/t_res/data/<new_file>.root
python3 analysis/plot_cfd_high_nphoton_timing_ROOT.py analysis/t_res/data/<new_file>.root --fraction 0.3 --mode single --top-events 3
```

## 분석

저장소 루트에서:

```bash
export LD_LIBRARY_PATH=$PWD/build/rootIO:$LD_LIBRARY_PATH
python3 analysis/plot_trigger_timing.py
python3 analysis/plot_trigger_timing.py noLG1_0.root --ref min
```

기본 입력은 `analysis/t_res/data/` → `build/CBDsim/stats/` → **`build/CBDsim/`** 순으로 `.root` 를 찾습니다 (`noLG1_0.root` 우선, 없으면 `An_test_1.root`, 그다음 첫 `.root`).  
기본 출력은 `analysis/t_res/figures/` (`timing_<stem>.png`).  
경로 바꿀 때: `CBDsim_STATS`, `CBDsim_FIGURES` 환경변수.

시뮬 예: `cd build/CBDsim && ./CBDsim run_ele_60gev.mac 0 noLG1` → `noLG1_0.root` (현재 기본 지오메트리는 Proto, LG 없음).

### LG vs no-LG Δt 비교 (타이밍 레졸루션 = σ/√2)

```bash
python3 analysis/compare_timing_resolution.py
```

기본 입력: `data/60GeV_e-_LG_0.root`, `data/60GeV_e-_noLG_0.root`.  
두 파일 트리 엔트리가 다르면 **작은 쪽(min)** 만큼만 각각 읽어 **동일 N**으로 비교합니다.

→ `figures/compare_delta_timing_LG_vs_noLG.png`, 같은 stem `.csv` (Gauss σ, σ/√2, `tree_events_used` 등).

### 한 이벤트 포톤 시간 분포 (`event_*_photon_time.png`)

**한 이벤트**만 골라, T1/T2 각각 **SiPM에 합쳐진 도착 시간**을 ROOT `HIST` 로 그립니다.  
**고정 폭 리빈(기본 20 ps)** 이며, **x 축은 포톤이 있는 시간대만**(여유 마진 포함) 자동 범위.

```bash
python3 analysis/plot_event_photon_time.py analysis/t_res/data/60GeV_e-_noLG_0.root --event 0
# python3 analysis/plot_event_photon_time.py ... --bin-width-ps 25
```

→ `figures/event_<stem>_ev<N>_photon_time.png`

### 이벤트 평균 포톤 시간 (`mean_photon_time_*.png`)

읽는 모든 이벤트에 대해 **빈마다 포톤 수의 산술 평균** TH1F (오버레이 곡선 없음).  
`--trig both` 이면 T1/T2 **같은 시간 축**(한쪽만 신호가 있어도 양쪽 패드에 그려짐).

```bash
python3 analysis/plot_overlay_photon_time.py analysis/t_res/data/60GeV_e-_noLG_0.root
```

→ `figures/mean_photon_time_<stem>.png`
