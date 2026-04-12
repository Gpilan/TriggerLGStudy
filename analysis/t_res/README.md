# t_res (타이밍 분석 베이스라인)

- **`data/`** — 시뮬레이션 산출 `.root` (깃에는 안 올림; 루트 `.gitignore`의 `*.root`)
- **`figures/`** — `plot_trigger_timing.py` 기본 PNG 출력

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
