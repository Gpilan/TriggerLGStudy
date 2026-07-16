# V5 overnight run log

> 자동 파이프라인: `run_v5_auto.sh`  
> 비교 스크립트: `../../compare_v5_optical_diag.py`  
> 진단 활성화: `export CBDsim_OPTICAL_DIAG=1`

## 목표 (5단계)

| # | 내용 | 산출물 |
|---|------|--------|
| ① | Optical stepping fate 집계 (100 ev) | `diag_4x4_baseline_100ev.txt` |
| ② | LG inlet kNPhi 64→256 | `diag_4x4_nphi256_100ev.txt` |
| ③ | kLGZGap 5 μm / 10 μm | `diag_4x4_nphi256_gap5um_100ev.txt`, `..._gap10um_...` |
| ④ | 1×1 vs 4×4 (diag + MPV) | `diag_1x1_nphi256_100ev.txt`, `v5_compare_summary.csv` |
| ⑤ | foil hollow shell (필요 시) | RUNLOG §⑤ + geometry patch |

## 도미넌트 판별 기준

- **A. scint↔world leak** (wedge/corner): `scint_to_world`, `world_to_scint` 높음
- **B. LG side leak**: `lg_to_world`, `world_to_lg` 높음
- **C. coupling loss**: `scint_to_lg` 낮음, path_scint 높음
- **D. SiPM/QE**: `sipm_detect` vs `sipm_qe_reject`
- **E. navigation/gap**: ③에서 gap 변경 시 boundary 카운트 변화

---

## ① baseline (kNPhi=64, kLGZGap=0, kHxWide=20 mm) — DONE

- **완료**: 2026-07-05 ~18:37 KST, wall ~21 min
- **파일**: `diag_4x4_baseline_100ev.txt`, `diag_4x4_baseline_100_0.root`

| metric | value / event |
|--------|----------------|
| scint→world | 505,932 |
| world→scint | 498,974 |
| scint→LG | 17,316 |
| LG→world | 79,763 |
| LG→SiPM | 5,727 |
| SiPM detect | 1,104 |
| path scint (mm/ev) | 4.32e6 |

**1차 결론**: **A (scint↔world wedge leak) 압도적**. scint↔world boundary가 scint→LG의 ~30×.

---

## ② kNPhi=256 — DONE

- **완료**: 2026-07-05 ~19:11 KST, wall ~33 min
- **파일**: `diag_4x4_nphi256_100ev.txt`, `diag_4x4_nphi256_100_0.root`

| metric | /event | vs ① delta |
|--------|--------|------------|
| scint→world | 526,124 | **+4.0%** (악화) |
| scint→LG | 18,452 | +6.6% |
| LG→SiPM | 5,835 | +1.9% |
| SiPM detect | 1,129 | +2.3% |

**결론**: kNPhi 64→256만으로는 **A(wedge leak) 해결 안 됨**. scint↔world 오히려 소폭 증가. ③ gap scan 계속 진행.

---

## ③ kLGZGap scan — DONE (gap=0 유지)

- **gap 5 μm** (~149 min): scint→LG **0**, SiPM detect **1/ev** → **coupling 붕괴** (navigation/분리)
- **gap 10 μm** (~145 min): 동일 — photon이 LG에 못 들어감
- **best 4×4**: `kLGZGap=0` (`diag_4x4_nphi256_100ev.txt`)

**결론**: gap은 **E(navigation) 실험**이었고, 5–10 μm는 **쓸 수 없음**. flush(0) 유지.

---

## ④ 1×1 vs 4×4 — DONE

| /event | 1×1 (hx=5 mm) | 4×4 nphi256 |
|--------|---------------|-------------|
| scint→world | 39.3M | 52.6M |
| scint→LG | 1.18M | 1.85M |
| SiPM detect | **203k** | 113k |
| path scint (mm) | 2.21M | 4.48M |

**결론**: 1×1은 면적 작아 scint↔world **절대값**은 낮지만, 4×4는 wedge 비율이 커서 leak 더 큼. SiPM 수집은 1×1이 유리 (기존 trend와 일치).

---

## ⑤ foil hollow shell — DEFERRED

- `scint_to_world / scint_to_lg = 28.5` (4×4 nphi256) → 여전히 **A 압도**
- kNPhi↑, gap μm **해결책 아님**
- **다음 작업**: scint-only wedge 메우기 / inlet 4-corner vertex / Al hollow shell (subtraction) — **낮에 geometry PR**

---

## 최종 요약

- **dominant mechanism**: **A — scint↔world wedge/corner leak** (전 단계 공통). ② nphi256 +4%, ③ gap은 coupling 파괴.
- **best geometry config**: `kNPhi=256`, `kLGZGap=0`, `kHxWide=20 mm` (V4 4×4 baseline 대비 nphi256 소폭 악화지만 gap 대비 유일 viable)
- **recommended next (Condor 3000 ev)**: geometry fix **전** large run 보류. 우선 wedge/corner patch prototype → 100 ev diag로 scint→world **>50% 감소** 확인 후 Condor.

---

## foil patch v5.1 (2026-07-06)

**변경** (`CBDsimDetectorConstructionProto.cc`):

1. **5개 분리 foil box** → **`G4SubtractionSolid` hollow shell** (측면 + beam +y, **-y 개방** LG coupling 유지)
2. **4 corner pad** (`protoFoilCorner*`) — scint -y face 모서리 wedge (LG polygon 밖 scint-only 구간) Al cover
3. `kFoilCornerInset` 제거 — side foil corner gap 제거

**검증 (100 ev diag, 수동)**:
```bash
cd Trigger && source envset.sh && cd build/CBDsim
export CBDsim_OPTICAL_DIAG=1
export CBDsim_OPTICAL_DIAG_OUT=$TRIGGER_ROOT/analysis/t_res/v5_optical_diag/diag_4x4_foilpatch_100ev.txt
./CBDsim run_diag_100ev.mac 0 $TRIGGER_ROOT/analysis/t_res/v5_optical_diag/diag_4x4_foilpatch_100
```
baseline(`diag_4x4_nphi256_100ev.txt`)과 `scint_to_world` 비교.

### foilpatch 100 ev 결과 (2026-07-06 ~12:30 KST)

| metric /event | nphi256 | foilpatch | delta |
|---------------|---------|-----------|-------|
| scint→world | 526,124 | 533,869 | **+1.5%** (거의 동일) |
| scint→LG | 18,452 | 18,720 | +1.5% |
| SiPM detect | 1,129 | 1,140 | +1.0% |
| path world (mm) | 124,198 | **38,460** | **−69%** |
| other_boundary | 18,064,000 | 18,524,000 | +2.5% |

**해석**: `scint_to_world` **boundary step 수는 안 줄었음** (통계 오차 수준). 다만 **world 안 path length 69% 감소** → ping-pong 깊이는 줄었을 가능성. scint→foil은 `other_boundary`로만 잡혀 `scint_to_world`에 반영 안 됨. wedge-only corner pad(1.1 mm)는 4×4 leak 대비 **부족**할 수 있음 → inlet 전체 ring foil 또는 corner vertex fix 검토.

### v5.2 tip ring foil (2026-07-06)

**변경**: `protoFoilTipRing` — `G4Tubs` annulus at LG outlet (`R_in ≈ kRguide·cos(π/kNPhi)`, `R_out = kRguide + gap + foil`), y at tip–SiPM junction, `AluminumSurf`.

**검증**:
```bash
export CBDsim_OPTICAL_DIAG=1
export CBDsim_OPTICAL_DIAG_OUT=$TRIGGER_ROOT/analysis/t_res/v5_optical_diag/diag_4x4_tipring_100ev.txt
./CBDsim run_diag_100ev.mac 0 .../diag_4x4_tipring_100
```
목표: `lg_to_world` ↓, `lg_to_sipm` / `sipm_detect` 유지.

### v5.3 SiPM coupling fix (2026-07-06)

**변경** (① tip ring 유지):
- `protoSipmWindowPhys` → **world에 flush** (`yWindowCenter = yLgTip − gap − windowHalf`)
- `protoSipmEnvPhys` → **wafer만** 감싸는 얇은 vacuum (full stack cylinder 제거)
- `LogicalBorderSurface`: **LG ↔ window** 만 `AirSurf` (LG ↔ env side **제거**)

**검증**: `diag_4x4_couplingfix_100ev.txt`

### v5.4 index-matching gel (2026-07-06)

**변경** (① tip ring + ② coupling 유지):
- `protoSipmGelPhys` — `Gelatin` n≈1.52, `kSiPMGelT=50 μm` between LG tip and window
- Stack: LG → gel → glass → wafer
- `Gelatin` bulk: RINDEX + ABSLENGTH 100 cm

**검증**: `diag_4x4_gel_100ev.txt`

### v5.5 LG–SiPM gap closure (2026-07-06)

**원인**: vis에서 큰 갭처럼 보이는 것 — (1) taper 옆 **vacuum mantle** (2) gel/window 얇은 원통 vs LG 원뿔 (3) tip ring 너무 얇음 (8 μm).

**변경** (open tip은 G4 tessellated 오류로 **cap 유지**):
- `kLGTipGelOverlap = 20 μm` — gel top이 LG tip plane에 **살짝 침투**
- `kSiPMGelT = 100 μm` (50→100)
- `kTipRingHalfY = 0.3 mm` — tip Al ring 축방향 확장
- ~~`kLgTipCylinderL`~~ **제거** — `y>=yCylStart` 버그로 **전체 LG가 R=7.5 원통**이 되어 tessellated taper 소멸. 원복함.

로그: `[Proto geometry] LG tip y=... gel top y=...`


**Compare ① vs ②** (`diag_4x4_baseline_100ev.txt` vs `diag_4x4_nphi256_100ev.txt`):
```
baseline: diag_4x4_baseline_100ev.txt
variant:  diag_4x4_nphi256_100ev.txt
key                                baseline        variant     delta%
----------------------------------------------------------------------
scint_to_world                    5.059e+07      5.261e+07       3.99
world_to_scint                     4.99e+07      5.194e+07       4.09
scint_to_lg                       1.732e+06      1.845e+06       6.56
lg_to_world                       7.976e+06      8.346e+06       4.63
lg_to_sipm                        5.727e+05      5.835e+05       1.90
sipm_detect                       1.104e+05      1.129e+05       2.27
path_scint_mm_per_event           4.316e+06      4.481e+06       3.81
path_lg_mm_per_event              7.472e+05      7.935e+05       6.19
```

### run 4x4_nphi256_gap5um
- finished: 2026-07-05 21:41:00 KST, wall 8931s
- file: `/u/user/rmsvlf000/Trigger/analysis/t_res/v5_optical_diag/diag_4x4_nphi256_gap5um_100ev.txt`
```
events 100
boundary_total 1532362017
scint_to_lg 0
lg_to_scint 0
scint_to_world 752765400
world_to_scint 751656320
lg_to_world 2042917
world_to_lg 2206834
lg_to_sipm 0
sipm_to_lg 0
other_boundary 23690546
bulk_absorb_scint 878390
bulk_absorb_lg 0
bulk_absorb_glass 0
bulk_absorb_si 2
sipm_detect 1
sipm_qe_reject 2
path_scint_mm_per_event 5.86515e+07
path_lg_mm_per_event 0
path_glass_mm_per_event 0.00669473
path_si_mm_per_event 0.000111687
path_world_mm_per_event 191486
```


**Compare ② vs ③ gap5um** (`diag_4x4_nphi256_100ev.txt` vs `diag_4x4_nphi256_gap5um_100ev.txt`):
```
baseline: diag_4x4_nphi256_100ev.txt
variant:  diag_4x4_nphi256_gap5um_100ev.txt
key                                baseline        variant     delta%
----------------------------------------------------------------------
scint_to_world                    5.261e+07      7.528e+08    1330.78
world_to_scint                    5.194e+07      7.517e+08    1347.15
scint_to_lg                       1.845e+06              0    -100.00
lg_to_world                       8.346e+06      2.043e+06     -75.52
lg_to_sipm                        5.835e+05              0    -100.00
sipm_detect                       1.129e+05              1    -100.00
path_scint_mm_per_event           4.481e+06      5.865e+07    1208.93
path_lg_mm_per_event              7.935e+05              0    -100.00
```

### run 4x4_nphi256_gap10um
- finished: 2026-07-06 00:06:26 KST, wall 8712s
- file: `/u/user/rmsvlf000/Trigger/analysis/t_res/v5_optical_diag/diag_4x4_nphi256_gap10um_100ev.txt`
```
events 100
boundary_total 1454870102
scint_to_lg 0
lg_to_scint 0
scint_to_world 714780051
world_to_scint 713730273
lg_to_world 1911146
world_to_lg 2064587
lg_to_sipm 0
sipm_to_lg 0
other_boundary 22384045
bulk_absorb_scint 834925
bulk_absorb_lg 0
bulk_absorb_glass 0
bulk_absorb_si 0
sipm_detect 0
sipm_qe_reject 0
path_scint_mm_per_event 5.56033e+07
path_lg_mm_per_event 0
path_glass_mm_per_event 0
path_si_mm_per_event 0
path_world_mm_per_event 199913
```


**Compare ② vs ③ gap10um** (`diag_4x4_nphi256_100ev.txt` vs `diag_4x4_nphi256_gap10um_100ev.txt`):
```
baseline: diag_4x4_nphi256_100ev.txt
variant:  diag_4x4_nphi256_gap10um_100ev.txt
key                                baseline        variant     delta%
----------------------------------------------------------------------
scint_to_world                    5.261e+07      7.148e+08    1258.58
world_to_scint                    5.194e+07      7.137e+08    1274.13
scint_to_lg                       1.845e+06              0    -100.00
lg_to_world                       8.346e+06      1.911e+06     -77.10
lg_to_sipm                        5.835e+05              0    -100.00
sipm_detect                       1.129e+05              0    -100.00
path_scint_mm_per_event           4.481e+06       5.56e+07    1140.90
path_lg_mm_per_event              7.935e+05              0    -100.00
```


## best 4x4 after ③
- file: `diag_4x4_nphi256_100ev.txt`

### run 1x1_nphi256
- finished: 2026-07-06 00:29:07 KST, wall 1348s
- file: `/u/user/rmsvlf000/Trigger/analysis/t_res/v5_optical_diag/diag_1x1_nphi256_100ev.txt`
```
events 100
boundary_total 103909436
scint_to_lg 1178071
lg_to_scint 544
scint_to_world 39316095
world_to_scint 38675198
lg_to_world 3060252
world_to_lg 2832558
lg_to_sipm 949364
sipm_to_lg 928
other_boundary 17896426
bulk_absorb_scint 79403
bulk_absorb_lg 1359
bulk_absorb_glass 67
bulk_absorb_si 142174
sipm_detect 203314
sipm_qe_reject 744967
path_scint_mm_per_event 2.21265e+06
path_lg_mm_per_event 391899
path_glass_mm_per_event 3094.74
path_si_mm_per_event 71.0366
path_world_mm_per_event 131631
```


**Compare 1x1 vs best 4x4** (`diag_1x1_nphi256_100ev.txt` vs `diag_4x4_nphi256_100ev.txt`):
```
baseline: diag_1x1_nphi256_100ev.txt
variant:  diag_4x4_nphi256_100ev.txt
key                                baseline        variant     delta%
----------------------------------------------------------------------
scint_to_world                    3.932e+07      5.261e+07      33.82
world_to_scint                    3.868e+07      5.194e+07      34.30
scint_to_lg                       1.178e+06      1.845e+06      56.63
lg_to_world                        3.06e+06      8.346e+06     172.72
lg_to_sipm                        9.494e+05      5.835e+05     -38.53
sipm_detect                       2.033e+05      1.129e+05     -44.46
path_scint_mm_per_event           2.213e+06      4.481e+06     102.51
path_lg_mm_per_event              3.919e+05      7.935e+05     102.47
```


## ⑤ foil decision
- scint_to_world/scint_to_lg = 28.5 on `diag_4x4_nphi256_100ev.txt`

- **ACTION**: scint↔world still dominant → foil hollow shell deferred to manual PR (geometry subtraction non-trivial overnight). Document as open item.


## pipeline finished 2026-07-06 00:29:09 KST
- summary: `v5_compare_summary.csv`
- full log: `pipeline.log`

---

## Path-length / fate scan (v5.6, 2026-07-06)

**목적**: trigger size 1×1…4×4에서 optical photon **total track length** 분포와 **fate** 비교 (Condor 전).

**코드**: `CBDsimOpticalDiagnostics` — `CBDsim_OPTICAL_DIAG_HIST=1`, `CBDsim_OPTICAL_DIAG_HIST_OUT`

**파이프라인**: `path_scan/run_path_scan.sh` (7 sizes × 100 ev, ~3h)

**결과**: `path_scan/PATH_SCAN_RESULTS.md`, `path_scan/path_vs_size_summary.csv`, `path_scan/figures/`

| size | mean path (mm) | detected % | kill_world % | kill_other % |
|------|----------------|------------|--------------|--------------|
| 1×1 | 130 | 10.4 | 34.7 | 12.2 |
| 2×2 | 158 | 8.8 | 33.9 | 20.2 |
| 4×4 | 253 | 5.5 | 35.3 | 31.6 |

**해석**: trigger ↑ → mean path **~2×**, detected fraction **절반**, `kill_other` **2.6×**. `kill_world` ~35% 유지.
