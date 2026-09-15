# SPE waveform: two-sided-gamma-v1

`plot_cfd_per_trigger_locallinear_peak.py`와 `plot_event_optical_mpv_matrix_root.py`는
`spe_response.py`의 동일한 SPE를 사용한다. 단위는 ns다.

2026-09-13 이전에는 `--response-fwhm-ns`가 kernel 길이만 바꾸었다.
rise=1.3, shape=1.915604733026일 때 FWHM=1/3/6 입력 모두 실제 폭 약3.900 ns (교차점 선형 보간 측정)였다.
따라서 과거 **“FWHM 변화가 timing에 영향 없음” 결론은 철회**한다.
기존 scan/파형은 옛 응답을 이용한 탐색 자료로만 남긴다. 새로운 응답으로 재분석해야 한다.

## 모델과 지원 범위

기본 gamma 함수는 `g(x) = exp(s*log(x/s) + s - x)` (x>0), 그 외0이다.
최대값은 x=s에서1. rising 10/50/90% root를 x10/x50/x90,
falling half root를 x50fall이라 쓰면:

- 상승 tau = 요청 rise / (x90-x10)
- peak 시각 = s × 상승 tau
- 최소 FWHM = (s-x50) × 상승 tau
- 하강 tau = (요청 FWHM - 최소 FWHM) / (x50fall-s)
- peak까지 x=t/상승 tau, 이후 x=s+(t-peak)/하강 tau

peak에서 값과 1차 미분이 연속이다. 2차 미분은 달라질 수 있다.
실측 R2076 fit이 아닌 **이상화 파형 모델 선택**이다. 기존 gamma의 상승부를 보존하면서
하강부를 독립 조정한다. 하강 tau=상승 tau일 때 기존 gamma 전체를 복원한다.

`--gamma-shape` 지원 범위는 0.1–100이며, FWHM은 위 최소값보다 커야 한다.
기본 shape/rise1.3의 최소 FWHM은1.41230813 ns. 따라서 1 ns는 거부한다.
예시 scan은 rise1.3 고정, FWHM2/3/6 ns. 경계에 가까운 짧은 하강부는 더 촘촘한
sampling 검증이 필요하다. 검증표의 수렴은 모든 가능한 극단 입력의 수렴 보증이 아니다.

`--gamma-tau-ns`는 상승 tau 호환 확인용이다. rise/shape와 불일치하면 오류이며
보통 생략한다. 예전처럼 rise를 무시하고 tau로 덮어쓰지 않는다.
기존 CLI를 재사용해도 기본 응답은 달라졌으므로 옛 분석과 구분한다.

## 커널과 정규화

`--kernel-tail-level` (기본1e-10)은 하강부 진폭/peak가 이 값에 도달하는 시점까지
커널을 유지한다. 지원0<level<=1e-4. FWHM과 별도 옵션이며 절단 시점은 해석적으로
찾고 샘플 격자까지 올림한다. 실행 로그에 모델명, 양쪽 tau, 절단 시각/기준을 출력한다.
SPE peak는1로 고정; 폭 변화에 따라 적분 면적은 바뀐다. 고정 전하/gain 모델이 아니다.
노이즈 없는 CFD는 공통 진폭 배율에 이상적으로 불변하지만 추출기 검증은 W02에서 다룬다.

기존 photon time bin 중심의 선형 분배, 선형 SPE 합성, 공통 transit offset은 유지한다.
SPE onset에 transit offset을 더하는 convention이며 실물 transit 정의의 검증은 별도다.
photon bin / waveform step / 최종 ΔT histogram bin은 서로 다르다.
이 수정만으로 DAQ200ps, TTS/noise 또는 절대 시간 분해능을 검증한 것은 아니다.

## 검증

```bash
python3 analysis/tests/test_spe_response.py
```

ROOT 없이 실제 두 frontend 함수를 실행한다. 단일 SPE의 요청/실측 rise·FWHM,
50/10/5ps 수렴, tail1e-10→1e-12, 다중 광자 합성 및 CFD, 불가능/충돌 입력,
기존 gamma 극한을 확인한다. 실측값은 `WIDTH_TABLE` JSON으로 출력된다.
수정 이전 ROOT의 광자 도착 bin은 출처가 검증되었다면 재분석 입력으로 쓸 수 있지만,
과거 geometry/QE 검증 문제는 이 수정으로 해결되지 않는다.

## 수치 합성과 CFD 품질 해석

기본 파형 격자는 10 ps다. 희소 입력은 이동·합산하고, 그 외 입력은 선형 합성 길이 이상으로
zero-padding한 real FFT를 사용한다. 음의 부동소수점 반올림 잔차만 0으로 제한한다.
실제 광자 도착 bin의 시간 정보보다 촘촘한 파형 격자가 도착 시각 정보를 새로 만들지는 않는다.

CFD 국소 직선 fit의 R²와 weighted residual은 적합 상태를 설명하는 수치다.
보정된 잡음 공분산 모델이 없으므로 해당 잔차를 통계적 χ²나 p-value로 해석하지 않는다.
실패 사건과 품질을 기록하고, 비교에서는 같은 fraction과 window를 사용한다.

`verified_timing.py`는 검증된 canonical-event JSON의 전체 완료 사건을 분모로 사용하며
T1−T2 중앙 68% 반폭과 같은 표본의 분산·공분산을 보고한다. √2 환산은 하지 않는다.
`photon_thinning.py`의 보조 비교는 원본 사건 단위 cluster 재표집을 사용한다.
