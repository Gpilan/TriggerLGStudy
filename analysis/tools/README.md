# 분석 도구 진입점

[공유 scan 결과](../results/scan_checkpoint_20260930/README.md) · 로컬 대시보드: `RESULTS.html`, 주제별 목록: `analysis/studies/`

| 목적 | 도구 | 설명 |
|---|---|---|
| 결과 포털 갱신 | `build_results.py --catalog <JSON>` | 명시적으로 선택한 최신 PNG를 주제별 상대 링크로 연결. 분석·시뮬 실행 없음 |
| timing 추출 | [verified_timing.py](../verified_timing.py) | 기존 CFD/SPE 분석 모듈 |
| SPE/CFD 및 Δt 분석 | [plot_cfd_per_trigger_locallinear_peak.py](../plot_cfd_per_trigger_locallinear_peak.py) | 기존 분석 함수; 연구 결과 재현에는 해당 frozen 버전 사용 |
| 광자 종료 성분 | [plot_optical_path_scan.py](../plot_optical_path_scan.py) | fate 분류 및 그림 |
| 폭scan 재현 | 로컬 `analysis/studies/tile_width/tools/` | 기존 source-relative 경로를 유지한 검증 스크립트 |
| 위치scan 재현 | 로컬 `analysis/studies/cross_position/tools/` | 최신 재피팅·트렌드·최종 감사 |

기존 import와 frozen manifest가 참조하는 파일은 이동하지 않았습니다. 이 폴더는 공용 도구의 진입점이며, 날짜별 스크립트를 복사해서 새 공용 구현으로 간주하지 않습니다. 기존 추출 도구는 새 출력 폴더를 전제로 하므로 결과 폴더에서 무작정 전체 재실행하지 마세요.

## 갱신 규칙

1. 새 분석은 독립된 기록 폴더에서 수행하고 ROOT·표본 수·fit·광자 수지 검증을 완료합니다.
2. 로컬 catalog에 해당 결과의 source, 날짜, 이벤트 수, 모델 차이, PNG 선택을 기록합니다.
3. `python3 analysis/tools/build_results.py --catalog <local-catalog.json>`로 최신 진입점을 갱신합니다.
4. `RESULTS.html`과 주제별 링크를 확인합니다. PNG는 링크이므로 그 경로로 덮어쓰면 원본도 변경됩니다.

주제별 목록과 루트 대시보드는 개인 로컬 기록으로 연결되므로 Git 제외 대상입니다. 공유·배포는 별도로 정식 산출물을 선택해 진행합니다. 기본 출력은 PNG만 사용합니다.
