# 정밀 광자 진단 (CBDsimPrecision)

일반 `CBDsim`에는 정밀 기록 코드가 컴파일되지 않는다. `CBDsimPrecision`은 같은
main, geometry, physics, actions를 사용하며 기존 광자 수지 코드에 읽기 전용 관찰자를
추가한다. 새 광자 생성, RNG 호출, 강제 kill, geometry 변경은 없다.

```bash
source envset.sh
cmake -S . -B build/precision -DCMAKE_BUILD_TYPE=Release -DCBDsim_BUILD_PRECISION=ON
cmake --build build/precision --target CBDsimPrecision -j4
python3 CBDsim/precision/run.py --binary build/precision/CBDsim/precision/CBDsimPrecision \
  --out /tmp/trigger-nolg-diagnostic-001 --mode nolg
```

출력 디렉토리는 매번 **새 이름**을 쓴다. 권장 실행 인터페이스는 위 runner다.
GUI나 여러 `beamOn`을 한 프로세스에서 실행하는 사용법은 지원하지 않는다.
단일 worker, 기본 e+ 60 GeV, (0,0,-100) mm/+z, seed42, 1 event다.
기존 사용자 사진의 z=-1000 mm와 다르다. 현재 world 안에서 시작하도록 정한
진단 기본값이며, 사진 조건 재현은 `--position-mm 0 0 -1000`으로 명시한다.
현재 world 밖 시작은 제어된 geometry 검증의 기준으로 사용하지 않는다.

## 기록 단계

- `--level off`: 정밀 관찰자 비활성화. 기존 aggregate 진단만 실행한다.
- `--level summary` (기본): 시작한 광자 최대 100,000개의 생성 및 최종 정보,
  world로 실제 투과한 횟수·첫/마지막 지점과 world에서 돌아온 횟수를 저장한다.
- `--level steps --event 0 --track 123`: 해당 event/track의 매 step을 추가 저장한다.
  track ID는 event 내에서만 유일하다. 먼저 summary에서 관심 ID를 찾아 같은
  binary, geometry, beam, seed로 재실행한다. `--track`을 생략하면 선택된 모든
  광자의 step을 저장하므로 파일이 커질 수 있다.

`--max-tracks 100000`, `--max-steps 100000`은 한 실행의 기록 개수 상한이다.
상한은 **기록만 제한하며 물리 추적을 중단하지 않는다.** track 상한은 먼저 추적된
광자부터 선택하므로 통계적 무작위 표본이 아니다. step 상한 이후의 경로는 저장되지
않으며 광자 요약은 계속 갱신한다. `--max-steps 0`은 step 기록 비활성화다.
필터에 제외된 광자도 정상 물리 추적과 기존 수지 계수에 포함된다.
정밀 관찰자의 메모리는 track 상한에 제한되지만 기존 Geant4/광자 수지 메모리는
별도이며, 전체 실행 비용까지 고정되는 것은 아니다.

## 산출물과 해석

- `photons-w0.tsv`: event/track/parent/creator, 생성 위치·시간·에너지·volume 경로,
  기존 수지의 정확한 종료 원인, 최종 위치·process·material·경계 상태,
  경로 길이, step 수, world 왕복 및 최종 touchable 경로·시간·방향.
- `steps-w0.tsv`: 선택된 광자의 step 전후 위치·방향·시간·에너지,
  volume/copy 계층, material, process, step/track/boundary 상태.
- `events-w0.tsv`: 이벤트 종료 및 광자/step 기록 한도에 따른 누락 개수.
- `diagnostics.txt`와 sidecars: 전체 광자 수지. 정밀 기록 cap과 무관하다.
- `manifest.json`: 실행 완료 여부, 조건, binary/source hash, 실행 시간·peak RSS,
  기록량. 실패/timeout/중단은 `failed_or_incomplete`로 남긴다.
- `--root`: ROOT 광자 검출수·도착시간 bin의 canonical digest 검증도 실행한다.

좌표 mm, 시간 ns, 에너지 eV. `first_to_world`/`last_to_world`는
`step|출발 touchable|world x,y,z|출발 volume local x,y,z` 형식이다.
`/volume[copy]` 계층으로 T1/T2와 nested wafer를 구분한다.
Boundary status는 해당 Geant4 enum 정수이며 현재 step이 `fGeomBoundary`가
아니면 -1이다. 반사는 투과 횟수에 포함하지 않는다.

**world 진입은 누설 확정이 아니다.** foil 공기 틈도 world이므로 돌아올 수 있다.
`fate=world_exit`가 world 경계 밖 최종 종료다. 마지막 world 진입 위치는 원인 후보를
찾는 단서이며 조립체 외부로 나간 최초 지점의 자동 판정은 아직 제공하지 않는다.
`unfinished`의 최종 필드에는 종료 의미가 없다. stack에서 추적 전 제거된 광자는
기존 수지와 마찬가지로 이 관찰자의 분모 밖이다.

`complete`는 실행 및 기록 수지 완주를 뜻하며 optical model의 타당성 판정이 아니다.
상한 도달은 manifest에서 별도로 확인한다. 정확한 누설 분포는 `summary`에서
track 누락이 없는 실행을 사용한다. 진단 TSV는 event 종료마다 flush하며, 비정상
종료 시 현재 event의 summary가 없을 수 있으므로 partial step 파일만으로 분모를 만들지 않는다.

## 디렉토리와 공유 범위

```text
CBDsim/
├── include/                  # 공용 헤더
├── src/                      # 공용 geometry·물성·actions·기존 수지
├── CBDsim.cc                 # 공용 main
├── CMakeLists.txt            # 일반 실행 및 정밀 target 선택
└── precision/
    ├── CMakeLists.txt        # 정밀 실행 target
    ├── src/PrecisionTrace.inc # 정밀 관찰자
    ├── run.py               # 정밀 실행·기록 검증
    └── README.md
```

공용 `.cc`/헤더 원본은 한 벌이다. 두 target이 공용 소스를 각각 컴파일하며,
정밀 target에만 `CBDsim_PRECISION_TRACE`를 정의한다. 공용 진단 소스의 작은
조건부 hook이 정밀 관찰자를 연결한다. 관찰자는 기존 내부 종료 원인 판정을
재사용하기 위해 `.inc`로 포함하며 별도의 geometry 복사본을 만들지 않는다.

CMake build root가 `<build>`이면 일반 실행 파일은 `<build>/CBDsim/CBDsim`,
정밀 실행 파일은 `<build>/CBDsim/precision/CBDsimPrecision`이다.
공용 소스를 수정하면 필요한 두 target을 다시 빌드한다.

## 흡수성 테이프 모델 (2026-09-16 확장)

`--gel-tape absorber --tape-extensions on`이 기본값이다. **noLG gel+window 네 옆면**,
**LG inlet gel 네 옆면**을 실제 sleeve로 덮는다. LG 출구 gel/window/ring은 바꾸지 않는다.
직접 접촉, 두께0.05 mm, PVC 대용 bulk, R=T=EFFICIENCY=0인 완전흡수 표면은
진단용 가정이다. 실물 물성 측정값이 아니다. gel 앞뒤 광학 결합면은 열려 있다.
기존 foil과 겹치는 부피는 foil에서 제거했다. noLG sleeve는 y=-30..-30.39 mm,
LG inlet sleeve는 y=-30..-30.02 mm다. 검출면과 앞뒤 광학 결합면은 덮지 않는다.
`--tape-extensions off` (`CBDsim_PROTO_TAPE_EXTENSIONS=0`)는 기존 noLG gel-only 조건이다.
이전 조건은 `--gel-tape off` (`CBDsim_PROTO_GEL_TAPE=0`)로 실행한다.
이 옵션은 공용 geometry에 적용되므로 일반 CBDsim에서도 환경변수로 선택한다.

전 광자 종료 원인/경계별 fraction:

```bash
python3 CBDsim/precision/summarize_terminations.py /path/to/completed-run --out /tmp/fates-new.json
```

분모는 실제 tracking에 진입한 광자다. 누락되거나 cap에 걸린 run은 거부한다.
`absorb_other_surface` 중 tape/foil 경계를 분리하고 pre/post volume과 material별로 집계한다.
검출·QE 탈락의 직접 SD kill과 Geant4의 물리적 흡수/세계 밖 종료를 구별해서 해석한다.
`--polarization x y z`는 optical probe의 편광을 지정한다(진행 방향에 수직으로 설정).

## Tile–foil 모서리 접합부 (2026-09-15)

`--corner-seal on`이 기본값이다 (`CBDsim_PROTO_CORNER_SEAL=1`). tile 주변 공기 틈의
열린 끝을 실제 Al rim으로 닫으며 기존 `AluminumSurf`의 파장별 반사율을 재사용한다.
Al rim은 LG tile 끝에서 멈춘다. inlet gel 옆면의 테이프는 별도 옵션이다. noLG는 기존 foil 끝까지 이어지고
이미 있는 gel tape/foil과 겹치는 부피를 제외한다. 이상적인 접촉 치수는 진단 모델이다.
2026-09-15 gel-tape 조건은 `--tape-extensions off --corner-seal off`,
그 이전 조건은 여기에 `--gel-tape off`도 지정한다.
분류기는 새 rim 흡수를 `corner_seal_absorption`으로 분리한다. 기존 `unknown=0`과
물리적인 `world_exit` 감소는 다른 검증 항목이다.

## LG 출구 sleeve 밀착 (2026-09-16)

`--tip-contact on`이 기본이며 일반 CBDsim은 `CBDsim_PROTO_TIP_CONTACT=1`에 해당한다.
LG sensor ring 내경 반지름만7.510→7.500mm로 조정하고 외경7.576mm 및 축 범위는 유지한다.
실물 밀착 설명을 반영한 optical contact 모델이며 scintillator foil 공기층은 유지한다.
`--tip-contact off` / `CBDsim_PROTO_TIP_CONTACT=0`으로 직전10um gap 모델을 재현한다.
과거 tape/corner 조건 재현에는 이 옵션도off로 지정한다. 반사율과 kill 로직은 변경하지 않는다.

## 원형 LG 끝단 시험 (2026-09-16)

`--round-tip on` / `CBDsim_PROTO_ROUND_TIP=1`이 기본 모델이다(2026-09-16 사용자 GUI 확인 후 승격).
LG 전체30mm와gel위치를유지하며길이0.2mm/반지름7.5mm 원통을taper와0.1mm겹쳐하나의LG solid로합친다.
Taper는29.9mm길이/끝반지름7.499mm이며끝단외형이변하므로측정된형상으로해석하지않는다.
`--round-tip off`/환경변수0은직전다각형LG를정확히재현한다.

정밀 photon TSV에 `shoulder_reflections`, `shoulder_absorptions`를 추가했다. 현재 원통 시작면(local y=-59.82mm)의 LG 내부 입사만 집계하며 0-length StepTooSmall은 제외한다. 기록 전용이고 일반 CBDsim에는 포함되지 않는다. 치수가 바뀌면 해당 진단 기준도 갱신해야 한다.
