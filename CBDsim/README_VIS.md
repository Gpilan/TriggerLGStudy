# CBDsim 시각화 / GUI (SSH X11에서 실행 시)

## 빌드 구조와 cmake 위치

이 프로젝트는 **한 단계 위**에서 cmake를 돌리는 구조입니다.

```
Trigger/                    ← 프로젝트 루트 (최상위 CMakeLists.txt)
├── CMakeLists.txt          ← rootIO, CBDsim, Reco, analysis 서브디렉터리 포함
├── build/                  ← 여기서 cmake 실행 (소스는 상위 디렉터리)
│   ├── CBDsim/             ← 실행 파일과 매크로(vis.mac 등)가 있는 곳
│   │   ├── CBDsim          ← 실행 파일
│   │   ├── vis.mac, init_vis.mac, gui.mac, run_ele.mac ...
│   │   └── ...
│   └── ...
└── CBDsim/                 ← 소스 (vis.mac 등 매크로 원본)
    ├── vis.mac, init_vis.mac, ...
    └── ...
```

**빌드 절차 (cms02 등 서버에서):**

```bash
cd ~/Trigger              # 프로젝트 루트로
source envset.sh           # LCG view 로드 (Geant4 등)

cd build                   # build 디렉터리로 (없으면 mkdir build 후 cd build)
cmake ..                   # 상위(Trigger/)를 소스로 cmake 실행 → 매크로가 build/CBDsim/로 복사됨
make

cd CBDsim
./CBDsim                   # GUI 모드 실행
```

- **cmake는 `Trigger/build/` 안에서 `cmake ..` 로 실행**합니다. (소스 = `Trigger/`)
- `vis.mac` 같은 매크로는 **cmake 할 때마다** `CBDsim/` → `build/CBDsim/` 으로 복사됩니다.  
  소스에서 `vis.mac`을 바꿨다면 **다시 `build`에서 `cmake ..` 한 번 돌린 뒤** `make` 하면 반영됩니다.

---

## GUI가 안 뜨는 경우

원격 서버(cms02 등)에 SSH로 접속한 뒤 `./CBDsim`만 실행하면 GUI가 안 뜨는 경우가 있습니다.

### 1. X11 포워딩으로 접속

- **반드시 X11 포워딩 옵션으로 접속**해야 합니다.
  - `ssh -X knu-x11` 또는 **`ssh -Y knu-x11`** (권장: Trusted X11)
- 접속 후 터미널에서 확인:
  ```bash
  echo $DISPLAY
  ```
  - `localhost:10.0` 같은 값이 나와야 합니다. 비어 있으면 GUI 불가.

### 2. 시각화 드라이버 (OGLX)

- `vis.mac`에서는 **OGLX** (X11용 OpenGL)를 사용하도록 설정해 두었습니다.
- SSH X11에서는 Qt OpenGL(OGL)이 잘 안 되는 경우가 많고, OGLX가 더 잘 동작합니다.
- 소스에서 `vis.mac`을 수정한 뒤에는 **다시 cmake/빌드**를 하거나, `build/CBDsim/vis.mac`을 직접 수정해 두어야 실행 시 반영됩니다.

### 3. 로컬 PC 측 (X11 클라이언트)

- **Linux**: 보통 별도 설정 없이 `ssh -Y` 만으로 동작합니다.
- **macOS**: XQuartz 설치 후, 터미널에서:
  ```bash
  defaults write org.macosforge.xquartz.X11 enable_iglx -bool true
  ```
  적용 후 XQuartz를 한 번 재시작한 다음 `ssh -Y`로 접속합니다.
- **Windows**: VcXsrv, Xming, MobaXterm 등 X 서버 실행 후 `DISPLAY`를 설정해 SSH 접속합니다.

### 4. GUI 없이 배치로만 실행

- 매크로 파일을 넘기면 GUI 없이 실행되고 바로 끝납니다.
  ```bash
  ./CBDsim run_ele.mac
  ```

정리: **`ssh -Y knu-x11`** 로 접속하고, `DISPLAY`가 설정된 상태에서 `./CBDsim`을 실행해 보세요. 여전히 창이 안 뜨면 로컬에서 indirect GLX(위 macOS 설정) 또는 방화벽/보안 정책을 확인하면 됩니다.
