Condor로 CBDsim 여러 시드 돌리기
===============================

파일
----
  run_cbdsim.sh   : job 인자 = 시드(보통 $(Process)). ROOT 는 build/CBDsim/stats/<OUT_SUBDIR>/job_<id>_<id>.root
  cbdsim_array.sub: submit 예제 (executable 경로 CHANGEME 수정)

한 번에 제출하기 전에
---------------------
  - build/CBDsim/CBDsim 빌드됨
  - 같은 노드에서 ./CBDsim/condor/run_cbdsim.sh 0 수동 실행으로 스모크 테스트
  - logs/ 디렉터리 생성: mkdir -p CBDsim/condor/logs (submit 위치에 맞춤)

출력 ROOT 경로
--------------
  RunAction 이 prefix_seed.root 로 쓰므로, 시드=5 이면
  stats/condor/job_5_5.root 형태 (OUT_SUBDIR=condor 기본)

환경변수
--------
  MACRO=run_muon.mac
  OUT_SUBDIR=mybatch
  (run_cbdsim.sh 안에서 읽음)

경북대 T3
---------
  https://t2-cms.knu.ac.kr/index.php/HTCondor
