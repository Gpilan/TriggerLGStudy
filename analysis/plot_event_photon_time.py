#!/usr/bin/env python3
"""
CBDsim .root 에서 **한 이벤트**의 트리거별 **광학 포톤 도착 시간** 분포를
**ROOT 기본 스타일(HIST)** 로 저장합니다.

원본은 SiPM SD 의 아주 잘게 잡힌 빈이며, 여기서는 **고정 폭(기본 20 ps)** 으로 리빈한 TH1F 를 그립니다.
**x 축 범위는 포톤이 있는 시간 구간만**(여유 마진 포함) 자동으로 잡습니다.

예:
  source envset.sh
  export LD_LIBRARY_PATH=$PWD/build/rootIO:$LD_LIBRARY_PATH
  python3 analysis/plot_event_photon_time.py analysis/t_res/data/60GeV_e-_noLG_0.root --event 0
"""

from __future__ import annotations

import argparse
import os
import sys

from photon_time_root import apply_photon_time_style, bin_width_ns_from_ps, th1_rebinned_from_merged
from plot_trigger_timing import (
    _figures_dir,
    _find_repo_root,
    _load_rootio,
    _merged_from_tower,
    _prepend_build_rootio_ld_path,
    _resolve_input_path,
    _input_root_dirs,
)


def _default_input(repo_root: str) -> str | None:
    for sd in _input_root_dirs(repo_root):
        for name in ("60GeV_e-_noLG_0.root", "noLG1_0.root", "An_test_1.root"):
            p = os.path.join(sd, name)
            if os.path.isfile(p):
                return p
    return None


def main() -> None:
    parser = argparse.ArgumentParser(
        description="한 이벤트 포톤 시간 분포 (ROOT HIST, 리빈)"
    )
    parser.add_argument("input", nargs="?", default=None, help=".root 경로")
    parser.add_argument("--event", "-e", type=int, default=0, help="트리 엔트리 번호")
    parser.add_argument("--trig", choices=("0", "1", "both"), default="both")
    parser.add_argument("-o", "--output", default=None)
    parser.add_argument(
        "--bin-width-ps",
        type=float,
        default=20.0,
        help="시간 축 빈 폭 (기본 20 ps)",
    )
    parser.add_argument("-l", "--rootio-lib", default=None)
    args = parser.parse_args()

    repo_root = _find_repo_root()
    if repo_root is None:
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    _prepend_build_rootio_ld_path(repo_root)

    if args.input:
        input_path = _resolve_input_path(args.input, repo_root)
        if not input_path:
            sys.stderr.write(f"입력 파일 없음: {args.input}\n")
            sys.exit(1)
    else:
        input_path = _default_input(repo_root)
        if not input_path:
            sys.stderr.write("입력 .root 를 찾지 못했습니다.\n")
            sys.exit(1)

    bw_ns = bin_width_ns_from_ps(args.bin_width_ps)

    _load_rootio(args.rootio_lib, repo_root)
    import ROOT

    ROOT.gROOT.SetBatch(True)
    apply_photon_time_style(ROOT)

    f = ROOT.TFile.Open(input_path)
    if not f or f.IsZombie():
        sys.stderr.write(f"파일 열기 실패: {input_path}\n")
        sys.exit(1)
    tree = f.Get("CBDsim")
    if not tree:
        sys.stderr.write("트리 CBDsim 없음\n")
        sys.exit(1)

    nent = int(tree.GetEntries())
    if args.event < 0 or args.event >= nent:
        sys.stderr.write(f"event 범위 오류: {args.event} (max {nent-1})\n")
        sys.exit(1)

    evt = ROOT.CBDsimInterface.CBDsimEventData()
    tree.SetBranchAddress("CBDsimEventData", evt)
    tree.GetEntry(args.event)

    def series(trig: int):
        lo, hi, cnt = _merged_from_tower(evt, trig)
        if not cnt:
            return [], [], 0
        centers = [0.5 * (float(lo[i]) + float(hi[i])) for i in range(len(cnt))]
        counts = [int(cnt[i]) for i in range(len(cnt))]
        return centers, counts, sum(counts)

    c0, y0, n0 = series(0)
    c1, y1, n1 = series(1)
    f.Close()

    bn = os.path.splitext(os.path.basename(input_path))[0]
    out = args.output or os.path.join(
        _figures_dir(repo_root),
        f"event_{bn}_ev{args.event}_photon_time.png",
    )
    out = os.path.abspath(out)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    c = ROOT.TCanvas("c_evt_photon", "photon time", 900, 820)
    if args.trig == "both":
        c.Divide(1, 2)

    def draw_pad(pad_num: int, trig: int, centers: list, counts: list, ntot: int, color: int) -> None:
        if args.trig == "both":
            c.cd(pad_num)
        else:
            c.cd(0)
        sub = f"T{trig+1}  #Sigma={ntot}"
        h = th1_rebinned_from_merged(
            ROOT,
            f"hEv{args.event}_tr{trig}",
            f"{sub};arrival time (ns);photons / {args.bin_width_ps:g} ps bin",
            centers,
            counts,
            bin_width_ns=bw_ns,
        )
        if h is None:
            return
        h.SetLineColor(color)
        h.Draw("HIST")

    if args.trig == "both":
        draw_pad(1, 0, c0, y0, n0, ROOT.kBlue + 1)
        draw_pad(2, 1, c1, y1, n1, ROOT.kOrange + 7)
    elif args.trig == "0":
        draw_pad(0, 0, c0, y0, n0, ROOT.kBlue + 1)
    else:
        draw_pad(0, 1, c1, y1, n1, ROOT.kOrange + 7)

    c.cd(0)
    c.SetTitle(f"{bn}  event {args.event}  |  rebin {args.bin_width_ps:g} ps")

    if not out.lower().endswith((".png", ".pdf", ".svg")):
        out += ".png"
    c.SaveAs(out)
    c.Close()
    print(f"저장: {out}")
    print(f"  T1 photons={n0}, T2 photons={n1}  (리빈 {args.bin_width_ps:g} ps)")


if __name__ == "__main__":
    main()
