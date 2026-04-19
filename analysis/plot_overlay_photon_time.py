#!/usr/bin/env python3
"""
읽는 모든 이벤트에 대해, **빈(시간)마다 포톤 수를 이벤트 평균**한 TH1F 하나를
트리거(T1/T2)별 패드에 ROOT `HIST` 로 그립니다.

- 오버레이 곡선 없음 — **평균 분포만** 표시.
- x 범위: **데이터가 있는 시간대만**(여유 마진). `--trig both` 이면 T1·T2 **공통 축**
  (한쪽에만 신호가 있어도 양쪽 패드에 같은 t 범위로 그려 T2가 비지 않게 함).
- 빈 폭: 기본 **20 ps** (`--bin-width-ps`).

예:
  python3 analysis/plot_overlay_photon_time.py analysis/t_res/data/60GeV_e-_noLG_0.root
"""

from __future__ import annotations

import argparse
import os
import sys

from photon_time_root import (
    apply_photon_time_style,
    bin_width_ns_from_ps,
    data_time_range_ns,
    th1_empty_same_bins,
    th1_rebinned_from_merged,
)
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


def _series_for_event(evt, trig: int):
    lo, hi, cnt = _merged_from_tower(evt, trig)
    if not cnt:
        return None, None
    centers = [0.5 * (float(lo[i]) + float(hi[i])) for i in range(len(cnt))]
    counts = [int(cnt[i]) for i in range(len(cnt))]
    return centers, counts


def _union_range_one_trigger(
    tree, evt, trig: int, n_read: int, stride: int, bw_ns: float
) -> tuple[float, float] | tuple[None, None]:
    g_lo = float("inf")
    g_hi = float("-inf")
    for i in range(0, n_read, stride):
        tree.GetEntry(i)
        centers, counts = _series_for_event(evt, trig)
        if not centers:
            continue
        for c, n in zip(centers, counts):
            if n > 0:
                cc = float(c)
                g_lo = min(g_lo, cc)
                g_hi = max(g_hi, cc)
    if g_lo == float("inf"):
        return None, None
    return data_time_range_ns([g_lo, g_hi], [1, 1], bw_ns)


def _union_range_both_triggers(
    tree, evt, n_read: int, stride: int, bw_ns: float
) -> tuple[float, float] | tuple[None, None]:
    """T1+T2 모두에 대해 n>0 인 빈 중심의 전역 min/max → T1·T2 패드 공통 축."""
    g_lo = float("inf")
    g_hi = float("-inf")
    for trig in (0, 1):
        for i in range(0, n_read, stride):
            tree.GetEntry(i)
            centers, counts = _series_for_event(evt, trig)
            if not centers:
                continue
            for c, n in zip(centers, counts):
                if n > 0:
                    cc = float(c)
                    g_lo = min(g_lo, cc)
                    g_hi = max(g_hi, cc)
    if g_lo == float("inf"):
        return None, None
    return data_time_range_ns([g_lo, g_hi], [1, 1], bw_ns)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="이벤트 평균 광학 포톤 도착 시간 (ROOT HIST)"
    )
    parser.add_argument("input", nargs="?", default=None)
    parser.add_argument("--trig", choices=("0", "1", "both"), default="both")
    parser.add_argument("--max-events", type=int, default=-1)
    parser.add_argument("--stride", type=int, default=1)
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
    n_read = nent if args.max_events < 0 else min(nent, args.max_events)
    stride = max(1, args.stride)
    n_loop = len(range(0, n_read, stride))

    evt = ROOT.CBDsimInterface.CBDsimEventData()
    tree.SetBranchAddress("CBDsimEventData", evt)

    triggers: list[tuple[int, str, int]] = []
    if args.trig == "both":
        triggers = [(0, "T1 (trigger 0): mean photons / bin", ROOT.kBlue + 1), (1, "T2 (trigger 1): mean photons / bin", ROOT.kOrange + 7)]
    elif args.trig == "0":
        triggers = [(0, "T1 (trigger 0): mean photons / bin", ROOT.kBlue + 1)]
    else:
        triggers = [(1, "T2 (trigger 1): mean photons / bin", ROOT.kOrange + 7)]

    # 공통 축: both 이면 T1+T2 모두에서 신호가 있는 구간; 아니면 해당 트리거만
    if args.trig == "both":
        t_lo, t_hi = _union_range_both_triggers(tree, evt, n_read, stride, bw_ns)
    else:
        t_lo, t_hi = _union_range_one_trigger(
            tree, evt, int(args.trig), n_read, stride, bw_ns
        )

    if t_lo is None or t_hi is None:
        sys.stderr.write("시간 범위를 정할 수 없습니다 (모든 이벤트에 해당 트리거 포톤이 없음?).\n")
        sys.exit(1)

    bn = os.path.splitext(os.path.basename(input_path))[0]
    out = args.output or os.path.join(
        _figures_dir(repo_root), f"mean_photon_time_{bn}.png"
    )
    out = os.path.abspath(out)
    if not out.lower().endswith((".png", ".pdf", ".svg")):
        out += ".png"
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    c = ROOT.TCanvas("c_mean_photon", "mean photon time", 1000, int(420 * len(triggers)))
    if len(triggers) > 1:
        c.Divide(1, len(triggers))

    for ip, (trig, subtitle, col) in enumerate(triggers, start=1):
        if len(triggers) > 1:
            c.cd(ip)
        else:
            c.cd(0)

        h_sum = None
        for i in range(0, n_read, stride):
            tree.GetEntry(i)
            centers, counts = _series_for_event(evt, trig)
            if centers is None or not counts:
                h = th1_empty_same_bins(
                    ROOT,
                    f"hz_{trig}_ev{i}",
                    f"{subtitle};arrival time (ns);photons / {args.bin_width_ps:g} ps bin",
                    t_lo,
                    t_hi,
                    bw_ns,
                )
            else:
                h = th1_rebinned_from_merged(
                    ROOT,
                    f"h_{trig}_ev{i}",
                    f"{subtitle};arrival time (ns);photons / {args.bin_width_ps:g} ps bin",
                    centers,
                    counts,
                    bin_width_ns=bw_ns,
                    t_min_ns=t_lo,
                    t_max_ns=t_hi,
                )
            if h is None:
                continue
            if h_sum is None:
                h_sum = h.Clone(f"hsum_{trig}")
                h_sum.SetDirectory(0)
            else:
                h_sum.Add(h)

        if h_sum is None:
            continue

        h_mean = h_sum.Clone(f"hmean_{trig}")
        h_mean.SetDirectory(0)
        h_mean.Scale(1.0 / float(n_loop))
        h_mean.SetTitle(
            f"{subtitle};arrival time (ns);#LT photons #GT / bin / event ({n_loop} ev.)"
        )
        h_mean.SetLineColor(col)
        h_mean.SetLineWidth(2)
        h_mean.Draw("HIST")
        ROOT.gPad.SetGrid(1, 1)
        ROOT.gPad.Modified()
        ROOT.gPad.Update()

    c.cd(0)
    c.SetTitle(
        f"{bn}  |  mean over {n_loop} events  |  rebin {args.bin_width_ps:g} ps  |  stride {stride}"
    )
    c.SaveAs(out)
    c.Close()
    f.Close()

    print(f"저장: {out}")
    print(f"  이벤트 수(평균 분모): {n_loop}, 공통 t 범위 [{t_lo:.6f}, {t_hi:.6f}] ns")


if __name__ == "__main__":
    main()
