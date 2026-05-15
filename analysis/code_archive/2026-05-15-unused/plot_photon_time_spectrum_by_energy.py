#!/usr/bin/env python3
"""
에너지·트리거별 **광자 도착시간 스펙트럼** (merged SiPM 시간 빈에 포톤 수를 가중).

- 이벤트마다 ``_merged_from_tower`` 의 각 빈 중심 시각에 ``counts`` 만큼 더한 뒤,
  기본(``--spectrum-mode mean-per-event``)은 **이벤트 수로 나눈 평균**
  (시간 빈마다 “이벤트당 평균 포톤 수”)으로 그린다.
- ``--spectrum-mode sum`` 이면 예전처럼 전체 합계 포톤 수.

LG / noLG 각각 2×4 패널 (위 T1, 아래 T2, 열=60·80·100·120 GeV).

입력·파일명 기본값은 ``plot_cfd_toa_vs_nphoton_by_energy.py`` 와 동일.

예:
  export LD_LIBRARY_PATH=$PWD/build/rootIO:$LD_LIBRARY_PATH
  python3 analysis/plot_photon_time_spectrum_by_energy.py

기본 시간 축은 **0~20 ns** 고정, 빈 **160** (약 0.125 ns/빈). 예전처럼 데이터에서 범위를
스캔하려면 ``--auto-time-range`` 를 쓴다.
"""

from __future__ import annotations

import argparse
import os
import sys

from plot_cfd_toa_vs_nphoton_by_energy import DEFAULT_LG, DEFAULT_NOLG
from plot_trigger_timing import (
    _figures_dir,
    _find_repo_root,
    _load_rootio,
    _merged_from_tower,
    _prepend_build_rootio_ld_path,
)

import ROOT


def _scan_time_range_merged(
    path: str,
    trig: int,
    *,
    max_events: int,
) -> tuple[float, float]:
    """merged 빈 중 count>0 인 구간의 최소/최대 (ns)."""
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise FileNotFoundError(path)
    tree = f.Get("CBDsim")
    if not tree:
        f.Close()
        raise RuntimeError(f"트리 CBDsim 없음: {path}")
    evt = ROOT.CBDsimInterface.CBDsimEventData()
    tree.SetBranchAddress("CBDsimEventData", evt)
    nmax = int(tree.GetEntries())
    if max_events >= 0:
        nmax = min(nmax, max_events)

    t_lo = float("inf")
    t_hi = float("-inf")
    for i in range(nmax):
        tree.GetEntry(i)
        lo, hi, cnt = _merged_from_tower(evt, trig)
        n = len(cnt)
        for j in range(n):
            if int(cnt[j]) <= 0:
                continue
            t_lo = min(t_lo, float(lo[j]))
            t_hi = max(t_hi, float(hi[j]))
    f.Close()
    if t_lo == float("inf"):
        return 0.0, 1.0
    return t_lo, t_hi


def _accumulate_time_spectrum(
    path: str,
    trig: int,
    *,
    max_events: int,
    nbins: int,
    t_lo: float,
    t_hi: float,
    edge_margin: bool,
):
    """(hist, edges) — 각 빈은 포톤 개수 가중 합.

    ``edge_margin`` 이 True 이면 [t_lo,t_hi]에 span의 3% 마진을 붙여 빈 경계를 잡는다.
    False 이면 정확히 t_lo~t_hi 구간을 nbins 로 나눈다. 범위 밖 중심값은 누적하지 않는다.
    """
    import numpy as np

    span = max(t_hi - t_lo, 1e-9)
    if edge_margin:
        margin = max(span * 0.03, 1e-9)
        edges = np.linspace(t_lo - margin, t_hi + margin, nbins + 1)
    else:
        edges = np.linspace(t_lo, t_hi, nbins + 1)
    h = np.zeros(nbins, dtype=float)

    f = ROOT.TFile.Open(path)
    tree = f.Get("CBDsim")
    evt = ROOT.CBDsimInterface.CBDsimEventData()
    tree.SetBranchAddress("CBDsimEventData", evt)
    nmax = int(tree.GetEntries())
    if max_events >= 0:
        nmax = min(nmax, max_events)

    for i in range(nmax):
        tree.GetEntry(i)
        lo, hi, cnt = _merged_from_tower(evt, trig)
        n = len(cnt)
        for j in range(n):
            w = float(cnt[j])
            if w <= 0.0:
                continue
            cen = 0.5 * (float(lo[j]) + float(hi[j]))
            idx = int(np.searchsorted(edges, cen, side="right")) - 1
            if idx < 0 or idx >= nbins:
                continue
            h[idx] += w
    f.Close()
    return h, edges, nmax


def _plot_series(
    *,
    data_dir: str,
    series: list[tuple[int, str]],
    tag: str,
    out_path: str,
    max_events: int,
    nbins: int,
    spectrum_mode: str,
    time_range: tuple[float, float] | None,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as e:
        print(f"matplotlib/numpy 없음: {e}", file=sys.stderr)
        sys.exit(1)

    color = "#c0392b" if tag == "LG" else "#2980b9"

    fig, axes = plt.subplots(2, 4, figsize=(18, 7.5))

    for col, (egev, fname) in enumerate(series):
        path = os.path.join(data_dir, fname)
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

        for row, trig, tlabel in ((0, 0, "T1"), (1, 1, "T2")):
            ax = axes[row, col]
            if time_range is None:
                t_lo, t_hi = _scan_time_range_merged(
                    path, trig, max_events=max_events
                )
                edge_margin = True
            else:
                t_lo, t_hi = time_range
                edge_margin = False
            h, edges, n_ev = _accumulate_time_spectrum(
                path,
                trig,
                max_events=max_events,
                nbins=nbins,
                t_lo=t_lo,
                t_hi=t_hi,
                edge_margin=edge_margin,
            )
            if spectrum_mode == "mean-per-event":
                h = h / max(float(n_ev), 1.0)
                ylabel = "mean photons / event / bin"
                stat_line = f"N_evt={n_ev}, sum/bin → ÷N_evt"
            else:
                ylabel = "photons (total in bin)"
                stat_line = f"N_evt={n_ev}, sum={float(np.sum(h)):.0f} photons"

            wbin = float(edges[1] - edges[0])
            centers = 0.5 * (edges[:-1] + edges[1:])
            ax.bar(
                centers,
                h,
                width=wbin,
                align="center",
                color=color,
                alpha=0.55,
                edgecolor="black",
                linewidth=0.2,
            )
            ax.set_xlabel("time (ns)")
            ax.set_ylabel(ylabel)
            ax.set_title(f"{egev} GeV  {tlabel}\n{fname}\n{stat_line}")
            ax.grid(True, alpha=0.28, axis="y")
            ax.set_xlim(float(edges[0]), float(edges[-1]))

    mode_title = (
        "mean per event (÷N_evt)"
        if spectrum_mode == "mean-per-event"
        else "total over events"
    )
    fig.suptitle(
        f"{tag}: photon arrival time ({mode_title}, merged SiPM) vs energy",
        fontsize=11,
    )
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"PNG: {out_path}")


def main() -> None:
    p = argparse.ArgumentParser(description="Photon time spectrum by energy (LG / noLG)")
    p.add_argument("--data-dir", default=None)
    p.add_argument("--max-events", type=int, default=-1)
    p.add_argument(
        "--t-min",
        type=float,
        default=0.0,
        help="고정 시간 축 하한 (ns); --auto-time-range 일 때 무시",
    )
    p.add_argument(
        "--t-max",
        type=float,
        default=20.0,
        help="고정 시간 축 상한 (ns); --auto-time-range 일 때 무시",
    )
    p.add_argument(
        "--auto-time-range",
        action="store_true",
        help="각 패널에서 merged 데이터로 min/max 스캔 (기본은 --t-min/--t-max 고정)",
    )
    p.add_argument(
        "--bins",
        type=int,
        default=160,
        help="시간 축 빈 개수 (기본 160 → 0~20 ns 에서 약 0.125 ns/빈)",
    )
    p.add_argument(
        "--spectrum-mode",
        choices=("mean-per-event", "sum"),
        default="mean-per-event",
        help="mean-per-event: 각 빈을 이벤트 수로 나눔(기본). sum: 전체 합.",
    )
    p.add_argument("-l", "--rootio-lib", default=None)
    p.add_argument("--out-lg", default=None)
    p.add_argument("--out-nolg", default=None)
    args = p.parse_args()

    if not args.auto_time_range and args.t_max <= args.t_min:
        print("--t-max 는 --t-min 보다 커야 합니다.", file=sys.stderr)
        sys.exit(1)

    repo_root = _find_repo_root()
    if repo_root is None:
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    _prepend_build_rootio_ld_path(repo_root)
    _load_rootio(args.rootio_lib, repo_root)
    ROOT.gROOT.SetBatch(True)

    data_dir = args.data_dir
    if data_dir is None:
        data_dir = os.path.join(repo_root, "analysis", "t_res", "data")
    else:
        data_dir = os.path.abspath(data_dir)

    fig_dir = _figures_dir(repo_root)
    suf = "" if args.spectrum_mode == "mean-per-event" else "_sum"
    out_lg = args.out_lg or os.path.join(
        fig_dir, f"photon_time_spectrum_LG_by_energy{suf}.png"
    )
    out_nolg = args.out_nolg or os.path.join(
        fig_dir, f"photon_time_spectrum_noLG_by_energy{suf}.png"
    )

    time_range: tuple[float, float] | None
    if args.auto_time_range:
        time_range = None
    else:
        time_range = (args.t_min, args.t_max)

    _plot_series(
        data_dir=data_dir,
        series=DEFAULT_LG,
        tag="LG",
        out_path=out_lg,
        max_events=args.max_events,
        nbins=args.bins,
        spectrum_mode=args.spectrum_mode,
        time_range=time_range,
    )
    _plot_series(
        data_dir=data_dir,
        series=DEFAULT_NOLG,
        tag="noLG",
        out_path=out_nolg,
        max_events=args.max_events,
        nbins=args.bins,
        spectrum_mode=args.spectrum_mode,
        time_range=time_range,
    )


if __name__ == "__main__":
    main()
