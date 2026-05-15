#!/usr/bin/env python3
"""
에너지별(60–120 GeV) **CFD(F×max bin) 시각 분포** (``plot_trigger_timing._cfd_absolute_time``,
기본 F=0.3). LG / noLG 각각 2×4 패널 (위 T1, 아래 T2).

스타일은 ``plot_cfd_per_trigger.py`` 와 맞춤: 기본은 **패널마다 데이터로 축 자동**
(``_auto_axis``, ``--bin-width`` 기본 0.02 ns), ``stairs`` + ``t=0`` 세로선 + 통계 박스.
``--fixed-range`` 를 주면 **고정 축** (기본 -5~5 ns, 200 빈, ``cfd_per_trigger --fixed-range`` 와 동일).

- ``--hist-mode mean-per-event``(기본): 빈 카운트를 **N_evt** 로 나눔.
- ``--hist-mode sum``: 이벤트 수(``stairs`` 높이) 그대로.

패널마다 **FWHM**(히스토그램 최고빈의 반고점 너비, 선형 보간)을 통계 박스에 표시.

예:
  export LD_LIBRARY_PATH=$PWD/build/rootIO:$LD_LIBRARY_PATH
  python3 analysis/plot_cfd_hist_by_energy.py
  python3 analysis/plot_cfd_hist_by_energy.py --fixed-range
"""

from __future__ import annotations

import argparse
import os
import sys

from plot_cfd_per_trigger import _auto_axis
from plot_cfd_toa_vs_nphoton_by_energy import DEFAULT_LG, DEFAULT_NOLG
from plot_trigger_timing import (
    _cfd_absolute_time,
    _figures_dir,
    _find_repo_root,
    _load_rootio,
    _merged_from_tower,
    _prepend_build_rootio_ld_path,
)

import ROOT


def _collect_cfd_times_with_nev(
    path: str,
    *,
    max_events: int,
    ns_per_unit: float,
    cfd_fraction: float,
) -> tuple[list[float], list[float], int, int, int]:
    """T1/T2 CFD 리스트, 처리 이벤트 수 n_ev, 트리거별 CFD 불가 스킵 수."""
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise FileNotFoundError(path)
    tree = f.Get("CBDsim")
    if not tree:
        f.Close()
        raise RuntimeError(f"트리 CBDsim 없음: {path}")

    evt = ROOT.CBDsimInterface.CBDsimEventData()
    tree.SetBranchAddress("CBDsimEventData", evt)

    n_ev = int(tree.GetEntries())
    if max_events >= 0:
        n_ev = min(n_ev, max_events)

    t1_vals: list[float] = []
    t2_vals: list[float] = []
    skipped_t1 = 0
    skipped_t2 = 0

    for i in range(n_ev):
        tree.GetEntry(i)
        lo1, hi1, c1 = _merged_from_tower(evt, 0)
        lo2, hi2, c2 = _merged_from_tower(evt, 1)
        v1 = _cfd_absolute_time(lo1, hi1, c1, cfd_fraction, ns_per_unit)
        v2 = _cfd_absolute_time(lo2, hi2, c2, cfd_fraction, ns_per_unit)
        if v1 is not None:
            t1_vals.append(v1)
        else:
            skipped_t1 += 1
        if v2 is not None:
            t2_vals.append(v2)
        else:
            skipped_t2 += 1

    f.Close()
    return t1_vals, t2_vals, n_ev, skipped_t1, skipped_t2


def _fwhm_ns_from_counts(counts, edges) -> float:
    """
    히스토그램 counts / bin edges 로부터 FWHM (ns).
    최고빈 높이의 절반에서 좌·우 경계를 bin 중심 선형 보간으로 구함 (단봉 가정).
    """
    import numpy as np

    c = np.asarray(counts, dtype=float)
    e = np.asarray(edges, dtype=float)
    if c.size < 2 or float(np.max(c)) <= 0.0:
        return float("nan")
    xc = 0.5 * (e[:-1] + e[1:])
    peak = float(np.max(c))
    hm = 0.5 * peak
    ip = int(np.argmax(c))

    def interp(i0: int, i1: int) -> float:
        ca, cb = float(c[i0]), float(c[i1])
        xa, xb = float(xc[i0]), float(xc[i1])
        if abs(cb - ca) < 1e-30:
            return 0.5 * (xa + xb)
        t = (hm - ca) / (cb - ca)
        t = max(0.0, min(1.0, t))
        return xa + t * (xb - xa)

    i = ip
    while i > 0 and c[i] >= hm:
        i -= 1
    if i >= 0 and c[i] < hm:
        x_left = interp(i, i + 1)
    else:
        x_left = float(e[0])

    k = ip
    while k < len(c) - 1 and c[k] >= hm:
        k += 1
    if k < len(c) and c[k] < hm:
        x_right = interp(k - 1, k)
    else:
        x_right = float(e[-1])

    return float(x_right - x_left)


def _plot_series(
    *,
    data_dir: str,
    series: list[tuple[int, str]],
    tag: str,
    out_path: str,
    max_events: int,
    nbins_fixed: int,
    hist_mode: str,
    fixed_range: bool,
    bin_width: float,
    t_fixed: tuple[float, float],
    cfd_fraction: float,
    ns_per_unit: float,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as e:
        print(f"matplotlib/numpy 없음: {e}", file=sys.stderr)
        sys.exit(1)

    fig, axes = plt.subplots(2, 4, figsize=(18, 7.5))

    for col, (egev, fname) in enumerate(series):
        path = os.path.join(data_dir, fname)
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

        t1_vals, t2_vals, n_ev, sk1, sk2 = _collect_cfd_times_with_nev(
            path,
            max_events=max_events,
            ns_per_unit=ns_per_unit,
            cfd_fraction=cfd_fraction,
        )

        for row, vals, tlabel, sk in (
            (0, t1_vals, "T1", sk1),
            (1, t2_vals, "T2", sk2),
        ):
            ax = axes[row, col]
            if not vals:
                ax.text(
                    0.5,
                    0.5,
                    "no data",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(f"{egev} GeV  {tlabel}\n{fname}")
                ax.set_xlabel(
                    f"CFD time (ns), F={cfd_fraction:g}×max bin photons"
                )
                continue

            if fixed_range:
                lo, hi = t_fixed[0], t_fixed[1]
                nb = max(1, int(nbins_fixed))
                mode_str = f"fixed [{lo:g}, {hi:g}], {nb} bins"
                counts, edges = np.histogram(vals, bins=nb, range=(lo, hi))
            else:
                lo, hi, nb = _auto_axis(vals, bin_width)
                mode_str = (
                    f"auto [{lo:.6f}, {hi:.6f}], {nb} bins (~{bin_width:g} ns/bin)"
                )
                counts, edges = np.histogram(vals, bins=nb, range=(lo, hi))

            fwhm = _fwhm_ns_from_counts(counts, edges)
            arr = np.asarray(vals, dtype=float)
            mu = float(np.mean(arr))
            sig = float(np.std(arr))

            y = counts.astype(float)
            if hist_mode == "mean-per-event":
                y = y / max(float(n_ev), 1.0)
                ylabel = "mean events / bin (÷N_evt)"
            else:
                ylabel = "events"

            ax.stairs(y, edges, color="black", linewidth=2.0, zorder=2)
            ax.axvline(0.0, color="k", ls="--", lw=0.8, alpha=0.45)
            ax.set_xlabel(
                f"CFD time (ns), F={cfd_fraction:g}×max bin photons\n({mode_str})"
            )
            ax.set_ylabel(ylabel)
            ax.set_title(f"{egev} GeV  {tlabel}\n{fname}")
            ax.grid(True, alpha=0.28, axis="y")
            ax.set_xlim(float(edges[0]), float(edges[-1]))

            fwhm_s = "nan" if np.isnan(fwhm) else f"{fwhm:.6f}"
            ax.text(
                0.97,
                0.97,
                rf"$N={len(vals)}$"
                + "\n"
                + rf"$\mu={mu:.6f}$, $\sigma={sig:.6f}$ ns"
                + "\n"
                + f"FWHM = {fwhm_s} ns",
                transform=ax.transAxes,
                fontsize=8,
                verticalalignment="top",
                horizontalalignment="right",
                bbox=dict(
                    boxstyle="round,pad=0.28",
                    facecolor="white",
                    edgecolor="0.75",
                    alpha=0.95,
                ),
            )

    mode_title = (
        "mean per event (÷N_evt)" if hist_mode == "mean-per-event" else "counts (stairs)"
    )
    axis_note = (
        f"fixed {t_fixed[0]:g}…{t_fixed[1]:g} ns"
        if fixed_range
        else f"auto axis, bin_width={bin_width:g} ns"
    )
    fig.suptitle(
        f"{tag}: CFD time F={cfd_fraction:g}×max bin ({mode_title}; {axis_note}) vs energy",
        fontsize=11,
    )
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"PNG: {out_path}")


def main() -> None:
    p = argparse.ArgumentParser(
        description="CFD time histogram by energy (LG / noLG); style aligned with plot_cfd_per_trigger.py"
    )
    p.add_argument("--data-dir", default=None)
    p.add_argument("--max-events", type=int, default=-1)
    p.add_argument("--cfd-fraction", type=float, default=0.3, dest="cfd_fraction")
    p.add_argument("--ns-per-unit", type=float, default=1.0)
    p.add_argument(
        "--fixed-range",
        action="store_true",
        help="고정 축: --t-min/--t-max/--bins (기본 -5~5 ns, 200빈; plot_cfd_per_trigger --fixed-range 와 동일)",
    )
    p.add_argument(
        "--t-min",
        type=float,
        default=-5.0,
        help="--fixed-range 일 때 CFD 시간 축 하한 (ns)",
    )
    p.add_argument(
        "--t-max",
        type=float,
        default=5.0,
        help="--fixed-range 일 때 CFD 시간 축 상한 (ns)",
    )
    p.add_argument(
        "--bin-width",
        type=float,
        default=0.02,
        dest="bin_width",
        help="자동 축일 때 빈 폭 (ns); plot_cfd_per_trigger 기본과 동일",
    )
    p.add_argument(
        "--bins",
        type=int,
        default=200,
        help="--fixed-range 일 때만 사용 (기본 200)",
    )
    p.add_argument(
        "--hist-mode",
        choices=("mean-per-event", "sum"),
        default="mean-per-event",
        help="mean-per-event: 빈 카운트÷N_evt(기본). sum: 이벤트 수 합.",
    )
    p.add_argument("-l", "--rootio-lib", default=None)
    p.add_argument("--out-lg", default=None)
    p.add_argument("--out-nolg", default=None)
    args = p.parse_args()

    if args.fixed_range and args.t_max <= args.t_min:
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
    frac_tag = f"{args.cfd_fraction:g}".replace(".", "p")
    suf = "" if args.hist_mode == "mean-per-event" else "_sum"
    out_lg = args.out_lg or os.path.join(
        fig_dir, f"cfd{frac_tag}_hist_LG_by_energy{suf}.png"
    )
    out_nolg = args.out_nolg or os.path.join(
        fig_dir, f"cfd{frac_tag}_hist_noLG_by_energy{suf}.png"
    )

    _plot_series(
        data_dir=data_dir,
        series=DEFAULT_LG,
        tag="LG",
        out_path=out_lg,
        max_events=args.max_events,
        nbins_fixed=args.bins,
        hist_mode=args.hist_mode,
        fixed_range=args.fixed_range,
        bin_width=args.bin_width,
        t_fixed=(args.t_min, args.t_max),
        cfd_fraction=args.cfd_fraction,
        ns_per_unit=args.ns_per_unit,
    )
    _plot_series(
        data_dir=data_dir,
        series=DEFAULT_NOLG,
        tag="noLG",
        out_path=out_nolg,
        max_events=args.max_events,
        nbins_fixed=args.bins,
        hist_mode=args.hist_mode,
        fixed_range=args.fixed_range,
        bin_width=args.bin_width,
        t_fixed=(args.t_min, args.t_max),
        cfd_fraction=args.cfd_fraction,
        ns_per_unit=args.ns_per_unit,
    )


if __name__ == "__main__":
    main()
