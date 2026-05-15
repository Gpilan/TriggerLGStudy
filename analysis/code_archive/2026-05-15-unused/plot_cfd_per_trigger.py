#!/usr/bin/env python3
"""
단일 CBDsim .root 에서 이벤트마다 **트리거별 CFD 시각** (``plot_trigger_timing`` 의
``_cfd_absolute_time``, 기본 fraction=0.3)을 구해 T1·T2 **각각** 히스토그램으로 저장합니다.
Δt(T1−T2)가 아니라, 들어온 신호(합산 시간 빈)에 대한 CFD 시각 분포만 봅니다.

**기본**: 입력 파일당 **PNG 한 개**(1×2: T1 | T2). 데이터 범위는 자동,
``--bin-width`` (기본 0.01 ns). 각 패널에 **N, μ, σ** 을 한 박스에 넣습니다.

필요: ROOT(PyROOT), ``build/rootIO/librootIO.so``, matplotlib, numpy

예:
  source envset.sh
  export LD_LIBRARY_PATH=$PWD/build/rootIO:$LD_LIBRARY_PATH
  python3 analysis/plot_cfd_per_trigger.py
  python3 analysis/plot_cfd_per_trigger.py --fixed-range --xmin -5 --xmax 5 --bins 400
"""

from __future__ import annotations

import argparse
import math
import os
import sys

from plot_trigger_timing import (
    _cfd_absolute_time,
    _figures_dir,
    _find_repo_root,
    _load_rootio,
    _merged_from_tower,
    _prepend_build_rootio_ld_path,
    _resolve_input_path,
)


def _collect_cfd_times(
    path: str,
    *,
    max_events: int,
    ns_per_unit: float,
    cfd_fraction: float,
) -> tuple[list[float], list[float], int, int, int, int]:
    """T1/T2 CFD 값 리스트와 스킵 카운트."""
    import ROOT

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

    t1_vals: list[float] = []
    t2_vals: list[float] = []
    skipped_t1 = 0
    skipped_t2 = 0

    for i in range(nmax):
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
    return (
        t1_vals,
        t2_vals,
        len(t1_vals),
        len(t2_vals),
        skipped_t1,
        skipped_t2,
    )


def _auto_axis(
    vals: list[float],
    bin_width: float,
    *,
    margin_frac: float = 0.08,
) -> tuple[float, float, int]:
    """
    데이터 범위 + 여백으로 (lo, hi, nbins).
    값이 거의 같으면 bin_width 기준으로 최소 폭을 갖게 펼침.
    """
    if not vals:
        return 0.0, 1.0, 1
    mn = min(vals)
    mx = max(vals)
    span = mx - mn
    if span <= 0.0:
        pad = max(bin_width * 5.0, 1e-9)
        lo, hi = mn - pad, mx + pad
    else:
        pad = max(span * margin_frac, bin_width * 2.0)
        lo, hi = mn - pad, mx + pad
    span_e = hi - lo
    nbins = max(1, int(math.ceil(span_e / max(bin_width, 1e-15))))
    nbins = min(nbins, 5000)
    return lo, hi, nbins


def _fixed_axis(xmin: float, xmax: float, bins: int) -> tuple[float, float, int]:
    if xmax <= xmin:
        raise SystemExit("--xmax must be > --xmin")
    return xmin, xmax, max(1, bins)


def _light_guide_title(stem: str) -> str:
    """파일 stem 으로 light guide 유무 문구 (플롯 제목용)."""
    if "nolg" in stem.lower():
        return "Without light guide"
    return "With light guide"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="단일 ROOT: T1/T2 각각 CFD 시각 분포 (plot_trigger_timing 과 동일 CFD 정의)"
    )
    parser.add_argument(
        "input",
        nargs="?",
        default=None,
        help="입력 .root (기본: analysis/t_res/data/60GeV_e-_noLG_1.root)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="출력 PNG (기본: figures/cfd_per_trigger_<stem>.png, T1·T2 동시 포함)",
    )
    parser.add_argument(
        "--fixed-range",
        action="store_true",
        help="고정 축: --xmin/--xmax/--bins 사용 (기본은 데이터 자동 범위 + --bin-width)",
    )
    parser.add_argument("--xmin", type=float, default=-5.0)
    parser.add_argument("--xmax", type=float, default=5.0)
    parser.add_argument(
        "--bins",
        type=int,
        default=400,
        help="--fixed-range 일 때만 사용 (기본 400, 더 촘촘한 빈)",
    )
    parser.add_argument(
        "--bin-width",
        type=float,
        default=0.01,
        dest="bin_width",
        help="자동 범위일 때 빈 폭 (ns). 작을수록 세밀 (기본 0.01 ns)",
    )
    parser.add_argument("--ns-per-unit", type=float, default=1.0)
    parser.add_argument(
        "--cfd-fraction",
        type=float,
        default=0.3,
        dest="cfd_fraction",
        help="임계 = F × (최대 빈 포톤 수) (기본 0.3)",
    )
    parser.add_argument("--max-events", type=int, default=-1)
    parser.add_argument("-l", "--rootio-lib", default=None)
    args = parser.parse_args()

    repo_root = _find_repo_root()
    if repo_root is None:
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    _prepend_build_rootio_ld_path(repo_root)

    def _resolve(p: str | None, default_rel: str) -> str:
        if p:
            r = _resolve_input_path(p, repo_root)
            if r:
                return r
            raise SystemExit(f"파일 없음: {p}")
        full = os.path.join(repo_root, default_rel)
        if os.path.isfile(full):
            return os.path.abspath(full)
        raise SystemExit(f"기본 경로에 파일 없음: {full}")

    in_path = _resolve(
        args.input,
        "analysis/t_res/data/60GeV_e-_noLG_1.root",
    )

    _load_rootio(args.rootio_lib, repo_root)
    import ROOT

    ROOT.gROOT.SetBatch(True)

    t1_vals, t2_vals, f1, f2, sk1, sk2 = _collect_cfd_times(
        in_path,
        max_events=args.max_events,
        ns_per_unit=args.ns_per_unit,
        cfd_fraction=args.cfd_fraction,
    )

    stem = os.path.splitext(os.path.basename(in_path))[0]
    fig_dir = _figures_dir(repo_root)
    if args.output:
        out = os.path.abspath(args.output)
        if not os.path.splitext(out)[1]:
            out += ".png"
    else:
        out = os.path.join(fig_dir, f"cfd_per_trigger_{stem}.png")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    print(
        f"입력: {in_path}\n"
        f"  T1: filled={f1}, skipped(empty CFD)={sk1}\n"
        f"  T2: filled={f2}, skipped(empty CFD)={sk2}"
    )

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        if t1_vals:
            print(
                f"  T1 CFD (ns): min={min(t1_vals):.6f}, max={max(t1_vals):.6f}, "
                f"std={float(np.std(t1_vals)):.6f}"
            )
        if t2_vals:
            print(
                f"  T2 CFD (ns): min={min(t2_vals):.6f}, max={max(t2_vals):.6f}, "
                f"std={float(np.std(t2_vals)):.6f}"
            )

        guide_line = _light_guide_title(stem)

        def _draw_panel(ax, vals: list[float], trigger_id: str) -> None:
            """단일 축에 히스토그램·통계(N, μ, σ) 표시."""
            ax.set_title(trigger_id, fontsize=12.5, fontweight="600", pad=8)

            if not vals:
                ax.text(
                    0.5,
                    0.5,
                    "no data",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=11,
                    color="0.45",
                )
                ax.set_xlabel("CFD time (ns)", fontsize=10)
                ax.set_ylabel("events", fontsize=10)
                return

            if args.fixed_range:
                lo, hi, nb = _fixed_axis(args.xmin, args.xmax, args.bins)
                mode = f"axis: [{lo:g}, {hi:g}] ns · {nb} bins"
            else:
                lo, hi, nb = _auto_axis(vals, args.bin_width)
                mode = (
                    f"axis: [{lo:.6f}, {hi:.6f}] ns · {nb} bins "
                    f"(Δ≈{args.bin_width:g} ns)"
                )

            counts, edges = np.histogram(vals, bins=nb, range=(lo, hi))
            ax.stairs(counts, edges, color="black", linewidth=2.0, zorder=2)

            arr = np.asarray(vals, dtype=float)
            mu = float(np.mean(arr))
            sig = float(np.std(arr))

            ax.axvline(0.0, color="k", ls="--", lw=0.85, alpha=0.4, zorder=1)
            ax.set_xlabel("CFD time (ns)", fontsize=10.5)
            ax.set_ylabel("events", fontsize=10.5)
            ax.tick_params(axis="both", labelsize=9)
            ax.grid(True, axis="y", alpha=0.2, linestyle="-", linewidth=0.55, zorder=0)

            stat_block = (
                r"$\mathbf{Statistics}$" + "\n\n"
                rf"$N$     $= {len(vals)}$" + "\n"
                rf"$\mu$    $= {mu:.4f}$ ns" + "\n"
                rf"$\sigma$ $= {sig:.4f}$ ns"
            )

            ax.text(
                0.98,
                0.98,
                stat_block,
                transform=ax.transAxes,
                fontsize=9,
                verticalalignment="top",
                horizontalalignment="right",
                linespacing=1.35,
                bbox=dict(
                    boxstyle="round,pad=0.55",
                    facecolor="#fffefb",
                    edgecolor="#b0b0b0",
                    linewidth=1.0,
                    alpha=0.98,
                ),
                zorder=8,
            )

            ax.text(
                0.5,
                -0.14,
                f"{mode} · F={args.cfd_fraction:g}× max bin photons",
                transform=ax.transAxes,
                ha="center",
                va="top",
                fontsize=7.5,
                color="0.42",
            )

        fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2), constrained_layout=False)
        fig.patch.set_facecolor("white")
        fig.suptitle(
            f"{guide_line}\n"
            f"Per-event CFD time (F={args.cfd_fraction:g}× max bin photons)",
            fontsize=11.5,
            fontweight="500",
            y=0.995,
        )
        _draw_panel(axes[0], t1_vals, "T1")
        _draw_panel(axes[1], t2_vals, "T2")
        fig.tight_layout(rect=[0, 0.02, 1, 0.86])
        fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"PNG 저장: {out}")
    except ImportError as e:
        print(f"matplotlib 없음 — PNG 생략 ({e})", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
