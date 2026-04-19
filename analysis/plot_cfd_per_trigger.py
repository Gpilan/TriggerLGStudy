#!/usr/bin/env python3
"""
단일 CBDsim .root 에서 이벤트마다 **트리거별 CFD 시각** (``plot_trigger_timing`` 의
``_cfd_absolute_time``, 기본 fraction=0.3)을 구해 T1·T2 **각각** 히스토그램으로 저장합니다.
Δt(T1−T2)가 아니라, 들어온 신호(합산 시간 빈)에 대한 CFD 시각 분포만 봅니다.

**기본**: 데이터에서 CFD 최솟값·최댓값을 읽어 **축 범위를 자동**으로 잡고,
``--bin-width`` (기본 0.02 ns)로 **빈 폭을 작게** 잡아 분포가 한 빈에 몰리지 않게 합니다.
(이전처럼 ±5 ns 에 10 빈이면 1 ns/빈이라, ns 단위로 거의 같은 CFD는 한 칸에만 들어갈 수 있음.)

필요: ROOT(PyROOT), ``build/rootIO/librootIO.so``, matplotlib, numpy

예:
  source envset.sh
  export LD_LIBRARY_PATH=$PWD/build/rootIO:$LD_LIBRARY_PATH
  python3 analysis/plot_cfd_per_trigger.py
  python3 analysis/plot_cfd_per_trigger.py --fixed-range --xmin -5 --xmax 5 --bins 200
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
        help="출력 PNG (기본: figures/cfd_per_trigger_<stem>.png)",
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
        default=200,
        help="--fixed-range 일 때만 사용 (기본 200)",
    )
    parser.add_argument(
        "--bin-width",
        type=float,
        default=0.02,
        dest="bin_width",
        help="자동 범위일 때 빈 폭 (ns). 작을수록 세밀 (기본 0.02 ns)",
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
    out = args.output
    if not out:
        out = os.path.join(_figures_dir(repo_root), f"cfd_per_trigger_{stem}.png")
    else:
        out = os.path.abspath(out)
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

        def _panel(ax, vals: list[float], tag: str, bn: str) -> None:
            if not vals:
                ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"{tag}\n{bn}")
                return
            if args.fixed_range:
                lo, hi, nb = _fixed_axis(args.xmin, args.xmax, args.bins)
                mode = f"fixed [{lo:g}, {hi:g}], {nb} bins"
            else:
                lo, hi, nb = _auto_axis(vals, args.bin_width)
                mode = f"auto [{lo:.6f}, {hi:.6f}], {nb} bins (~{args.bin_width:g} ns/bin)"
            counts, edges = np.histogram(vals, bins=nb, range=(lo, hi))
            ax.stairs(counts, edges, color="black", linewidth=2.0, zorder=2)
            ax.set_xlabel(
                f"CFD time (ns), F={args.cfd_fraction:g}×max bin photons\n({mode})"
            )
            ax.set_ylabel("events")
            ax.set_title(f"{tag}\n{bn}")
            ax.axvline(0.0, color="k", ls="--", lw=0.8, alpha=0.45)
            arr = np.asarray(vals, dtype=float)
            mu = float(np.mean(arr))
            sig = float(np.std(arr))
            ax.text(
                0.97,
                0.97,
                rf"$N={len(vals)}$"
                + "\n"
                + rf"$\mu={mu:.6f}$, $\sigma={sig:.6f}$ ns",
                transform=ax.transAxes,
                fontsize=9,
                verticalalignment="top",
                horizontalalignment="right",
                bbox=dict(
                    boxstyle="round,pad=0.3",
                    facecolor="white",
                    edgecolor="0.75",
                    alpha=0.95,
                ),
            )

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
        bn = os.path.basename(in_path)
        _panel(axes[0], t1_vals, "T1 (trigger 0)", bn)
        _panel(axes[1], t2_vals, "T2 (trigger 1)", bn)

        fig.suptitle(
            "Per-event CFD absolute time (same definition as plot_trigger_timing.py)",
            fontsize=10,
        )
        fig.tight_layout()
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"PNG 저장: {out}")
    except ImportError as e:
        print(f"matplotlib 없음 — PNG 생략 ({e})", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
