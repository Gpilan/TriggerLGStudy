#!/usr/bin/env python3
"""
동일 이벤트에서 T1 CFD 시각 vs T2 CFD 시각 (X, Y)으로 **상관**을 봅니다.

- 이벤트마다 ``(t_CFD^T1, t_CFD^T2)`` 한 점 (둘 다 유효할 때만).
- 산점도 + 2D 히스토그램, Pearson r / Spearman ρ 출력.

필요: ROOT, librootIO, matplotlib, numpy, scipy(선택: Spearman)

예:
  python3 analysis/plot_cfd_t1_t2_correlation.py analysis/t_res/data/60GeV_e-_noLG_1.root
"""

from __future__ import annotations

import argparse
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


def _collect_pairs(
    path: str,
    *,
    max_events: int,
    ns_per_unit: float,
    cfd_fraction: float,
) -> tuple[list[float], list[float], int]:
    """같은 이벤트 인덱스의 (T1, T2) CFD 쌍. 둘 중 하나만 없으면 해당 이벤트는 제외."""
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

    x_t1: list[float] = []
    y_t2: list[float] = []
    skipped = 0

    for i in range(nmax):
        tree.GetEntry(i)
        lo1, hi1, c1 = _merged_from_tower(evt, 0)
        lo2, hi2, c2 = _merged_from_tower(evt, 1)
        v1 = _cfd_absolute_time(lo1, hi1, c1, cfd_fraction, ns_per_unit)
        v2 = _cfd_absolute_time(lo2, hi2, c2, cfd_fraction, ns_per_unit)
        if v1 is not None and v2 is not None:
            x_t1.append(v1)
            y_t2.append(v2)
        else:
            skipped += 1

    f.Close()
    return x_t1, y_t2, skipped


def main() -> None:
    parser = argparse.ArgumentParser(
        description="T1 CFD vs T2 CFD (동일 이벤트) 상관 플롯"
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
        help="출력 PNG (기본: figures/cfd_T1_vs_T2_<stem>.png)",
    )
    parser.add_argument(
        "--cfd-fraction",
        type=float,
        default=0.3,
        dest="cfd_fraction",
    )
    parser.add_argument("--ns-per-unit", type=float, default=1.0)
    parser.add_argument("--max-events", type=int, default=-1)
    parser.add_argument(
        "--bins2d",
        type=int,
        default=40,
        help="2D 히스토그램 빈 수 (한 축)",
    )
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

    x, y, skipped = _collect_pairs(
        in_path,
        max_events=args.max_events,
        ns_per_unit=args.ns_per_unit,
        cfd_fraction=args.cfd_fraction,
    )

    stem = os.path.splitext(os.path.basename(in_path))[0]
    out = args.output
    if not out:
        out = os.path.join(_figures_dir(repo_root), f"cfd_T1_vs_T2_{stem}.png")
    else:
        out = os.path.abspath(out)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    import numpy as np

    if len(x) < 2:
        print(f"쌍 부족 (N={len(x)}), skipped_events={skipped}", file=sys.stderr)
        sys.exit(1)

    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    r_pearson = float(np.corrcoef(xa, ya)[0, 1])

    rho = None
    try:
        from scipy import stats

        rho = float(stats.spearmanr(xa, ya).statistic)
    except Exception:
        pass

    print(
        f"입력: {in_path}\n"
        f"  유효 쌍 N={len(x)}, T1/T2 둘 중 하나라도 CFD 불가였던 이벤트={skipped}\n"
        f"  Pearson r = {r_pearson:.6f}"
        + (f"\n  Spearman ρ = {rho:.6f}" if rho is not None else "")
    )

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

        ax = axes[0]
        ax.scatter(xa, ya, s=8, alpha=0.35, c="#1f4e79", edgecolors="none")
        ax.set_xlabel(r"T1 CFD time (ns)")
        ax.set_ylabel(r"T2 CFD time (ns)")
        ax.set_aspect("equal", adjustable="box")
        lim_lo = float(min(xa.min(), ya.min()))
        lim_hi = float(max(xa.max(), ya.max()))
        pad = (lim_hi - lim_lo) * 0.05 + 1e-6
        ax.set_xlim(lim_lo - pad, lim_hi + pad)
        ax.set_ylim(lim_lo - pad, lim_hi + pad)
        ax.plot([lim_lo - pad, lim_hi + pad], [lim_lo - pad, lim_hi + pad], "k--", lw=0.8, alpha=0.4)
        ax.set_title("Scatter (same event)")
        ax.grid(True, alpha=0.3)

        ax2 = axes[1]
        h2d, xe, ye, im = ax2.hist2d(
            xa,
            ya,
            bins=args.bins2d,
            cmap="Blues",
        )
        ax2.set_xlabel(r"T1 CFD time (ns)")
        ax2.set_ylabel(r"T2 CFD time (ns)")
        ax2.set_title("2D histogram")
        plt.colorbar(im, ax=ax2, label="events")

        stat = rf"$N={len(x)}$, Pearson $r={r_pearson:.4f}$"
        if rho is not None:
            stat += rf", Spearman $\rho={rho:.4f}$"
        fig.suptitle(
            f"CFD(F={args.cfd_fraction:g}) T1 vs T2 — {os.path.basename(in_path)}\n{stat}",
            fontsize=10,
        )
        fig.tight_layout()
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"PNG 저장: {out}")
    except ImportError as e:
        print(f"matplotlib 없음 ({e})", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
