#!/usr/bin/env python3
"""
에너지별(60–120 GeV)로 **CFD TOA vs 트리거별 총 포톤 수 N** 산점도.

- `plot_cfd_toa_vs_nphoton.py` 와 동일한 `_collect_toa_n` 정의.
- LG 시리즈 / noLG 시리즈 각각 **한 개 PNG** (2행×4열: 위 T1, 아래 T2).

기본 입력 (``analysis/t_res/data/``):

  LG:   LG_60GeV_e-_LG_1, LG_80GeV_0, LG_100GeV_0, LG_120GeV_0
  noLG: NoLG_60GeV_e-_noLG_1, NoLG_stats_{80,100,120}GeV_0

예:
  export LD_LIBRARY_PATH=$PWD/build/rootIO:$LD_LIBRARY_PATH
  python3 analysis/plot_cfd_toa_vs_nphoton_by_energy.py
"""

from __future__ import annotations

import argparse
import os
import sys

from plot_cfd_toa_vs_nphoton import _collect_toa_n
from plot_trigger_timing import (
    _figures_dir,
    _find_repo_root,
    _load_rootio,
    _prepend_build_rootio_ld_path,
)

# (GeV, 파일명) — analysis/t_res/data/ 기준
DEFAULT_LG = [
    (60, "LG_60GeV_e-_LG_1.root"),
    (80, "LG_80GeV_0.root"),
    (100, "LG_100GeV_0.root"),
    (120, "LG_120GeV_0.root"),
]

DEFAULT_NOLG = [
    (60, "NoLG_60GeV_e-_noLG_1.root"),
    (80, "NoLG_stats_80GeV_0.root"),
    (100, "NoLG_stats_100GeV_0.root"),
    (120, "NoLG_stats_120GeV_0.root"),
]


def _resolve_data_dir(repo_root: str, user_dir: str | None) -> str:
    if user_dir:
        return os.path.abspath(user_dir)
    return os.path.join(repo_root, "analysis", "t_res", "data")


def _plot_one_series(
    *,
    repo_root: str,
    data_dir: str,
    series: list[tuple[int, str]],
    tag: str,
    out_path: str,
    cfd_fraction: float,
    ns_per_unit: float,
    max_events: int,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        print(f"matplotlib 없음: {e}", file=sys.stderr)
        sys.exit(1)

    fig, axes = plt.subplots(2, 4, figsize=(18, 7.5))
    color = "#c0392b" if tag == "LG" else "#2980b9"

    for col, (egev, fname) in enumerate(series):
        path = os.path.join(data_dir, fname)
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
        t1, n1, t2, n2, _sk = _collect_toa_n(
            path,
            max_events=max_events,
            ns_per_unit=ns_per_unit,
            cfd_fraction=cfd_fraction,
        )
        for ax, t, n, trig in (
            (axes[0, col], t1, n1, "T1"),
            (axes[1, col], t2, n2, "T2"),
        ):
            ax.scatter(t, n, s=10, alpha=0.35, c=color, edgecolors="none")
            ax.set_xlabel("CFD TOA (ns)")
            ax.set_ylabel("N photons")
            ax.set_title(f"{egev} GeV  {trig}\n{fname}")
            ax.grid(True, alpha=0.25)

    fig.suptitle(
        f"{tag}: CFD TOA vs merged photon count (F={cfd_fraction:g}×max bin) — per energy",
        fontsize=11,
    )
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"PNG: {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="CFD vs N by energy (LG / noLG 각 1 PNG)")
    parser.add_argument(
        "--data-dir",
        default=None,
        help="ROOT 모음 디렉터리 (기본: analysis/t_res/data)",
    )
    parser.add_argument("--cfd-fraction", type=float, default=0.3, dest="cfd_fraction")
    parser.add_argument("--ns-per-unit", type=float, default=1.0)
    parser.add_argument("--max-events", type=int, default=-1)
    parser.add_argument("-l", "--rootio-lib", default=None)
    parser.add_argument(
        "--out-lg",
        default=None,
        help="LG 출력 PNG 경로",
    )
    parser.add_argument(
        "--out-nolg",
        default=None,
        help="noLG 출력 PNG 경로",
    )
    args = parser.parse_args()

    repo_root = _find_repo_root()
    if repo_root is None:
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    _prepend_build_rootio_ld_path(repo_root)
    _load_rootio(args.rootio_lib, repo_root)

    import ROOT

    ROOT.gROOT.SetBatch(True)

    data_dir = _resolve_data_dir(repo_root, args.data_dir)
    fig_dir = _figures_dir(repo_root)

    out_lg = args.out_lg or os.path.join(fig_dir, "cfd_toa_vs_N_LG_by_energy.png")
    out_nolg = args.out_nolg or os.path.join(fig_dir, "cfd_toa_vs_N_noLG_by_energy.png")

    _plot_one_series(
        repo_root=repo_root,
        data_dir=data_dir,
        series=DEFAULT_LG,
        tag="LG",
        out_path=out_lg,
        cfd_fraction=args.cfd_fraction,
        ns_per_unit=args.ns_per_unit,
        max_events=args.max_events,
    )
    _plot_one_series(
        repo_root=repo_root,
        data_dir=data_dir,
        series=DEFAULT_NOLG,
        tag="noLG",
        out_path=out_nolg,
        cfd_fraction=args.cfd_fraction,
        ns_per_unit=args.ns_per_unit,
        max_events=args.max_events,
    )


if __name__ == "__main__":
    main()
