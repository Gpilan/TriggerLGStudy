#!/usr/bin/env python3
"""
에너지별로 트리거당 **이벤트별 총 포톤 수 N** (merged SiPM sum) 히스토그램.

LG / noLG 각각 2×4 패널 (위 T1, 아래 T2, 열=60·80·100·120 GeV).

입력 ROOT·파일명 기본값은 ``plot_cfd_toa_vs_nphoton_by_energy.py`` 와 동일.

예:
  export LD_LIBRARY_PATH=$PWD/build/rootIO:$LD_LIBRARY_PATH
  python3 analysis/plot_nphoton_hist_by_energy.py
"""

from __future__ import annotations

import argparse
import os
import sys

from plot_cfd_toa_vs_nphoton import _nphoton
from plot_cfd_toa_vs_nphoton_by_energy import DEFAULT_LG, DEFAULT_NOLG
from plot_trigger_timing import (
    _figures_dir,
    _find_repo_root,
    _load_rootio,
    _merged_from_tower,
    _prepend_build_rootio_ld_path,
)

import ROOT


def _collect_nphotons_per_trigger(
    path: str,
    *,
    max_events: int,
) -> tuple[list[int], list[int]]:
    """각 이벤트의 (T1 총 포톤, T2 총 포톤)."""
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

    n1: list[int] = []
    n2: list[int] = []
    for i in range(nmax):
        tree.GetEntry(i)
        lo1, hi1, c1 = _merged_from_tower(evt, 0)
        lo2, hi2, c2 = _merged_from_tower(evt, 1)
        n1.append(_nphoton(lo1, hi1, c1))
        n2.append(_nphoton(lo2, hi2, c2))
    f.Close()
    return n1, n2


def _bins_integer(arr: list[int]):
    import numpy as np

    if not arr:
        return np.array([-0.5, 0.5, 1.5])
    lo, hi = min(arr), max(arr)
    if lo == hi:
        return np.arange(lo - 0.5, lo + 2.5, 1.0)
    return np.arange(lo - 0.5, hi + 1.5, 1.0)


def _plot_series(
    *,
    data_dir: str,
    series: list[tuple[int, str]],
    tag: str,
    out_path: str,
    max_events: int,
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
        arr1, arr2 = _collect_nphotons_per_trigger(path, max_events=max_events)
        for ax, arr, trig in (
            (axes[0, col], arr1, "T1"),
            (axes[1, col], arr2, "T2"),
        ):
            bins = _bins_integer(arr)
            ax.hist(
                arr,
                bins=bins,
                color=color,
                edgecolor="black",
                linewidth=0.35,
                alpha=0.85,
            )
            mu = float(np.mean(arr)) if arr else 0.0
            sig = float(np.std(arr)) if arr else 0.0
            ax.set_xlabel("N photons (merged)")
            ax.set_ylabel("events")
            ax.set_title(
                f"{egev} GeV  {trig}\n{fname}\n"
                rf"$\mu={mu:.1f}$, $\sigma={sig:.1f}$"
            )
            ax.grid(True, alpha=0.28, axis="y")

    fig.suptitle(
        f"{tag}: photon multiplicity per event (merged SiPM), by beam energy",
        fontsize=11,
    )
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"PNG: {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="N photons/event histogram by energy")
    parser.add_argument("--data-dir", default=None)
    parser.add_argument("--max-events", type=int, default=-1)
    parser.add_argument("-l", "--rootio-lib", default=None)
    parser.add_argument("--out-lg", default=None)
    parser.add_argument("--out-nolg", default=None)
    args = parser.parse_args()

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
    out_lg = args.out_lg or os.path.join(fig_dir, "nphoton_per_event_LG_by_energy.png")
    out_nolg = args.out_nolg or os.path.join(fig_dir, "nphoton_per_event_noLG_by_energy.png")

    _plot_series(
        data_dir=data_dir,
        series=DEFAULT_LG,
        tag="LG",
        out_path=out_lg,
        max_events=args.max_events,
    )
    _plot_series(
        data_dir=data_dir,
        series=DEFAULT_NOLG,
        tag="noLG",
        out_path=out_nolg,
        max_events=args.max_events,
    )


if __name__ == "__main__":
    main()
