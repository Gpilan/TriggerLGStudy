#!/usr/bin/env python3
"""
이벤트마다 **CFD 시각(TOA)** 과 **해당 트리거의 총 광자 수 N** 을 (X, Y)로 산점도.

- CFD: ``plot_trigger_timing._cfd_absolute_time`` (기본 F=0.3), 트리거별 ``_merged_from_tower`` 합산 빈.
- N: 같은 merged 빈의 ``sum(counts)`` (``siPMPhotonSumTrig*`` 와 동일해야 함).

두 샘플(예: LG / noLG)을 색으로 겹쳐 비교 (T1 패널, T2 패널).

예:
  export LD_LIBRARY_PATH=$PWD/build/rootIO:$LD_LIBRARY_PATH
  python3 analysis/plot_cfd_toa_vs_nphoton.py
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


def _nphoton(lo, hi, cnt) -> int:
    if not cnt:
        return 0
    return int(sum(int(c) for c in cnt))


def _collect_toa_n(
    path: str,
    *,
    max_events: int,
    ns_per_unit: float,
    cfd_fraction: float,
) -> tuple[list[float], list[int], list[float], list[int], int]:
    """
    (t1, n1, t2, n2, skipped) — CFD 불가인 트리거는 해당 축에서 제외(쌍 불일치 허용:
    한 이벤트에서 T1만 있고 T2 없을 수 있음 → T1만 점 하나).
    """
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

    t1: list[float] = []
    n1: list[int] = []
    t2: list[float] = []
    n2: list[int] = []
    skipped = 0

    for i in range(nmax):
        tree.GetEntry(i)
        lo1, hi1, c1 = _merged_from_tower(evt, 0)
        lo2, hi2, c2 = _merged_from_tower(evt, 1)
        v1 = _cfd_absolute_time(lo1, hi1, c1, cfd_fraction, ns_per_unit)
        v2 = _cfd_absolute_time(lo2, hi2, c2, cfd_fraction, ns_per_unit)
        np1 = _nphoton(lo1, hi1, c1)
        np2 = _nphoton(lo2, hi2, c2)

        if v1 is not None:
            t1.append(v1)
            n1.append(np1)
        else:
            skipped += 1
        if v2 is not None:
            t2.append(v2)
            n2.append(np2)
        else:
            skipped += 1

    f.Close()
    return t1, n1, t2, n2, skipped


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CFD TOA vs N photons (per trigger), LG vs noLG scatter"
    )
    parser.add_argument(
        "--lg",
        default=None,
        help="LG .root (기본: analysis/t_res/data/60GeV_e-_LG_1.root)",
    )
    parser.add_argument(
        "--no-lg",
        dest="no_lg",
        default=None,
        help="no-LG .root (기본: analysis/t_res/data/60GeV_e-_noLG_1.root)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="출력 PNG (기본: figures/cfd_toa_vs_N_<lg_stem>_vs_<nlg_stem>.png)",
    )
    parser.add_argument("--cfd-fraction", type=float, default=0.3, dest="cfd_fraction")
    parser.add_argument("--ns-per-unit", type=float, default=1.0)
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

    path_lg = _resolve(args.lg, "analysis/t_res/data/60GeV_e-_LG_1.root")
    path_nlg = _resolve(args.no_lg, "analysis/t_res/data/60GeV_e-_noLG_1.root")

    _load_rootio(args.rootio_lib, repo_root)
    import ROOT

    ROOT.gROOT.SetBatch(True)

    t1_lg, n1_lg, t2_lg, n2_lg, sk_lg = _collect_toa_n(
        path_lg,
        max_events=args.max_events,
        ns_per_unit=args.ns_per_unit,
        cfd_fraction=args.cfd_fraction,
    )
    t1_n, n1_n, t2_n, n2_n, sk_n = _collect_toa_n(
        path_nlg,
        max_events=args.max_events,
        ns_per_unit=args.ns_per_unit,
        cfd_fraction=args.cfd_fraction,
    )

    stem_lg = os.path.splitext(os.path.basename(path_lg))[0]
    stem_nlg = os.path.splitext(os.path.basename(path_nlg))[0]
    out = args.output
    if not out:
        out = os.path.join(
            _figures_dir(repo_root),
            f"cfd_toa_vs_N_{stem_lg}_vs_{stem_nlg}.png",
        )
    else:
        out = os.path.abspath(out)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    print(
        f"LG:     {path_lg}\n"
        f"  T1 points={len(t1_lg)}, T2 points={len(t2_lg)}, skipped(trig CFD None)={sk_lg}\n"
        f"noLG:   {path_nlg}\n"
        f"  T1 points={len(t1_n)}, T2 points={len(t2_n)}, skipped={sk_n}"
    )

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

        def _scatter(ax, xlg, ylg, xn, yn, title: str) -> None:
            ax.scatter(
                xlg,
                ylg,
                s=12,
                alpha=0.35,
                c="#c0392b",
                edgecolors="none",
                label=f"LG ({stem_lg})",
            )
            ax.scatter(
                xn,
                yn,
                s=12,
                alpha=0.35,
                c="#2980b9",
                edgecolors="none",
                label=f"no LG ({stem_nlg})",
            )
            ax.set_xlabel(r"CFD TOA (ns), $F=" + f"{args.cfd_fraction:g}" + r"\times N_{\mathrm{bin}}^{\max}$")
            ax.set_ylabel("N photons (merged SiPM sum)")
            ax.set_title(title)
            ax.legend(loc="upper right", fontsize=8)
            ax.grid(True, alpha=0.3)

        _scatter(
            axes[0],
            t1_lg,
            n1_lg,
            t1_n,
            n1_n,
            "T1 (trigger 0)",
        )
        _scatter(
            axes[1],
            t2_lg,
            n2_lg,
            t2_n,
            n2_n,
            "T2 (trigger 1)",
        )

        fig.suptitle(
            "CFD time-of-arrival vs total photon count (per trigger, per event)",
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
