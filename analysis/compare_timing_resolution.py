#!/usr/bin/env python3
"""
두 CBDsim .root 를 읽어, `plot_trigger_timing.py` 의 **세 번째 패널(Δt)** 과 같은
히스토그램을 각각 그리고 가우시안 σ를 구합니다. 기본은 이벤트별 가중 평균 시간 차
(⟨t⟩_T1 − ⟨t⟩_T2); ``--ref cfd`` 이면 CFD 시각 차입니다.

``--ref cfd`` 일 때: merged SiPM 시간 빈을 **먼저 균일 리빈(기본 10 ps)**,
**선형 보간**으로 빈 중심을 촘촘히 만든 뒤 그 시퀀스로 CFD 시각을 계산합니다.
Δt 히스토그램 기본 범위는 **−2 ~ 2 ns** (``--xmin`` / ``--xmax`` / ``--bins`` 로 조절).

산출:
  - `analysis/t_res/figures/` 에 비교 PNG (1×2 패널, Δt 히스토그램만)
  - 동일 stem 의 `.csv` (엔트리 수, 스킵 등; Gauss 피팅 없음)

필요: ROOT(PyROOT), `build/rootIO/librootIO.so`, numpy, matplotlib

예:
  source envset.sh
  export LD_LIBRARY_PATH=$PWD/build/rootIO:$LD_LIBRARY_PATH
  python3 analysis/compare_timing_resolution.py

기본 입력은 `analysis/t_res/data/60GeV_e-_LG_0.root` 와 `60GeV_e-_noLG_0.root` 이며,
두 파일의 `CBDsim` 트리 엔트리 수가 다르면 **작은 쪽(min)** 만큼만 각각 읽어 동일 N으로 비교합니다.
`--max-events N` 을 주면 추가로 `min(N1, N2, N)` 으로 제한합니다.
`--no-match-entries` 를 주면 각 파일 전체 엔트리를 사용합니다 (N 불일치 가능).
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys

# plot_trigger_timing 과 동일 디렉터리에서 헬퍼 재사용
from plot_trigger_timing import (
    _figures_dir,
    _find_repo_root,
    _load_rootio,
    _merged_from_tower,
    _prepend_build_rootio_ld_path,
    _resolve_input_path,
    _weighted_mean_absolute,
)


def _cfd_absolute_time_float(
    lows: list[float],
    highs: list[float],
    counts: list[float],
    fraction: float,
    scale: float,
) -> float | None:
    """``_cfd_absolute_time`` 과 동일 정의; 빈 카운트는 실수 가능."""
    if not counts:
        return None
    n_max = max(counts)
    if n_max <= 0:
        return None
    thr = float(fraction) * float(n_max)
    if thr <= 0:
        return None
    idx = list(range(len(counts)))
    idx.sort(key=lambda i: 0.5 * (float(lows[i]) + float(highs[i])))
    cum = 0.0
    for i in idx:
        c = float(counts[i])
        if c <= 0:
            continue
        t_lo = float(lows[i])
        t_hi = float(highs[i])
        cum_prev = cum
        cum += c
        if cum >= thr:
            need = thr - cum_prev
            frac_in_bin = need / c if c > 0 else 0.0
            frac_in_bin = max(0.0, min(1.0, frac_in_bin))
            t_abs = t_lo + frac_in_bin * (t_hi - t_lo)
            return float(t_abs * scale)
    return None


def _rebin_merged_to_width(
    lo: list[float],
    hi: list[float],
    counts: list[float],
    bin_width_ns: float,
) -> tuple[list[float], list[float], list[float]]:
    """원래 빈을 겹침 비율로 나누어 균일 폭(예: 10 ps) 빈으로 합침."""
    import numpy as np

    if not counts:
        return [], [], []
    active = [i for i, c in enumerate(counts) if float(c) > 0.0]
    if not active:
        return lo, hi, [float(c) for c in counts]
    t_min = min(float(lo[i]) for i in active)
    t_max = max(float(hi[i]) for i in active)
    t_min -= bin_width_ns
    t_max += bin_width_ns
    if t_max <= t_min or bin_width_ns <= 0:
        return lo, hi, [float(c) for c in counts]
    edges = np.arange(t_min, t_max + bin_width_ns, bin_width_ns, dtype=float)
    if edges[-1] < t_max:
        edges = np.append(edges, t_max)
    nbin = len(edges) - 1
    new_counts = np.zeros(nbin, dtype=float)
    for blo, bhi, c in zip(lo, hi, counts):
        cf = float(c)
        if cf <= 0:
            continue
        bw = float(bhi) - float(blo)
        if bw <= 0:
            continue
        blo, bhi = float(blo), float(bhi)
        j0 = max(0, int(np.searchsorted(edges, blo, side="right") - 1))
        j1 = min(nbin, int(np.searchsorted(edges, bhi, side="left")))
        for j in range(j0, j1):
            e0, e1 = float(edges[j]), float(edges[j + 1])
            olo = max(blo, e0)
            ohi = min(bhi, e1)
            ov = max(0.0, ohi - olo)
            if ov > 0:
                new_counts[j] += cf * ov / bw
    new_lo = [float(edges[i]) for i in range(nbin)]
    new_hi = [float(edges[i + 1]) for i in range(nbin)]
    return new_lo, new_hi, new_counts.tolist()


def _interpolate_counts_linear(
    lo: list[float],
    hi: list[float],
    counts: list[float],
) -> tuple[list[float], list[float], list[float]]:
    """빈 중심에서 선형 보간해 빈 개수를 약 2배로 촘촘히 (CFD 경계 보간용)."""
    import numpy as np

    if len(counts) < 2:
        return lo, hi, [float(x) for x in counts]
    centers = 0.5 * (np.array(lo, dtype=float) + np.array(hi, dtype=float))
    c = np.array(counts, dtype=float)
    n_out = 2 * len(centers) - 1
    new_centers = np.linspace(float(centers[0]), float(centers[-1]), n_out)
    new_c = np.interp(new_centers, centers, c)
    dc = float(new_centers[1] - new_centers[0])
    new_lo = (new_centers - 0.5 * dc).tolist()
    new_hi = (new_centers + 0.5 * dc).tolist()
    return new_lo, new_hi, new_c.tolist()


def _fill_delta_hist(
    path: str,
    *,
    max_events: int,
    d_xmin: float,
    d_xmax: float,
    d_bins: int,
    ns_per_unit: float,
    hist_name: str,
    hist_title: str,
    ref_mode: str = "mean",
    cfd_fraction: float = 0.3,
    cfd_input_bin_ns: float = 0.01,
    cfd_interp: bool = True,
) -> tuple[object, int, int]:
    """ROOT TH1F(Δt) 반환. (histogram, n_used_events, skipped_delta).

    ref_mode ``mean``: 가중 평균 절대시간 차.
    ``cfd``: merged 빈을 ``cfd_input_bin_ns`` 로 리빈 후(기본 10 ps),
    옵션 선형 보간 뒤 ``_cfd_absolute_time_float`` 로 CFD.
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

    h_delta = ROOT.TH1F(
        hist_name,
        hist_title,
        d_bins,
        d_xmin,
        d_xmax,
    )
    h_delta.Sumw2()
    h_delta.SetDirectory(0)

    skipped_delta = 0
    for i in range(nmax):
        tree.GetEntry(i)
        lo1, hi1, cnt1 = _merged_from_tower(evt, 0)
        lo2, hi2, cnt2 = _merged_from_tower(evt, 1)
        if ref_mode == "cfd":
            lo1, hi1, cnt1 = _rebin_merged_to_width(lo1, hi1, cnt1, cfd_input_bin_ns)
            lo2, hi2, cnt2 = _rebin_merged_to_width(lo2, hi2, cnt2, cfd_input_bin_ns)
            if cfd_interp:
                lo1, hi1, cnt1 = _interpolate_counts_linear(lo1, hi1, cnt1)
                lo2, hi2, cnt2 = _interpolate_counts_linear(lo2, hi2, cnt2)
            mu1 = _cfd_absolute_time_float(lo1, hi1, cnt1, cfd_fraction, ns_per_unit)
            mu2 = _cfd_absolute_time_float(lo2, hi2, cnt2, cfd_fraction, ns_per_unit)
        else:
            mu1 = _weighted_mean_absolute(lo1, hi1, cnt1, ns_per_unit)
            mu2 = _weighted_mean_absolute(lo2, hi2, cnt2, ns_per_unit)
        if mu1 is not None and mu2 is not None:
            h_delta.Fill(mu1 - mu2)
        else:
            skipped_delta += 1

    f.Close()
    n_dt = int(h_delta.GetEntries())
    return h_delta, n_dt, skipped_delta


def _cbdsim_tree_entries(path: str) -> int:
    """ROOT 파일의 CBDsim 트리 엔트리 수 (파일 열었다 닫음)."""
    import ROOT

    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise FileNotFoundError(path)
    t = f.Get("CBDsim")
    if not t:
        f.Close()
        raise RuntimeError(f"트리 CBDsim 없음: {path}")
    n = int(t.GetEntries())
    f.Close()
    return n


def _plot_pair_mpl(
    axes,
    h_a,
    label_a: str,
    h_b,
    label_b: str,
    ns_per_unit: float,
    *,
    ref_mode: str = "mean",
) -> None:
    import numpy as np

    def hist_xy(h):
        cx = np.array([h.GetBinCenter(b) for b in range(1, h.GetNbinsX() + 1)])
        cy = np.array([h.GetBinContent(b) for b in range(1, h.GetNbinsX() + 1)])
        return cx, cy

    for ax, h, title in (
        (axes[0], h_a, label_a),
        (axes[1], h_b, label_b),
    ):
        x, y = hist_xy(h)
        w = (x[1] - x[0]) if len(x) > 1 else 0.1
        edges = np.linspace(float(x[0]) - w / 2, float(x[-1]) + w / 2, len(x) + 1)
        ax.stairs(y, edges, color="#1f77b4", linewidth=2.0, zorder=2)
        if ref_mode == "cfd":
            xl = (
                r"$\mathrm{CFD}_{T1} - \mathrm{CFD}_{T2}$ (ns)"
                + (f" [×{ns_per_unit}]" if ns_per_unit != 1 else "")
            )
        else:
            xl = (
                r"$\langle t\rangle_{\mathrm{T1}} - \langle t\rangle_{\mathrm{T2}}$ (ns)"
                + (f" [×{ns_per_unit}]" if ns_per_unit != 1 else "")
            )
        ax.set_xlabel(xl)
        ax.set_ylabel("events")
        ax.set_title(title)
        ax.axvline(0.0, color="k", ls="--", lw=0.8, alpha=0.5)
        mu_h = float(h.GetMean())
        sigma_h = float(h.GetRMS())
        lines = [
            rf"Entries = {int(h.GetEntries())}",
            rf"Mean = {mu_h:.4f} ns",
            rf"Std Dev = {sigma_h:.4f} ns",
        ]
        ax.text(
            0.97,
            0.97,
            "\n".join(lines),
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment="top",
            horizontalalignment="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.75", alpha=0.95),
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="두 ROOT 의 Δt 히스토그램 비교 (CFD: 10 ps 리빈·보간 후 CFD; Gauss 피팅 없음)"
    )
    parser.add_argument(
        "with_lg",
        nargs="?",
        default=None,
        help="LG 케이스 .root (기본: analysis/t_res/data/60GeV_e-_LG_0.root)",
    )
    parser.add_argument(
        "no_lg",
        nargs="?",
        default=None,
        help="no-LG 케이스 .root (기본: analysis/t_res/data/60GeV_e-_noLG_0.root)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="PNG 경로 (기본: figures/compare_delta_timing_<stem>.png)",
    )
    parser.add_argument(
        "--csv",
        dest="csv_path",
        default=None,
        help="CSV 경로 (기본: PNG 와 같은 stem .csv)",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=-1,
        help="각 파일에서 읽을 최대 이벤트 (-1=제한 없음; 기본은 두 파일 min 엔트리와 함께 사용)",
    )
    parser.add_argument(
        "--no-match-entries",
        action="store_true",
        help="트리 엔트리 수를 맞추지 않고 각 파일 전체 사용",
    )
    parser.add_argument(
        "--xmin",
        type=float,
        default=-2.0,
        help="Δt 히스토그램 하한 (ns, 기본 -2)",
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=2.0,
        help="Δt 히스토그램 상한 (ns, 기본 2)",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=400,
        help="Δt 히스토그램 빈 수 (기본 400, ±2 ns → 약 10 ps/빈)",
    )
    parser.add_argument("--delta-xmin", type=float, default=None)
    parser.add_argument("--delta-xmax", type=float, default=None)
    parser.add_argument("--delta-bins", type=int, default=None)
    parser.add_argument("--ns-per-unit", type=float, default=1.0)
    parser.add_argument(
        "--ref",
        choices=("mean", "cfd"),
        default="mean",
        help="Δt 정의: mean=가중평균 시간 차, cfd=plot_trigger_timing 과 동일 CFD",
    )
    parser.add_argument(
        "--cfd-fraction",
        type=float,
        default=0.3,
        dest="cfd_fraction",
        help="--ref cfd 일 때 임계 = F × (최대 빈 포톤 수) (기본 0.3)",
    )
    parser.add_argument(
        "--cfd-input-bin-ps",
        type=float,
        default=10.0,
        dest="cfd_input_bin_ps",
        help="--ref cfd 일 때 merged 시간축 리빈 폭 (ps, 기본 10)",
    )
    parser.add_argument(
        "--no-cfd-interp",
        action="store_true",
        help="--ref cfd 일 때 리빈 후 선형 보간 생략",
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

    path_lg = _resolve(
        args.with_lg,
        "analysis/t_res/data/60GeV_e-_LG_0.root",
    )
    path_nlg = _resolve(
        args.no_lg,
        "analysis/t_res/data/60GeV_e-_noLG_0.root",
    )

    d_xmin = args.delta_xmin if args.delta_xmin is not None else args.xmin
    d_xmax = args.delta_xmax if args.delta_xmax is not None else args.xmax
    d_bins = args.delta_bins if args.delta_bins is not None else args.bins

    _load_rootio(args.rootio_lib, repo_root)
    import ROOT

    ROOT.gROOT.SetBatch(True)

    n_lg_tree = _cbdsim_tree_entries(path_lg)
    n_nlg_tree = _cbdsim_tree_entries(path_nlg)
    if args.no_match_entries:
        max_ev_lg = n_lg_tree if args.max_events < 0 else min(n_lg_tree, args.max_events)
        max_ev_nlg = n_nlg_tree if args.max_events < 0 else min(n_nlg_tree, args.max_events)
        print(
            f"트리 엔트리: LG={n_lg_tree}, noLG={n_nlg_tree} "
            f"— 맞춤 없이 각각 {max_ev_lg} / {max_ev_nlg} 이벤트 사용"
        )
    else:
        max_ev = min(n_lg_tree, n_nlg_tree)
        if args.max_events >= 0:
            max_ev = min(max_ev, args.max_events)
        max_ev_lg = max_ev_nlg = max_ev
        print(
            f"트리 엔트리: LG={n_lg_tree}, noLG={n_nlg_tree} "
            f"→ 각각 처음 {max_ev} 이벤트로 비교 (동일 N)"
        )

    if args.ref == "cfd":
        title_lg = "LG: CFD(t)_{T1}-CFD(t)_{T2};ns;events"
        title_nlg = "no LG: CFD(t)_{T1}-CFD(t)_{T2};ns;events"
    else:
        title_lg = "LG: #LT t#GT_{T1}-#LT t#GT_{T2};ns;events"
        title_nlg = "no LG: #LT t#GT_{T1}-#LT t#GT_{T2};ns;events"

    cfd_input_bin_ns = args.cfd_input_bin_ps / 1000.0

    h_lg, n_lg, sk_lg = _fill_delta_hist(
        path_lg,
        max_events=max_ev_lg,
        d_xmin=d_xmin,
        d_xmax=d_xmax,
        d_bins=d_bins,
        ns_per_unit=args.ns_per_unit,
        hist_name="hDeltaLG",
        hist_title=title_lg,
        ref_mode=args.ref,
        cfd_fraction=args.cfd_fraction,
        cfd_input_bin_ns=cfd_input_bin_ns,
        cfd_interp=not args.no_cfd_interp,
    )
    h_nlg, n_nlg, sk_nlg = _fill_delta_hist(
        path_nlg,
        max_events=max_ev_nlg,
        d_xmin=d_xmin,
        d_xmax=d_xmax,
        d_bins=d_bins,
        ns_per_unit=args.ns_per_unit,
        hist_name="hDeltaNoLG",
        hist_title=title_nlg,
        ref_mode=args.ref,
        cfd_fraction=args.cfd_fraction,
        cfd_input_bin_ns=cfd_input_bin_ns,
        cfd_interp=not args.no_cfd_interp,
    )

    stem = "compare_delta_timing_LG_vs_noLG"
    if args.ref == "cfd":
        stem = f"{stem}_cfd"
    out_png = args.output
    if not out_png:
        out_png = os.path.join(_figures_dir(repo_root), f"{stem}.png")
    else:
        out_png = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)

    csv_path = args.csv_path
    if not csv_path:
        base, _ = os.path.splitext(out_png)
        csv_path = base + ".csv"
    else:
        csv_path = os.path.abspath(args.csv_path)

    rows = []
    meta = (
        ("LG", path_lg, h_lg, n_lg, sk_lg, n_lg_tree, max_ev_lg),
        ("no_LG", path_nlg, h_nlg, n_nlg, sk_nlg, n_nlg_tree, max_ev_nlg),
    )
    for tag, path, h_hist, n_ent, sk, n_tree, n_used in meta:
        row: dict[str, str | float | int] = {
            "sample": tag,
            "path": path,
            "tree_entries_in_file": n_tree,
            "tree_events_used": n_used,
            "matched_N": (not args.no_match_entries),
            "delta_hist_entries": n_ent,
            "skipped_delta": sk,
            "hist_mu_ns": float(h_hist.GetMean()),
            "hist_sigma_ns": float(h_hist.GetRMS()),
        }
        rows.append(row)

    fieldnames = [
        "sample",
        "path",
        "tree_entries_in_file",
        "tree_events_used",
        "matched_N",
        "hist_mu_ns",
        "hist_sigma_ns",
        "delta_hist_entries",
        "skipped_delta",
    ]
    with open(csv_path, "w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)
    print(f"CSV 저장: {csv_path}")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
        _plot_pair_mpl(
            axes,
            h_lg,
            r"w\ Light guide (Δt)",
            h_nlg,
            r"w\o Light guide (Δt)",
            args.ns_per_unit,
            ref_mode=args.ref,
        )
        fig.tight_layout()
        fig.savefig(out_png, dpi=150)
        print(f"PNG 저장: {out_png}")
        plt.close(fig)
    except ImportError as e:
        print(f"matplotlib 없음 — PNG 생략 ({e})", file=sys.stderr)


if __name__ == "__main__":
    main()
