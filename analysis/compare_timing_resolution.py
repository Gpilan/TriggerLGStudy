#!/usr/bin/env python3
"""
두 CBDsim .root 를 읽어, `plot_trigger_timing.py` 의 **세 번째 패널(Δt)** 과 동일한
히스토그램(이벤트별 ⟨t⟩_T1 − ⟨t⟩_T2)을 각각 그리고 가우시안 σ를 구합니다.

**타이밍 레졸루션** (정의): σ / √2  (ns)

산출:
  - `analysis/t_res/figures/` 에 비교 PNG (1×2 패널, Δt + Gauss 피트)
  - 동일 stem 의 `.csv` (σ, σ_err, μ, μ_err, timing_resolution_ns 등)

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

SQRT2 = math.sqrt(2.0)


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
) -> tuple[object, int, int]:
    """ROOT TH1F(Δt) 반환. (histogram, n_used_events, skipped_delta)."""
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


def _fit_gaus_delta(h, ROOT, tag: str) -> dict | None:
    """ROOT TF1 gaus 피팅. 실패 시 None."""
    xmin = float(h.GetXaxis().GetXmin())
    xmax = float(h.GetXaxis().GetXmax())
    if h.GetSumOfWeights() <= 0:
        return None
    fitfn = ROOT.TF1(f"fit_delta_gaus_{tag}", "gaus", xmin, xmax)
    fitfn.SetNpx(500)
    fitfn.SetParameter(0, float(h.GetMaximum()))
    fitfn.SetParameter(1, float(h.GetMean()))
    fitfn.SetParameter(2, max(float(h.GetRMS()), 1e-6))
    h.Fit(fitfn, "QN", "", xmin, xmax)
    return {
        "mu": float(fitfn.GetParameter(1)),
        "mu_err": float(fitfn.GetParError(1)),
        "sigma": float(fitfn.GetParameter(2)),
        "sigma_err": float(fitfn.GetParError(2)),
        "fitfn": fitfn,
    }


def _plot_pair_mpl(
    axes,
    h_a,
    fit_a: dict | None,
    label_a: str,
    h_b,
    fit_b: dict | None,
    label_b: str,
    ns_per_unit: float,
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    def hist_xy(h):
        cx = np.array([h.GetBinCenter(b) for b in range(1, h.GetNbinsX() + 1)])
        cy = np.array([h.GetBinContent(b) for b in range(1, h.GetNbinsX() + 1)])
        return cx, cy

    def _gauss(x, a, mu, sigma):
        return a * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

    for ax, h, fit, title in (
        (axes[0], h_a, fit_a, label_a),
        (axes[1], h_b, fit_b, label_b),
    ):
        x, y = hist_xy(h)
        w = (x[1] - x[0]) if len(x) > 1 else 0.1
        edges = np.linspace(float(x[0]) - w / 2, float(x[-1]) + w / 2, len(x) + 1)
        ax.stairs(y, edges, color="black", linewidth=2.0, zorder=2)
        if fit is not None:
            xf = np.linspace(float(x.min()), float(x.max()), 300)
            p = fit["fitfn"]
            ax.plot(
                xf,
                _gauss(xf, p.GetParameter(0), p.GetParameter(1), p.GetParameter(2)),
                color="red",
                ls="--",
                lw=2.0,
                zorder=3,
            )
        sig = fit["sigma"] if fit else float(h.GetRMS())
        sig_e = fit["sigma_err"] if fit else 0.0
        tres = sig / SQRT2
        tres_e = sig_e / SQRT2
        ax.set_xlabel(
            r"$\langle t\rangle_{\mathrm{T1}} - \langle t\rangle_{\mathrm{T2}}$ (ns)"
            + (f" [×{ns_per_unit}]" if ns_per_unit != 1 else "")
        )
        ax.set_ylabel("events")
        ax.set_title(title)
        ax.axvline(0.0, color="k", ls="--", lw=0.8, alpha=0.5)
        lines = [
            rf"$N={int(h.GetEntries())}$",
            rf"Gauss $\sigma={sig:.4f}\pm{sig_e:.4f}$ ns",
            rf"timing res. $=\sigma/\sqrt{{2}}={tres:.4f}\pm{tres_e:.4f}$ ns",
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
        description="두 ROOT 의 Δt(패널 3 동일) 비교 + σ/√2 타이밍 레졸루션"
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
    parser.add_argument("--xmin", type=float, default=-15.0)
    parser.add_argument("--xmax", type=float, default=15.0)
    parser.add_argument("--bins", type=int, default=120)
    parser.add_argument("--delta-xmin", type=float, default=None)
    parser.add_argument("--delta-xmax", type=float, default=None)
    parser.add_argument("--delta-bins", type=int, default=None)
    parser.add_argument("--ns-per-unit", type=float, default=1.0)
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

    h_lg, n_lg, sk_lg = _fill_delta_hist(
        path_lg,
        max_events=max_ev_lg,
        d_xmin=d_xmin,
        d_xmax=d_xmax,
        d_bins=d_bins,
        ns_per_unit=args.ns_per_unit,
        hist_name="hDeltaLG",
        hist_title="LG: #LT t#GT_{T1}-#LT t#GT_{T2};ns;events",
    )
    h_nlg, n_nlg, sk_nlg = _fill_delta_hist(
        path_nlg,
        max_events=max_ev_nlg,
        d_xmin=d_xmin,
        d_xmax=d_xmax,
        d_bins=d_bins,
        ns_per_unit=args.ns_per_unit,
        hist_name="hDeltaNoLG",
        hist_title="no LG: #LT t#GT_{T1}-#LT t#GT_{T2};ns;events",
    )

    fit_lg = _fit_gaus_delta(h_lg, ROOT, "lg")
    fit_nlg = _fit_gaus_delta(h_nlg, ROOT, "nlg")

    if fit_lg:
        print(
            "  [LG]     Gauss: "
            f"μ = {fit_lg['mu']:.6f} ± {fit_lg['mu_err']:.6f} ns, "
            f"σ = {fit_lg['sigma']:.6f} ± {fit_lg['sigma_err']:.6f} ns"
        )
    if fit_nlg:
        print(
            "  [no LG]  Gauss: "
            f"μ = {fit_nlg['mu']:.6f} ± {fit_nlg['mu_err']:.6f} ns, "
            f"σ = {fit_nlg['sigma']:.6f} ± {fit_nlg['sigma_err']:.6f} ns"
        )

    stem = "compare_delta_timing_LG_vs_noLG"
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
        ("LG", path_lg, fit_lg, n_lg, sk_lg, n_lg_tree, max_ev_lg),
        ("no_LG", path_nlg, fit_nlg, n_nlg, sk_nlg, n_nlg_tree, max_ev_nlg),
    )
    for tag, path, fit, n_ent, sk, n_tree, n_used in meta:
        row: dict[str, str | float | int] = {
            "sample": tag,
            "path": path,
            "tree_entries_in_file": n_tree,
            "tree_events_used": n_used,
            "matched_N": (not args.no_match_entries),
            "delta_hist_entries": n_ent,
            "skipped_delta": sk,
        }
        if fit:
            row["mu_ns"] = fit["mu"]
            row["mu_err_ns"] = fit["mu_err"]
            row["sigma_ns"] = fit["sigma"]
            row["sigma_err_ns"] = fit["sigma_err"]
            row["timing_resolution_ns"] = fit["sigma"] / SQRT2
            row["timing_resolution_err_ns"] = fit["sigma_err"] / SQRT2
        else:
            row["mu_ns"] = ""
            row["mu_err_ns"] = ""
            row["sigma_ns"] = ""
            row["sigma_err_ns"] = ""
            row["timing_resolution_ns"] = ""
            row["timing_resolution_err_ns"] = ""
        rows.append(row)

    fieldnames = [
        "sample",
        "path",
        "tree_entries_in_file",
        "tree_events_used",
        "matched_N",
        "mu_ns",
        "mu_err_ns",
        "sigma_ns",
        "sigma_err_ns",
        "timing_resolution_ns",
        "timing_resolution_err_ns",
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
        bn_lg = os.path.basename(path_lg)
        bn_nlg = os.path.basename(path_nlg)
        _plot_pair_mpl(
            axes,
            h_lg,
            fit_lg,
            f"LG (Δt)\n{bn_lg}",
            h_nlg,
            fit_nlg,
            f"no LG (Δt)\n{bn_nlg}",
            args.ns_per_unit,
        )
        fig.suptitle(
            r"Per-event $\langle t\rangle_{\mathrm{T1}}-\langle t\rangle_{\mathrm{T2}}$ "
            r"(timing res. $=\sigma/\sqrt{2}$)",
            fontsize=11,
        )
        fig.tight_layout()
        fig.savefig(out_png, dpi=150)
        print(f"PNG 저장: {out_png}")
        plt.close(fig)
    except ImportError as e:
        print(f"matplotlib 없음 — PNG 생략 ({e})", file=sys.stderr)


if __name__ == "__main__":
    main()
