#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass

import matplotlib.pyplot as plt

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
CODE_ARCHIVE_DIR = os.path.join(THIS_DIR, "code_archive")
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)
if CODE_ARCHIVE_DIR not in sys.path:
    sys.path.insert(0, CODE_ARCHIVE_DIR)

from plot_cfd_per_trigger_locallinear_peak import (
    _auto_axis,
    _collect_toa,
    _fit_gaus,
    _make_hist,
)
from plot_event_optical_mpv_matrix_root import (
    _build_nph_hist_and_landau_fit,
    _scan_events,
)
from code_archive.plot_trigger_timing import (
    _figures_dir,
    _find_repo_root,
    _load_rootio,
    _prepend_build_rootio_ld_path,
)


@dataclass
class TrendPoint:
    size: float
    root_path: str
    mpv_t1: float | None
    mpv_t2: float | None
    mpv_mean: float | None
    sigma_ns: float | None
    sigma_ps: float | None
    n_dt: int


def _default_inputs(repo_root: str) -> dict[float, str]:
    return {
        1.0: os.path.join(repo_root, "analysis/t_res/data/v4_60GeV_e-_LG_1x1_3000_0.root"),
        1.5: os.path.join(repo_root, "analysis/t_res/data/v4_60GeV_e-_LG_1p5x1p5_3000_0.root"),
        2.0: os.path.join(repo_root, "analysis/t_res/data/v4_60GeV_e-_LG_2x2_3000_0.root"),
        2.5: os.path.join(repo_root, "analysis/t_res/data/v4_60GeV_e-_LG_2p5x2p5_3000_0.root"),
        3.0: os.path.join(repo_root, "analysis/t_res/data/v4_60GeV_e-_LG_3x3_3000_0.root"),
        3.5: os.path.join(repo_root, "analysis/t_res/data/v4_60GeV_e-_LG_3p5x3p5_3000_0.root"),
        4.0: os.path.join(repo_root, "analysis/t_res/data/v3_60GeV_e-_LG_3000_0.root"),
    }


def _fit_mpv(ROOT, root_path: str, max_events: int) -> tuple[float | None, float | None, float | None]:
    ev_t1 = _scan_events(root_path, trig=0, max_events=max_events)
    ev_t2 = _scan_events(root_path, trig=1, max_events=max_events)
    vals_t1 = [ev.nphoton for ev in ev_t1]
    vals_t2 = [ev.nphoton for ev in ev_t2]
    h1, f1, mpv1 = _build_nph_hist_and_landau_fit(ROOT, "h_mpv_t1_tmp", vals_t1, "tmp")
    h2, f2, mpv2 = _build_nph_hist_and_landau_fit(ROOT, "h_mpv_t2_tmp", vals_t2, "tmp")
    _ = h1, h2, f1, f2
    mpv_mean = None
    if mpv1 is not None and mpv2 is not None:
        mpv_mean = 0.5 * (float(mpv1) + float(mpv2))
    elif mpv1 is not None:
        mpv_mean = float(mpv1)
    elif mpv2 is not None:
        mpv_mean = float(mpv2)
    return mpv1, mpv2, mpv_mean


def _fit_timing_sigma_ns(
    ROOT,
    root_path: str,
    max_events: int,
    peak_frac: float,
    fit_half_window: int,
    fit_r2_min: float,
    plot_bin_ps: float,
    sample_step_ns: float,
    response_rise_ns: float,
    response_fwhm_ns: float,
    transit_time_ns: float,
    gamma_shape: float,
    gamma_tau_ns: float | None,
) -> tuple[float | None, int]:
    _, _, dt_vals, _, _, _, _ = _collect_toa(
        root_path,
        max_events=max_events,
        peak_frac=peak_frac,
        fit_half_window=fit_half_window,
        fit_r2_min=fit_r2_min,
        sample_step_ns=sample_step_ns,
        response_rise_ns=response_rise_ns,
        response_fwhm_ns=response_fwhm_ns,
        transit_time_ns=transit_time_ns,
        gamma_shape=gamma_shape,
        gamma_tau_ns=gamma_tau_ns,
    )
    if not dt_vals:
        return None, 0
    bin_ns = float(plot_bin_ps) / 1000.0
    lo, hi, nb = _auto_axis(dt_vals, bin_ns)
    h = _make_hist(ROOT, "h_dt_tmp", dt_vals, lo, hi, nb, "tmp")
    fit = _fit_gaus(ROOT, h, "fit_dt_tmp")
    sigma = abs(float(fit[1])) if fit is not None else None
    return sigma, len(dt_vals)


def main() -> None:
    p = argparse.ArgumentParser(description="LG size trend (MPV + timing resolution) for v4 data")
    p.add_argument("--max-events", type=int, default=3000)
    p.add_argument("--peak-frac", type=float, default=0.3)
    p.add_argument("--fit-half-window", type=int, default=2)
    p.add_argument("--fit-r2-min", type=float, default=0.90)
    p.add_argument("--plot-bin-ps", type=float, default=5.0)
    p.add_argument("--sample-step-ns", type=float, default=0.05)
    p.add_argument("--response-rise-ns", type=float, default=1.3)
    p.add_argument("--response-fwhm-ns", type=float, default=3.0)
    p.add_argument("--transit-time-ns", type=float, default=14.0)
    p.add_argument("--gamma-shape", type=float, default=1.915604733026)
    p.add_argument("--gamma-tau-ns", type=float, default=None)
    p.add_argument("-o", "--output", default=None, help="output PNG path")
    p.add_argument("--csv", default=None, help="output CSV path")
    p.add_argument("-l", "--rootio-lib", default=None)
    args = p.parse_args()

    repo_root = _find_repo_root() or os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    _prepend_build_rootio_ld_path(repo_root)
    _load_rootio(args.rootio_lib, repo_root)
    import ROOT

    ROOT.gROOT.SetBatch(True)

    files = _default_inputs(repo_root)
    points: list[TrendPoint] = []
    for size in sorted(files):
        root_path = files[size]
        if not os.path.isfile(root_path):
            continue
        mpv_t1, mpv_t2, mpv_mean = _fit_mpv(ROOT, root_path, args.max_events)
        sigma_ns, n_dt = _fit_timing_sigma_ns(
            ROOT,
            root_path,
            args.max_events,
            args.peak_frac,
            args.fit_half_window,
            args.fit_r2_min,
            args.plot_bin_ps,
            args.sample_step_ns,
            args.response_rise_ns,
            args.response_fwhm_ns,
            args.transit_time_ns,
            args.gamma_shape,
            args.gamma_tau_ns,
        )
        sigma_ps = None if sigma_ns is None else sigma_ns * 1000.0
        points.append(
            TrendPoint(
                size=size,
                root_path=root_path,
                mpv_t1=mpv_t1,
                mpv_t2=mpv_t2,
                mpv_mean=mpv_mean,
                sigma_ns=sigma_ns,
                sigma_ps=sigma_ps,
                n_dt=n_dt,
            )
        )

    if not points:
        raise SystemExit("유효한 입력 ROOT가 없습니다.")

    out_dir = os.path.join(_figures_dir(repo_root), "version4_lg_scan_trend")
    os.makedirs(out_dir, exist_ok=True)
    out_png = (
        os.path.abspath(args.output)
        if args.output
        else os.path.join(out_dir, "lg_size_trend_mpv_sigma_v4plusv3_from1_to4_xmax4.png")
    )
    out_csv = (
        os.path.abspath(args.csv)
        if args.csv
        else os.path.join(out_dir, "lg_size_trend_mpv_sigma_v4plusv3_from1_to4_xmax4.csv")
    )

    xs = [pt.size for pt in points]
    y_mpv_t1 = [pt.mpv_t1 for pt in points]
    y_mpv_t2 = [pt.mpv_t2 for pt in points]
    y_mpv_mean = [pt.mpv_mean for pt in points]
    y_sigma_ps = [pt.sigma_ps for pt in points]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.8, 7.2), sharex=True)
    ax1.plot(xs, y_mpv_t1, marker="o", linewidth=1.5, label="T1 MPV")
    ax1.plot(xs, y_mpv_t2, marker="s", linewidth=1.5, label="T2 MPV")
    ax1.plot(xs, y_mpv_mean, marker="D", linewidth=2.0, label="Mean MPV")
    ax1.set_ylabel("Optical Photon MPV")
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="best")
    ax1.set_title("Trigger Size Trend")

    ax2.plot(xs, y_sigma_ps, marker="o", linewidth=2.0, color="tab:red", label="Timing resolution sigma")
    ax2.set_ylabel("Timing Resolution (ps)")
    ax2.set_xlabel("Trigger overlap size (x = y, cm)")
    ax2.grid(True, alpha=0.25)
    ax2.legend(loc="best")

    ax2.set_xlim(0.0, 5.0)
    ax2.set_xticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    reserved_start = max(xs) if xs else 0.0
    if reserved_start < 5.0:
        ax2.axvspan(reserved_start, 5.0, color="gray", alpha=0.06)
        ax2.text(reserved_start + 0.05, ax2.get_ylim()[1] * 0.92, "reserved (no data yet)", fontsize=9, alpha=0.7)

    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["size_cm", "root_path", "mpv_t1", "mpv_t2", "mpv_mean", "sigma_ns", "sigma_ps", "n_dt"])
        for pt in points:
            w.writerow(
                [
                    f"{pt.size:.2f}",
                    pt.root_path,
                    "" if pt.mpv_t1 is None else f"{pt.mpv_t1:.6f}",
                    "" if pt.mpv_t2 is None else f"{pt.mpv_t2:.6f}",
                    "" if pt.mpv_mean is None else f"{pt.mpv_mean:.6f}",
                    "" if pt.sigma_ns is None else f"{pt.sigma_ns:.6f}",
                    "" if pt.sigma_ps is None else f"{pt.sigma_ps:.3f}",
                    pt.n_dt,
                ]
            )

    print(f"PNG 저장: {out_png}")
    print(f"CSV 저장: {out_csv}")
    for pt in points:
        print(
            f"size={pt.size:.1f} cm | MPV(T1/T2/mean)=({pt.mpv_t1:.2f},{pt.mpv_t2:.2f},{pt.mpv_mean:.2f}) "
            f"| sigma={pt.sigma_ps:.2f} ps | n_dt={pt.n_dt}"
        )


if __name__ == "__main__":
    main()
