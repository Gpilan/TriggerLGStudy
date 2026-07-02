#!/usr/bin/env python3
"""
이벤트별 타이밍 분포에서 peak의 일정 비율(기본 0.3) 교차점을
상승부 국소 선형 피팅으로 TOA 정의해 T1/T2/Δ(T1-T2) 1x3 히스토그램 생성.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from typing import Any, Optional

from code_archive.plot_trigger_timing import (
    _figures_dir,
    _find_repo_root,
    _load_rootio,
    _merged_from_tower,
    _prepend_build_rootio_ld_path,
    _resolve_input_path,
)


def _auto_axis(vals: list[float], bin_width: float, *, margin_frac: float = 0.08):
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
    nbins = max(1, int(math.ceil((hi - lo) / max(bin_width, 1e-15))))
    nbins = min(nbins, 5000)
    return lo, hi, nbins


def _weighted_linear_fit(x: list[float], y: list[float], w: list[float]):
    sw = sum(w)
    if sw <= 0.0:
        return None
    sx = sum(wi * xi for wi, xi in zip(w, x))
    sy = sum(wi * yi for wi, yi in zip(w, y))
    sxx = sum(wi * xi * xi for wi, xi in zip(w, x))
    sxy = sum(wi * xi * yi for wi, xi, yi in zip(w, x, y))
    den = sw * sxx - sx * sx
    if abs(den) <= 1e-18:
        return None
    a = (sw * sxy - sx * sy) / den
    b = (sy - a * sx) / sw
    return a, b


def _weighted_r2(x: list[float], y: list[float], w: list[float], a: float, b: float) -> float | None:
    sw = sum(w)
    if sw <= 0.0:
        return None
    ybar = sum(wi * yi for wi, yi in zip(w, y)) / sw
    ss_tot = sum(wi * (yi - ybar) ** 2 for wi, yi in zip(w, y))
    if ss_tot <= 1e-18:
        return None
    ss_res = sum(wi * (yi - (a * xi + b)) ** 2 for wi, xi, yi in zip(w, x, y))
    return float(1.0 - ss_res / ss_tot)


def _percentile(vals: list[float], q: float) -> float | None:
    if not vals:
        return None
    arr = sorted(vals)
    if len(arr) == 1:
        return float(arr[0])
    qq = max(0.0, min(1.0, float(q)))
    pos = qq * (len(arr) - 1)
    i0 = int(math.floor(pos))
    i1 = int(math.ceil(pos))
    if i0 == i1:
        return float(arr[i0])
    t = pos - i0
    return float(arr[i0] * (1.0 - t) + arr[i1] * t)


def _toa_local_linear_from_xy(
    x: list[float],
    y: list[float],
    peak_frac: float,
    fit_half_window: int,
):
    def _fallback_interp() -> float | None:
        x0, y0 = x[i_cross - 1], y[i_cross - 1]
        x1, y1 = x[i_cross], y[i_cross]
        if y1 == y0:
            return 0.5 * (x0 + x1)
        t = (thr - y0) / (y1 - y0)
        return x0 + t * (x1 - x0)

    if not y:
        return None, {"method": "no_data", "r2": None, "nfit": 0}
    n = len(y)
    if n < 2 or len(x) != len(y):
        return None, {"method": "too_few_bins", "r2": None, "nfit": n}

    i_peak = max(range(n), key=lambda i: y[i])
    peak = y[i_peak]
    if peak <= 0.0:
        return None, {"method": "nonpositive_peak", "r2": None, "nfit": 0}
    thr = float(peak_frac) * peak

    i_cross = None
    for i in range(1, i_peak + 1):
        if y[i - 1] <= thr <= y[i]:
            i_cross = i
            break
    if i_cross is None:
        return None, {"method": "no_crossing", "r2": None, "nfit": 0}

    i0 = max(0, i_cross - fit_half_window)
    i1 = min(i_peak, i_cross + fit_half_window)
    xs = x[i0 : i1 + 1]
    ys = y[i0 : i1 + 1]
    if len(xs) < 2:
        return _fallback_interp(), {"method": "fallback_short_window", "r2": None, "nfit": len(xs)}

    ws = [1.0 / max(ys_i**0.5, 1.0) for ys_i in ys]
    fit = _weighted_linear_fit(xs, ys, ws)
    if fit is None:
        return _fallback_interp(), {"method": "fallback_singular_fit", "r2": None, "nfit": len(xs)}

    a, b = fit
    r2 = _weighted_r2(xs, ys, ws, a, b)
    if a <= 1e-15:
        return _fallback_interp(), {"method": "fallback_nonpositive_slope", "r2": r2, "nfit": len(xs)}

    t_cross = (thr - b) / a
    xmin = min(xs)
    xmax = max(xs)
    if t_cross < xmin or t_cross > xmax:
        return _fallback_interp(), {"method": "fallback_out_of_window", "r2": r2, "nfit": len(xs)}

    yfit = [a * xx + b for xx in xs]
    chi2 = float(sum(wi * (yy - ff) ** 2 for wi, yy, ff in zip(ws, ys, yfit)))
    ndf = int(max(len(xs) - 2, 0))
    pval = float(math.exp(-0.5 * chi2)) if ndf > 0 else 1.0
    return float(t_cross), {
        "method": "local_linear",
        "r2": r2,
        "nfit": len(xs),
        "chi2": chi2,
        "ndf": ndf,
        "pvalue": pval,
    }


def _gamma_response(shape: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    return float((x / shape) ** shape * math.exp(shape - x))


def _gamma_root(shape: float, level: float, low: float, high: float, rising: bool) -> float:
    lo = float(low)
    hi = float(high)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if (_gamma_response(shape, mid) < level) == rising:
            lo = mid
        else:
            hi = mid
    return float(0.5 * (lo + hi))


def _gamma_dimensionless_rise(shape: float) -> float:
    t10 = _gamma_root(shape, 0.10, 0.0, shape, True)
    t90 = _gamma_root(shape, 0.90, 0.0, shape, True)
    return float(max(t90 - t10, 1e-9))


def _gamma_tau_from_rise(response_rise_ns: float, gamma_shape: float) -> float:
    return float(response_rise_ns) / _gamma_dimensionless_rise(float(gamma_shape))


def _load_r2076_time_response(spec_csv: str) -> tuple[Optional[float], Optional[float]]:
    if not os.path.isfile(spec_csv):
        return None, None
    rise_ns: Optional[float] = None
    transit_ns: Optional[float] = None
    with open(spec_csv, "r", encoding="utf-8-sig", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if str(row.get("category", "")).strip() != "time_response":
                continue
            parameter = str(row.get("parameter", "")).strip()
            qualifier = str(row.get("qualifier", "")).strip()
            val = str(row.get("value", "")).strip()
            try:
                fv = float(val)
            except ValueError:
                continue
            if parameter == "rise_time" and qualifier in ("typ", ""):
                rise_ns = fv
            elif parameter == "transit_time" and qualifier in ("typ", ""):
                transit_ns = fv
    return rise_ns, transit_ns


def _build_reconstructed_waveform(
    lo: list[float],
    hi: list[float],
    cnt: list[int],
    *,
    sample_step_ns: float,
    response_rise_ns: float,
    response_fwhm_ns: float,
    gamma_shape: float,
    gamma_tau_ns: Optional[float],
    transit_time_ns: float,
):
    centers = [0.5 * (float(loi) + float(hii)) for loi, hii in zip(lo, hi)]
    valid = [(t, int(c)) for t, c in zip(centers, cnt) if int(c) > 0]
    if not valid:
        return [], []

    first_hit = min(t for t, _ in valid)
    last_hit = max(t for t, _ in valid)
    start = math.floor((first_hit - 2.0 * sample_step_ns) / sample_step_ns) * sample_step_ns
    stop = last_hit + 12.0 * response_fwhm_ns
    n_samples = int(math.ceil((stop - start) / sample_step_ns)) + 1
    n_samples = max(n_samples, 2)

    impulse = [0.0] * n_samples
    for hit_time, count in valid:
        pos = (float(hit_time) - start) / sample_step_ns
        i0 = int(math.floor(pos))
        frac = float(pos - i0)
        amp = float(count)
        if 0 <= i0 < n_samples:
            impulse[i0] += amp * (1.0 - frac)
        i1 = i0 + 1
        if 0 <= i1 < n_samples:
            impulse[i1] += amp * frac

    time = [start + i * sample_step_ns + float(transit_time_ns) for i in range(n_samples)]
    tau_ns = float(gamma_tau_ns) if gamma_tau_ns is not None else _gamma_tau_from_rise(response_rise_ns, gamma_shape)
    kernel_samples = int(math.ceil(12.0 * response_fwhm_ns / sample_step_ns)) + 1
    kernel = [0.0] * kernel_samples
    for k in range(1, kernel_samples):
        x = (float(k) * sample_step_ns) / tau_ns
        kernel[k] = float((x / gamma_shape) ** gamma_shape * math.exp(gamma_shape - x)) if x > 0.0 else 0.0

    voltage = [0.0] * n_samples
    for i, amp in enumerate(impulse):
        if amp <= 0.0:
            continue
        available = min(kernel_samples, n_samples - i)
        for k in range(1, available):
            voltage[i + k] += amp * kernel[k]
    return time, voltage


def _collect_toa(
    path: str,
    *,
    max_events: int,
    peak_frac: float,
    fit_half_window: int,
    fit_r2_min: float,
    sample_step_ns: float,
    response_rise_ns: float,
    response_fwhm_ns: float,
    transit_time_ns: float,
    gamma_shape: float,
    gamma_tau_ns: Optional[float],
):
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
    dt_vals: list[float] = []
    sk1 = sk2 = skd = 0
    fit_stats: dict[str, Any] = {
        "events": nmax,
        "t1_linear": 0,
        "t1_fallback": 0,
        "t1_no_toa": 0,
        "t1_bad_r2": 0,
        "t1_r2": [],
        "t1_chi2ndf": [],
        "t1_pvalue": [],
        "t2_linear": 0,
        "t2_fallback": 0,
        "t2_no_toa": 0,
        "t2_bad_r2": 0,
        "t2_r2": [],
        "t2_chi2ndf": [],
        "t2_pvalue": [],
    }

    for i in range(nmax):
        tree.GetEntry(i)
        lo1, hi1, c1 = _merged_from_tower(evt, 0)
        lo2, hi2, c2 = _merged_from_tower(evt, 1)
        wx1, wy1 = _build_reconstructed_waveform(
            lo1,
            hi1,
            c1,
            sample_step_ns=sample_step_ns,
            response_rise_ns=response_rise_ns,
            response_fwhm_ns=response_fwhm_ns,
            gamma_shape=gamma_shape,
            gamma_tau_ns=gamma_tau_ns,
            transit_time_ns=transit_time_ns,
        )
        wx2, wy2 = _build_reconstructed_waveform(
            lo2,
            hi2,
            c2,
            sample_step_ns=sample_step_ns,
            response_rise_ns=response_rise_ns,
            response_fwhm_ns=response_fwhm_ns,
            gamma_shape=gamma_shape,
            gamma_tau_ns=gamma_tau_ns,
            transit_time_ns=transit_time_ns,
        )

        v1, m1 = _toa_local_linear_from_xy(wx1, wy1, peak_frac, fit_half_window)
        v2, m2 = _toa_local_linear_from_xy(wx2, wy2, peak_frac, fit_half_window)

        for prefix, v, meta in (("t1", v1, m1), ("t2", v2, m2)):
            method = str(meta.get("method", "unknown"))
            r2 = meta.get("r2", None)
            if method == "local_linear":
                fit_stats[f"{prefix}_linear"] += 1
            elif v is not None:
                fit_stats[f"{prefix}_fallback"] += 1
            else:
                fit_stats[f"{prefix}_no_toa"] += 1

            if r2 is not None:
                rr = float(r2)
                fit_stats[f"{prefix}_r2"].append(rr)
                if rr < fit_r2_min:
                    fit_stats[f"{prefix}_bad_r2"] += 1
            chi2 = meta.get("chi2", None)
            ndf = int(meta.get("ndf", 0) or 0)
            if chi2 is not None and ndf > 0:
                fit_stats[f"{prefix}_chi2ndf"].append(float(chi2) / float(ndf))
            pval = meta.get("pvalue", None)
            if pval is not None:
                fit_stats[f"{prefix}_pvalue"].append(float(pval))

        if v1 is not None:
            t1_vals.append(v1)
        else:
            sk1 += 1
        if v2 is not None:
            t2_vals.append(v2)
        else:
            sk2 += 1
        if v1 is not None and v2 is not None:
            dt_vals.append(float(v1 - v2))
        else:
            skd += 1

    f.Close()
    return t1_vals, t2_vals, dt_vals, sk1, sk2, skd, fit_stats


def _format_fit_quality(fit_stats: dict[str, Any], fit_r2_min: float) -> str:
    lines = [f"  Fit quality check (waveform local-linear, R^2 threshold={fit_r2_min:.3f})"]
    for trig in ("t1", "t2"):
        linear = int(fit_stats[f"{trig}_linear"])
        fallback = int(fit_stats[f"{trig}_fallback"])
        no_toa = int(fit_stats[f"{trig}_no_toa"])
        r2_vals = [float(v) for v in fit_stats[f"{trig}_r2"]]
        bad_r2 = int(fit_stats[f"{trig}_bad_r2"])
        total_valid = linear + fallback
        fallback_pct = (100.0 * fallback / total_valid) if total_valid > 0 else 0.0

        tag = "T1" if trig == "t1" else "T2"
        lines.append(
            f"    {tag}: linear={linear}, fallback={fallback} ({fallback_pct:.1f}%), no_toa={no_toa}"
        )
        if r2_vals:
            mean_r2 = sum(r2_vals) / len(r2_vals)
            p05 = _percentile(r2_vals, 0.05)
            p50 = _percentile(r2_vals, 0.50)
            p95 = _percentile(r2_vals, 0.95)
            lines.append(
                f"        R^2: mean={mean_r2:.4f}, p05={p05:.4f}, p50={p50:.4f}, p95={p95:.4f}, below_thr={bad_r2}"
            )
            chi2ndf_vals = [float(v) for v in fit_stats[f"{trig}_chi2ndf"]]
            pval_vals = [float(v) for v in fit_stats[f"{trig}_pvalue"]]
            if chi2ndf_vals:
                c50 = _percentile(chi2ndf_vals, 0.50)
                c95 = _percentile(chi2ndf_vals, 0.95)
                lines.append(f"        chi2/ndf: p50={c50:.4f}, p95={c95:.4f}")
            if pval_vals:
                p50v = _percentile(pval_vals, 0.50)
                p05v = _percentile(pval_vals, 0.05)
                lines.append(f"        p-value: p50={p50v:.4f}, p05={p05v:.4f}")
        else:
            lines.append("        R^2: n/a (local-linear fit not used)")
    return "\n".join(lines)


def _apply_style(ROOT):
    s = ROOT.gStyle
    s.SetOptStat(111111)
    s.SetOptFit(1111)
    s.SetOptTitle(1)
    s.SetHistLineColor(ROOT.kBlack)
    s.SetHistLineWidth(2)
    s.SetTitleFont(42, "XYZ")
    s.SetLabelFont(42, "XYZ")
    s.SetStatFont(42)


def _make_hist(ROOT, name: str, vals: list[float], lo: float, hi: float, nb: int, title: str):
    h = ROOT.TH1F(name, title, nb, lo, hi)
    h.Sumw2()
    h.SetDirectory(0)
    h.SetLineColor(ROOT.kBlack)
    h.SetLineWidth(2)
    for v in vals:
        h.Fill(float(v))
    return h


def _fit_gaus(ROOT, h, fit_name: str):
    if h.GetEntries() <= 0:
        return None
    xmin = float(h.GetXaxis().GetXmin())
    xmax = float(h.GetXaxis().GetXmax())
    fitfn = ROOT.TF1(fit_name, "gaus", xmin, xmax)
    # Increase sampling density so the fit curve appears visually smooth.
    fitfn.SetNpx(2000)
    fitfn.SetParameter(0, max(float(h.GetMaximum()), 1e-9))
    fitfn.SetParameter(1, float(h.GetMean()))
    fitfn.SetParameter(2, max(float(h.GetRMS()), 1e-6))
    status = int(h.Fit(fitfn, "QNR", "", xmin, xmax))
    if status != 0:
        return None
    fitfn.SetLineColor(ROOT.kRed)
    fitfn.SetLineWidth(2)
    fitfn.SetLineStyle(1)
    sigma = float(fitfn.GetParameter(2))
    sigma_err = float(fitfn.GetParError(2))
    chi2 = float(fitfn.GetChisquare())
    ndf = int(fitfn.GetNDF())
    pval = float(fitfn.GetProb())
    return fitfn, sigma, sigma_err, chi2, ndf, pval


def _draw_panel(ROOT, h, vals, fit_name: str, bin_ps: float, show_no_data_x: float = 0.40):
    h.Draw("HIST")
    if not vals:
        t = ROOT.TLatex()
        t.SetNDC()
        t.SetTextFont(42)
        t.SetTextSize(0.04)
        t.DrawLatex(show_no_data_x, 0.50, "no data")
        return
    fit = _fit_gaus(ROOT, h, fit_name)
    if fit is not None:
        fitfn, sigma, sigma_err, chi2, ndf, pval = fit
        fitfn.Draw("same")
        tx = ROOT.TLatex()
        tx.SetNDC()
        tx.SetTextFont(42)
        tx.SetTextSize(0.042)
        tx.SetTextColor(ROOT.kRed + 1)
        tx.DrawLatex(0.14, 0.86, f"#sigma = {sigma:.4f} #pm {sigma_err:.4f} ns")
        if ndf > 0:
            tx.DrawLatex(0.14, 0.80, f"#chi^{{2}}/ndf = {chi2:.2f}/{ndf} = {chi2/ndf:.3f}")
        else:
            tx.DrawLatex(0.14, 0.80, f"#chi^{{2}}/ndf = {chi2:.2f}/0")
        tx.DrawLatex(0.14, 0.74, f"p-value = {pval:.3g}")
    bw = ROOT.TLatex()
    bw.SetNDC()
    bw.SetTextFont(42)
    bw.SetTextSize(0.028)
    bw.DrawLatex(0.14, 0.08, f"hist bin = {bin_ps:.1f} ps")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Waveform reconstruction + local-linear CFD(peak fraction) 기반 T1/T2/Δ 3패널"
    )
    parser.add_argument("input", nargs="?", default=None)
    parser.add_argument("-o", "--output", default=None)
    parser.add_argument("--max-events", type=int, default=3000)
    parser.add_argument("--peak-frac", type=float, default=0.3)
    parser.add_argument(
        "--rebin-ps",
        type=float,
        default=None,
        help="deprecated: waveform CFD 경로에서는 사용하지 않음",
    )
    parser.add_argument("--plot-bin-ps", type=float, default=1.0)
    parser.add_argument("--fit-half-window", type=int, default=2, help="교차점 주변 반쪽 윈도우(bin)")
    parser.add_argument(
        "--fit-r2-min",
        type=float,
        default=0.90,
        help="피팅 품질 경고 기준 R^2 (요약 통계용, 기본 0.90)",
    )
    parser.add_argument("--fixed-range", action="store_true")
    parser.add_argument("--xmin", type=float, default=-5.0)
    parser.add_argument("--xmax", type=float, default=5.0)
    parser.add_argument(
        "--pmt-spec-csv",
        default="analysis/reference/pmt_r2076/R2076_pmt_spec_reference.csv",
        help="PMT 스펙 CSV 경로 (rise/transit typ 자동 로드)",
    )
    parser.add_argument("--response-rise-ns", type=float, default=None, help="파형 10-90 rise(ns), 미지정시 spec typ")
    parser.add_argument("--response-fwhm-ns", type=float, default=None, help="파형 FWHM(ns), 미지정시 3.0ns")
    parser.add_argument("--transit-time-ns", type=float, default=None, help="파형 transit(ns), 미지정시 spec typ")
    parser.add_argument(
        "--sample-step-ns",
        type=float,
        default=0.05,
        help="waveform 샘플 간격(ns). 작게 줄수록 파형이 더 부드러워짐",
    )
    parser.add_argument("--gamma-shape", type=float, default=1.915604733026, help="gamma pulse shape")
    parser.add_argument(
        "--gamma-tau-ns",
        type=float,
        default=None,
        help="gamma pulse tau(ns), 미지정시 rise+shape로 자동 보정",
    )
    parser.add_argument("-l", "--rootio-lib", default=None)
    args = parser.parse_args()

    if not (0.0 < args.peak_frac < 1.0):
        raise SystemExit("--peak-frac는 0~1 사이여야 합니다.")
    if args.fit_half_window < 1:
        raise SystemExit("--fit-half-window는 1 이상이어야 합니다.")
    if args.plot_bin_ps <= 0.0:
        raise SystemExit("--plot-bin-ps는 0보다 커야 합니다.")
    if args.sample_step_ns <= 0.0 or args.gamma_shape <= 0.0:
        raise SystemExit("--sample-step-ns, --gamma-shape는 0보다 커야 합니다.")
    if args.response_rise_ns is not None and args.response_rise_ns <= 0.0:
        raise SystemExit("--response-rise-ns는 0보다 커야 합니다.")
    if args.response_fwhm_ns is not None and args.response_fwhm_ns <= 0.0:
        raise SystemExit("--response-fwhm-ns는 0보다 커야 합니다.")
    if args.transit_time_ns is not None and args.transit_time_ns < 0.0:
        raise SystemExit("--transit-time-ns는 0 이상이어야 합니다.")
    if args.gamma_tau_ns is not None and args.gamma_tau_ns <= 0.0:
        raise SystemExit("--gamma-tau-ns는 0보다 커야 합니다.")
    if args.rebin_ps is not None:
        print(f"[note] --rebin-ps={args.rebin_ps:g} 는 waveform CFD 모드에서 무시됩니다.")

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

    in_path = _resolve(args.input, "analysis/t_res/data/v2_60GeV_e-_LG_3000_0.root")
    spec_csv_path = _resolve(args.pmt_spec_csv, "analysis/reference/pmt_r2076/R2076_pmt_spec_reference.csv")
    spec_rise_ns, spec_transit_ns = _load_r2076_time_response(spec_csv_path)
    response_rise_ns = float(args.response_rise_ns) if args.response_rise_ns is not None else float(spec_rise_ns if spec_rise_ns is not None else 1.0)
    response_fwhm_ns = float(args.response_fwhm_ns) if args.response_fwhm_ns is not None else 3.0
    transit_time_ns = float(args.transit_time_ns) if args.transit_time_ns is not None else float(spec_transit_ns if spec_transit_ns is not None else 0.0)

    _load_rootio(args.rootio_lib, repo_root)
    import ROOT

    ROOT.gROOT.SetBatch(True)
    _apply_style(ROOT)

    t1_vals, t2_vals, dt_vals, sk1, sk2, skd, fit_stats = _collect_toa(
        in_path,
        max_events=args.max_events,
        peak_frac=args.peak_frac,
        fit_half_window=args.fit_half_window,
        fit_r2_min=args.fit_r2_min,
        sample_step_ns=args.sample_step_ns,
        response_rise_ns=response_rise_ns,
        response_fwhm_ns=response_fwhm_ns,
        transit_time_ns=transit_time_ns,
        gamma_shape=args.gamma_shape,
        gamma_tau_ns=args.gamma_tau_ns,
    )

    plot_bin_ns = args.plot_bin_ps / 1000.0
    if args.fixed_range and args.xmax <= args.xmin:
        raise SystemExit("--xmax must be > --xmin")
    if args.fixed_range:
        nb_fixed = max(1, int(math.ceil((args.xmax - args.xmin) / max(plot_bin_ns, 1e-15))))
        lo1, hi1, nb1 = args.xmin, args.xmax, nb_fixed
        lo2, hi2, nb2 = args.xmin, args.xmax, nb_fixed
        lo3, hi3, nb3 = args.xmin, args.xmax, nb_fixed
    else:
        lo1, hi1, nb1 = _auto_axis(t1_vals, plot_bin_ns)
        lo2, hi2, nb2 = _auto_axis(t2_vals, plot_bin_ns)
        lo3, hi3, nb3 = _auto_axis(dt_vals, plot_bin_ns)

    h1 = _make_hist(ROOT, "hToaT1", t1_vals, lo1, hi1, nb1, "T1;TOA waveform CFD (ns);events")
    h2 = _make_hist(ROOT, "hToaT2", t2_vals, lo2, hi2, nb2, "T2;TOA waveform CFD (ns);events")
    h3 = _make_hist(ROOT, "hToaDiffT1T2", dt_vals, lo3, hi3, nb3, "T1 - T2;#DeltaTOA waveform CFD (ns);events")

    stem = os.path.splitext(os.path.basename(in_path))[0]
    if args.output:
        out = os.path.abspath(args.output)
        if not os.path.splitext(out)[1]:
            out += ".png"
    else:
        pbin = f"{args.plot_bin_ps:g}".replace(".", "p")
        out = os.path.join(_figures_dir(repo_root), f"cfd_waveform_recon_per_trigger_{stem}_1x3_to{pbin}ps.png")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    c = ROOT.TCanvas("cToaPerTrigWaveform", "waveform recon CFD", 1860, 560)
    c.Divide(3, 1)
    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)
    _draw_panel(ROOT, h1, t1_vals, "fit_t1", ((hi1 - lo1) / float(nb1)) * 1000.0)
    c.cd(2)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)
    _draw_panel(ROOT, h2, t2_vals, "fit_t2", ((hi2 - lo2) / float(nb2)) * 1000.0)
    c.cd(3)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)
    _draw_panel(ROOT, h3, dt_vals, "fit_t1_minus_t2", ((hi3 - lo3) / float(nb3)) * 1000.0, show_no_data_x=0.35)
    c.SaveAs(out)
    c.Close()

    print(
        "[PMT spec] "
        f"csv={spec_csv_path}, "
        f"rise(ns)={'cli' if args.response_rise_ns is not None else ('spec' if spec_rise_ns is not None else 'default')}:{response_rise_ns:.3f}, "
        f"transit(ns)={'cli' if args.transit_time_ns is not None else ('spec' if spec_transit_ns is not None else 'default')}:{transit_time_ns:.3f}, "
        f"fwhm(ns)={'cli' if args.response_fwhm_ns is not None else 'default'}:{response_fwhm_ns:.3f}\n"
        f"입력: {in_path}\n"
        f"  T1: filled={len(t1_vals)}, skipped={sk1}\n"
        f"  T2: filled={len(t2_vals)}, skipped={sk2}\n"
        f"  T1-T2: filled={len(dt_vals)}, skipped={skd}\n"
        f"  TOA 정의: waveform reconstruction + local-linear rising edge at {args.peak_frac:g} * peak\n"
        f"  waveform: dt={args.sample_step_ns:.3f} ns, rise={response_rise_ns:.3f} ns, fwhm={response_fwhm_ns:.3f} ns, transit={transit_time_ns:.3f} ns, shape={args.gamma_shape:.6f}, tau={float(args.gamma_tau_ns) if args.gamma_tau_ns is not None else _gamma_tau_from_rise(response_rise_ns, args.gamma_shape):.6f} ns\n"
        f"PNG 저장: {out}\n"
        f"{_format_fit_quality(fit_stats, args.fit_r2_min)}"
    )


if __name__ == "__main__":
    main()
