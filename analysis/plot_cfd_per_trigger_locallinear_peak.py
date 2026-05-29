#!/usr/bin/env python3
"""
이벤트별 타이밍 분포에서 peak의 일정 비율(기본 0.3) 교차점을
상승부 국소 선형 피팅으로 TOA 정의해 T1/T2/Δ(T1-T2) 1x3 히스토그램 생성.
"""

from __future__ import annotations

import argparse
import math
import os

from compare_timing_resolution import _rebin_merged_to_width
from plot_trigger_timing import (
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


def _toa_local_linear_from_peak_fraction(
    lo: list[float],
    hi: list[float],
    cnt: list[int],
    peak_frac: float,
    fit_half_window: int,
) -> float | None:
    if not cnt:
        return None
    n = len(cnt)
    if n < 2:
        return None

    y = [float(v) for v in cnt]
    x = [0.5 * (float(lo[i]) + float(hi[i])) for i in range(n)]

    i_peak = max(range(n), key=lambda i: y[i])
    peak = y[i_peak]
    if peak <= 0.0:
        return None
    thr = float(peak_frac) * peak

    i_cross = None
    for i in range(1, i_peak + 1):
        if y[i - 1] <= thr <= y[i]:
            i_cross = i
            break
    if i_cross is None:
        return None

    i0 = max(0, i_cross - fit_half_window)
    i1 = min(i_peak, i_cross + fit_half_window)
    xs = x[i0 : i1 + 1]
    ys = y[i0 : i1 + 1]
    if len(xs) < 2:
        x0, y0 = x[i_cross - 1], y[i_cross - 1]
        x1, y1 = x[i_cross], y[i_cross]
        if y1 == y0:
            return 0.5 * (x0 + x1)
        t = (thr - y0) / (y1 - y0)
        return x0 + t * (x1 - x0)

    ws = [1.0 / max(ys_i**0.5, 1.0) for ys_i in ys]
    fit = _weighted_linear_fit(xs, ys, ws)
    if fit is None:
        x0, y0 = x[i_cross - 1], y[i_cross - 1]
        x1, y1 = x[i_cross], y[i_cross]
        if y1 == y0:
            return 0.5 * (x0 + x1)
        t = (thr - y0) / (y1 - y0)
        return x0 + t * (x1 - x0)

    a, b = fit
    if a <= 1e-15:
        x0, y0 = x[i_cross - 1], y[i_cross - 1]
        x1, y1 = x[i_cross], y[i_cross]
        if y1 == y0:
            return 0.5 * (x0 + x1)
        t = (thr - y0) / (y1 - y0)
        return x0 + t * (x1 - x0)

    t_cross = (thr - b) / a
    xmin = min(xs)
    xmax = max(xs)
    if t_cross < xmin or t_cross > xmax:
        x0, y0 = x[i_cross - 1], y[i_cross - 1]
        x1, y1 = x[i_cross], y[i_cross]
        if y1 == y0:
            return 0.5 * (x0 + x1)
        t = (thr - y0) / (y1 - y0)
        return x0 + t * (x1 - x0)
    return float(t_cross)


def _collect_toa(
    path: str,
    *,
    max_events: int,
    peak_frac: float,
    rebin_ns: float,
    fit_half_window: int,
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

    for i in range(nmax):
        tree.GetEntry(i)
        lo1, hi1, c1 = _merged_from_tower(evt, 0)
        lo2, hi2, c2 = _merged_from_tower(evt, 1)
        lo1, hi1, c1 = _rebin_merged_to_width(lo1, hi1, c1, rebin_ns)
        lo2, hi2, c2 = _rebin_merged_to_width(lo2, hi2, c2, rebin_ns)

        v1 = _toa_local_linear_from_peak_fraction(lo1, hi1, c1, peak_frac, fit_half_window)
        v2 = _toa_local_linear_from_peak_fraction(lo2, hi2, c2, peak_frac, fit_half_window)

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
    return t1_vals, t2_vals, dt_vals, sk1, sk2, skd


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
    fitfn.SetNpx(500)
    fitfn.SetParameter(0, max(float(h.GetMaximum()), 1e-9))
    fitfn.SetParameter(1, float(h.GetMean()))
    fitfn.SetParameter(2, max(float(h.GetRMS()), 1e-6))
    status = int(h.Fit(fitfn, "QNR", "", xmin, xmax))
    if status != 0:
        return None
    fitfn.SetLineColor(ROOT.kRed)
    fitfn.SetLineWidth(2)
    fitfn.SetLineStyle(2)
    sigma = float(fitfn.GetParameter(2))
    sigma_err = float(fitfn.GetParError(2))
    return fitfn, sigma, sigma_err


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
        fitfn, sigma, sigma_err = fit
        fitfn.Draw("same")
        tx = ROOT.TLatex()
        tx.SetNDC()
        tx.SetTextFont(42)
        tx.SetTextSize(0.050)
        tx.SetTextColor(ROOT.kRed + 1)
        tx.DrawLatex(0.14, 0.86, f"#sigma = {sigma:.4f} #pm {sigma_err:.4f} ns")
    bw = ROOT.TLatex()
    bw.SetNDC()
    bw.SetTextFont(42)
    bw.SetTextSize(0.028)
    bw.DrawLatex(0.14, 0.08, f"hist bin = {bin_ps:.1f} ps")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Local-linear rising-edge CFD(peak fraction) 기반 T1/T2/Δ 3패널"
    )
    parser.add_argument("input", nargs="?", default=None)
    parser.add_argument("-o", "--output", default=None)
    parser.add_argument("--max-events", type=int, default=3000)
    parser.add_argument("--peak-frac", type=float, default=0.3)
    parser.add_argument("--rebin-ps", type=float, default=100.0)
    parser.add_argument("--plot-bin-ps", type=float, default=1.0)
    parser.add_argument("--fit-half-window", type=int, default=2, help="교차점 주변 반쪽 윈도우(bin)")
    parser.add_argument("--fixed-range", action="store_true")
    parser.add_argument("--xmin", type=float, default=-5.0)
    parser.add_argument("--xmax", type=float, default=5.0)
    parser.add_argument("-l", "--rootio-lib", default=None)
    args = parser.parse_args()

    if not (0.0 < args.peak_frac < 1.0):
        raise SystemExit("--peak-frac는 0~1 사이여야 합니다.")
    if args.fit_half_window < 1:
        raise SystemExit("--fit-half-window는 1 이상이어야 합니다.")

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
    _load_rootio(args.rootio_lib, repo_root)
    import ROOT

    ROOT.gROOT.SetBatch(True)
    _apply_style(ROOT)

    t1_vals, t2_vals, dt_vals, sk1, sk2, skd = _collect_toa(
        in_path,
        max_events=args.max_events,
        peak_frac=args.peak_frac,
        rebin_ns=args.rebin_ps / 1000.0,
        fit_half_window=args.fit_half_window,
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

    h1 = _make_hist(ROOT, "hToaT1", t1_vals, lo1, hi1, nb1, "T1;TOA local-linear CFD (ns);events")
    h2 = _make_hist(ROOT, "hToaT2", t2_vals, lo2, hi2, nb2, "T2;TOA local-linear CFD (ns);events")
    h3 = _make_hist(ROOT, "hToaDiffT1T2", dt_vals, lo3, hi3, nb3, "T1 - T2;#DeltaTOA (ns);events")

    stem = os.path.splitext(os.path.basename(in_path))[0]
    if args.output:
        out = os.path.abspath(args.output)
        if not os.path.splitext(out)[1]:
            out += ".png"
    else:
        out = os.path.join(_figures_dir(repo_root), f"cfd_locallinear_per_trigger_{stem}_1x3.png")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    c = ROOT.TCanvas("cToaPerTrigLocalLinear", "local linear CFD", 1860, 560)
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
        f"입력: {in_path}\n"
        f"  T1: filled={len(t1_vals)}, skipped={sk1}\n"
        f"  T2: filled={len(t2_vals)}, skipped={sk2}\n"
        f"  T1-T2: filled={len(dt_vals)}, skipped={skd}\n"
        f"  TOA 정의: local-linear rising edge at {args.peak_frac:g} * peak, rebin={args.rebin_ps:g} ps\n"
        f"PNG 저장: {out}"
    )


if __name__ == "__main__":
    main()
