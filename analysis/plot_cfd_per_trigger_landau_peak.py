#!/usr/bin/env python3
"""
이벤트별 타이밍 분포를 Landau로 피팅하고, peak의 일정 비율(기본 0.3) 교차점을 TOA로 정의해
T1/T2/Δ(T1-T2) 1x3 히스토그램을 그립니다.
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


def _auto_axis(
    vals: list[float],
    bin_width: float,
    *,
    margin_frac: float = 0.08,
) -> tuple[float, float, int]:
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


def _fit_landau_toa(
    ROOT,
    lo: list[float],
    hi: list[float],
    cnt: list[int],
    peak_frac: float,
    *,
    fit_name: str,
) -> float | None:
    if not cnt:
        return None
    nb = len(cnt)
    xlo = float(lo[0])
    xhi = float(hi[-1])
    if xhi <= xlo:
        return None

    htmp = ROOT.TH1F(f"{fit_name}_h", "", nb, xlo, xhi)
    htmp.SetDirectory(0)
    for i, v in enumerate(cnt, start=1):
        htmp.SetBinContent(i, float(v))

    fit = ROOT.TF1(fit_name, "landau", xlo, xhi)
    fit.SetNpx(700)
    fit.SetParameter(0, max(float(htmp.GetMaximum()), 1.0))
    fit.SetParameter(1, float(htmp.GetXaxis().GetBinCenter(htmp.GetMaximumBin())))
    fit.SetParameter(2, max(float(htmp.GetRMS()) * 0.35, 1e-4))
    status = int(htmp.Fit(fit, "QNR", "", xlo, xhi))
    if status != 0:
        return None

    mpv_x = float(fit.GetParameter(1))
    mpv_y = float(fit.Eval(mpv_x))
    thr_y = float(peak_frac) * mpv_y
    if not (xlo < mpv_x <= xhi):
        return None

    nscan = 2000
    xs = [xlo + (mpv_x - xlo) * i / float(nscan) for i in range(nscan + 1)]
    ys = [float(fit.Eval(x)) for x in xs]
    for i in range(1, len(xs)):
        y0 = ys[i - 1]
        y1 = ys[i]
        if y0 <= thr_y <= y1:
            x0 = xs[i - 1]
            x1 = xs[i]
            if y1 == y0:
                return 0.5 * (x0 + x1)
            t = (thr_y - y0) / (y1 - y0)
            return x0 + t * (x1 - x0)
    return None


def _collect_toa_landau_peak(
    path: str,
    *,
    max_events: int,
    peak_frac: float,
    cfd_input_bin_ns: float,
) -> tuple[list[float], list[float], list[float], int, int, int]:
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
    skipped_t1 = 0
    skipped_t2 = 0
    skipped_dt = 0

    for i in range(nmax):
        tree.GetEntry(i)
        lo1, hi1, c1 = _merged_from_tower(evt, 0)
        lo2, hi2, c2 = _merged_from_tower(evt, 1)
        lo1, hi1, c1 = _rebin_merged_to_width(lo1, hi1, c1, cfd_input_bin_ns)
        lo2, hi2, c2 = _rebin_merged_to_width(lo2, hi2, c2, cfd_input_bin_ns)

        v1 = _fit_landau_toa(ROOT, lo1, hi1, c1, peak_frac, fit_name=f"fit_e{i}_t1")
        v2 = _fit_landau_toa(ROOT, lo2, hi2, c2, peak_frac, fit_name=f"fit_e{i}_t2")

        if v1 is not None:
            t1_vals.append(v1)
        else:
            skipped_t1 += 1
        if v2 is not None:
            t2_vals.append(v2)
        else:
            skipped_t2 += 1
        if v1 is not None and v2 is not None:
            dt_vals.append(float(v1 - v2))
        else:
            skipped_dt += 1

    f.Close()
    return t1_vals, t2_vals, dt_vals, skipped_t1, skipped_t2, skipped_dt


def _apply_root_style(ROOT) -> None:
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
    return fitfn


def _draw_panel(ROOT, h, vals: list[float], fit_name: str, *, bin_ps: float, show_no_data_x: float = 0.40):
    h.Draw("HIST")
    if not vals:
        t = ROOT.TLatex()
        t.SetNDC()
        t.SetTextFont(42)
        t.SetTextSize(0.04)
        t.DrawLatex(show_no_data_x, 0.50, "no data")
        return
    fitfn = _fit_gaus(ROOT, h, fit_name)
    if fitfn is not None:
        fitfn.Draw("same")
    bw = ROOT.TLatex()
    bw.SetNDC()
    bw.SetTextFont(42)
    bw.SetTextSize(0.028)
    bw.DrawLatex(0.14, 0.08, f"hist bin = {bin_ps:.1f} ps")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Landau fit peak-fraction TOA 기반 T1/T2/Δ 3패널"
    )
    parser.add_argument("input", nargs="?", default=None)
    parser.add_argument("-o", "--output", default=None)
    parser.add_argument("--max-events", type=int, default=3000)
    parser.add_argument("--peak-frac", type=float, default=0.3, help="Landau peak 높이 대비 비율")
    parser.add_argument("--rebin-ps", type=float, default=100.0, help="TOA 계산용 rebin 폭 (ps)")
    parser.add_argument("--plot-bin-ps", type=float, default=1.0, help="최종 히스토그램 bin 폭 (ps)")
    parser.add_argument("--fixed-range", action="store_true")
    parser.add_argument("--xmin", type=float, default=-5.0)
    parser.add_argument("--xmax", type=float, default=5.0)
    parser.add_argument("-l", "--rootio-lib", default=None)
    args = parser.parse_args()

    if not (0.0 < args.peak_frac < 1.0):
        raise SystemExit("--peak-frac는 0~1 사이여야 합니다.")

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
    _apply_root_style(ROOT)

    t1_vals, t2_vals, dt_vals, sk1, sk2, skd = _collect_toa_landau_peak(
        in_path,
        max_events=args.max_events,
        peak_frac=args.peak_frac,
        cfd_input_bin_ns=args.rebin_ps / 1000.0,
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

    h1 = _make_hist(ROOT, "hToaT1", t1_vals, lo1, hi1, nb1, "T1;TOA from Landau-peak CFD (ns);events")
    h2 = _make_hist(ROOT, "hToaT2", t2_vals, lo2, hi2, nb2, "T2;TOA from Landau-peak CFD (ns);events")
    h3 = _make_hist(ROOT, "hToaDiffT1T2", dt_vals, lo3, hi3, nb3, "T1 - T2;#DeltaTOA (ns);events")

    stem = os.path.splitext(os.path.basename(in_path))[0]
    if args.output:
        out = os.path.abspath(args.output)
        if not os.path.splitext(out)[1]:
            out += ".png"
    else:
        out = os.path.join(_figures_dir(repo_root), f"cfd_landau_peak_per_trigger_{stem}_1x3.png")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    c = ROOT.TCanvas("cToaPerTrigLandauPeak", "landau peak CFD", 1860, 560)
    c.Divide(3, 1)
    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)
    _draw_panel(ROOT, h1, t1_vals, "fit_t1", bin_ps=((hi1 - lo1) / float(nb1)) * 1000.0)
    c.cd(2)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)
    _draw_panel(ROOT, h2, t2_vals, "fit_t2", bin_ps=((hi2 - lo2) / float(nb2)) * 1000.0)
    c.cd(3)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)
    _draw_panel(
        ROOT,
        h3,
        dt_vals,
        "fit_t1_minus_t2",
        bin_ps=((hi3 - lo3) / float(nb3)) * 1000.0,
        show_no_data_x=0.35,
    )
    c.SaveAs(out)
    c.Close()

    print(
        f"입력: {in_path}\n"
        f"  T1: filled={len(t1_vals)}, skipped={sk1}\n"
        f"  T2: filled={len(t2_vals)}, skipped={sk2}\n"
        f"  T1-T2: filled={len(dt_vals)}, skipped={skd}\n"
        f"  TOA 정의: Landau fit rising edge at {args.peak_frac:g} * peak, rebin={args.rebin_ps:g} ps\n"
        f"PNG 저장: {out}"
    )


if __name__ == "__main__":
    main()
