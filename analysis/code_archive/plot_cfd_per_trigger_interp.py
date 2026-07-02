#!/usr/bin/env python3
"""
단일 CBDsim .root에서 이벤트별 TOA를 계산해 T1/T2 히스토그램으로 그립니다.

TOA 정의:
1) merged SiPM 시간 빈을 균일 폭으로 리빈(기본 100 ps)
2) 리빈된 빈 중심 카운트를 선형 보간(기본 on)
3) CFD 시각 계산 (기본 fraction=0.3)

출력:
- 1x3 패널 PNG (T1 | T2 | T1-T2)
- x축: event-level TOA (ns)
- y축: events
"""

from __future__ import annotations

import argparse
import math
import os

from compare_timing_resolution import (
    _cfd_absolute_time_float,
    _interpolate_counts_linear,
    _rebin_merged_to_width,
)
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


def _collect_cfd_times_interp(
    path: str,
    *,
    max_events: int,
    ns_per_unit: float,
    cfd_fraction: float,
    cfd_input_bin_ns: float,
    cfd_interp: bool,
) -> tuple[list[float], list[float], list[float], int, int, int]:
    """T1/T2/Δ(T1-T2) TOA 리스트와 스킵 수."""
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

        if cfd_interp:
            lo1, hi1, c1 = _interpolate_counts_linear(lo1, hi1, c1)
            lo2, hi2, c2 = _interpolate_counts_linear(lo2, hi2, c2)

        v1 = _cfd_absolute_time_float(lo1, hi1, c1, cfd_fraction, ns_per_unit)
        v2 = _cfd_absolute_time_float(lo2, hi2, c2, cfd_fraction, ns_per_unit)

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
    """ROOT 기본 스타일에 가깝게 설정."""
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
    """가우시안 피팅 결과 (A, Aerr, mu, mu_err, sigma, sigma_err, chi2, ndf, prob, fitfn)."""
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
    amp = float(fitfn.GetParameter(0))
    amp_err = float(fitfn.GetParError(0))
    mu = float(fitfn.GetParameter(1))
    mu_err = float(fitfn.GetParError(1))
    sigma = float(fitfn.GetParameter(2))
    sigma_err = float(fitfn.GetParError(2))
    chi2 = float(fitfn.GetChisquare())
    ndf = int(fitfn.GetNDF())
    prob = float(fitfn.GetProb())
    fitfn.SetLineColor(ROOT.kRed)
    fitfn.SetLineWidth(2)
    fitfn.SetLineStyle(2)
    return amp, amp_err, mu, mu_err, sigma, sigma_err, chi2, ndf, prob, fitfn


def _draw_panel(ROOT, h, vals: list[float], fit_name: str, *, bin_ps: float, show_no_data_x: float = 0.40):
    h.Draw("HIST")
    if not vals:
        t = ROOT.TLatex()
        t.SetNDC()
        t.SetTextFont(42)
        t.SetTextSize(0.04)
        t.DrawLatex(show_no_data_x, 0.50, "no data")
        return None

    fit = _fit_gaus(ROOT, h, fit_name)
    if fit is not None:
        amp, amp_err, mu, mu_err, sigma, sigma_err, chi2, ndf, prob, fitfn = fit
        fitfn.Draw("same")

        # 통계 박스를 우상단으로 고정
        ROOT.gPad.Update()
        st = h.FindObject("stats")
        if st:
            st.SetX1NDC(0.67)
            st.SetX2NDC(0.97)
            st.SetY1NDC(0.62)
            st.SetY2NDC(0.94)
            ROOT.gPad.Modified()
            ROOT.gPad.Update()

        # sigma를 크게 표시
        tx = ROOT.TLatex()
        tx.SetNDC()
        tx.SetTextFont(42)
        tx.SetTextSize(0.055)
        tx.SetTextColor(ROOT.kRed + 1)
        tx.DrawLatex(0.14, 0.86, f"#sigma = {sigma:.4f} #pm {sigma_err:.4f} ns")

        # 피팅 품질/파라미터 요약
        info = ROOT.TLatex()
        info.SetNDC()
        info.SetTextFont(42)
        info.SetTextSize(0.028)
        red_chi2 = (chi2 / float(ndf)) if ndf > 0 else 0.0
        info.DrawLatex(0.14, 0.79, f"A = {amp:.2f} #pm {amp_err:.2f}")
        info.DrawLatex(0.14, 0.73, f"#mu = {mu:.4f} #pm {mu_err:.4f} ns")
        info.DrawLatex(0.14, 0.67, f"#chi^{{2}}/NDF = {red_chi2:.3f} ({chi2:.1f}/{ndf})")
        info.DrawLatex(0.14, 0.61, f"p-value = {prob:.3g}")

        # 실제 히스토그램 bin 폭을 명시 (ps)
        bw = ROOT.TLatex()
        bw.SetNDC()
        bw.SetTextFont(42)
        bw.SetTextSize(0.028)
        bw.DrawLatex(0.14, 0.08, f"hist bin = {bin_ps:.1f} ps")
        return amp, amp_err, mu, mu_err, sigma, sigma_err, chi2, ndf, prob
    return None


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Per-event TOA histogram (CFD 0.3 + rebin/interpolation), T1/T2/Δ 3패널"
    )
    parser.add_argument(
        "input",
        nargs="?",
        default=None,
        help="입력 .root (기본: analysis/t_res/data/NoLG_60GeV_e-_noLG_1.root)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="출력 PNG (기본: figures/cfd0p3_interp_per_trigger_<stem>_1x3.png)",
    )
    parser.add_argument(
        "--fixed-range",
        action="store_true",
        help="고정 x축 사용 (--xmin/--xmax/--bins)",
    )
    parser.add_argument("--xmin", type=float, default=-5.0)
    parser.add_argument("--xmax", type=float, default=5.0)
    parser.add_argument(
        "--bin-width",
        type=float,
        default=0.01,
        help="자동 축에서 히스토그램 bin 폭 (ns)",
    )
    parser.add_argument("--ns-per-unit", type=float, default=1.0)
    parser.add_argument("--cfd-fraction", type=float, default=0.3, dest="cfd_fraction")
    parser.add_argument(
        "--cfd-input-bin-ps",
        type=float,
        default=100.0,
        dest="cfd_input_bin_ps",
        help="CFD 전 리빈 폭 (ps, 기본 100)",
    )
    parser.add_argument(
        "--plot-bin-ps",
        type=float,
        default=None,
        dest="plot_bin_ps",
        help="플롯 히스토그램 bin 폭 (ps). 미지정 시 CFD 리빈 폭과 동일",
    )
    parser.add_argument(
        "--no-cfd-interp",
        action="store_true",
        help="리빈 후 선형 interpolation 비활성화",
    )
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

    in_path = _resolve(
        args.input,
        "analysis/t_res/data/NoLG_60GeV_e-_noLG_1.root",
    )

    _load_rootio(args.rootio_lib, repo_root)
    import ROOT

    ROOT.gROOT.SetBatch(True)
    _apply_root_style(ROOT)

    cfd_input_bin_ns = args.cfd_input_bin_ps / 1000.0
    plot_bin_ns = (args.plot_bin_ps / 1000.0) if args.plot_bin_ps is not None else cfd_input_bin_ns
    t1_vals, t2_vals, dt_vals, sk1, sk2, skd = _collect_cfd_times_interp(
        in_path,
        max_events=args.max_events,
        ns_per_unit=args.ns_per_unit,
        cfd_fraction=args.cfd_fraction,
        cfd_input_bin_ns=cfd_input_bin_ns,
        cfd_interp=not args.no_cfd_interp,
    )

    stem = os.path.splitext(os.path.basename(in_path))[0]
    fig_dir = _figures_dir(repo_root)
    if args.output:
        out = os.path.abspath(args.output)
        if not os.path.splitext(out)[1]:
            out += ".png"
    else:
        out = os.path.join(fig_dir, f"cfd0p3_interp_per_trigger_{stem}_1x3.png")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    print(
        f"입력: {in_path}\n"
        f"  T1: filled={len(t1_vals)}, skipped(empty CFD)={sk1}\n"
        f"  T2: filled={len(t2_vals)}, skipped(empty CFD)={sk2}\n"
        f"  T1-T2: filled={len(dt_vals)}, skipped(pair missing)={skd}\n"
        f"  CFD fraction={args.cfd_fraction}, rebin={args.cfd_input_bin_ps} ps, "
        f"interp={'off' if args.no_cfd_interp else 'on'}, plot-bin={plot_bin_ns*1000.0:g} ps"
    )

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

    bw1_ps = ((hi1 - lo1) / float(nb1)) * 1000.0
    bw2_ps = ((hi2 - lo2) / float(nb2)) * 1000.0
    bw3_ps = ((hi3 - lo3) / float(nb3)) * 1000.0
    print(
        f"  Hist bin widths: T1={bw1_ps:.3f} ps, "
        f"T2={bw2_ps:.3f} ps, T1-T2={bw3_ps:.3f} ps"
    )

    h1 = _make_hist(
        ROOT,
        "hToaT1",
        t1_vals,
        lo1,
        hi1,
        nb1,
        "T1;TOA from CFD (ns);events",
    )
    h2 = _make_hist(
        ROOT,
        "hToaT2",
        t2_vals,
        lo2,
        hi2,
        nb2,
        "T2;TOA from CFD (ns);events",
    )
    h3 = _make_hist(
        ROOT,
        "hToaDiffT1T2",
        dt_vals,
        lo3,
        hi3,
        nb3,
        "T1 - T2;#DeltaTOA (ns);events",
    )

    canvas_title = (
        f"Per-event TOA (CFD={args.cfd_fraction:g}, "
        f"rebin={args.cfd_input_bin_ps:g} ps, interp={'off' if args.no_cfd_interp else 'on'})"
    )
    c = ROOT.TCanvas("cToaPerTrig", canvas_title, 1860, 560)
    c.Divide(3, 1)

    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)
    fit1 = _draw_panel(ROOT, h1, t1_vals, "fit_t1", bin_ps=bw1_ps)
    if fit1 is not None:
        amp, amp_err, mu, mu_err, sigma, sigma_err, chi2, ndf, prob = fit1
        print(
            f"  T1 Gauss fit: A={amp:.3f}±{amp_err:.3f}, mu={mu:.6f}±{mu_err:.6f} ns, "
            f"sigma={sigma:.6f}±{sigma_err:.6f} ns, chi2/ndf="
            f"{(chi2/ndf if ndf>0 else 0.0):.4f}, p={prob:.3g}"
        )

    c.cd(2)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)
    fit2 = _draw_panel(ROOT, h2, t2_vals, "fit_t2", bin_ps=bw2_ps)
    if fit2 is not None:
        amp, amp_err, mu, mu_err, sigma, sigma_err, chi2, ndf, prob = fit2
        print(
            f"  T2 Gauss fit: A={amp:.3f}±{amp_err:.3f}, mu={mu:.6f}±{mu_err:.6f} ns, "
            f"sigma={sigma:.6f}±{sigma_err:.6f} ns, chi2/ndf="
            f"{(chi2/ndf if ndf>0 else 0.0):.4f}, p={prob:.3g}"
        )

    c.cd(3)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)
    fit3 = _draw_panel(ROOT, h3, dt_vals, "fit_t1_minus_t2", bin_ps=bw3_ps, show_no_data_x=0.35)
    if fit3 is not None:
        amp, amp_err, mu, mu_err, sigma, sigma_err, chi2, ndf, prob = fit3
        print(
            f"  T1-T2 Gauss fit: A={amp:.3f}±{amp_err:.3f}, mu={mu:.6f}±{mu_err:.6f} ns, "
            f"sigma={sigma:.6f}±{sigma_err:.6f} ns, chi2/ndf="
            f"{(chi2/ndf if ndf>0 else 0.0):.4f}, p={prob:.3g}"
        )

    c.SaveAs(out)
    c.Close()
    print(f"PNG 저장: {out}")


if __name__ == "__main__":
    main()
