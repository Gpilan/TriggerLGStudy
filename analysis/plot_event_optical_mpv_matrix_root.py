#!/usr/bin/env python3
"""
MPV(란다우 피팅) 기반 이벤트 선택으로 3x2 매트릭스 플롯을 그립니다 (PyROOT).

레이아웃 (2열 x 3행):
  1행: T1/T2 이벤트별 총 포톤수 히스토그램 + Landau 피팅(MPV 표시)
  2행: MPV에 가장 가까운 이벤트의 100 ps 타이밍 히스토그램 + Landau 피팅
  3행: 같은 이벤트의 200 ps 타이밍 히스토그램 + Landau 피팅
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass

from compare_timing_resolution import _rebin_merged_to_width
from plot_trigger_timing import (
    _figures_dir,
    _find_repo_root,
    _load_rootio,
    _merged_from_tower,
    _prepend_build_rootio_ld_path,
    _resolve_input_path,
)


@dataclass
class EventSummary:
    event_id: int
    nphoton: int
    lo: list[float]
    hi: list[float]
    counts: list[int]


def _sum_counts(counts: list[int]) -> int:
    return int(sum(int(c) for c in counts)) if counts else 0


def _scan_events(path: str, trig: int, max_events: int) -> list[EventSummary]:
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

    out: list[EventSummary] = []
    for i in range(nmax):
        tree.GetEntry(i)
        lo, hi, counts = _merged_from_tower(evt, trig)
        nph = _sum_counts(counts)
        if nph <= 0:
            continue
        out.append(
            EventSummary(
                event_id=i,
                nphoton=nph,
                lo=list(lo),
                hi=list(hi),
                counts=[int(x) for x in counts],
            )
        )

    f.Close()
    return out


def _pick_nearest_to_target(events: list[EventSummary], target: float) -> EventSummary:
    if not events:
        raise RuntimeError("유효 이벤트가 없습니다.")
    return min(events, key=lambda ev: abs(float(ev.nphoton) - float(target)))


def _make_time_hist(ROOT, name: str, ev: EventSummary, bin_ps: float, title: str):
    bin_ns = float(bin_ps) / 1000.0
    lo, hi, cnt = _rebin_merged_to_width(ev.lo, ev.hi, ev.counts, bin_ns)
    if not cnt:
        h = ROOT.TH1F(name, title, 1, 0.0, 1.0)
        h.SetDirectory(0)
        h.SetEntries(0.0)
        return h, 0.0

    nb = len(cnt)
    xlo = float(lo[0])
    xhi = float(hi[-1])
    if xhi <= xlo:
        xhi = xlo + bin_ns
    h = ROOT.TH1F(name, title, nb, xlo, xhi)
    h.Sumw2()
    h.SetDirectory(0)
    h.SetLineColor(ROOT.kBlack)
    h.SetLineWidth(2)
    for i, v in enumerate(cnt, start=1):
        y = float(v)
        h.SetBinContent(i, y)
        h.SetBinError(i, y**0.5 if y > 0 else 0.0)
    total_photons = float(sum(float(v) for v in cnt))
    h.SetEntries(total_photons)
    return h, total_photons


def _fit_landau_on_hist(ROOT, h, fit_name: str):
    if h.GetEntries() <= 0:
        return None, None
    xmin = float(h.GetXaxis().GetXmin())
    xmax = float(h.GetXaxis().GetXmax())
    fitfn = ROOT.TF1(fit_name, "landau", xmin, xmax)
    fitfn.SetNpx(700)
    fitfn.SetParameter(0, max(float(h.GetMaximum()), 1.0))
    fitfn.SetParameter(1, float(h.GetXaxis().GetBinCenter(h.GetMaximumBin())))
    fitfn.SetParameter(2, max(float(h.GetRMS()) * 0.35, 1e-4))
    status = int(h.Fit(fitfn, "QNR", "", xmin, xmax))
    fitfn.SetLineColor(ROOT.kRed + 1)
    fitfn.SetLineWidth(2)
    fitfn.SetLineStyle(1)
    if status != 0:
        return None, None
    return fitfn, float(fitfn.GetParameter(1))


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


def _fit_local_linear_threshold_on_hist(ROOT, h, peak_frac: float, half_window: int, fit_name: str):
    nb = int(h.GetNbinsX())
    if nb < 2:
        return None, None
    y = [float(h.GetBinContent(i)) for i in range(1, nb + 1)]
    x = [float(h.GetXaxis().GetBinCenter(i)) for i in range(1, nb + 1)]
    i_peak = max(range(nb), key=lambda i: y[i])
    peak = y[i_peak]
    if peak <= 0.0:
        return None, None
    thr = float(peak_frac) * peak

    i_cross = None
    for i in range(1, i_peak + 1):
        if y[i - 1] <= thr <= y[i]:
            i_cross = i
            break
    if i_cross is None:
        return None, None

    i0 = max(0, i_cross - half_window)
    i1 = min(i_peak, i_cross + half_window)
    xs = x[i0 : i1 + 1]
    ys = y[i0 : i1 + 1]
    if len(xs) < 2:
        return None, None
    ws = [1.0 / max(ys_i**0.5, 1.0) for ys_i in ys]
    fit = _weighted_linear_fit(xs, ys, ws)
    if fit is None:
        return None, None
    a, b = fit
    if a <= 1e-15:
        return None, None

    toa = (thr - b) / a
    xmin = min(xs)
    xmax = max(xs)
    if toa < xmin or toa > xmax:
        return None, None

    fline = ROOT.TF1(fit_name, "[0]*x+[1]", xmin, xmax)
    fline.SetParameter(0, a)
    fline.SetParameter(1, b)
    fline.SetLineColor(ROOT.kBlue + 2)
    fline.SetLineWidth(2)
    fline.SetLineStyle(1)
    return fline, float(toa)


def _build_nph_hist_and_landau_fit(ROOT, name: str, vals: list[int], title: str):
    if not vals:
        h = ROOT.TH1F(name, title, 1, 0.0, 1.0)
        h.SetDirectory(0)
        return h, None, None

    xmin = float(min(vals))
    xmax = float(max(vals))
    if xmax <= xmin:
        xmax = xmin + 1.0
    span = xmax - xmin
    nb = int(max(40, min(200, math.ceil(span / 15.0))))

    h = ROOT.TH1F(name, title, nb, xmin, xmax)
    h.Sumw2()
    h.SetDirectory(0)
    h.SetLineColor(ROOT.kBlack)
    h.SetLineWidth(2)
    for v in vals:
        h.Fill(float(v))

    mpv_guess = float(h.GetXaxis().GetBinCenter(h.GetMaximumBin()))
    sigma_guess = max(float(h.GetRMS()) * 0.4, 1.0)
    fitfn = ROOT.TF1(f"{name}_landau", "landau", xmin, xmax)
    fitfn.SetNpx(600)
    fitfn.SetParameter(0, max(float(h.GetMaximum()), 1.0))
    fitfn.SetParameter(1, mpv_guess)
    fitfn.SetParameter(2, sigma_guess)
    status = int(h.Fit(fitfn, "QNR", "", xmin, xmax))

    fitfn.SetLineColor(ROOT.kRed + 1)
    fitfn.SetLineWidth(2)
    fitfn.SetLineStyle(1)

    if status == 0:
        mpv = float(fitfn.GetParameter(1))
    else:
        fitfn = None
        mpv = mpv_guess
    return h, fitfn, mpv


def _draw_top_panel(ROOT, h, fitfn, mpv: float | None, ev: EventSummary, trig_label: str):
    h.Draw("HIST")
    if fitfn is not None:
        fitfn.Draw("same")

    tx = ROOT.TLatex()
    tx.SetNDC()
    tx.SetTextFont(42)
    tx.SetTextSize(0.033)
    if mpv is not None:
        tx.DrawLatex(0.14, 0.90, f"{trig_label}: MPV = {mpv:.1f}")
    tx.DrawLatex(0.14, 0.84, f"selected event = {ev.event_id}")
    tx.DrawLatex(0.14, 0.78, f"selected N = {ev.nphoton}")


def _draw_timing_panel(ROOT, h, trig_label: str, ev: EventSummary, peak_frac: float, half_window: int):
    h.Draw("HIST")
    fline, toa = _fit_local_linear_threshold_on_hist(
        ROOT, h, peak_frac, half_window, f"fit_local_{h.GetName()}"
    )
    if fline is not None:
        fline.Draw("same")
    tx = ROOT.TLatex()
    tx.SetNDC()
    tx.SetTextFont(42)
    tx.SetTextSize(0.033)
    tx.DrawLatex(0.14, 0.93, f"{trig_label} event={ev.event_id}, N={ev.nphoton}")
    if toa is not None:
        tx.DrawLatex(0.14, 0.87, f"Local-linear TOA ({peak_frac:.2f}*peak) = {toa:.4f} ns")


def draw_one(
    input_path: str,
    output_path: str,
    max_events: int,
    xmax: float,
    fine_bin_ps: float,
    coarse_bin_ps: float,
    peak_frac: float,
    fit_half_window: int,
) -> None:
    repo_root = _find_repo_root()
    if repo_root is None:
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    _prepend_build_rootio_ld_path(repo_root)
    _load_rootio(None, repo_root)

    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(1110)
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetHistLineWidth(2)
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")

    evs_t1 = _scan_events(input_path, trig=0, max_events=max_events)
    evs_t2 = _scan_events(input_path, trig=1, max_events=max_events)
    if not evs_t1 or not evs_t2:
        raise RuntimeError("T1/T2 유효 이벤트가 부족합니다.")

    nvals_t1 = [ev.nphoton for ev in evs_t1]
    nvals_t2 = [ev.nphoton for ev in evs_t2]

    h_n1, f_n1, mpv1 = _build_nph_hist_and_landau_fit(
        ROOT, "hNphT1", nvals_t1, "T1 N photons / event;N photons;events"
    )
    h_n2, f_n2, mpv2 = _build_nph_hist_and_landau_fit(
        ROOT, "hNphT2", nvals_t2, "T2 N photons / event;N photons;events"
    )

    ev_t1 = _pick_nearest_to_target(evs_t1, mpv1 if mpv1 is not None else 0.0)
    ev_t2 = _pick_nearest_to_target(evs_t2, mpv2 if mpv2 is not None else 0.0)

    h_t1_fine, _ = _make_time_hist(
        ROOT,
        "h_t1_fine_mpv",
        ev_t1,
        fine_bin_ps,
        f"T1 selected event ({fine_bin_ps:g} ps);time (ns);photons/bin",
    )
    h_t2_fine, _ = _make_time_hist(
        ROOT,
        "h_t2_fine_mpv",
        ev_t2,
        fine_bin_ps,
        f"T2 selected event ({fine_bin_ps:g} ps);time (ns);photons/bin",
    )
    h_t1_coarse, _ = _make_time_hist(
        ROOT,
        "h_t1_coarse_mpv",
        ev_t1,
        coarse_bin_ps,
        f"T1 selected event ({coarse_bin_ps:g} ps);time (ns);photons/bin",
    )
    h_t2_coarse, _ = _make_time_hist(
        ROOT,
        "h_t2_coarse_mpv",
        ev_t2,
        coarse_bin_ps,
        f"T2 selected event ({coarse_bin_ps:g} ps);time (ns);photons/bin",
    )

    c = ROOT.TCanvas("c_evt_mpv_matrix", "event optical MPV matrix", 1400, 1350)
    c.Divide(2, 3)

    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)
    _draw_top_panel(ROOT, h_n1, f_n1, mpv1, ev_t1, "T1")

    c.cd(2)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)
    _draw_top_panel(ROOT, h_n2, f_n2, mpv2, ev_t2, "T2")

    for pad_idx, h, trig_label, ev in (
        (3, h_t1_fine, "T1", ev_t1),
        (4, h_t2_fine, "T2", ev_t2),
        (5, h_t1_coarse, "T1", ev_t1),
        (6, h_t2_coarse, "T2", ev_t2),
    ):
        c.cd(pad_idx)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetBottomMargin(0.12)
        h.GetXaxis().SetRangeUser(0.0, float(xmax))
        _draw_timing_panel(ROOT, h, trig_label, ev, peak_frac, fit_half_window)

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    c.SaveAs(output_path)
    c.Close()

    print(
        f"입력: {input_path}\n"
        f"  T1 MPV={mpv1:.3f}, selected event={ev_t1.event_id}, N={ev_t1.nphoton}\n"
        f"  T2 MPV={mpv2:.3f}, selected event={ev_t2.event_id}, N={ev_t2.nphoton}\n"
        f"PNG 저장: {output_path}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Landau MPV 기반 이벤트 선택으로 T1/T2 3x2 매트릭스 생성"
    )
    parser.add_argument(
        "input",
        nargs="?",
        default=None,
        help="입력 .root (기본: analysis/t_res/data/v2_60GeV_e-_LG_3000_0.root)",
    )
    parser.add_argument("--max-events", type=int, default=-1)
    parser.add_argument("--xmax", type=float, default=20.0, help="타이밍 패널 x축 최대(ns)")
    parser.add_argument("--fine-bin-ps", type=float, default=100.0, help="2행 타이밍 bin (ps)")
    parser.add_argument("--coarse-bin-ps", type=float, default=200.0, help="3행 타이밍 bin (ps)")
    parser.add_argument("--peak-frac", type=float, default=0.3, help="국소 선형 TOA 기준 peak 비율")
    parser.add_argument("--fit-half-window", type=int, default=2, help="교차점 주변 반쪽 윈도우(bin)")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="출력 PNG (기본: figures/event_optical_mpv_matrix_<stem>.png)",
    )
    parser.add_argument("-l", "--rootio-lib", default=None)
    args = parser.parse_args()

    if args.xmax <= 0.0:
        raise SystemExit("--xmax는 0보다 커야 합니다.")
    if args.fine_bin_ps <= 0.0 or args.coarse_bin_ps <= 0.0:
        raise SystemExit("--fine-bin-ps, --coarse-bin-ps는 0보다 커야 합니다.")
    if not (0.0 < args.peak_frac < 1.0):
        raise SystemExit("--peak-frac는 0~1 사이여야 합니다.")
    if args.fit_half_window < 1:
        raise SystemExit("--fit-half-window는 1 이상이어야 합니다.")

    repo_root = _find_repo_root()
    if repo_root is None:
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

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
    stem = os.path.splitext(os.path.basename(in_path))[0]
    out = (
        os.path.abspath(args.output)
        if args.output
        else os.path.join(_figures_dir(repo_root), f"event_optical_mpv_matrix_{stem}.png")
    )
    if not os.path.splitext(out)[1]:
        out += ".png"

    draw_one(
        in_path,
        out,
        max_events=args.max_events,
        xmax=args.xmax,
        fine_bin_ps=args.fine_bin_ps,
        coarse_bin_ps=args.coarse_bin_ps,
        peak_frac=args.peak_frac,
        fit_half_window=args.fit_half_window,
    )


if __name__ == "__main__":
    main()
