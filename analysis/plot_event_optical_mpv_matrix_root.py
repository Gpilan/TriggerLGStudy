#!/usr/bin/env python3
"""
MPV(anchor) 기반 이벤트 선택으로 3x2 매트릭스 플롯을 그립니다 (PyROOT).

레이아웃 (2열 x 3행):
  1행: T1/T2 이벤트별 총 포톤수 히스토그램 + Landau 피팅(MPV 표시)
  2행: 선택 이벤트의 200 ps photon-arrival 분포
  3행: 같은 이벤트의 reconstructed waveform + local-linear CFD
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from dataclasses import dataclass
from typing import Optional

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
CODE_ARCHIVE_DIR = os.path.join(THIS_DIR, "code_archive")
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)
if CODE_ARCHIVE_DIR not in sys.path:
    sys.path.insert(0, CODE_ARCHIVE_DIR)

from code_archive.compare_timing_resolution import _rebin_merged_to_width
from code_archive.plot_trigger_timing import (
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


@dataclass
class ReconstructedPulse:
    time: list[float]
    voltage: list[float]
    impulse: list[float]
    peak: float
    t10: Optional[float]
    t_frac: Optional[float]
    t90: Optional[float]
    photon_t_frac: Optional[float]


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
    pos = [ev for ev in events if ev.nphoton > 0]
    pool = pos if pos else events
    return min(pool, key=lambda ev: abs(float(ev.nphoton) - float(target)))


def _event_map(events: list[EventSummary]) -> dict[int, EventSummary]:
    return {int(ev.event_id): ev for ev in events}


def _must_get_event(events_by_id: dict[int, EventSummary], event_id: int, label: str) -> EventSummary:
    ev = events_by_id.get(int(event_id))
    if ev is None:
        raise RuntimeError(f"{label}에서 event_id={event_id}를 찾지 못했습니다.")
    return ev


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


def _fit_local_linear_threshold_on_hist(ROOT, h, peak_frac: float, half_window: int, fit_name: str):
    nb = int(h.GetNbinsX())
    if nb < 2:
        return None, None, None
    y = [float(h.GetBinContent(i)) for i in range(1, nb + 1)]
    x = [float(h.GetXaxis().GetBinCenter(i)) for i in range(1, nb + 1)]
    i_peak = max(range(nb), key=lambda i: y[i])
    peak = y[i_peak]
    if peak <= 0.0:
        return None, None, None
    thr = float(peak_frac) * peak

    i_cross = None
    for i in range(1, i_peak + 1):
        if y[i - 1] <= thr <= y[i]:
            i_cross = i
            break
    if i_cross is None:
        return None, None, None

    i0 = max(0, i_cross - half_window)
    i1 = min(i_peak, i_cross + half_window)
    xs = x[i0 : i1 + 1]
    ys = y[i0 : i1 + 1]
    if len(xs) < 2:
        return None, None, None
    ws = [1.0 / max(ys_i**0.5, 1.0) for ys_i in ys]
    fit = _weighted_linear_fit(xs, ys, ws)
    if fit is None:
        return None, None, None
    a, b = fit
    if a <= 1e-15:
        return None, None, None

    toa = (thr - b) / a
    xmin = min(xs)
    xmax = max(xs)
    if toa < xmin or toa > xmax:
        return None, None, None

    r2 = _weighted_r2(xs, ys, ws, a, b)
    n = len(xs)
    yfit = [a * xx + b for xx in xs]
    # Weighted proxy chi2 for local linear fit quality.
    chi2 = float(sum(wi * (yy - ff) ** 2 for wi, yy, ff in zip(ws, ys, yfit)))
    ndf = int(max(n - 2, 0))
    # Approximate p-value for weighted residual metric.
    pval = float(ROOT.TMath.Prob(chi2, ndf)) if ndf > 0 else 1.0

    fline = ROOT.TF1(fit_name, "[0]*x+[1]", xmin, xmax)
    fline.SetParameter(0, a)
    fline.SetParameter(1, b)
    fline.SetLineColor(ROOT.kBlue + 2)
    fline.SetLineWidth(2)
    fline.SetLineStyle(1)
    meta = {
        "threshold": float(thr),
        "peak": float(peak),
        "r2": r2,
        "chi2": chi2,
        "ndf": ndf,
        "pvalue": pval,
        "xmin": float(xmin),
        "xmax": float(xmax),
    }
    return fline, float(toa), meta


def _crossing(time: list[float], values: list[float], peak_index: int, level: float) -> Optional[float]:
    if not time or not values or len(time) != len(values):
        return None
    i_max = int(min(max(peak_index, 0), len(values) - 1))
    for i in range(1, i_max + 1):
        y0 = float(values[i - 1])
        y1 = float(values[i])
        if y0 <= level <= y1:
            dy = y1 - y0
            if abs(dy) <= 1e-18:
                return float(time[i])
            x0 = float(time[i - 1])
            x1 = float(time[i])
            return float(x0 + (level - y0) * (x1 - x0) / dy)
    return None


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
            value_txt = str(row.get("value", "")).strip()
            try:
                value = float(value_txt)
            except ValueError:
                continue
            if parameter == "rise_time" and qualifier in ("typ", ""):
                rise_ns = value
            elif parameter == "transit_time" and qualifier in ("typ", ""):
                transit_ns = value
    return rise_ns, transit_ns


def _build_reconstructed_waveform(
    ev: EventSummary,
    sample_step_ns: float,
    response_rise_ns: float,
    response_fwhm_ns: float,
    gamma_shape: float,
    gamma_tau_ns: Optional[float],
    transit_time_ns: float,
    peak_frac: float,
) -> ReconstructedPulse:
    centers = [0.5 * (float(lo) + float(hi)) for lo, hi in zip(ev.lo, ev.hi)]
    valid = [(t, int(c)) for t, c in zip(centers, ev.counts) if int(c) > 0]
    if not valid:
        return ReconstructedPulse([], [], [], 0.0, None, None, None, None)

    first_hit = min(t for t, _ in valid)
    last_hit = max(t for t, _ in valid)
    start = math.floor((first_hit - 2.0 * sample_step_ns) / sample_step_ns) * sample_step_ns
    stop = last_hit + 12.0 * response_fwhm_ns
    n_samples = int(math.ceil((stop - start) / sample_step_ns)) + 1
    if n_samples < 2:
        n_samples = 2

    impulse = [0.0] * n_samples
    for hit_time, count in valid:
        pos = (float(hit_time) - start) / sample_step_ns
        i0 = int(math.floor(pos))
        frac = float(pos - i0)
        amp = float(count)
        # Spread each bin content to neighboring samples (linear interpolation)
        # to avoid stair-step quantization from nearest-sample assignment.
        if 0 <= i0 < n_samples:
            impulse[i0] += amp * (1.0 - frac)
        i1 = i0 + 1
        if 0 <= i1 < n_samples:
            impulse[i1] += amp * frac
    time_offset_ns = float(transit_time_ns)
    time = [start + i * sample_step_ns + time_offset_ns for i in range(n_samples)]

    photon_peak = max(range(n_samples), key=lambda i: impulse[i])
    photon_thr = float(peak_frac) * float(impulse[photon_peak])
    photon_t_frac = _crossing(time, impulse, photon_peak, photon_thr)

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

    peak_idx = max(range(n_samples), key=lambda i: voltage[i])
    peak = float(voltage[peak_idx])
    if peak <= 0.0:
        return ReconstructedPulse(time, voltage, impulse, 0.0, None, None, None, photon_t_frac)

    t10 = _crossing(time, voltage, peak_idx, 0.10 * peak)
    t_frac = _crossing(time, voltage, peak_idx, float(peak_frac) * peak)
    t90 = _crossing(time, voltage, peak_idx, 0.90 * peak)
    return ReconstructedPulse(time, voltage, impulse, peak, t10, t_frac, t90, photon_t_frac)


def _make_hist_from_samples(ROOT, name: str, xs: list[float], ys: list[float], title: str):
    if not xs or not ys or len(xs) != len(ys):
        h = ROOT.TH1F(name, title, 1, 0.0, 1.0)
        h.SetDirectory(0)
        h.SetEntries(0.0)
        return h

    if len(xs) >= 2:
        step = float(xs[1] - xs[0])
    else:
        step = 0.2
    if step <= 0.0:
        step = 0.2
    xlo = float(xs[0] - 0.5 * step)
    xhi = float(xs[-1] + 0.5 * step)
    if xhi <= xlo:
        xhi = xlo + step

    h = ROOT.TH1F(name, title, len(xs), xlo, xhi)
    h.Sumw2()
    h.SetDirectory(0)
    h.SetLineColor(ROOT.kBlue + 1)
    h.SetLineWidth(3)
    for i, y in enumerate(ys, start=1):
        yy = float(y)
        h.SetBinContent(i, yy)
        h.SetBinError(i, max(yy, 0.0) ** 0.5 if yy > 0.0 else 0.0)
    h.SetEntries(float(sum(max(float(v), 0.0) for v in ys)))
    return h


def _waveform_xlim(pulse: ReconstructedPulse, fallback_xmax: float) -> tuple[float, float]:
    if pulse.time:
        xlo = max(0.0, float(min(pulse.time)) - 0.5)
        xhi = float(max(pulse.time)) + 0.5
        if xhi > xlo:
            return xlo, xhi
    return 0.0, float(fallback_xmax)


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


def _draw_top_panel(
    ROOT,
    h,
    fitfn,
    mpv: float | None,
    ev: EventSummary,
    trig_label: str,
    anchor_mpv: float | None,
):
    h.Draw("HIST")
    if fitfn is not None:
        fitfn.Draw("same")

    tx = ROOT.TLatex()
    tx.SetNDC()
    tx.SetTextFont(42)
    tx.SetTextSize(0.033)
    if mpv is not None:
        tx.DrawLatex(0.14, 0.90, f"{trig_label}: MPV = {mpv:.1f}")
    if anchor_mpv is not None:
        tx.DrawLatex(0.14, 0.84, f"anchor MPV(LG T1) = {anchor_mpv:.1f}")
        tx.DrawLatex(0.14, 0.78, f"selected event = {ev.event_id}, N = {ev.nphoton}")
    else:
        tx.DrawLatex(0.14, 0.84, f"selected event = {ev.event_id}")
        tx.DrawLatex(0.14, 0.78, f"selected N = {ev.nphoton}")


def _draw_arrival_panel(ROOT, h, trig_label: str, ev: EventSummary):
    y_upper = max(float(h.GetMaximum()) * 1.08, 1.0)
    h.GetYaxis().SetRangeUser(0.0, y_upper)
    h.Draw("HIST")
    tx = ROOT.TLatex()
    tx.SetNDC()
    tx.SetTextFont(42)
    tx.SetTextSize(0.029)
    tx.DrawLatex(0.14, 0.88, f"{trig_label} photon-arrival (event={ev.event_id}, N={ev.nphoton})")
    h._overlay_keepalive = [tx]


def _draw_waveform_panel(
    ROOT,
    h,
    trig_label: str,
    ev: EventSummary,
    pulse: ReconstructedPulse,
    peak_frac: float,
    half_window: int,
):
    keepalive = []
    fline, toa, meta = _fit_local_linear_threshold_on_hist(
        ROOT, h, peak_frac, half_window, f"fit_local_{h.GetName()}"
    )
    y_upper = max(float(h.GetMaximum()) * 1.05, 1.0)
    h.GetYaxis().SetRangeUser(0.0, y_upper)
    h.Draw("HIST")
    if fline is not None:
        fline.Draw("same")
        keepalive.append(fline)
    # Highlight the local-linear fit window in red for quick visual inspection.
    if meta is not None:
        xfit_lo = float(meta["xmin"])
        xfit_hi = float(meta["xmax"])
        y_top = float(h.GetYaxis().GetXmax())
        l_lo = ROOT.TLine(xfit_lo, 0.0, xfit_lo, y_top)
        l_lo.SetLineColor(ROOT.kRed + 1)
        l_lo.SetLineStyle(2)
        l_lo.SetLineWidth(2)
        l_lo.Draw("same")
        keepalive.append(l_lo)
        l_hi = ROOT.TLine(xfit_hi, 0.0, xfit_hi, y_top)
        l_hi.SetLineColor(ROOT.kRed + 1)
        l_hi.SetLineStyle(2)
        l_hi.SetLineWidth(2)
        l_hi.Draw("same")
        keepalive.append(l_hi)
        l_span = ROOT.TLine(xfit_lo, y_top * 0.96, xfit_hi, y_top * 0.96)
        l_span.SetLineColor(ROOT.kRed + 1)
        l_span.SetLineStyle(1)
        l_span.SetLineWidth(3)
        l_span.Draw("same")
        keepalive.append(l_span)
    # Draw CFD threshold and selected CFD point for visual validation.
    if meta is not None:
        thr = float(meta["threshold"])
        xmin = float(h.GetXaxis().GetXmin())
        xmax = float(h.GetXaxis().GetXmax())
        l_thr = ROOT.TLine(xmin, thr, xmax, thr)
        l_thr.SetLineColor(ROOT.kGreen + 2)
        l_thr.SetLineStyle(2)
        l_thr.SetLineWidth(2)
        l_thr.Draw("same")
        keepalive.append(l_thr)
    if toa is not None and meta is not None:
        p_cfd = ROOT.TMarker(float(toa), float(meta["threshold"]), 20)
        p_cfd.SetMarkerColor(ROOT.kMagenta + 2)
        p_cfd.SetMarkerSize(1.2)
        p_cfd.Draw("same")
        keepalive.append(p_cfd)
    tx = ROOT.TLatex()
    tx.SetNDC()
    tx.SetTextFont(42)
    tx.SetTextSize(0.029)
    tx.SetTextAlign(31)  # right-aligned text block
    x_text = 0.88
    tx.DrawLatex(x_text, 0.72, f"{trig_label} waveform event={ev.event_id}, N={ev.nphoton}")
    if toa is not None:
        tx.DrawLatex(x_text, 0.66, f"Local-linear TOA ({peak_frac:.2f}*peak) = {toa:.4f} ns")
    if pulse.t10 is not None and pulse.t90 is not None:
        tx.DrawLatex(x_text, 0.60, f"rise t90-t10 = {pulse.t90 - pulse.t10:.4f} ns")
    if meta is not None:
        r2 = meta.get("r2", None)
        chi2 = float(meta["chi2"])
        ndf = int(meta["ndf"])
        pval = float(meta["pvalue"])
        if r2 is not None:
            tx.DrawLatex(x_text, 0.54, f"R^{{2}} = {float(r2):.4f}")
        if ndf > 0:
            tx.DrawLatex(x_text, 0.48, f"#chi^{{2}}/ndf = {chi2:.2f}/{ndf} = {chi2/ndf:.3f}")
        else:
            tx.DrawLatex(x_text, 0.48, f"#chi^{{2}}/ndf = {chi2:.2f}/0")
        tx.DrawLatex(x_text, 0.42, f"p-value = {pval:.3g}")
    keepalive.append(tx)
    # PyROOT 객체 소멸로 overlay가 사라지는 문제 방지.
    h._overlay_keepalive = keepalive


def draw_one(
    input_path: str,
    output_path: str,
    anchor_path: str | None,
    max_events: int,
    xmax: float,
    arrival_bin_ps: float,
    peak_frac: float,
    fit_half_window: int,
    response_rise_ns: float,
    response_fwhm_ns: float,
    transit_time_ns: float,
    sample_step_ns: float,
    gamma_shape: float,
    gamma_tau_ns: Optional[float],
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
    evs_t1_map = _event_map(evs_t1)
    evs_t2_map = _event_map(evs_t2)

    nvals_t1 = [ev.nphoton for ev in evs_t1 if ev.nphoton > 0]
    nvals_t2 = [ev.nphoton for ev in evs_t2 if ev.nphoton > 0]

    h_n1, f_n1, mpv1 = _build_nph_hist_and_landau_fit(
        ROOT, "hNphT1", nvals_t1, "T1 N photons / event;N photons;events"
    )
    h_n2, f_n2, mpv2 = _build_nph_hist_and_landau_fit(
        ROOT, "hNphT2", nvals_t2, "T2 N photons / event;N photons;events"
    )

    anchor_input = anchor_path if anchor_path else input_path
    anchor_evs_t1 = _scan_events(anchor_input, trig=0, max_events=max_events)
    anchor_nvals_t1 = [ev.nphoton for ev in anchor_evs_t1 if ev.nphoton > 0]
    _, _, anchor_mpv = _build_nph_hist_and_landau_fit(
        ROOT, "hNphAnchorT1_tmp", anchor_nvals_t1, "anchor"
    )
    if anchor_mpv is None:
        anchor_mpv = mpv1 if mpv1 is not None else 0.0

    ev_t1 = _pick_nearest_to_target(evs_t1, float(anchor_mpv))
    ev_t2 = _must_get_event(evs_t2_map, ev_t1.event_id, "T2")
    # T1/T2 event pairing rule for this dataset.
    ev_t1 = _must_get_event(evs_t1_map, ev_t1.event_id, "T1")

    h_t1_arrival, _ = _make_time_hist(
        ROOT,
        "h_t1_arrival_mpv",
        ev_t1,
        arrival_bin_ps,
        f"T1 selected event ({arrival_bin_ps:g} ps);time (ns);photons/bin",
    )
    h_t2_arrival, _ = _make_time_hist(
        ROOT,
        "h_t2_arrival_mpv",
        ev_t2,
        arrival_bin_ps,
        f"T2 selected event ({arrival_bin_ps:g} ps);time (ns);photons/bin",
    )
    pulse_t1 = _build_reconstructed_waveform(
        ev_t1,
        sample_step_ns=sample_step_ns,
        response_rise_ns=response_rise_ns,
        response_fwhm_ns=response_fwhm_ns,
        gamma_shape=gamma_shape,
        gamma_tau_ns=gamma_tau_ns,
        transit_time_ns=transit_time_ns,
        peak_frac=peak_frac,
    )
    pulse_t2 = _build_reconstructed_waveform(
        ev_t2,
        sample_step_ns=sample_step_ns,
        response_rise_ns=response_rise_ns,
        response_fwhm_ns=response_fwhm_ns,
        gamma_shape=gamma_shape,
        gamma_tau_ns=gamma_tau_ns,
        transit_time_ns=transit_time_ns,
        peak_frac=peak_frac,
    )
    h_wf_t1 = _make_hist_from_samples(
        ROOT,
        "h_t1_waveform_mpv",
        pulse_t1.time,
        pulse_t1.voltage,
        f"T1 reconstructed waveform ({sample_step_ns:.3f} ns);time (ns);a.u.",
    )
    h_wf_t2 = _make_hist_from_samples(
        ROOT,
        "h_t2_waveform_mpv",
        pulse_t2.time,
        pulse_t2.voltage,
        f"T2 reconstructed waveform ({sample_step_ns:.3f} ns);time (ns);a.u.",
    )

    c = ROOT.TCanvas("c_evt_mpv_matrix", "event optical MPV matrix", 1400, 1350)
    c.Divide(2, 3)

    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)
    _draw_top_panel(ROOT, h_n1, f_n1, mpv1, ev_t1, "T1", float(anchor_mpv))

    c.cd(2)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)
    _draw_top_panel(ROOT, h_n2, f_n2, mpv2, ev_t2, "T2", float(anchor_mpv))

    for pad_idx, h, trig_label, ev in (
        (3, h_t1_arrival, "T1", ev_t1),
        (4, h_t2_arrival, "T2", ev_t2),
    ):
        c.cd(pad_idx)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetBottomMargin(0.12)
        h.GetXaxis().SetRangeUser(0.0, float(xmax))
        _draw_arrival_panel(ROOT, h, trig_label, ev)

    for pad_idx, h, trig_label, ev, pulse in (
        (5, h_wf_t1, "T1", ev_t1, pulse_t1),
        (6, h_wf_t2, "T2", ev_t2, pulse_t2),
    ):
        c.cd(pad_idx)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetBottomMargin(0.12)
        wx0, wx1 = _waveform_xlim(pulse, xmax)
        h.GetXaxis().SetRangeUser(float(wx0), float(wx1))
        _draw_waveform_panel(ROOT, h, trig_label, ev, pulse, peak_frac, fit_half_window)

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    c.SaveAs(output_path)
    c.Close()

    tau_effective = float(gamma_tau_ns) if gamma_tau_ns is not None else _gamma_tau_from_rise(response_rise_ns, gamma_shape)
    print(
        f"입력: {input_path}\n"
        f"  anchor(LG T1) path={anchor_input}\n"
        f"  anchor(LG T1) MPV={float(anchor_mpv):.3f}\n"
        f"  T1 MPV={float(mpv1) if mpv1 is not None else float('nan'):.3f}, selected event={ev_t1.event_id}, N={ev_t1.nphoton}\n"
        f"  T2 MPV={float(mpv2) if mpv2 is not None else float('nan'):.3f}, paired event={ev_t2.event_id}, N={ev_t2.nphoton}\n"
        f"  waveform: rise={response_rise_ns:.3f} ns, fwhm={response_fwhm_ns:.3f} ns, transit={transit_time_ns:.3f} ns, dt={sample_step_ns:.3f} ns, shape={gamma_shape:.6f}, tau={tau_effective:.6f} ns\n"
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
    parser.add_argument("--arrival-bin-ps", type=float, default=200.0, help="2행 photon-arrival bin (ps)")
    parser.add_argument("--peak-frac", type=float, default=0.3, help="국소 선형 TOA 기준 peak 비율")
    parser.add_argument("--fit-half-window", type=int, default=2, help="교차점 주변 반쪽 윈도우(bin)")
    parser.add_argument(
        "--anchor-root",
        default=None,
        help="선택 기준 anchor MPV를 계산할 ROOT 파일 (기본: 입력 파일 자체)",
    )
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
    if args.arrival_bin_ps <= 0.0:
        raise SystemExit("--arrival-bin-ps는 0보다 커야 합니다.")
    if not (0.0 < args.peak_frac < 1.0):
        raise SystemExit("--peak-frac는 0~1 사이여야 합니다.")
    if args.fit_half_window < 1:
        raise SystemExit("--fit-half-window는 1 이상이어야 합니다.")
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
    anchor_path = _resolve(args.anchor_root, "analysis/t_res/data/v2_60GeV_e-_LG_3000_0.root") if args.anchor_root else None
    spec_csv_path = _resolve(args.pmt_spec_csv, "analysis/reference/pmt_r2076/R2076_pmt_spec_reference.csv")
    spec_rise_ns, spec_transit_ns = _load_r2076_time_response(spec_csv_path)

    response_rise_ns = float(args.response_rise_ns) if args.response_rise_ns is not None else float(spec_rise_ns if spec_rise_ns is not None else 1.0)
    response_fwhm_ns = float(args.response_fwhm_ns) if args.response_fwhm_ns is not None else 3.0
    transit_time_ns = float(args.transit_time_ns) if args.transit_time_ns is not None else float(spec_transit_ns if spec_transit_ns is not None else 0.0)

    print(
        "[PMT spec] "
        f"csv={spec_csv_path}, "
        f"rise(ns)={'cli' if args.response_rise_ns is not None else ('spec' if spec_rise_ns is not None else 'default')}:{response_rise_ns:.3f}, "
        f"transit(ns)={'cli' if args.transit_time_ns is not None else ('spec' if spec_transit_ns is not None else 'default')}:{transit_time_ns:.3f}, "
        f"fwhm(ns)={'cli' if args.response_fwhm_ns is not None else 'default'}:{response_fwhm_ns:.3f}"
    )

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
        anchor_path=anchor_path,
        max_events=args.max_events,
        xmax=args.xmax,
        arrival_bin_ps=args.arrival_bin_ps,
        peak_frac=args.peak_frac,
        fit_half_window=args.fit_half_window,
        response_rise_ns=response_rise_ns,
        response_fwhm_ns=response_fwhm_ns,
        transit_time_ns=transit_time_ns,
        sample_step_ns=args.sample_step_ns,
        gamma_shape=args.gamma_shape,
        gamma_tau_ns=args.gamma_tau_ns,
    )


if __name__ == "__main__":
    main()
