#!/usr/bin/env python3
"""
고광자 수 이벤트만 골라, merged SiPM **시간 빈**을 **균일 시간축으로 리빈**한 뒤 합산하고 ROOT로 그립니다.
(ROOT에 기록된 merged 빈은 이벤트마다 ``비어 있지 않은 빈``만 담겨 길이가 다름 → 직접 합산 불가.)

수직선: ``plot_trigger_timing._cfd_absolute_time`` (기본 F=0.3).

- 선별: ``--combine sum`` 이면 트리거별 상위 ``--top-k`` **합산**; ``--combine single``(기본)이면
  ``N_{T1}+N_{T2}`` 가 큰 순으로 ``--n-events`` 개(기본 **1**) **각각** 표시(N열×2행: T1 위·T2 아래).
- 시간축: **top-k** 이벤트의 시간 구간 + 여백 후 균일 분할. 빈 폭은 ``--target-nbins`` 로
  자동(nice step)하거나 ``--bin-width-ns`` 로 고정.
- 리빈: sparse 빈 ``[lo,hi]`` 를 균일 축에 **길이 비율로 분배**(중심만 넣는 것보다 매끈함).
- 스타일: 연한 채움 + 외곽선, 그리드, 여백·폰트 약간 키움.

필요: PyROOT, numpy, ``build/rootIO/librootIO.so``

예:
  source envset.sh
  export LD_LIBRARY_PATH=$PWD/build/rootIO:$LD_LIBRARY_PATH
  python3 analysis/plot_cfd_high_nphoton_timing_ROOT.py
  python3 analysis/plot_cfd_high_nphoton_timing_ROOT.py --combine sum --top-k 20
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

from plot_trigger_timing import (
    _cfd_absolute_time,
    _figures_dir,
    _find_repo_root,
    _load_rootio,
    _merged_from_tower,
    _prepend_build_rootio_ld_path,
    _resolve_input_path,
)


def _nphoton(lo: list[float], hi: list[float], cnt: list[int]) -> int:
    if not cnt:
        return 0
    return int(sum(int(c) for c in cnt))


def _collect_ranked(
    path: str,
    trig: int,
    *,
    max_events: int,
) -> list[tuple[int, int, list[float], list[float], list[int]]]:
    """(N, evt_idx, lo, hi, c) 리스트, N 내림차순으로 정렬해 반환."""
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

    rows: list[tuple[int, int, list[float], list[float], list[int]]] = []

    for i in range(nmax):
        tree.GetEntry(i)
        lo, hi, c = _merged_from_tower(evt, trig)
        if not c or not lo or len(lo) != len(hi) or len(lo) != len(c):
            continue
        nph = _nphoton(lo, hi, c)
        rows.append((nph, i, list(lo), list(hi), list(int(x) for x in c)))

    f.Close()
    rows.sort(key=lambda r: r[0], reverse=True)
    return rows


def _time_span_topk(
    ranked: list[tuple[int, int, list[float], list[float], list[int]]],
    top_k: int,
    *,
    margin_frac: float,
) -> tuple[float, float]:
    k = min(max(top_k, 1), len(ranked))
    t_lo = min(float(ranked[i][2][0]) for i in range(k))
    t_hi = max(float(ranked[i][3][-1]) for i in range(k))
    span = max(t_hi - t_lo, 1e-9)
    m = margin_frac * span
    return t_lo - m, t_hi + m


def _nice_bin_width(span_ns: float, target_nbins: int) -> float:
    """목표 빈 수에 가깝도록 1–2–5×10^n 형태의 빈 폭(ns)."""
    n = max(int(target_nbins), 8)
    g = max(span_ns / float(n), 1e-15)
    exp = math.floor(math.log10(g))
    base = 10.0**exp
    for f in (1.0, 2.0, 5.0, 10.0):
        bw = f * base
        if span_ns / bw <= n * 1.15:
            return bw
    return 10.0 * base


def _rebin_sparse_to_edges(
    lo: list[float],
    hi: list[float],
    c: list[int],
    edges: np.ndarray,
) -> np.ndarray:
    """단일 sparse 스펙트럼을 균일 edges에 겹침 비율로 분배."""
    nb = int(edges.size) - 1
    acc = np.zeros(nb, dtype=np.float64)
    if not c or not lo:
        return acc
    e_lo = edges[:-1]
    e_hi = edges[1:]
    for j in range(len(c)):
        w = float(c[j])
        if w <= 0.0:
            continue
        t0, t1 = float(lo[j]), float(hi[j])
        width = max(t1 - t0, 1e-15)
        for i in range(nb):
            ov = min(t1, float(e_hi[i])) - max(t0, float(e_lo[i]))
            if ov > 0.0:
                acc[i] += w * (ov / width)
    return acc


def _rebin_sparse_to_edges_integer(
    lo: list[float],
    hi: list[float],
    c: list[int],
    edges: np.ndarray,
) -> np.ndarray:
    """단일 이벤트용: 각 sparse 빈 카운트를 중심 빈에 통째로 넣어 정수 카운트 유지."""
    nb = int(edges.size) - 1
    acc = np.zeros(nb, dtype=np.int64)
    if not c or not lo:
        return acc.astype(np.float64)
    for j in range(len(c)):
        w = int(c[j])
        if w <= 0:
            continue
        cen = 0.5 * (float(lo[j]) + float(hi[j]))
        idx = int(np.searchsorted(edges, cen, side="right") - 1)
        if 0 <= idx < nb:
            acc[idx] += w
    return acc.astype(np.float64)


def _rebin_sum_topk(
    ranked: list[tuple[int, int, list[float], list[float], list[int]]],
    top_k: int,
    edges: np.ndarray,
) -> np.ndarray:
    """상위 top_k 이벤트의 스펙트럼 합."""
    acc = np.zeros(int(edges.size) - 1, dtype=np.float64)
    k = min(top_k, len(ranked))
    for row in ranked[:k]:
        lo, hi, c = row[2], row[3], row[4]
        acc += _rebin_sparse_to_edges(lo, hi, c, edges)
    return acc


def _scan_events_by_total_nphoton(
    path: str,
    *,
    max_events: int,
) -> list[tuple[int, int, int, int]]:
    """(N_T1+N_T2, evt_idx, N_T1, N_T2) 내림차순."""
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

    rows: list[tuple[int, int, int, int]] = []

    for i in range(nmax):
        tree.GetEntry(i)
        lo1, hi1, c1 = _merged_from_tower(evt, 0)
        lo2, hi2, c2 = _merged_from_tower(evt, 1)
        n1 = _nphoton(lo1, hi1, c1) if c1 else 0
        n2 = _nphoton(lo2, hi2, c2) if c2 else 0
        rows.append((n1 + n2, i, n1, n2))

    f.Close()
    rows.sort(key=lambda r: r[0], reverse=True)
    return rows


def _time_span_one_event(
    lo: list[float],
    hi: list[float],
    *,
    margin_frac: float,
) -> tuple[float, float]:
    if not lo or not hi:
        return 0.0, 1.0
    t_lo = float(lo[0])
    t_hi = float(hi[-1])
    span = max(t_hi - t_lo, 1e-9)
    m = margin_frac * span
    return t_lo - m, t_hi + m


def _y_range_integer_with_headroom(
    h: "ROOT.TH1D", *, overflow_frac: float = 0.1, overflow_min: int = 3
) -> float:
    """
    y 최댓값을 정수로 올리고, 위에 overflow(여유) 몇 단위를 더해 CFD 선·막대가 안 잘리게 함.
    반환: 설정한 y_max.
    """
    ymax = float(max(h.GetMaximum(), 0.0))
    base = int(math.ceil(ymax))
    extra = max(overflow_min, int(math.ceil(overflow_frac * max(base, 1))))
    n_major = 6
    step = max(1, int(math.ceil((base + extra) / n_major)))
    y_max_int = int(step * n_major)
    y_max = float(y_max_int)
    h.SetMinimum(0.0)
    h.SetMaximum(y_max)
    ya = h.GetYaxis()
    try:
        ya.SetDecimals(0)
    except Exception:
        pass
    # 가능한 경우(작은 범위) y축 단위를 1로 고정.
    # ROOT의 ndivisions 포맷: 100*n2 + n1 (n1: major ticks)
    if y_max_int <= 60:
        ya.SetNdivisions(y_max_int, 0, 0, False)
    else:
        # 큰 범위에서는 과도한 눈금 방지
        ya.SetNdivisions(100 * n_major + 5, 0)
    return y_max


def _make_th1d_uniform(
    name: str,
    title: str,
    t_lo: float,
    t_hi: float,
    counts: np.ndarray,
) -> "ROOT.TH1D":
    import ROOT

    nb = int(counts.size)
    if nb < 1 or t_hi <= t_lo:
        raise ValueError("invalid uniform histogram")
    h = ROOT.TH1D(name, title, nb, float(t_lo), float(t_hi))
    h.SetDirectory(0)
    for j in range(nb):
        h.SetBinContent(j + 1, float(counts[j]))
    return h


def _apply_plot_style(h: "ROOT.TH1D", *, pad: "ROOT.TPad", line_color: int) -> None:
    import ROOT

    h.SetLineColor(line_color)
    h.SetLineWidth(2)
    h.SetFillColorAlpha(line_color, 0.22)
    h.GetXaxis().SetTitleSize(0.048)
    h.GetYaxis().SetTitleSize(0.048)
    h.GetYaxis().SetTitle("")  # y축 제목 숨김 (틱 숫자는 정수)
    h.GetXaxis().SetLabelSize(0.042)
    h.GetYaxis().SetLabelSize(0.042)
    h.GetYaxis().SetTitleOffset(1.25)
    h.GetXaxis().SetTitleOffset(1.05)
    h.GetXaxis().SetNdivisions(505, 0)
    pad.SetGrid(1, 1)
    pad.SetLeftMargin(0.12)
    pad.SetRightMargin(0.05)
    pad.SetBottomMargin(0.14)
    pad.SetTopMargin(0.07)


def _draw_one_pad(
    pad: "ROOT.TPad",
    h: "ROOT.TH1D",
    cfd_t: float | None,
    cfd_fraction: float,
    *,
    line_color: int,
    lines: list[str],
    x_max_override: float | None = None,
    y_max_override: float | None = None,
) -> None:
    import ROOT

    pad.cd()
    _apply_plot_style(h, pad=pad, line_color=line_color)
    if x_max_override is not None:
        x_lo = h.GetXaxis().GetXmin()
        x_hi = h.GetXaxis().GetXmax()
        if x_max_override > x_lo:
            h.GetXaxis().SetRangeUser(x_lo, min(x_max_override, x_hi))
    if y_max_override is not None and y_max_override > 0:
        h.SetMinimum(0.0)
        h.SetMaximum(float(y_max_override))
        y_max = float(y_max_override)
    else:
        y_max = _y_range_integer_with_headroom(h)
    h.Draw("HIST F")  # 연한 채움 + 막대

    ymax = float(y_max)
    if ymax <= 0.0:
        ymax = 1.0
    if cfd_t is not None:
        ln = ROOT.TLine(float(cfd_t), 0.0, float(cfd_t), ymax)
        ln.SetLineColor(ROOT.kRed + 1)
        ln.SetLineWidth(2)
        ln.SetLineStyle(2)
        ln.Draw()

    if cfd_t is not None:
        lines = list(lines) + [
            f"CFD {cfd_fraction:g} #times max(bin): {cfd_t:.4f} ns"
        ]
    else:
        lines = list(lines) + ["CFD: undefined (empty / zero max)"]

    pt = ROOT.TPaveText(0.46, 0.58, 0.985, 0.98, "NDC BR")
    pt.SetFillColorAlpha(ROOT.kWhite, 0.92)
    pt.SetBorderSize(1)
    pt.SetCornerRadius(0.02)
    pt.SetTextFont(42)
    pt.SetTextSize(0.032)
    for line in lines:
        pt.AddText(line)
    pt.Draw()


def _edges_from_span(
    t_lo: float,
    t_hi: float,
    args: argparse.Namespace,
) -> tuple[np.ndarray, float, int, str]:
    span = max(t_hi - t_lo, 1e-9)
    if args.bin_width_ns is not None:
        bw = max(float(args.bin_width_ns), 1e-12)
        mode = "fixed width"
    else:
        bw = _nice_bin_width(span, max(args.target_nbins, 8))
        mode = f"auto (~{args.target_nbins} bins, nice step)"
    nbins = max(1, min(50000, int(math.ceil(span / bw))))
    t_hi = t_lo + nbins * bw
    edges = np.linspace(t_lo, t_hi, nbins + 1, dtype=np.float64)
    return edges, bw, nbins, mode


def _plot_combine_sum(
    args: argparse.Namespace,
    in_path: str,
    out: str,
    ROOT: object,
) -> None:
    colors = (ROOT.kAzure + 2, ROOT.kOrange + 7)
    # 1×2 패드가 지나치게 가로로 길지 않게 (가로≈세로에 가깝게)
    c = ROOT.TCanvas("c", "high-N photon timing + CFD (sum)", 980, 520)
    c.Divide(2, 1)

    for pad_idx, (trig, tlabel) in enumerate(((0, "T1"), (1, "T2")), start=1):
        ranked = _collect_ranked(in_path, trig, max_events=args.max_events)
        if not ranked:
            c.cd(pad_idx)
            t = ROOT.TLatex(0.15, 0.5, f"{tlabel}: no merged time data")
            t.SetNDC()
            t.SetTextSize(0.06)
            t.Draw()
            continue

        k_use = min(args.top_k, len(ranked))
        t_lo, t_hi = _time_span_topk(ranked, k_use, margin_frac=args.time_margin_frac)
        edges, bw, nbins, mode = _edges_from_span(t_lo, t_hi, args)
        acc = _rebin_sum_topk(ranked, k_use, edges)
        used = [ranked[i][1] for i in range(k_use)]

        nb = int(edges.size) - 1
        lows = [float(edges[i]) for i in range(nb)]
        highs = [float(edges[i + 1]) for i in range(nb)]
        counts = [int(round(x)) for x in acc.tolist()]
        cfd_t = _cfd_absolute_time(
            lows, highs, counts, args.cfd_fraction, args.ns_per_unit
        )
        n_tot = float(np.sum(acc))

        rebin_note = f"{mode}: {nbins} bins × {bw:.4g} ns"
        h = _make_th1d_uniform(
            f"h_{tlabel}",
            f"{tlabel};time (ns);optical photons (sum, top-{k_use} evt)",
            float(edges[0]),
            float(edges[-1]),
            acc,
        )

        c.cd(pad_idx)
        lines = [
            f"{tlabel}  sum of top-{k_use} by N_opt",
            rebin_note,
            f"#Sigma photons #approx {n_tot:.0f}",
            f"events: {used[:8]}{'...' if len(used) > 8 else ''}",
        ]
        _draw_one_pad(
            ROOT.gPad,
            h,
            cfd_t,
            args.cfd_fraction,
            line_color=colors[pad_idx - 1],
            lines=lines,
        )

    c.SaveAs(out)


def _plot_combine_single(
    args: argparse.Namespace,
    in_path: str,
    out: str,
    ROOT: object,
) -> None:
    ranked_totals_all = _scan_events_by_total_nphoton(
        in_path, max_events=args.max_events
    )
    if not ranked_totals_all:
        raise SystemExit("이벤트가 없습니다.")

    # 발표용 single 이벤트는 두 트리거(T1/T2)가 모두 보이는 케이스를 우선 사용.
    ranked_totals = [r for r in ranked_totals_all if r[2] > 0 and r[3] > 0]
    if not ranked_totals:
        ranked_totals = ranked_totals_all

    nev = min(max(args.n_events, 1), len(ranked_totals))
    ncols = nev
    cw = 420 * ncols
    if ncols == 1:
        cw = max(cw, 760)
    cw = min(int(cw), 3600)
    # 2행 패드가 한 화면에 함께 보이도록 세로 비율 완화
    ch = int(max(700, min(1050, 1.2 * cw)))
    c = ROOT.TCanvas(
        "c",
        "per-event high-N timing + CFD",
        cw,
        ch,
    )
    c.Divide(ncols, 2)

    palette = (ROOT.kAzure + 2, ROOT.kSpring + 5, ROOT.kOrange + 7, ROOT.kViolet + 1)

    f = ROOT.TFile.Open(in_path)
    if not f or f.IsZombie():
        raise FileNotFoundError(in_path)
    tree = f.Get("CBDsim")
    if not tree:
        f.Close()
        raise RuntimeError(f"트리 CBDsim 없음: {in_path}")
    evt = ROOT.CBDsimInterface.CBDsimEventData()
    tree.SetBranchAddress("CBDsimEventData", evt)

    selected_info: list[str] = []
    for col in range(ncols):
        _ncombo, evt_idx, n1, n2 = ranked_totals[col]
        selected_info.append(f"evt={evt_idx} (N_T1={n1}, N_T2={n2}, total={_ncombo})")
        for row, (trig, tlabel) in enumerate(((0, "T1"), (1, "T2"))):
            pad_num = col + 1 + row * ncols
            c.cd(pad_num)
            tree.GetEntry(evt_idx)
            lo, hi, cc = _merged_from_tower(evt, trig)
            n_this = n1 if trig == 0 else n2

            if not cc or not lo:
                lab = ROOT.TLatex(
                    0.12,
                    0.55,
                    f"{tlabel} evt {evt_idx}: no merged hits",
                )
                lab.SetNDC()
                lab.SetTextSize(0.06)
                lab.Draw()
                continue

            t_lo, t_hi = _time_span_one_event(
                lo, hi, margin_frac=args.time_margin_frac
            )
            edges, bw, nbins, mode = _edges_from_span(t_lo, t_hi, args)
            acc = _rebin_sparse_to_edges_integer(lo, hi, cc, edges)

            nb = int(edges.size) - 1
            lows = [float(edges[i]) for i in range(nb)]
            highs = [float(edges[i + 1]) for i in range(nb)]
            counts = [int(round(x)) for x in acc.tolist()]
            cfd_t = _cfd_absolute_time(
                lows, highs, counts, args.cfd_fraction, args.ns_per_unit
            )

            rebin_note = f"{mode}: {nbins} bins × {bw:.4g} ns"
            h = _make_th1d_uniform(
                f"h_{tlabel}_ev{evt_idx}",
                f"{tlabel};time (ns);optical photons / bin",
                float(edges[0]),
                float(edges[-1]),
                acc,
            )

            lines = [
                f"{tlabel}  evt {evt_idx}  (#{col + 1} by N_{{T1}}+N_{{T2}})",
                f"N_opt({tlabel}) = {n_this},  N_T1+N_T2 = {_ncombo}",
                rebin_note,
            ]
            _draw_one_pad(
                ROOT.gPad,
                h,
                cfd_t,
                args.cfd_fraction,
                line_color=palette[col % len(palette)],
                lines=lines,
                x_max_override=args.x_max_single,
                y_max_override=args.y_max_single,
            )

    f.Close()
    c.SaveAs(out)
    if selected_info:
        print("single mode selected events:")
        for s in selected_info:
            print("  " + s)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="고 N_opt 이벤트 → 균일 리빈 합산 스펙트럼 + CFD (ROOT)"
    )
    parser.add_argument(
        "input",
        nargs="?",
        default=None,
        help="입력 .root (기본: analysis/t_res/data/LG_60GeV_e-_LG_1.root)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="출력 PNG (기본: figures/cfd_highN_timing_<stem>_ROOT.png)",
    )
    parser.add_argument(
        "--combine",
        choices=("single", "sum"),
        default="single",
        help="single: N_T1+N_T2 상위 --n-events 개를 이벤트별로; sum: 트리거별 top-k 합산",
    )
    parser.add_argument(
        "--n-events",
        type=int,
        default=1,
        dest="n_events",
        help="--combine single 일 때 표시할 이벤트 개수 (기본 1)",
    )
    parser.add_argument("--top-k", type=int, default=20, dest="top_k", help="--combine sum: 트리거별 상위 K")
    parser.add_argument(
        "--bin-width-ns",
        type=float,
        default=None,
        dest="bin_width_ns",
        help="균일 리빈 빈 폭 (ns). 생략 시 --target-nbins 로 자동(nice step)",
    )
    parser.add_argument(
        "--target-nbins",
        type=int,
        default=280,
        dest="target_nbins",
        help="자동 빈 폭일 때 목표 빈 개수 (기본 280)",
    )
    parser.add_argument(
        "--time-margin-frac",
        type=float,
        default=0.05,
        dest="time_margin_frac",
        help="top-k 이벤트 시간 구간에 붙일 상대 여백 (기본 5%%)",
    )
    parser.add_argument(
        "--cfd-fraction",
        type=float,
        default=0.3,
        dest="cfd_fraction",
        help="F in F×max(bin) (기본 0.3)",
    )
    parser.add_argument("--ns-per-unit", type=float, default=1.0, dest="ns_per_unit")
    parser.add_argument(
        "--x-max-single",
        type=float,
        default=15.0,
        dest="x_max_single",
        help="single 모드 x축 상한 (기본 15 ns)",
    )
    parser.add_argument(
        "--y-max-single",
        type=float,
        default=12.0,
        dest="y_max_single",
        help="single 모드 y축 상한 (기본 12)",
    )
    parser.add_argument("--max-events", type=int, default=-1)
    parser.add_argument("-l", "--rootio-lib", default=None)
    args = parser.parse_args()

    repo_root = _find_repo_root()
    if repo_root is None:
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    _prepend_build_rootio_ld_path(repo_root)

    default_rel = os.path.join("analysis", "t_res", "data", "LG_60GeV_e-_LG_1.root")
    if args.input:
        in_path = _resolve_input_path(args.input, repo_root)
        if not in_path:
            raise SystemExit(f"파일 없음: {args.input}")
    else:
        full = os.path.join(repo_root, default_rel)
        if not os.path.isfile(full):
            raise SystemExit(f"기본 경로에 파일 없음: {full}")
        in_path = os.path.abspath(full)

    _load_rootio(args.rootio_lib, repo_root)
    import ROOT

    ROOT.gROOT.SetBatch(True)

    stem = os.path.splitext(os.path.basename(in_path))[0]
    fig_dir = _figures_dir(repo_root)
    if args.output:
        out = os.path.abspath(args.output)
        if not os.path.splitext(out)[1]:
            out += ".png"
    elif args.combine == "single":
        out = os.path.join(
            fig_dir, f"cfd_highN_timing_{stem}_perEvt{args.n_events}_ROOT.png"
        )
    else:
        out = os.path.join(fig_dir, f"cfd_highN_timing_{stem}_ROOT.png")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    if args.combine == "sum":
        _plot_combine_sum(args, in_path, out, ROOT)
    else:
        _plot_combine_single(args, in_path, out, ROOT)
    print(f"입력: {in_path}")
    print(f"저장: {out}")


if __name__ == "__main__":
    main()
