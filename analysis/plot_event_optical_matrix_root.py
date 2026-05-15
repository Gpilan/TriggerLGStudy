#!/usr/bin/env python3
"""
단일 이벤트 광자 시간분포를 2x2 매트릭스로 그립니다 (PyROOT).

행/열 구성:
  [T1(10 ps),  T2(10 ps)]
  [T1(100 ps), T2(100 ps)]

이벤트 선택:
  - T1: 목표 포톤수(--target-t1, 기본 2200)에 가장 가까운 이벤트
  - T2: 목표 포톤수(--target-t2, 기본 2100)에 가장 가까운 이벤트
"""

from __future__ import annotations

import argparse
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
class SelectedEvent:
    event_id: int
    nphoton: int
    lo: list[float]
    hi: list[float]
    counts: list[int]


def _sum_counts(counts: list[int]) -> int:
    return int(sum(int(c) for c in counts)) if counts else 0


def _pick_event_near_target(path: str, trig: int, target_n: int, max_events: int) -> SelectedEvent:
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

    best: SelectedEvent | None = None
    best_diff = 10**18

    for i in range(nmax):
        tree.GetEntry(i)
        lo, hi, counts = _merged_from_tower(evt, trig)
        nph = _sum_counts(counts)
        if nph <= 0:
            continue
        diff = abs(nph - target_n)
        if diff < best_diff:
            best_diff = diff
            best = SelectedEvent(
                event_id=i,
                nphoton=nph,
                lo=list(lo),
                hi=list(hi),
                counts=[int(x) for x in counts],
            )

    f.Close()
    if best is None:
        raise RuntimeError(f"trigger {trig}에서 유효 이벤트를 찾지 못함")
    return best


def _make_hist_from_rebin(ROOT, name: str, ev: SelectedEvent, bin_ps: float, title: str):
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
    # ROOT stat box "Entries"를 실제 포톤 합으로 보이게 맞춤.
    h.SetEntries(total_photons)
    return h, total_photons


def _pad_hist_for_rebin(ROOT, h, name: str, factor: int, title: str):
    """Rebin factor로 정확히 나눠지도록 오른쪽 0-bin 패딩."""
    nb = int(h.GetNbinsX())
    if factor <= 1 or (nb % factor == 0):
        return h

    nb_pad = ((nb + factor - 1) // factor) * factor
    xmin = float(h.GetXaxis().GetXmin())
    xmax = float(h.GetXaxis().GetXmax())
    bw = (xmax - xmin) / float(nb)
    xmax_pad = xmin + bw * float(nb_pad)

    hp = ROOT.TH1F(name, title, nb_pad, xmin, xmax_pad)
    hp.Sumw2()
    hp.SetDirectory(0)
    hp.SetLineColor(ROOT.kBlack)
    hp.SetLineWidth(2)

    for i in range(1, nb + 1):
        hp.SetBinContent(i, float(h.GetBinContent(i)))
        hp.SetBinError(i, float(h.GetBinError(i)))
    hp.SetEntries(float(h.GetEntries()))
    return hp


def main() -> None:
    parser = argparse.ArgumentParser(
        description="단일 이벤트 옵티컬 포톤 시간분포 2x2 (10ps/100ps, T1/T2)"
    )
    parser.add_argument(
        "input",
        nargs="?",
        default=None,
        help="입력 .root (기본: analysis/t_res/data/v2_60GeV_e-_LG_3000_0.root)",
    )
    parser.add_argument("--target-t1", type=int, default=2200, help="T1 목표 포톤수")
    parser.add_argument("--target-t2", type=int, default=2100, help="T2 목표 포톤수")
    parser.add_argument("--max-events", type=int, default=-1)
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="출력 PNG (기본: figures/event_optical_matrix_<stem>.png)",
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

    in_path = _resolve(args.input, "analysis/t_res/data/v2_60GeV_e-_LG_3000_0.root")

    _load_rootio(args.rootio_lib, repo_root)
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(1110)
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetHistLineWidth(2)
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")

    ev_t1 = _pick_event_near_target(
        in_path, trig=0, target_n=args.target_t1, max_events=args.max_events
    )
    ev_t2 = _pick_event_near_target(
        in_path, trig=1, target_n=args.target_t2, max_events=args.max_events
    )

    stem = os.path.splitext(os.path.basename(in_path))[0]
    if args.output:
        out = os.path.abspath(args.output)
        if not os.path.splitext(out)[1]:
            out += ".png"
    else:
        out = os.path.join(_figures_dir(repo_root), f"event_optical_matrix_{stem}.png")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    h_t1_10, t1_total = _make_hist_from_rebin(
        ROOT,
        "h_t1_10",
        ev_t1,
        10.0,
        f"T1 (10 ps);time (ns);photons/bin",
    )
    h_t2_10, t2_total = _make_hist_from_rebin(
        ROOT,
        "h_t2_10",
        ev_t2,
        10.0,
        f"T2 (10 ps);time (ns);photons/bin",
    )

    # 100 ps 히스토는 10 ps 히스토를 ROOT Rebin(10)으로 생성:
    # "10개(10 ps) -> 1개(100 ps)"를 코드 경로에서 명시적으로 보장.
    h_t1_10_pad = _pad_hist_for_rebin(
        ROOT, h_t1_10, "h_t1_10_pad", 10, "T1 (10 ps, padded for Rebin);time (ns);photons/bin"
    )
    h_t2_10_pad = _pad_hist_for_rebin(
        ROOT, h_t2_10, "h_t2_10_pad", 10, "T2 (10 ps, padded for Rebin);time (ns);photons/bin"
    )

    h_t1_100 = h_t1_10_pad.Rebin(10, "h_t1_100")
    h_t2_100 = h_t2_10_pad.Rebin(10, "h_t2_100")
    h_t1_100.SetDirectory(0)
    h_t2_100.SetDirectory(0)
    h_t1_100.SetTitle("T1 (100 ps, from 10 ps Rebin);time (ns);photons/bin")
    h_t2_100.SetTitle("T2 (100 ps, from 10 ps Rebin);time (ns);photons/bin")
    h_t1_100.SetLineColor(ROOT.kBlack)
    h_t2_100.SetLineColor(ROOT.kBlack)
    h_t1_100.SetLineWidth(2)
    h_t2_100.SetLineWidth(2)
    # Rebin 이후에도 Entries(포톤 합) 일관성 유지
    h_t1_100.SetEntries(t1_total)
    h_t2_100.SetEntries(t2_total)

    c = ROOT.TCanvas("c_evt_matrix", "event optical matrix", 1400, 900)
    c.Divide(2, 2)

    pad_specs = [
        (1, h_t1_10, f"T1 event={ev_t1.event_id}, N={ev_t1.nphoton}"),
        (2, h_t2_10, f"T2 event={ev_t2.event_id}, N={ev_t2.nphoton}"),
        (3, h_t1_100, f"T1 event={ev_t1.event_id}, N={ev_t1.nphoton}"),
        (4, h_t2_100, f"T2 event={ev_t2.event_id}, N={ev_t2.nphoton}"),
    ]
    for pad_idx, h, caption in pad_specs:
        c.cd(pad_idx)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetBottomMargin(0.12)
        h.Draw("HIST")
        tx = ROOT.TLatex()
        tx.SetNDC()
        tx.SetTextFont(42)
        tx.SetTextSize(0.034)
        tx.DrawLatex(0.14, 0.94, caption)

    c.SaveAs(out)
    c.Close()

    print(
        f"입력: {in_path}\n"
        f"  T1 선택: event={ev_t1.event_id}, N={ev_t1.nphoton} (target={args.target_t1})\n"
        f"  T2 선택: event={ev_t2.event_id}, N={ev_t2.nphoton} (target={args.target_t2})\n"
        f"PNG 저장: {out}"
    )


if __name__ == "__main__":
    main()
