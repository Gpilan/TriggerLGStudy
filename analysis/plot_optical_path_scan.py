#!/usr/bin/env python3
"""Parse CBDsim optical path histogram files and plot trigger-size trends."""

from __future__ import annotations

import argparse
import csv
import os
import re
from dataclasses import dataclass, field

import matplotlib.pyplot as plt
import numpy as np

FATE_ORDER = [
    "detected",
    "qe_reject",
    "absorb_scint",
    "absorb_lg",
    "absorb_glass",
    "absorb_si",
    "absorb_other",
    "kill_world",
    "kill_other",
]

FATE_COLORS = {
    "detected": "#2ca02c",
    "qe_reject": "#98df8a",
    "absorb_scint": "#ff7f0e",
    "absorb_lg": "#1f77b4",
    "absorb_glass": "#aec7e8",
    "absorb_si": "#9467bd",
    "absorb_other": "#c5b0d5",
    "kill_world": "#d62728",
    "kill_other": "#8c564b",
}

SIZE_FROM_TAG = {
    "1x1": 1.0,
    "1p5x1p5": 1.5,
    "2x2": 2.0,
    "2p5x2p5": 2.5,
    "3x3": 3.0,
    "3p5x3p5": 3.5,
    "4x4": 4.0,
    "4x4_nolg": 4.0,
}

TAG_DISPLAY = {
    "1x1": "1×1",
    "1p5x1p5": "1.5×1.5",
    "2x2": "2×2",
    "2p5x2p5": "2.5×2.5",
    "3x3": "3×3",
    "3p5x3p5": "3.5×3.5",
    "4x4": "4×4",
    "4x4_nolg": "4×4 no-LG",
}

LEGEND_LABELS = {
    "4x4": "4x4 LG",
    "4x4_nolg": "4x4 no-LG",
}

TRIGGER_XLABEL = "trigger size (scintillator x half-width)"


def _tag_label_root(tag: str) -> str:
    """ASCII size label for PyROOT (Unicode multiplication sign renders broken)."""
    if tag in LEGEND_LABELS:
        return LEGEND_LABELS[tag]
    return tag.replace("p", ".")


def _trigger_axis_labels(hists: list[PathHistFile]) -> tuple[list[float], list[str]]:
    hists = sorted(hists, key=lambda h: h.trigger_size)
    sizes = [h.trigger_size for h in hists]
    tick_labels = [
        f"{TAG_DISPLAY.get(h.tag, h.tag)}\n({h.trigger_size * 5.0:g} mm)" for h in hists
    ]
    return sizes, tick_labels


def _style_trigger_axis(ax, hists: list[PathHistFile]) -> list[float]:
    sizes, tick_labels = _trigger_axis_labels(hists)
    ax.set_xticks(sizes)
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel(TRIGGER_XLABEL)
    return sizes


@dataclass
class PathHistFile:
    tag: str
    trigger_size: float
    nbins: int
    xmax_mm: float
    total_tracks: int
    mean_path_mm_all: float
    total_counts: np.ndarray
    fates: dict[str, dict] = field(default_factory=dict)


def _is_size_scan_tag(tag: str) -> bool:
    return tag in SIZE_FROM_TAG and not tag.endswith("_nolg")


def _parse_tag(path: str) -> str:
    base = os.path.basename(path)
    m = re.match(r"path_(.+)_100ev\.path_hist\.txt$", base)
    if not m:
        raise ValueError(f"unexpected hist filename: {path}")
    return m.group(1)


def load_path_hist(path: str) -> PathHistFile:
    tag = _parse_tag(path)
    meta: dict[str, str] = {}
    total_counts: list[int] | None = None
    fates: dict[str, dict] = {}
    pending_fate: str | None = None
    expect_total_counts = False

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("nbins "):
                meta["nbins"] = line.split()[1]
            elif line.startswith("xmax_mm "):
                meta["xmax_mm"] = line.split()[1]
            elif line.startswith("total_tracks "):
                meta["total_tracks"] = line.split()[1]
            elif line.startswith("mean_path_mm_all "):
                meta["mean_path_mm_all"] = line.split()[1]
            elif line.startswith("hist total "):
                expect_total_counts = True
            elif line.startswith("fate "):
                parts = line.split()
                pending_fate = parts[1]
                fates[pending_fate] = {
                    "count": int(parts[3]),
                    "frac": float(parts[5]),
                    "mean_path_mm": float(parts[7]),
                    "under": int(parts[9]),
                    "over": int(parts[11]),
                    "counts": None,
                }
                expect_total_counts = False
            elif line.startswith("counts"):
                vals = [int(x) for x in line.split()[1:]]
                if expect_total_counts:
                    total_counts = vals
                    expect_total_counts = False
                elif pending_fate:
                    fates[pending_fate]["counts"] = np.array(vals, dtype=np.int64)
                    pending_fate = None

    if total_counts is None:
        raise ValueError(f"missing total histogram in {path}")

    return PathHistFile(
        tag=tag,
        trigger_size=SIZE_FROM_TAG.get(tag, float("nan")),
        nbins=int(meta["nbins"]),
        xmax_mm=float(meta["xmax_mm"]),
        total_tracks=int(meta["total_tracks"]),
        mean_path_mm_all=float(meta.get("mean_path_mm_all", "0")),
        total_counts=np.array(total_counts, dtype=np.int64),
        fates=fates,
    )


def _hist_mean(counts: np.ndarray, xmax: float, under: int = 0, over: int = 0) -> float:
    n = int(counts.sum()) + under + over
    if n <= 0:
        return float("nan")
    w = xmax / len(counts)
    xs = (np.arange(len(counts)) + 0.5) * w
    return float((counts * xs).sum() + over * (xmax + 0.5 * w) + under * (0.5 * w)) / n


def write_summary_csv(hists: list[PathHistFile], out_csv: str) -> None:
    hists = sorted(hists, key=lambda h: h.trigger_size)
    fields = [
        "tag",
        "trigger_size",
        "total_tracks",
        "mean_path_mm_all",
        "sipm_detect",
    ]
    for fate in FATE_ORDER:
        fields.append(f"frac_{fate}")
        fields.append(f"mean_path_mm_{fate}")

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for h in hists:
            row = {
                "tag": h.tag,
                "trigger_size": h.trigger_size,
                "total_tracks": h.total_tracks,
                "mean_path_mm_all": h.mean_path_mm_all,
                "sipm_detect": h.fates.get("detected", {}).get("count", 0),
            }
            for fate in FATE_ORDER:
                fd = h.fates.get(fate, {})
                row[f"frac_{fate}"] = fd.get("frac", 0.0)
                row[f"mean_path_mm_{fate}"] = fd.get("mean_path_mm", float("nan"))
            w.writerow(row)


def _add_fate_stacked_legend(ax) -> None:
    """Place fate legend below axes so stacked bars stay visible."""
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=5,
        fontsize=8,
        frameon=True,
    )


def plot_fate_fractions(hists: list[PathHistFile], out_png: str) -> None:
    hists = sorted(hists, key=lambda h: h.trigger_size)
    fig, ax = plt.subplots(figsize=(9, 5.4))
    sizes = _style_trigger_axis(ax, hists)
    bottom = np.zeros(len(hists))
    for fate in FATE_ORDER:
        vals = np.array([h.fates.get(fate, {}).get("frac", 0.0) for h in hists])
        ax.bar(sizes, vals, bottom=bottom, width=0.35, label=fate, color=FATE_COLORS.get(fate))
        bottom += vals
    ax.set_ylabel("fraction of optical photon tracks")
    ax.set_title("Optical photon fate vs trigger size (100 ev diag, v5 geometry)")
    ax.set_ylim(0, 1.05)
    _add_fate_stacked_legend(ax)
    fig.subplots_adjust(bottom=0.28)
    fig.savefig(out_png, dpi=150, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def plot_fate_counts(hists: list[PathHistFile], out_png: str) -> None:
    """Stacked bar by absolute fate counts; total bar height = generated photons (100 ev)."""
    hists = sorted(hists, key=lambda h: h.trigger_size)
    fig, ax = plt.subplots(figsize=(9, 5.4))
    sizes = _style_trigger_axis(ax, hists)
    bottom = np.zeros(len(hists))
    for fate in FATE_ORDER:
        vals = np.array([h.fates.get(fate, {}).get("count", 0) for h in hists], dtype=float)
        ax.bar(sizes, vals, bottom=bottom, width=0.35, label=fate, color=FATE_COLORS.get(fate))
        bottom += vals
    totals = np.array([h.total_tracks for h in hists], dtype=float)
    for x, tot in zip(sizes, totals):
        ax.text(x, tot, f"{tot / 100.0:.0f}/ev", ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("optical photon tracks (100 events)")
    ax.set_title("Optical photon fate vs trigger size — absolute counts (v5, 100 ev)")
    ax.set_ylim(0, totals.max() * 1.1)
    _add_fate_stacked_legend(ax)
    fig.subplots_adjust(bottom=0.28)
    fig.savefig(out_png, dpi=150, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def plot_mean_path(hists: list[PathHistFile], out_png: str) -> None:
    hists = sorted(hists, key=lambda h: h.trigger_size)
    fig, ax = plt.subplots(figsize=(8, 5.2))
    sizes = _style_trigger_axis(ax, hists)
    for fate, style in [
        ("detected", "o-"),
        ("kill_world", "s--"),
        ("absorb_scint", "^:"),
        ("absorb_lg", "v-."),
    ]:
        ys = [h.fates.get(fate, {}).get("mean_path_mm", float("nan")) for h in hists]
        ax.plot(sizes, ys, style, label=fate, color=FATE_COLORS.get(fate))
    ax.plot(sizes, [h.mean_path_mm_all for h in hists], "k*-", label="all")
    ax.set_ylabel("mean total track length (mm)")
    ax.set_title("Mean optical photon path length vs trigger size")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.subplots_adjust(bottom=0.18)
    fig.savefig(out_png, dpi=150, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)


def plot_total_overlay(
    hists: list[PathHistFile],
    out_png: str,
    tags: list[str],
    *,
    xmax_plot: float = 1000.0,
    log_y: bool = False,
) -> None:
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetLegendFont(42)
    ROOT.gStyle.SetLegendTextSize(0.035)
    ROOT.gStyle.SetHistLineWidth(2)

    color_by_tag = {
        "1x1": ROOT.kBlue + 1,
        "2x2": ROOT.kRed + 1,
        "4x4": ROOT.kGreen + 2,
    }

    keepalive: list = []
    canvas = ROOT.TCanvas("c_path_overlay", "c_path_overlay", 900, 560)
    keepalive.append(canvas)
    if log_y:
        canvas.SetLogy()
    legend = ROOT.TLegend(0.52, 0.62, 0.88, 0.88)
    keepalive.append(legend)

    first = True
    ymax = 0.0
    ymin_pos = float("inf")
    lead_hist = None

    for h in sorted(hists, key=lambda x: x.trigger_size):
        if h.tag not in tags:
            continue
        bin_w = h.xmax_mm / h.nbins
        hist = ROOT.TH1F(
            f"h_total_{h.tag}_{'log' if log_y else 'lin'}",
            "Total path length distribution;total optical photon track length (mm);normalized density",
            h.nbins,
            0.0,
            h.xmax_mm,
        )
        keepalive.append(hist)
        for i, count in enumerate(h.total_counts, start=1):
            hist.SetBinContent(i, float(count))
        integral = hist.Integral()
        if integral > 0:
            hist.Scale(1.0 / (integral * bin_w))
        ymax = max(ymax, hist.GetMaximum())
        for i in range(1, h.nbins + 1):
            x_center = (i - 0.5) * bin_w
            if x_center > xmax_plot:
                continue
            dens = hist.GetBinContent(i)
            if dens > 0:
                ymin_pos = min(ymin_pos, dens)
        hist.SetLineColor(color_by_tag.get(h.tag, ROOT.kBlack))
        hist.SetLineWidth(2)
        hist.Draw("HIST" if first else "HIST SAME")
        legend.AddEntry(
            hist,
            f"{_tag_label_root(h.tag)} ({h.trigger_size * 5.0:g} mm, mean={h.mean_path_mm_all:.0f} mm)",
            "l",
        )
        lead_hist = hist
        first = False

    if lead_hist is None:
        raise RuntimeError(f"no histograms for tags={tags}")

    lead_hist.GetXaxis().SetRangeUser(0.0, xmax_plot)
    if ymax > 0:
        if log_y:
            ymin = (ymin_pos * 0.5) if ymin_pos < float("inf") else ymax * 1e-4
            ymin = max(ymin, ymax * 1e-5)
            lead_hist.GetYaxis().SetRangeUser(ymin, ymax * 5.0)
        else:
            lead_hist.GetYaxis().SetRangeUser(0.0, ymax * 1.25)
    legend.Draw()
    canvas.SetGridy()
    canvas.SaveAs(out_png)


def plot_detected_overlay(
    hists: list[PathHistFile],
    out_png: str,
    tags: list[str],
    *,
    xmax_plot: float = 1000.0,
    log_y: bool = False,
    y_mode: str = "density",
    color_by_tag: dict[str, int] | None = None,
) -> None:
    import ROOT

    if y_mode not in ("density", "counts"):
        raise ValueError(f"y_mode must be 'density' or 'counts', got {y_mode!r}")

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetLegendFont(42)
    ROOT.gStyle.SetLegendTextSize(0.035)
    ROOT.gStyle.SetHistLineWidth(2)

    color_by_tag = color_by_tag or {
        "1x1": ROOT.kBlue + 1,
        "2x2": ROOT.kRed + 1,
        "4x4": ROOT.kGreen + 2,
        "4x4_nolg": ROOT.kMagenta + 1,
    }
    tag_order = {t: i for i, t in enumerate(tags)}

    if y_mode == "density":
        title = "Detected photon path length;total track length (mm);normalized density"
    else:
        title = "Detected photon path length (raw counts);total track length (mm);detected photons (100 events)"

    keepalive: list = []
    canvas = ROOT.TCanvas(
        f"c_detected_overlay_{y_mode}_{'log' if log_y else 'lin'}",
        f"c_detected_overlay_{y_mode}",
        900,
        560,
    )
    keepalive.append(canvas)
    if log_y:
        canvas.SetLogy()
    legend = ROOT.TLegend(0.48, 0.62, 0.88, 0.88)
    keepalive.append(legend)

    first = True
    ymax = 0.0
    ymin_pos = float("inf")
    lead_hist = None

    for h in sorted(hists, key=lambda x: tag_order.get(x.tag, 99)):
        if h.tag not in tags:
            continue
        fd = h.fates.get("detected", {})
        counts = fd.get("counts")
        n_detect = fd.get("count", 0)
        if counts is None or n_detect <= 0:
            continue
        mean_mm = fd.get("mean_path_mm", float("nan"))
        bin_w = h.xmax_mm / h.nbins
        hist = ROOT.TH1F(
            f"h_detected_{h.tag}_{y_mode}_{'log' if log_y else 'lin'}",
            title,
            h.nbins,
            0.0,
            h.xmax_mm,
        )
        keepalive.append(hist)
        for i, count in enumerate(counts, start=1):
            hist.SetBinContent(i, float(count))
        if y_mode == "density":
            hist.Scale(1.0 / (float(n_detect) * bin_w))
        ymax = max(ymax, hist.GetMaximum())
        for i in range(1, h.nbins + 1):
            x_center = (i - 0.5) * bin_w
            if x_center > xmax_plot:
                continue
            val = hist.GetBinContent(i)
            if val > 0:
                ymin_pos = min(ymin_pos, val)
        hist.SetLineColor(color_by_tag.get(h.tag, ROOT.kBlack))
        hist.SetLineWidth(2)
        hist.Draw("HIST" if first else "HIST SAME")
        if y_mode == "density":
            legend_label = (
                f"{_tag_label_root(h.tag)} ({h.trigger_size * 5.0:g} mm, mean={mean_mm:.0f} mm)"
            )
        else:
            legend_label = (
                f"{_tag_label_root(h.tag)} ({h.trigger_size * 5.0:g} mm, "
                f"N={n_detect / 100.0:.0f}/ev, mean={mean_mm:.0f} mm)"
            )
        legend.AddEntry(hist, legend_label, "l")
        lead_hist = hist
        first = False

    if lead_hist is None:
        raise RuntimeError(f"no detected histograms for tags={tags}")

    lead_hist.GetXaxis().SetRangeUser(0.0, xmax_plot)
    if ymax > 0:
        if log_y:
            ymin = (ymin_pos * 0.5) if ymin_pos < float("inf") else ymax * 1e-4
            ymin = max(ymin, ymax * 1e-5)
            lead_hist.GetYaxis().SetRangeUser(ymin, ymax * 5.0)
        else:
            lead_hist.GetYaxis().SetRangeUser(0.0, ymax * 1.25)
    legend.Draw()
    canvas.SetGridy()
    canvas.SaveAs(out_png)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot optical path scan results")
    parser.add_argument("--scan-dir", required=True, help="path_scan output directory")
    args = parser.parse_args()

    scan_dir = os.path.abspath(args.scan_dir)
    fig_dir = os.path.join(scan_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    hist_paths = sorted(
        os.path.join(scan_dir, f)
        for f in os.listdir(scan_dir)
        if f.endswith("_100ev.path_hist.txt")
    )
    if not hist_paths:
        raise SystemExit(f"no path_hist files in {scan_dir}")

    hists = [load_path_hist(p) for p in hist_paths]
    scan_hists = [h for h in hists if _is_size_scan_tag(h.tag)]
    lg_nolg_tags = ["4x4", "4x4_nolg"]
    lg_nolg_hists = [h for h in hists if h.tag in lg_nolg_tags]

    summary_csv = os.path.join(scan_dir, "path_vs_size_summary.csv")
    write_summary_csv(scan_hists, summary_csv)
    print(f"wrote {summary_csv}")

    plot_fate_fractions(scan_hists, os.path.join(fig_dir, "fate_fraction_vs_size.png"))
    plot_fate_counts(scan_hists, os.path.join(fig_dir, "fate_counts_vs_size.png"))
    plot_mean_path(scan_hists, os.path.join(fig_dir, "mean_path_vs_size.png"))
    overlay_tags = ["1x1", "2x2", "4x4"]
    xmax_hist = max(h.xmax_mm for h in scan_hists)
    plot_total_overlay(
        scan_hists,
        os.path.join(fig_dir, "total_path_overlay_1x1_2x2_4x4.png"),
        overlay_tags,
        xmax_plot=1000.0,
        log_y=False,
    )
    plot_total_overlay(
        scan_hists,
        os.path.join(fig_dir, "total_path_overlay_1x1_2x2_4x4_logy.png"),
        overlay_tags,
        xmax_plot=xmax_hist,
        log_y=True,
    )
    plot_detected_overlay(
        scan_hists,
        os.path.join(fig_dir, "detected_path_overlay_1x1_2x2_4x4.png"),
        overlay_tags,
        xmax_plot=1000.0,
        log_y=False,
    )
    plot_detected_overlay(
        scan_hists,
        os.path.join(fig_dir, "detected_path_overlay_1x1_2x2_4x4_logy.png"),
        overlay_tags,
        xmax_plot=xmax_hist,
        log_y=True,
    )
    plot_detected_overlay(
        scan_hists,
        os.path.join(fig_dir, "detected_path_counts_overlay_1x1_2x2_4x4.png"),
        overlay_tags,
        xmax_plot=1000.0,
        log_y=False,
        y_mode="counts",
    )
    plot_detected_overlay(
        scan_hists,
        os.path.join(fig_dir, "detected_path_counts_overlay_1x1_2x2_4x4_logy.png"),
        overlay_tags,
        xmax_plot=xmax_hist,
        log_y=True,
        y_mode="counts",
    )

    if len(lg_nolg_hists) >= 2:
        import ROOT

        lg_nolg_colors = {
            "4x4": ROOT.kGreen + 2,
            "4x4_nolg": ROOT.kMagenta + 1,
        }
        plot_detected_overlay(
            lg_nolg_hists,
            os.path.join(fig_dir, "detected_path_counts_overlay_4x4_LG_vs_noLG.png"),
            lg_nolg_tags,
            xmax_plot=1000.0,
            log_y=False,
            y_mode="counts",
            color_by_tag=lg_nolg_colors,
        )
        plot_detected_overlay(
            lg_nolg_hists,
            os.path.join(fig_dir, "detected_path_counts_overlay_4x4_LG_vs_noLG_logy.png"),
            lg_nolg_tags,
            xmax_plot=xmax_hist,
            log_y=True,
            y_mode="counts",
            color_by_tag=lg_nolg_colors,
        )
        plot_detected_overlay(
            lg_nolg_hists,
            os.path.join(fig_dir, "detected_path_overlay_4x4_LG_vs_noLG.png"),
            lg_nolg_tags,
            xmax_plot=1000.0,
            log_y=False,
            y_mode="density",
            color_by_tag=lg_nolg_colors,
        )
        plot_detected_overlay(
            lg_nolg_hists,
            os.path.join(fig_dir, "detected_path_overlay_4x4_LG_vs_noLG_logy.png"),
            lg_nolg_tags,
            xmax_plot=xmax_hist,
            log_y=True,
            y_mode="density",
            color_by_tag=lg_nolg_colors,
        )
    else:
        print("skip 4x4 LG vs no-LG overlays (need path_4x4 and path_4x4_nolg hist files)")
    print(f"figures in {fig_dir}")


if __name__ == "__main__":
    main()
