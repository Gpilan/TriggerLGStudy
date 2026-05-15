"""SiPM 병합 시간 벡터 → ROOT TH1F (고정 폭 리빈, 기본 20 ps)."""

from __future__ import annotations

import math


def bin_width_ns_from_ps(width_ps: float) -> float:
    """1 ps = 1e-3 ns."""
    return float(width_ps) * 1e-3


def data_time_range_ns(
    centers: list[float],
    counts: list[int],
    bin_width_ns: float,
    *,
    margin_ns: float | None = None,
) -> tuple[float, float]:
    """
    포톤이 하나라도 있는 빈만 보고 [min, max] 중심을 잡고,
    여유(margin)를 더한 뒤 `bin_width_ns` 경계에 맞춤.
    """
    w = float(bin_width_ns)
    if w <= 0:
        raise ValueError("bin_width_ns must be > 0")
    mins: list[float] = []
    maxs: list[float] = []
    for c, n in zip(centers, counts):
        if n > 0:
            mins.append(float(c))
            maxs.append(float(c))
    if not mins:
        return 0.0, 240.0

    t_min = min(mins)
    t_max = max(maxs)
    span = t_max - t_min
    if margin_ns is None:
        margin_ns = max(5.0 * w, 0.03 * max(span, w))

    lo = t_min - margin_ns
    hi = t_max + margin_ns
    lo = math.floor(lo / w) * w
    hi = math.ceil(hi / w) * w
    if hi <= lo:
        hi = lo + w
    return lo, hi


def th1_empty_same_bins(
    ROOT,
    name: str,
    title: str,
    t_min_ns: float,
    t_max_ns: float,
    bin_width_ns: float,
):
    """동일 빈 구조의 0만 있는 TH1F (이벤트 평균 시 빈 이벤트 누적용)."""
    w = float(bin_width_ns)
    nbins = max(1, int(round((t_max_ns - t_min_ns) / w)))
    h = ROOT.TH1F(name, title, nbins, t_min_ns, t_min_ns + nbins * w)
    h.Sumw2()
    h.SetDirectory(0)
    return h


def th1_rebinned_from_merged(
    ROOT,
    name: str,
    title: str,
    centers: list[float],
    counts: list[int],
    *,
    bin_width_ns: float,
    t_min_ns: float | None = None,
    t_max_ns: float | None = None,
):
    """
    `bin_width_ns` 폭의 빈으로 `Fill(center, weight)` 누적.
    `t_min_ns` / `t_max_ns` 가 None 이면 **데이터가 있는 시간 구간만** (여유 포함).
    """
    if not centers or not counts:
        return None
    w = float(bin_width_ns)
    if w <= 0:
        raise ValueError("bin_width_ns must be > 0")

    if t_min_ns is None or t_max_ns is None:
        t_min_ns, t_max_ns = data_time_range_ns(centers, counts, w)

    nbins = max(1, int(round((t_max_ns - t_min_ns) / w)))
    h = ROOT.TH1F(name, title, nbins, t_min_ns, t_min_ns + nbins * w)
    h.Sumw2()
    h.SetDirectory(0)
    for c, n in zip(centers, counts):
        if n:
            h.Fill(float(c), float(n))
    return h


def apply_photon_time_style(ROOT) -> None:
    s = ROOT.gStyle
    s.SetOptStat(1111)
    s.SetHistLineWidth(2)
    s.SetTitleFont(42, "XYZ")
    s.SetLabelFont(42, "XYZ")
