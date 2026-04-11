#!/usr/bin/env python3
"""
CBDsim ROOT에서 트리거별(T1/T2) 옵티컬 포톤 도착 시간 분포를 그립니다.

- 우선 이벤트 단위 합산 `timeMergedCountsTrig0/1` (+ edge 벡터) 사용.
- 비어 있으면 `towerT1`/`towerT2`의 SiPM `timeBin*`를 C++과 같이 채널 합산.
- 이벤트마다 기준 시간 t_ref를 빼서 분포를 가운데(0 근처)에 맞춥니다 (기본: 가중 평균).
- 모든 이벤트를 합친 히스토그램에 대해 가중 mean·std(RMS)·총 광자 수 N을 플롯·콘솔에 표시합니다.

필요: ROOT(PyROOT), 빌드 산출물 `build/rootIO/librootIO.so`.

기본 동작 (인자 생략 시):
  저장소 루트를 스크립트 위치에서 추정하고, `build/rootIO` 를 LD_LIBRARY_PATH 앞에 붙인 뒤
  입력: `build/CBDsim/stats/` 안의 `An_test_1.root` (없으면 그 디렉터리의 첫 .root)
  출력: 같은 디렉터리에 `timing_<입력파일stem>.png`

예:
  python3 plot_trigger_timing.py
  python3 plot_trigger_timing.py my.root --ref min
"""

from __future__ import annotations

import argparse
import glob
import os
import sys


def _find_repo_root() -> str | None:
    """`build/rootIO/librootIO.so` 가 나올 때까지 상위 디렉터리 탐색."""
    p = os.path.dirname(os.path.abspath(__file__))
    for _ in range(10):
        cand = os.path.join(p, "build", "rootIO", "librootIO.so")
        if os.path.isfile(cand):
            return p
        parent = os.path.dirname(p)
        if parent == p:
            break
        p = parent
    return None


def _stats_dir(repo_root: str) -> str:
    env = os.environ.get("CBDsim_STATS")
    if env and os.path.isdir(env):
        return os.path.abspath(env)
    return os.path.join(repo_root, "build", "CBDsim", "stats")


def _prepend_build_rootio_ld_path(repo_root: str) -> None:
    libdir = os.path.join(repo_root, "build", "rootIO")
    if not os.path.isdir(libdir):
        return
    ld = os.environ.get("LD_LIBRARY_PATH", "")
    parts = [x for x in ld.split(":") if x]
    if libdir in parts:
        return
    os.environ["LD_LIBRARY_PATH"] = libdir + (":" + ld if ld else "")


def _default_input_root(repo_root: str) -> str | None:
    sd = _stats_dir(repo_root)
    preferred = os.path.join(sd, "An_test_1.root")
    if os.path.isfile(preferred):
        return preferred
    roots = sorted(glob.glob(os.path.join(sd, "*.root")))
    return roots[0] if roots else None


def _resolve_input_path(user_arg: str | None, repo_root: str) -> str | None:
    if user_arg:
        if os.path.isfile(user_arg):
            return os.path.abspath(user_arg)
        base = os.path.basename(user_arg)
        in_stats = os.path.join(_stats_dir(repo_root), base)
        if os.path.isfile(in_stats):
            return os.path.abspath(in_stats)
        under = os.path.join(repo_root, user_arg)
        if os.path.isfile(under):
            return os.path.abspath(under)
        return None
    return _default_input_root(repo_root)


def _resolve_output_path(
    user_out: str | None, input_path: str, repo_root: str
) -> str:
    if user_out:
        return os.path.abspath(user_out)
    stem = os.path.splitext(os.path.basename(input_path))[0]
    sd = _stats_dir(repo_root)
    return os.path.join(sd, f"timing_{stem}.png")


def _load_rootio(lib: str | None, repo_root: str | None) -> None:
    import ROOT

    candidates = []
    if lib:
        candidates.append(lib)
    env = os.environ.get("ROOTIO_LIB")
    if env:
        candidates.append(env)
    here = os.path.dirname(os.path.abspath(__file__))
    candidates.extend(
        [
            os.path.join(here, "..", "build", "rootIO", "librootIO.so"),
            os.path.join(here, "..", "rootIO", "librootIO.so"),
        ]
    )
    if repo_root:
        candidates.append(os.path.join(repo_root, "build", "rootIO", "librootIO.so"))
    for path in candidates:
        ap = os.path.abspath(path)
        if os.path.isfile(ap) and ROOT.gSystem.Load(ap) == 0:
            return
    sys.stderr.write(
        "librootIO.so 를 찾지 못했습니다. Trigger 에서 cmake 빌드 후 "
        "build/rootIO/librootIO.so 가 있는지 확인하거나 -l / ROOTIO_LIB 로 지정하세요.\n"
    )
    sys.exit(1)


def _hist_mean_std_rms(h) -> tuple[float | None, float | None, float | None]:
    """
    합산 히스토그램(t - t_ref 분포)의 가중 평균, 표준편차, 총 가중치.
    ROOT TH1::GetMean / GetRMS 는 빈 중심에 대한 가중 통계와 일치한다.
    """
    import ROOT

    w = float(h.GetSumOfWeights())
    if w <= 0:
        return None, None, None
    return float(h.GetMean()), float(h.GetRMS()), w


def _vec_to_lists(lo, hi, cnt):
    """PyROOT vector -> Python lists."""
    n = cnt.size()
    if n == 0:
        return [], [], []
    lows = [float(lo[i]) for i in range(n)]
    highs = [float(hi[i]) for i in range(n)]
    counts = [int(cnt[i]) for i in range(n)]
    return lows, highs, counts


def _merge_sipm_side(acc, sipm_vec) -> None:
    """Merge one SiPM vector into acc = (lows, highs, counts) or None."""
    n = sipm_vec.size()
    for i in range(n):
        s = sipm_vec[i]
        if s.timeBinCounts.size() == 0:
            continue
        lo, hi, c = _vec_to_lists(s.timeBinEdgeLow, s.timeBinEdgeHigh, s.timeBinCounts)
        if not c:
            continue
        if acc[0] is None:
            acc[0] = lo
            acc[1] = hi
            acc[2] = list(c)
        else:
            m = min(len(acc[2]), len(c))
            for j in range(m):
                acc[2][j] += c[j]


def _merged_from_tower(evt, trig: int):
    """(lows, highs, counts) for trigger 0=T1, 1=T2."""
    if trig == 0:
        lo, hi, cnt = (
            evt.timeMergedEdgeLowTrig0,
            evt.timeMergedEdgeHighTrig0,
            evt.timeMergedCountsTrig0,
        )
    else:
        lo, hi, cnt = (
            evt.timeMergedEdgeLowTrig1,
            evt.timeMergedEdgeHighTrig1,
            evt.timeMergedCountsTrig1,
        )
    if cnt.size() > 0:
        return _vec_to_lists(lo, hi, cnt)

    tw = evt.towerT1 if trig == 0 else evt.towerT2
    acc = [None, None, None]
    _merge_sipm_side(acc, tw.SiPMs)
    if acc[0] is None:
        return [], [], []
    return acc[0], acc[1], acc[2]


def _ref_time(lows, highs, counts, ref_mode: str) -> float | None:
    """기준 시간 t_ref (ns). 분포를 이 값 기준으로 0 근처에 맞춤."""
    try:
        import numpy as np
    except ImportError:
        np = None

    total = sum(counts)
    if total <= 0:
        return None

    if np is not None:
        lo = np.asarray(lows, dtype=float)
        hi = np.asarray(highs, dtype=float)
        c = np.asarray(counts, dtype=float)
        centers = 0.5 * (lo + hi)
    else:
        centers = [0.5 * (lows[i] + highs[i]) for i in range(len(counts))]

    if ref_mode == "mean":
        if np is not None:
            return float(np.average(centers, weights=c))
        s = sum(centers[i] * counts[i] for i in range(len(counts)))
        return s / total

    if ref_mode == "min":
        # 신호가 있는 가장 이른 빈의 중심
        for i in range(len(counts)):
            if counts[i] > 0:
                return centers[i] if np is None else float(centers[i])
        return None

    if ref_mode == "mode":
        # 광자 수가 가장 많은 빈의 중심
        imax = max(range(len(counts)), key=lambda i: counts[i])
        return centers[imax] if np is None else float(centers[imax])

    raise ValueError(f"unknown ref mode: {ref_mode}")


def _accumulate_np(
    acc,
    lows,
    highs,
    counts,
    t_ref: float,
    scale: float,
    edges,
) -> None:
    """t' = (t_center - t_ref) * scale 로 가중치 counts 만큼 edges 빈에 누적 (numpy)."""
    import numpy as np

    lo = np.asarray(lows, dtype=np.float64)
    hi = np.asarray(highs, dtype=np.float64)
    c = np.asarray(counts, dtype=np.float64)
    centers = 0.5 * (lo + hi)
    x = (centers - t_ref) * scale
    m = c > 0
    if not np.any(m):
        return
    hist, _ = np.histogram(x[m], bins=edges, weights=c[m])
    acc += hist


def _accumulate_root(global_hist, lows, highs, counts, t_ref: float, scale: float) -> None:
    """느린 경로: PyROOT TH1F.Fill (numpy 없을 때만)."""
    for i in range(len(counts)):
        if counts[i] == 0:
            continue
        t_c = 0.5 * (lows[i] + highs[i])
        x = (t_c - t_ref) * scale
        global_hist.Fill(x, float(counts[i]))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="T1/T2 옵티컬 포톤 시간 구조 (빈 합산) — 기준시간 보정 후 히스토그램"
    )
    parser.add_argument(
        "input",
        nargs="?",
        default=None,
        help="CBDsim .root (생략 시 build/CBDsim/stats/An_test_1.root 또는 그 디렉터리의 첫 .root)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="출력 이미지 (생략 시 입력과 같은 디렉터리에 timing_<이름>.png)",
    )
    parser.add_argument(
        "-l",
        "--rootio-lib",
        default=None,
        help="librootIO.so 절대 경로",
    )
    parser.add_argument(
        "--ref",
        choices=("mean", "min", "mode"),
        default="mean",
        help="이벤트별 기준 시간: mean=가중평균, min=가장 이른 빈, mode=최다빈 (기본 mean)",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=-1,
        help="처리할 최대 이벤트 수 (-1=전부)",
    )
    parser.add_argument(
        "--xmin",
        type=float,
        default=-15.0,
        help="보정 후 시간 축 하한 (ns)",
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=15.0,
        help="보정 후 시간 축 상한 (ns)",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=120,
        help="히스토그램 빈 수",
    )
    parser.add_argument(
        "--ns-per-unit",
        type=float,
        default=1.0,
        help="축 스케일 (1.0이면 단위 ns)",
    )
    args = parser.parse_args()

    repo_root = _find_repo_root()
    if repo_root is None:
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    _prepend_build_rootio_ld_path(repo_root)

    input_path = _resolve_input_path(args.input, repo_root)
    if not input_path:
        sys.stderr.write(
            "입력 .root 를 찾을 수 없습니다. 파일 경로를 주거나 "
            f"{_stats_dir(repo_root)}/ 에 .root 를 두세요.\n"
        )
        sys.exit(1)
    output_path = _resolve_output_path(args.output, input_path, repo_root)

    # matplotlib 는 첫 import 시 폰트 캐시 등으로 수 분 걸릴 수 있음 → ROOT 루프 전에 로드
    _mpl_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".mplconfig")
    os.makedirs(_mpl_dir, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", _mpl_dir)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot  # noqa: F401
    except ImportError:
        pass

    _load_rootio(args.rootio_lib, repo_root)
    import ROOT

    ROOT.gROOT.SetBatch(True)

    f = ROOT.TFile.Open(input_path)
    if not f or f.IsZombie():
        sys.stderr.write(f"파일을 열 수 없습니다: {input_path}\n")
        sys.exit(1)
    tree = f.Get("CBDsim")
    if not tree:
        sys.stderr.write("트리 CBDsim 없음\n")
        sys.exit(1)
    total_entries = int(tree.GetEntries())

    evt = ROOT.CBDsimInterface.CBDsimEventData()
    tree.SetBranchAddress("CBDsimEventData", evt)

    nmax = tree.GetEntries()
    if args.max_events >= 0:
        nmax = min(nmax, args.max_events)

    try:
        import numpy as np
    except ImportError:
        np = None

    edges = None
    acc1 = acc2 = None
    if np is not None:
        edges = np.linspace(args.xmin, args.xmax, args.bins + 1)
        acc1 = np.zeros(args.bins, dtype=np.float64)
        acc2 = np.zeros(args.bins, dtype=np.float64)

    h1 = ROOT.TH1F(
        "hT1",
        "Trigger 1 (T1);t - t_{ref} (ns);optical photons",
        args.bins,
        args.xmin,
        args.xmax,
    )
    h2 = ROOT.TH1F(
        "hT2",
        "Trigger 2 (T2);t - t_{ref} (ns);optical photons",
        args.bins,
        args.xmin,
        args.xmax,
    )
    h1.Sumw2()
    h2.Sumw2()

    skipped = [0, 0]
    for i in range(int(nmax)):
        tree.GetEntry(i)
        for trig, sk, acc in ((0, 0, acc1), (1, 1, acc2)):
            lows, highs, counts = _merged_from_tower(evt, trig)
            if not counts or sum(counts) == 0:
                skipped[sk] += 1
                continue
            t_ref = _ref_time(lows, highs, counts, args.ref)
            if t_ref is None:
                skipped[sk] += 1
                continue
            if np is not None and edges is not None and acc is not None:
                _accumulate_np(acc, lows, highs, counts, t_ref, args.ns_per_unit, edges)
            else:
                h = h1 if trig == 0 else h2
                _accumulate_root(h, lows, highs, counts, t_ref, args.ns_per_unit)

    if np is not None and acc1 is not None and acc2 is not None:
        for b in range(1, args.bins + 1):
            h1.SetBinContent(b, acc1[b - 1])
            h2.SetBinContent(b, acc2[b - 1])
            h1.SetBinError(b, np.sqrt(max(acc1[b - 1], 0.0)))
            h2.SetBinError(b, np.sqrt(max(acc2[b - 1], 0.0)))

    # matplotlib 플롯 (위에서 Agg 이미 설정됨)
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        sys.stderr.write("matplotlib/numpy 없음 — ROOT Canvas 로 저장합니다.\n")
        c = ROOT.TCanvas("c", "", 900, 400)
        c.Divide(2, 1)
        c.cd(1)
        h1.Draw("HIST")
        c.cd(2)
        h2.Draw("HIST")
        out = output_path
        if not out.lower().endswith((".png", ".pdf", ".svg")):
            out += ".png"
        c.SaveAs(out)
        m1, s1, w1 = _hist_mean_std_rms(h1)
        m2, s2, w2 = _hist_mean_std_rms(h2)
        f.Close()
        print(
            f"저장: {out} (skipped T1={skipped[0]} T2={skipped[1]} empty events, "
            f"n={int(nmax)}/{total_entries})"
        )
        for lab, m, s, w in (
            ("T1", m1, s1, w1),
            ("T2", m2, s2, w2),
        ):
            if m is None:
                print(f"  {lab}: 합산 분포 없음 (mean/std N/A)")
            else:
                print(f"  {lab}: mean={m:.4f} ns, std(RMS)={s:.4f} ns, N={w:.0f}")
        return

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)

    def root_hist_to_xy(h):
        cx = np.array([h.GetBinCenter(b) for b in range(1, h.GetNbinsX() + 1)])
        cy = np.array([h.GetBinContent(b) for b in range(1, h.GetNbinsX() + 1)])
        ce = np.array([h.GetBinError(b) for b in range(1, h.GetNbinsX() + 1)])
        return cx, cy, ce

    for ax, h, title in zip(
        axes,
        (h1, h2),
        ("T1 (trigger 0)", "T2 (trigger 1)"),
    ):
        x, y, err = root_hist_to_xy(h)
        ax.bar(x, y, width=(x[1] - x[0]) if len(x) > 1 else 0.1, align="center", alpha=0.85)
        ax.set_xlabel(r"$t - t_{\mathrm{ref}}$ (ns)" + (f" [×{args.ns_per_unit}]" if args.ns_per_unit != 1 else ""))
        ax.set_ylabel("optical photons (summed)")
        m, s, w = _hist_mean_std_rms(h)
        stat_line = (
            f"mean={m:.3f} ns, std={s:.3f} ns, N={w:.0f}"
            if m is not None
            else "합산 분포 없음"
        )
        ax.set_title(f"{title}\n(ref={args.ref} per event)\n{stat_line}", fontsize=9)
        ax.axvline(0.0, color="k", ls="--", lw=0.8, alpha=0.5)

    fig.suptitle(os.path.basename(input_path))
    fig.tight_layout()
    out = output_path
    if not out.lower().endswith((".png", ".pdf", ".svg")):
        out += ".png"
    fig.savefig(out, dpi=150)
    m1, s1, w1 = _hist_mean_std_rms(h1)
    m2, s2, w2 = _hist_mean_std_rms(h2)
    f.Close()
    print(f"저장: {out}")
    print(
        f"  이벤트 사용: {int(nmax)} / {total_entries} "
        f"(빈 이벤트 스킵 T1={skipped[0]}, T2={skipped[1]})"
    )
    for lab, m, s, w in (
        ("T1", m1, s1, w1),
        ("T2", m2, s2, w2),
    ):
        if m is None:
            print(f"  {lab}: 합산 분포 없음 (mean/std N/A)")
        else:
            print(f"  {lab}: mean={m:.4f} ns, std(RMS)={s:.4f} ns, N={w:.0f}")


if __name__ == "__main__":
    main()
