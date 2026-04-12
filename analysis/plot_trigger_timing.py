#!/usr/bin/env python3
"""
CBDsim ROOT에서 트리거별(T1/T2) 옵티컬 포톤 도착 시간 분포를 그립니다.

- 막대: **에러바 없음** (HIST만). Δt 패널은 **가우시안 피팅**으로 σ 추출(ROOT: `Fit("gaus")`, mpl: `scipy.optimize.curve_fit`, 없으면 RMS만 표시).
- 패널 우상단(matplotlib): N, μ±σ_μ, σ(RMS 또는 Gauss σ) (ns).
- **세 번째 패널**: 이벤트마다 빈 중심 가중 평균 도착시각 ⟨t⟩(절대 ns)을 T1/T2 각각 구한 뒤 **⟨t⟩_T1 − ⟨t⟩_T2** 를 한 점으로 두고 **이벤트 수**를 쌓은 히스토그램 (y ≥ 0).
- **`--style`**: `matplotlib` | `root` | `both`(기본) — ROOT는 통상 스타일(OptStat, HIST E1)로 3패널 저장.

- 우선 이벤트 단위 합산 `timeMergedCountsTrig0/1` (+ edge 벡터) 사용.
- 비어 있으면 `towerT1`/`towerT2`의 SiPM `timeBin*`를 C++과 같이 채널 합산.
- 이벤트마다 기준 시간 t_ref를 빼서 분포를 가운데(0 근처)에 맞춥니다 (기본: 가중 평균).

필요: ROOT(PyROOT), 빌드 산출물 `build/rootIO/librootIO.so`.

기본 동작 (인자 생략 시):
  저장소 루트를 스크립트 위치에서 추정하고, `build/rootIO` 를 LD_LIBRARY_PATH 앞에 붙인 뒤
  입력: `analysis/t_res/data/` → `build/CBDsim/stats/` → `build/CBDsim/` 순으로 검색.
         우선 파일명 `noLG1_0.root`, 없으면 `An_test_1.root`, 없으면 해당 디렉터리의 첫 `.root`
  출력: `analysis/t_res/figures/` 에 `timing_<입력파일stem>.png`
  환경변수 `CBDsim_STATS` → **첫 번째** 입력 디렉터리만 덮어쓰기 (나머지 폴백은 그대로), `CBDsim_FIGURES` → PNG 출력

예 (저장소 루트에서, Proto 지오메트리·LG 없음 가정):
  python3 analysis/plot_trigger_timing.py
  python3 analysis/plot_trigger_timing.py noLG1_0.root --ref min
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


def _data_dir(repo_root: str) -> str:
    """기본 입력 .root 디렉터리 (`analysis/t_res/data`). `CBDsim_STATS`로 재정의."""
    env = os.environ.get("CBDsim_STATS")
    if env and os.path.isdir(env):
        return os.path.abspath(env)
    return os.path.join(repo_root, "analysis", "t_res", "data")


def _figures_dir(repo_root: str) -> str:
    """기본 PNG 출력 디렉터리 (`analysis/t_res/figures`). `CBDsim_FIGURES`로 재정의."""
    env = os.environ.get("CBDsim_FIGURES")
    if env and os.path.isdir(env):
        return os.path.abspath(env)
    return os.path.join(repo_root, "analysis", "t_res", "figures")


def _legacy_stats_dir(repo_root: str) -> str:
    """`build/CBDsim/stats` (cmake가 만드는 통계 디렉터리)."""
    return os.path.join(repo_root, "build", "CBDsim", "stats")


def _cbdsim_build_dir(repo_root: str) -> str:
    """실행 파일 cwd에 생기는 .root (`build/CBDsim/`)."""
    return os.path.join(repo_root, "build", "CBDsim")


def _input_root_dirs(repo_root: str) -> list[str]:
    """입력 .root 검색 순서 (중복 제거). `CBDsim_STATS`가 있으면 첫 항목만 대체."""
    seen: set[str] = set()
    out: list[str] = []
    for d in (_data_dir(repo_root), _legacy_stats_dir(repo_root), _cbdsim_build_dir(repo_root)):
        ad = os.path.abspath(d)
        if os.path.isdir(ad) and ad not in seen:
            seen.add(ad)
            out.append(ad)
    return out


# Proto(LG 없음) 이후 권장 접두어 산출명과 구버전 예시
_PREFERRED_ROOT_NAMES: tuple[str, ...] = ("noLG1_0.root", "An_test_1.root")


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
    for sd in _input_root_dirs(repo_root):
        for name in _PREFERRED_ROOT_NAMES:
            cand = os.path.join(sd, name)
            if os.path.isfile(cand):
                return cand
        roots = sorted(glob.glob(os.path.join(sd, "*.root")))
        if roots:
            return roots[0]
    return None


def _resolve_input_path(user_arg: str | None, repo_root: str) -> str | None:
    if user_arg:
        if os.path.isfile(user_arg):
            return os.path.abspath(user_arg)
        base = os.path.basename(user_arg)
        for sd in _input_root_dirs(repo_root):
            cand = os.path.join(sd, base)
            if os.path.isfile(cand):
                return os.path.abspath(cand)
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
    return os.path.join(_figures_dir(repo_root), f"timing_{stem}.png")


def _path_with_tag(path: str, tag: str) -> str:
    """foo/bar.png -> foo/bar_tag.png"""
    d, f = os.path.split(path)
    base, ext = os.path.splitext(f)
    return os.path.join(d, f"{base}_{tag}{ext}")


def _apply_root_plot_style(ROOT) -> None:
    """통상 ROOT 배치 플롯에 가깝게 (투박)."""
    s = ROOT.gStyle
    s.SetOptStat(1111)
    s.SetOptTitle(1)
    s.SetHistLineWidth(2)
    s.SetFrameLineWidth(1)
    s.SetTitleFont(42, "XYZ")
    s.SetLabelFont(42, "XYZ")
    s.SetStatFont(42)


def _save_root_canvas(
    ROOT,
    h1,
    h2,
    h_delta,
    out_path: str,
    input_basename: str,
) -> None:
    """T1, T2, per-event ⟨t⟩ 차이 세 패널 — 막대만(HIST), Δt 패널은 가우시안 피팅 곡선."""
    _apply_root_plot_style(ROOT)
    c = ROOT.TCanvas("c_timing", "timing", 1350, 420)
    if input_basename:
        c.SetTitle(input_basename)
    c.Divide(3, 1)
    hs = (
        (h1, "T1 (trigger 0)", "t - t_{ref} (ns)", "optical photons"),
        (h2, "T2 (trigger 1)", "t - t_{ref} (ns)", "optical photons"),
        (
            h_delta,
            "Per-event #LT t#GT_{T1}-#LT t#GT_{T2}",
            "#LT t#GT_{T1}-#LT t#GT_{T2} (ns)",
            "events",
        ),
    )
    for i, (h, subtitle, xax, yax) in enumerate(hs, start=1):
        c.cd(i)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetBottomMargin(0.15)
        h.SetLineWidth(2)
        if i < 3:
            h.SetLineColor(ROOT.kBlue + 1)
            h.SetLineStyle(1)
            h.SetTitle(f"{subtitle};{xax};{yax}")
            h.Draw("HIST")
        else:
            # 데이터: 검은 실선(스텝), 피트: 빨간 대시 — 피트 먼저 그린 뒤 히스토를 SAME으로 덮어 둘 다 보이게
            h.SetLineColor(ROOT.kBlack)
            h.SetLineStyle(1)  # solid
            h.SetLineWidth(2)
            h.SetFillStyle(0)
            h.SetMarkerStyle(0)
            h.SetTitle(f"{subtitle};{xax};{yax}")
            xmin = float(h.GetXaxis().GetXmin())
            xmax = float(h.GetXaxis().GetXmax())
            fitfn = None
            if h.GetSumOfWeights() > 0:
                fitfn = ROOT.TF1("fit_delta_gaus", "gaus", xmin, xmax)
                fitfn.SetNpx(500)
                fitfn.SetParameter(0, float(h.GetMaximum()))
                fitfn.SetParameter(1, float(h.GetMean()))
                fitfn.SetParameter(2, max(float(h.GetRMS()), 1e-6))
                h.Fit(fitfn, "QN", "", xmin, xmax)
            h.Draw("HIST")
            if fitfn is not None:
                fitfn.SetLineColor(ROOT.kRed)
                fitfn.SetLineWidth(2)
                fitfn.SetLineStyle(1)  
                fitfn.Draw("same")
                pad = ROOT.gPad
                pad.Modified()
                pad.Update()
                mu = fitfn.GetParameter(1)
                mu_e = fitfn.GetParError(1)
                sig = fitfn.GetParameter(2)
                sig_e = fitfn.GetParError(2)
                print(
                    "  ROOT Gauss fit (Δt): "
                    f"μ = {mu:.6f} ± {mu_e:.6f} ns, "
                    f"σ = {sig:.6f} ± {sig_e:.6f} ns"
                )
                lat = ROOT.TLatex()
                lat.SetNDC()
                lat.SetTextFont(42)
                lat.SetTextSize(0.034)
                lat.DrawLatex(0.14, 0.84, f"Gauss #mu = {mu:.3f} #pm {mu_e:.3f} ns")
                lat.DrawLatex(0.14, 0.77, f"Gauss #sigma = {sig:.3f} #pm {sig_e:.3f} ns")
    out = out_path
    if not out.lower().endswith((".png", ".pdf", ".svg")):
        out += ".png"
    c.SaveAs(out)
    c.Close()


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
    if hasattr(tw, "SiPMFronts"):
        _merge_sipm_side(acc, tw.SiPMFronts)
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


def _weighted_mean_absolute(
    lows: list[float],
    highs: list[float],
    counts: list[int],
    scale: float,
) -> float | None:
    """빈 중심의 포톤 수 가중 평균 시각 (절대 ns, 트리거별 t_ref 보정 없음)."""
    if not counts:
        return None
    tot = sum(counts)
    if tot <= 0:
        return None
    s = 0.0
    for i in range(len(counts)):
        if counts[i] == 0:
            continue
        t_c = 0.5 * (lows[i] + highs[i])
        s += t_c * counts[i]
    return (s / float(tot)) * scale


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
        description=(
            "T1/T2 옵티컬 포톤 시간 구조 (빈 합산, 기준시간 보정) + "
            "이벤트별 ⟨t⟩_T1−⟨t⟩_T2 (절대 ns)"
        )
    )
    parser.add_argument(
        "input",
        nargs="?",
        default=None,
        help="CBDsim .root (생략 시 noLG1_0.root / An_test_1.root 우선, data·stats·build/CBDsim 순 검색)",
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
    parser.add_argument(
        "--style",
        choices=("matplotlib", "root", "both"),
        default="both",
        help="matplotlib: mpl 3패널(T1,T2,이벤트별⟨t⟩차이), root: ROOT 3패널, both: 둘 다 저장",
    )
    parser.add_argument(
        "--delta-xmin",
        type=float,
        default=None,
        help="세 번째 패널(이벤트별 Δt) x 하한 (기본: --xmin 과 동일)",
    )
    parser.add_argument(
        "--delta-xmax",
        type=float,
        default=None,
        help="세 번째 패널 x 상한 (기본: --xmax 와 동일)",
    )
    parser.add_argument(
        "--delta-bins",
        type=int,
        default=None,
        help="세 번째 패널 빈 수 (기본: --bins 와 동일)",
    )
    args = parser.parse_args()

    repo_root = _find_repo_root()
    if repo_root is None:
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    _prepend_build_rootio_ld_path(repo_root)

    input_path = _resolve_input_path(args.input, repo_root)
    if not input_path:
        dirs = ", ".join(_input_root_dirs(repo_root))
        sys.stderr.write(
            "입력 .root 를 찾을 수 없습니다. 경로를 직접 주거나 다음 중 한 곳에 .root 를 두세요: "
            f"{dirs}\n"
        )
        sys.exit(1)
    output_path = _resolve_output_path(args.output, input_path, repo_root)
    _od = os.path.dirname(output_path)
    if _od:
        os.makedirs(_od, exist_ok=True)

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

    d_xmin = args.delta_xmin if args.delta_xmin is not None else args.xmin
    d_xmax = args.delta_xmax if args.delta_xmax is not None else args.xmax
    d_bins = args.delta_bins if args.delta_bins is not None else args.bins

    h_delta = ROOT.TH1F(
        "hDeltaMeanT",
        "Per-event mean time difference;#LT t#GT_{T1}-#LT t#GT_{T2} (ns);events",
        d_bins,
        d_xmin,
        d_xmax,
    )
    h_delta.Sumw2()
    # TFile 이 열려 있을 때 gDirectory 가 파일을 가리키면 TH1 이 파일 소유가 되어 Close() 시 삭제됨
    h1.SetDirectory(0)
    h2.SetDirectory(0)
    h_delta.SetDirectory(0)

    skipped = [0, 0]
    skipped_delta = 0
    for i in range(int(nmax)):
        tree.GetEntry(i)
        lo1, hi1, cnt1 = _merged_from_tower(evt, 0)
        lo2, hi2, cnt2 = _merged_from_tower(evt, 1)
        mu1 = _weighted_mean_absolute(lo1, hi1, cnt1, args.ns_per_unit)
        mu2 = _weighted_mean_absolute(lo2, hi2, cnt2, args.ns_per_unit)
        if mu1 is not None and mu2 is not None:
            h_delta.Fill(mu1 - mu2)
        else:
            skipped_delta += 1

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

    # T1/T2: OptStat「Entries」= 해당 트리거에 포톤이 있어 누적에 포함된 이벤트 수 (SetBinContent만 쓸 때 120=빈 개수로 잘못 보이던 것 수정)
    evt_t1 = int(nmax) - skipped[0]
    evt_t2 = int(nmax) - skipped[1]
    h1.SetEntries(float(evt_t1))
    h2.SetEntries(float(evt_t2))

    bn = os.path.basename(input_path)
    style = args.style

    root_out = output_path if style == "root" else _path_with_tag(output_path, "root")
    if style in ("root", "both"):
        _save_root_canvas(ROOT, h1, h2, h_delta, root_out, bn)
        print(f"저장 (ROOT): {root_out}")

    if style in ("matplotlib", "both"):
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            if style == "matplotlib":
                sys.stderr.write("matplotlib/numpy 없음 — 종료.\n")
                f.Close()
                sys.exit(1)
            sys.stderr.write("matplotlib/numpy 없음 — ROOT 출력만 사용.\n")
        else:
            fig, axes = plt.subplots(1, 3, figsize=(14, 4))

            def root_hist_to_xy(h):
                cx = np.array([h.GetBinCenter(b) for b in range(1, h.GetNbinsX() + 1)])
                cy = np.array([h.GetBinContent(b) for b in range(1, h.GetNbinsX() + 1)])
                return cx, cy

            def _gauss(x, a, mu, sigma):
                return a * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

            for ax, h, title, is_delta in zip(
                axes,
                (h1, h2, h_delta),
                (
                    "T1 (trigger 0)",
                    "T2 (trigger 1)",
                    r"$\langle t\rangle_{\mathrm{T1}}-\langle t\rangle_{\mathrm{T2}}$ (per event)",
                ),
                (False, False, True),
            ):
                x, y = root_hist_to_xy(h)
                w = (x[1] - x[0]) if len(x) > 1 else 0.1
                if is_delta:
                    # 데이터: 검은 실선 스텝 / 피트는 아래에서 빨간 대시로 그림
                    edges = np.linspace(
                        float(x[0]) - w / 2,
                        float(x[-1]) + w / 2,
                        len(x) + 1,
                    )
                    ax.stairs(
                        y,
                        edges,
                        color="black",
                        linewidth=2.2,
                        linestyle="-",
                        zorder=2,
                    )
                else:
                    ax.bar(
                        x,
                        y,
                        width=w,
                        align="center",
                        alpha=0.85,
                        color="#4477aa",
                    )
                if is_delta:
                    ax.set_xlabel(
                        r"$\langle t\rangle_{\mathrm{T1}} - \langle t\rangle_{\mathrm{T2}}$ (ns)"
                        + (f" [×{args.ns_per_unit}]" if args.ns_per_unit != 1 else "")
                    )
                    ax.set_ylabel("events")
                    ax.set_title("Absolute weighted mean time\n(both triggers)")
                else:
                    ax.set_xlabel(
                        r"$t - t_{\mathrm{ref}}$ (ns)"
                        + (f" [×{args.ns_per_unit}]" if args.ns_per_unit != 1 else "")
                    )
                    ax.set_ylabel("optical photons (summed)")
                    ax.set_title(f"{title}\n(ref={args.ref} per event)")
                ax.axvline(0.0, color="k", ls="--", lw=0.8, alpha=0.5)
                if is_delta:
                    ntot = float(h.GetEntries())
                else:
                    ntot = float(evt_t1 if h.GetName() == "hT1" else evt_t2)
                mu = float(h.GetMean())
                sigma = float(h.GetRMS())
                mu_err = float(h.GetMeanError()) if h.GetSumOfWeights() > 0 else 0.0
                try:
                    sig_err = float(h.GetRMSError())
                except Exception:
                    sig_err = 0.0

                stat_lines = [
                    rf"$N={ntot:.0f}$",
                    rf"$\mu={mu:.3f}\pm{mu_err:.3f}$ ns",
                ]
                fit_drawn = False
                if is_delta and np.sum(y) > 0 and len(x) >= 3:
                    sw = float(np.sum(y))
                    mu0 = float(np.sum(x * y) / sw)
                    var0 = float(np.sum(y * (x - mu0) ** 2) / sw)
                    sig0 = max(np.sqrt(max(var0, 1e-18)), 1e-6)
                    p0 = [max(float(np.max(y)), 1e-9), mu0, sig0]
                    fit_ok = False
                    try:
                        from scipy.optimize import curve_fit

                        popt, pcov = curve_fit(
                            _gauss,
                            x,
                            y,
                            p0=p0,
                            sigma=np.sqrt(np.maximum(y, 1.0)),
                            absolute_sigma=True,
                            maxfev=10000,
                        )
                        perr = np.sqrt(np.diag(pcov))
                        xf = np.linspace(float(x.min()), float(x.max()), 200)
                        ax.plot(
                            xf,
                            _gauss(xf, *popt),
                            color="red",
                            ls="--",
                            lw=2.0,
                            zorder=3,
                            label="Gauss fit",
                        )
                        fit_ok = True
                        fit_drawn = True
                        print(
                            "  mpl Gauss fit (Δt): "
                            f"μ = {popt[1]:.6f} ± {perr[1]:.6f} ns, "
                            f"σ = {popt[2]:.6f} ± {perr[2]:.6f} ns"
                        )
                        stat_lines.append(
                            rf"Gauss $\sigma={popt[2]:.3f}\pm{perr[2]:.3f}$ ns"
                        )
                    except Exception:
                        pass
                    if not fit_ok:
                        if sig_err > 0.0:
                            stat_lines.append(rf"RMS $\sigma={sigma:.3f}\pm{sig_err:.3f}$ ns")
                        else:
                            stat_lines.append(rf"RMS $\sigma={sigma:.3f}$ ns")
                elif not is_delta:
                    if sig_err > 0.0:
                        stat_lines.append(rf"$\sigma={sigma:.3f}\pm{sig_err:.3f}$ ns")
                    else:
                        stat_lines.append(rf"$\sigma={sigma:.3f}$ ns")
                else:
                    if sig_err > 0.0:
                        stat_lines.append(rf"RMS $\sigma={sigma:.3f}\pm{sig_err:.3f}$ ns")
                    else:
                        stat_lines.append(rf"RMS $\sigma={sigma:.3f}$ ns")

                stat_txt = "\n".join(stat_lines)
                ax.text(
                    0.97,
                    0.97,
                    stat_txt,
                    transform=ax.transAxes,
                    fontsize=8,
                    verticalalignment="top",
                    horizontalalignment="right",
                    bbox=dict(
                        boxstyle="round,pad=0.25",
                        facecolor="white",
                        edgecolor="0.7",
                        alpha=0.92,
                    ),
                )
                if fit_drawn:
                    ax.legend(loc="upper left", fontsize=7)

            fig.suptitle(bn)
            fig.tight_layout()
            mpl_out = output_path
            if not mpl_out.lower().endswith((".png", ".pdf", ".svg")):
                mpl_out += ".png"
            fig.savefig(mpl_out, dpi=150)
            print(f"저장 (matplotlib): {mpl_out}")

    n_dt = int(h_delta.GetEntries())
    f.Close()
    print(
        f"  이벤트 사용: {int(nmax)} / {total_entries} "
        f"(빈 이벤트 스킵 T1={skipped[0]}, T2={skipped[1]})"
    )
    print(
        f"  Δt 히스토그램: {n_dt} entries (한쪽 트리거만/무포톤으로 Δt 제외: {skipped_delta})"
    )


if __name__ == "__main__":
    main()
