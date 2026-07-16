#!/usr/bin/env python3
"""Compare V5 optical diagnostic summary files (key=value per line)."""
from __future__ import annotations

import argparse
from pathlib import Path


def load(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        key, _, val = line.partition(" ")
        try:
            out[key] = float(val)
        except ValueError:
            out[key] = val  # type: ignore[assignment]
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("baseline", type=Path)
    p.add_argument("variant", type=Path)
    p.add_argument(
        "--keys",
        nargs="*",
        default=[
            "scint_to_world",
            "world_to_scint",
            "scint_to_lg",
            "lg_to_world",
            "lg_to_sipm",
            "sipm_detect",
            "path_scint_mm_per_event",
            "path_lg_mm_per_event",
        ],
    )
    args = p.parse_args()
    a = load(args.baseline)
    b = load(args.variant)
    print(f"baseline: {args.baseline.name}")
    print(f"variant:  {args.variant.name}")
    print(f"{'key':<28} {'baseline':>14} {'variant':>14} {'delta%':>10}")
    print("-" * 70)
    for k in args.keys:
        va = float(a.get(k, 0))
        vb = float(b.get(k, 0))
        pct = (100.0 * (vb - va) / va) if va else float("nan")
        print(f"{k:<28} {va:14.4g} {vb:14.4g} {pct:10.2f}")


if __name__ == "__main__":
    main()
