#!/usr/bin/env python3
"""Plot audited, existing ROOT fit and photon statistics; performs no new fit."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PANELS = (
    ('width', 'Width', 'Width W [mm]', 'L = 60 mm, T = 5 mm | 60 GeV $e^-$', [10,15,20,25,30,35,40]),
    ('length', 'Length', 'Length L [mm]', 'W = 10 mm, T = 5 mm | 60 GeV $e^+$', [30,40,50,60,80]),
    ('thickness', 'Thickness', 'Thickness T [mm]', 'W = 10 mm, L = 60 mm | 60 GeV $e^+$', [2,3,5,7,10]),
)
METRICS = (
    ('sigma_ps', 'sigma_error_ps', r'Gaussian fit $\sigma(\Delta t)$ [ps]'),
    ('detected_mean', 'detected_sem', 'Detected photons / event (T1 + T2)'),
)

def groups(rows, scan):
    selected = [r for r in rows if r['scan'] == scan]
    width = None if scan == 'width' else 10
    specs = [('nolg', width, 'noLG', '#0072B2', 'o', '-'), ('lg', width, 'LG', '#009E73', 's', '--')]
    for mode, width, label, color, marker, line in specs:
        points = sorted((r for r in selected if r['mode'] == mode and (width is None or float(r['width_mm']) == width)), key=lambda r: float(r['scan_mm']))
        assert len(points) == (7 if scan == 'width' else 5)
        yield points, label, color, marker, line

def draw(rows, output, metrics, name, title):
    fig, axes = plt.subplots(len(metrics), 3, figsize=(19, 5.1*len(metrics)+1.3), squeeze=False)
    for ri, (value, error, ylabel) in enumerate(metrics):
        lo = min(float(r[value])-float(r[error]) for r in rows)
        hi = max(float(r[value])+float(r[error]) for r in rows)
        ylim = (max(0, lo-(hi-lo)*.1), hi+(hi-lo)*.3)
        for ci, (scan, heading, xlabel, fixed, ticks) in enumerate(PANELS):
            ax = axes[ri, ci]
            for points, label, color, marker, line in groups(rows, scan):
                ax.errorbar([float(r['scan_mm']) for r in points], [float(r[value]) for r in points], yerr=[float(r[error]) for r in points], label=label, color=color, marker=marker, linestyle=line, linewidth=2.1, markersize=7, capsize=4, elinewidth=1.5)
            ax.set(xlabel=xlabel, ylabel=ylabel if ci == 0 else '', xticks=ticks, ylim=ylim)
            ax.axvline({'width':10, 'length':60, 'thickness':5}[scan], color='#889099', linestyle=':', linewidth=1.3, zorder=0)
            ax.grid(alpha=.2)
            ax.tick_params(direction='in', top=True, right=True)
            if ri == 0:
                ax.set_title(heading+'\n'+fixed, fontsize=14, pad=13)
                ax.legend(loc='upper left', fontsize=11, framealpha=.95, ncol=2)
    fig.suptitle(title, fontsize=23, y=.975)
    fig.text(.07, .062, 'Reference geometry: (W, L, T) = (10, 60, 5) mm | vary one dimension at a time', fontsize=12)
    fig.text(.07, .037, '3000 events / condition | error bars: fit SE (timing), SEM (photons)', fontsize=12)
    fig.text(.07, .012, 'Separate production campaigns: beam species and model revisions differ.', fontsize=12, color='#4b5563')
    fig.subplots_adjust(left=.07, right=.985, bottom=.17 if len(metrics)==1 else .125, top=.79 if len(metrics)==1 else .86, wspace=.2, hspace=.3)
    fig.savefig(output/name, dpi=180, facecolor='white')
    plt.close(fig)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args=parser.parse_args()
    rows=list(csv.DictReader(args.data.open()))
    assert len(rows)==44 and len({r['dataset'] for r in rows})==41
    rows = [r for r in rows if r['scan'] == 'width' or float(r['width_mm']) == 10]
    assert len(rows) == 34 and len({r['dataset'] for r in rows}) == 32
    for r in rows:
        fixed = {'width': {'length_mm':60, 'thickness_mm':5}, 'length': {'width_mm':10, 'thickness_mm':5}, 'thickness': {'width_mm':10, 'length_mm':60}}[r['scan']]
        assert all(float(r[k]) == v for k, v in fixed.items())
    args.out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':13, 'axes.labelsize':14, 'xtick.labelsize':12, 'ytick.labelsize':12, 'axes.spines.top':True, 'axes.spines.right':True})
    draw(rows,args.out,METRICS,'tile_dimensions_overview.png','Scintillator dimensions | timing and detected light')
    draw(rows,args.out,METRICS[:1],'timing_sigma_vs_dimensions.png','Timing resolution vs scintillator dimensions')
    draw(rows,args.out,METRICS[1:],'detected_photons_vs_dimensions.png','Detected light vs scintillator dimensions')
    print(json.dumps({'figures':[p.name for p in sorted(args.out.glob('*.png'))], 'input_sha256':hashlib.sha256(args.data.read_bytes()).hexdigest()}))
if __name__=='__main__':
    main()
