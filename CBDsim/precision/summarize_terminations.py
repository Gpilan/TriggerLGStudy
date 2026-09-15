#!/usr/bin/env python3
"""Classify every recorded optical terminal point; fractions use tracked photons, not steps."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path


def classify(row):
    fate = row['fate']
    boundary = row['final_pre_volume'] + ' ' + row['final_post_volume']
    if fate == 'absorb_other_surface':
        if 'protoGelTapePhys' in boundary:
            return 'tape_surface_absorption'
        if 'protoFoilCornerSealPhys' in boundary:
            return 'corner_seal_absorption'
        if 'Foil' in boundary:
            return 'foil_surface_absorption'
        return 'other_surface_absorption'
    return fate


def summarize(run):
    manifest = json.loads((run/'manifest.json').read_text())
    if manifest['status'] != 'complete' or manifest['precision']['tracks_omitted_by_cap']:
        raise ValueError(f'Incomplete or capped photon record: {run}')
    categories, fates, states, points = Counter(), Counter(), Counter(), Counter()
    ids = set()
    for file in run.glob('photons-w*.tsv'):
        with file.open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                key = (file.stem, row['event'], row['track'])
                if key in ids:
                    raise ValueError(f'Duplicate photon: {key}')
                ids.add(key)
                category = classify(row)
                categories[category] += 1
                fates[row['fate']] += 1
                states[row['track_status']] += 1
                points[(category, row['final_pre_volume'], row['final_post_volume'], row['final_material'])] += 1
    n = int(manifest['budget']['budget_started'])
    if len(ids) != n or sum(categories.values()) != n:
        raise ValueError('Photon accounting mismatch')
    for fate, count in fates.items():
        if count != int(manifest['budget']['budget_fate_'+fate]):
            raise ValueError(f'Fate accounting mismatch: {fate}')
    # These two SD branches explicitly call SetTrackStatus(fStopAndKill).
    # A generic Geant4 terminal status also includes physical absorption/world exit.
    explicit_sensor = fates['detected'] + fates['qe_reject']
    return dict(run=str(run), tracked=n, categories=dict(categories),
                fractions={k:v/n if n else 0. for k,v in categories.items()},
                track_status_counts=dict(states), sensor_direct_kill_count=explicit_sensor,
                sensor_direct_kill_fraction=explicit_sensor/n if n else 0.,
                terminal_boundaries=[dict(category=k[0],pre_volume=k[1],post_volume=k[2],material=k[3],count=v,
                                          fraction=v/n if n else 0.) for k,v in sorted(points.items())])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs',nargs='+',type=Path)
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args()
    result=[summarize(run.resolve()) for run in args.runs]
    with args.out.open('x') as f:
        json.dump(result,f,indent=2)
        f.write('\n')


if __name__=='__main__':
    main()
