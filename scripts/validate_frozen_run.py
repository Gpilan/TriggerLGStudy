#!/usr/bin/env python3
"""Validate a small frozen CBDsim run; produce a canonical event digest, not timing claims."""
import hashlib
import json
import math
import struct
from pathlib import Path
import sys


def primary_matches(actual, requested):
    """The ROOT schema stores float32. Directions also undergo gun normalization."""
    if len(actual)!=7 or len(requested)!=7 or not all(math.isfinite(v) for v in (*actual,*requested)):
        return False
    expected=[struct.unpack('f',struct.pack('f',v))[0] for v in requested]
    return (actual[:4]==expected[:4] and
            all(math.isclose(a,b,rel_tol=2**-23,abs_tol=2**-149)
                for a,b in zip(actual[4:],expected[4:])))


def validate(run, events, seed):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    if ROOT.gSystem.Load(str(run / 'lib/librootIO.so')) < 0:
        raise RuntimeError('Cannot load the frozen ROOT dictionary')
    f = ROOT.TFile.Open(str(run / f'events_{seed}.root'), 'READ')
    if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError('Missing, corrupt or recovered ROOT output')
    tree = f.Get('CBDsim')
    if not tree or tree.GetEntries() != events:
        raise RuntimeError('ROOT event count mismatch')
    diag = {k: float(v) for k, v in (x.split() for x in (run / 'diagnostics.txt').read_text().splitlines())}
    for k in ['fate_unknown', 'fate_no_rindex', 'fate_unfinished', 'balance_delta', 'missing_event_records']:
        if diag.get('budget_' + k) != 0:
            raise RuntimeError('Nonzero or absent diagnostic: ' + k)
    if diag.get('budget_closed') != 1 or diag.get('events') != events:
        raise RuntimeError('Unclosed diagnostic budget')
    budgets = {}
    for line in (run / 'diagnostics.txt.events.txt').read_text().splitlines():
        if line.startswith('#'):
            continue
        parts = line.split()
        eid = int(parts[1])
        if eid in budgets:
            raise RuntimeError('Duplicate event ID')
        budgets[eid] = {k: float(v) for k, v in (x.split('=') for x in parts[2:])}
    if set(budgets) != set(range(events)):
        raise RuntimeError('Missing event budget IDs')
    if len((run / 'diagnostics.txt.exceptions.txt').read_text().splitlines()) != 1:
        raise RuntimeError('Unexpected optical exceptions')
    request=json.loads((run/'manifest.json').read_text()).get('primary_request',[60000.,0.,0.,0.,0.,0.,1.])
    records = []
    for i in range(events):
        tree.GetEntry(i)
        ev = tree.CBDsimEventData
        eid = int(ev.event_number)
        kin = [float(getattr(ev, k)) for k in ['primaryEkin', 'primaryVx', 'primaryVy', 'primaryVz', 'primaryDirX', 'primaryDirY', 'primaryDirZ']]
        if not primary_matches(kin,request):
            raise RuntimeError('Primary kinematics differ from the frozen beam request')
        counts = [int(ev.siPMPhotonSumTrig0), int(ev.siPMPhotonSumTrig1)]
        b = budgets[eid]
        if sum(counts) != b['fate_detected']:
            raise RuntimeError('ROOT/diagnostic detected count mismatch')
        if sum(v for k, v in b.items() if k.startswith('source_')) != b['started']:
            raise RuntimeError('Source budget mismatch')
        if sum(v for k, v in b.items() if k.startswith('fate_')) != b['started'] or b['started'] != b['terminal']:
            raise RuntimeError('Fate budget mismatch')
        hist = []
        for trig in (0, 1):
            c = getattr(ev, f'timeMergedCountsTrig{trig}')
            lo = getattr(ev, f'timeMergedEdgeLowTrig{trig}')
            hi = getattr(ev, f'timeMergedEdgeHighTrig{trig}')
            if len(c) != len(lo) or len(c) != len(hi):
                raise RuntimeError('Invalid ROOT time arrays')
            hist.append([[j, int(n), float(lo[j]), float(hi[j])] for j, n in enumerate(c) if n])
        records.append({'event': eid, 'primary': kin, 'detected': counts, 'time_bins': hist, 'budget': b})
    f.Close()
    records.sort(key=lambda x: x['event'])
    encoded = json.dumps(records, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
    result = {'events': events, 'event_digest_sha256': hashlib.sha256(encoded).hexdigest(),
              'tracked': int(diag['budget_started']), 'detected': int(diag['budget_fate_detected']),
              'purpose': 'reproducibility smoke; no timing or efficiency precision claim'}
    (run / 'canonical-events.json').write_bytes(encoded)
    (run / 'validation.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    print(json.dumps(validate(Path(sys.argv[1]).resolve(), int(sys.argv[2]), int(sys.argv[3]))))
