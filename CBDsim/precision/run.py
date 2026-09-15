#!/usr/bin/env python3
"""Bounded, single-worker photon diagnostics; every invocation needs a NEW output dir."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

REPO = Path(__file__).resolve().parents[2]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def rows(path):
    with path.open() as stream:
        yield from csv.DictReader(stream, delimiter='\t')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--mode', choices=['lg', 'nolg'], required=True)
    p.add_argument('--round-tip', choices=['on','off'], default='on', help='default 0.2 mm cylindrical LG tip; off restores polygon outlet')
    p.add_argument('--tip-contact', choices=['on','off'], default='on',
                   help='LG sensor sleeve contact by default; off restores 10 um air gap')
    p.add_argument('--corner-seal', choices=['on','off'], default='on')
    p.add_argument('--tape-extensions', choices=['on','off'], default='on',
                   help='extend noLG sleeve over window; add LG inlet gel sleeve')
    p.add_argument('--gel-tape', choices=['absorber','off'], default='absorber',
                   help='absorbing sleeves; ideal R=T=0, 0.05 mm direct contact')
    p.add_argument('--level', choices=['off', 'summary', 'steps'], default='summary')
    p.add_argument('--events', type=int, default=1)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--particle', choices=['e+', 'e-', 'opticalphoton'], default='e+')
    p.add_argument('--energy', type=float, default=60.)
    p.add_argument('--energy-unit', choices=['GeV', 'MeV', 'eV'], default='GeV')
    p.add_argument('--position-mm', nargs=3, type=float, default=[0., 0., -100.])
    p.add_argument('--direction', nargs=3, type=float, default=[0., 0., 1.])
    p.add_argument('--polarization', nargs=3, type=float, default=[0., 0., 1.])
    p.add_argument('--max-tracks', type=int, default=100000)
    p.add_argument('--max-steps', type=int, default=100000)
    p.add_argument('--event', type=int, default=-1, help='record one event; -1 selects all')
    p.add_argument('--track', type=int, default=-1, help='step rows for this ID; -1 selects all recorded tracks')
    p.add_argument('--timeout', type=float, default=180.)
    p.add_argument('--root', action='store_true', help='also validate canonical ROOT timing/count payload')
    a = p.parse_args()
    import math
    if (a.events < 1 or a.max_tracks < 0 or a.max_steps < 0 or a.event < -1 or a.track < -1
            or not math.isfinite(a.timeout) or a.timeout <= 0
            or not math.isfinite(a.energy) or a.energy <= 0
            or not all(math.isfinite(v) for v in a.position_mm + a.direction + a.polarization)
            or not any(a.direction)):
        p.error('Invalid event/cap/filter/beam/timeout value')
    binary = a.binary.resolve(strict=True)
    rootio = next((parent/'rootIO' for parent in binary.parents
                   if (parent/'rootIO/librootIO.so').is_file()), None)
    if rootio is None:
        p.error('Cannot find rootIO/librootIO.so in the binary build ancestors')
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    env = {k: v for k, v in os.environ.items() if not k.startswith('CBDsim_')}
    env.update(CBDsim_PROTO_ROUND_TIP=str(int(a.round_tip == 'on')),
               CBDsim_PROTO_TIP_CONTACT=str(int(a.tip_contact == 'on')),
               CBDsim_PROTO_TAPE_EXTENSIONS=str(int(a.tape_extensions == 'on')),
               CBDsim_PROTO_CORNER_SEAL=str(int(a.corner_seal == 'on')),
               CBDsim_PROTO_GEL_TAPE=str(int(a.gel_tape == 'absorber')),
               CBDsim_PROTO_NO_LG=str(int(a.mode == 'nolg')), CBDsim_PROTO_SCAN_S_MM='0',
               CBDsim_OPTICAL_DIAG='1', CBDsim_OPTICAL_DIAG_OUT=str(out/'diagnostics.txt'))
    if a.level != 'off':
        env.update(CBDsim_PRECISION_OUT=str(out), CBDsim_PRECISION_MAX_TRACKS=str(a.max_tracks),
                   CBDsim_PRECISION_MAX_STEPS=str(a.max_steps if a.level == 'steps' else 0),
                   CBDsim_PRECISION_EVENT=str(a.event), CBDsim_PRECISION_TRACK=str(a.track))
    macro = (f'/run/numberOfThreads 1\n/run/initialize\n/gun/particle {a.particle}\n'
             f'/gun/energy {a.energy} {a.energy_unit}\n'
             f'/gun/position {" ".join(map(str, a.position_mm))} mm\n'
             f'/gun/direction {" ".join(map(str, a.direction))}\n'
             f'/gun/polarization {" ".join(map(str, a.polarization))}\n/tracking/verbose 0\n'
             f'/run/beamOn {a.events}\n')
    (out/'run.mac').write_text(macro)
    cmd = [str(binary), str(out/'run.mac'), str(a.seed)]
    if a.root:
        cmd.append(str(out/'events'))
        (out/'lib').mkdir()
        for name in ['librootIO.so', 'librootIO_rdict.pcm', 'librootIO.rootmap']:
            shutil.copy2(rootio/name, out/'lib'/name)
    norm = math.sqrt(sum(v*v for v in a.direction))
    manifest = dict(status='running', command=cmd, cwd=str(REPO),
                    settings={k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()},
                    binary_sha256=digest(binary), runner_sha256=digest(Path(__file__)),
                    rootio_sha256=digest(rootio/'librootIO.so'),
                    diagnostic_environment={k:v for k,v in env.items() if k.startswith('CBDsim_')},
                    primary_request=[a.energy*{'GeV':1000., 'MeV':1., 'eV':1.e-6}[a.energy_unit],
                                     *a.position_mm, *[v/norm for v in a.direction]],
                    source_sha256={str(f.relative_to(REPO)): digest(f) for folder in ['CBDsim','rootIO']
                                   for f in (REPO/folder).rglob('*')
                                   if f.is_file() and f.suffix in ['.cc','.hh','.h','.inc']})
    def save():
        (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    save()
    start = time.monotonic()
    try:
        with (out/'run.log').open('w') as log:
            proc = subprocess.Popen(cmd,
                                    cwd=REPO, env=env, stdout=log, stderr=subprocess.STDOUT,
                                    start_new_session=True)
            try:
                while True:
                    pid, status, usage = os.wait4(proc.pid, os.WNOHANG)
                    if pid:
                        code = os.waitstatus_to_exitcode(status)
                        proc.returncode = code
                        break
                    if time.monotonic() - start >= a.timeout:
                        raise subprocess.TimeoutExpired(cmd, a.timeout)
                    time.sleep(0.1)
            except (subprocess.TimeoutExpired, KeyboardInterrupt):
                import signal
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
                raise
        manifest['returncode'] = code
        if code:
            raise RuntimeError(f'Process failed: {code}')
        manifest.update(process_wall_seconds=time.monotonic()-start, peak_rss_kib=usage.ru_maxrss,
                        user_cpu_seconds=usage.ru_utime, system_cpu_seconds=usage.ru_stime)
        diag = dict(line.split() for line in (out/'diagnostics.txt').read_text().splitlines())
        if (int(diag['events']) != a.events or int(diag['budget_events']) != a.events
                or int(diag['budget_aborted_events']) != 0 or int(diag['budget_closed']) != 1):
            raise RuntimeError('Missing/aborted/unclosed events')
        manifest['budget'] = {k: v for k,v in diag.items() if k.startswith('budget_')}
        if a.level != 'off':
            records = [r for f in out.glob('events-w*.tsv') for r in rows(f)]
            if (len(records) != a.events or {int(r['event']) for r in records} != set(range(a.events))
                    or any(int(r['aborted']) for r in records)
                    or sum(int(r['seen']) for r in records) != int(diag['budget_started'])):
                raise RuntimeError('Precision event accounting failed or wrong executable')
            written = sum(int(r['tracks_written']) for r in records)
            selected = sum(int(r['selected']) for r in records)
            omitted = sum(int(r['omitted_by_track_cap']) for r in records)
            if written + omitted != selected or sum(1 for f in out.glob('photons-w*.tsv') for _ in rows(f)) != written:
                raise RuntimeError('Precision photon accounting failed')
            step_count = sum(1 for f in out.glob('steps-w*.tsv') for _ in rows(f))
            if step_count != max(int(r['steps_written_total']) for r in records):
                raise RuntimeError('Precision step accounting failed')
            manifest['precision'] = dict(tracks_written=written, tracks_omitted_by_cap=omitted,
                steps_written=step_count, steps_omitted_by_cap=max(int(r['steps_omitted_by_cap_total']) for r in records))
        if a.root:
            import ROOT
            ROOT.gSystem.AddDynamicPath(str(out/'lib'))
            sys.path.insert(0, str(REPO/'scripts'))
            from validate_frozen_run import validate
            manifest['root_validation'] = validate(out, a.events, a.seed)
        manifest['status'] = 'complete'
    except BaseException as exc:
        manifest.update(status='failed_or_incomplete', error=repr(exc))
        raise
    finally:
        manifest['elapsed_seconds'] = time.monotonic()-start
        manifest['trace_bytes'] = sum(f.stat().st_size for f in out.glob('*-w*.tsv'))
        save()
    print(json.dumps({k: manifest[k] for k in ['status','process_wall_seconds','peak_rss_kib','trace_bytes']}, indent=2))


if __name__ == '__main__':
    main()
