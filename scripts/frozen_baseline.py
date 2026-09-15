#!/usr/bin/env python3
"""Build from a private source snapshot; run immutable-input, collision-safe central-beam smokes.
Use after `source envset.sh`. No git mutation, shared-build rebuild, or batch submission.
"""
import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
QE = Path('analysis/reference/pmt_r2076/R2076_quantum_efficiency_percent.csv')
ENV_KEYS = ('PATH', 'LD_LIBRARY_PATH', 'ROOTSYS', 'LCG_VERSION', 'CMAKE_PREFIX_PATH')


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, obj):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n')
    tmp.replace(path)


def check_files(base, hashes):
    for name, expected in hashes.items():
        if sha(base / name) != expected:
            raise RuntimeError('Frozen input changed: ' + name)


def run_logged(argv, cwd, log, env=None):
    with log.open('w') as f:
        subprocess.run(argv, cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT, check=True)


def libraries(binary, env):
    raw = subprocess.check_output(['ldd', str(binary)], env=env, text=True)
    if 'not found' in raw:
        raise RuntimeError('Unresolved dynamic libraries: ' + raw)
    paths = set(re.findall(r'(?:=>\s+|^\s*)(/\S+)\s+\(', raw, re.MULTILINE))
    return raw, {str(Path(x).resolve()): sha(x) for x in sorted(paths)}


def prepare(bundle, build_type=""):
    bundle.mkdir(parents=True, exist_ok=False)
    m = {'schema': 1, 'state': 'preparing', 'created_utc': now(), 'repository': str(REPO),
         'head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip(),
         'scope': 'nominal G01 geometry; 60 GeV central e-; one worker; R03 technical baseline'}
    save(bundle / 'manifest.json', m)
    try:
        for name, args in [('dirty.patch', ['diff', 'HEAD', '--binary']), ('index.patch', ['diff', '--cached', '--binary'])]:
            (bundle / name).write_bytes(subprocess.check_output(['git'] + args, cwd=REPO))
        source = bundle / 'source'
        paths = {Path('CMakeLists.txt'), Path('envset.sh'), Path('analysis/CMakeLists.txt'), QE,
                 Path('scripts/frozen_baseline.py'), Path('scripts/validate_frozen_run.py')}
        for directory in ['CBDsim', 'rootIO', 'Reco']:
            for f in (REPO / directory).rglob('*'):
                if f.is_file() and (f.name == 'CMakeLists.txt' or f.suffix in {'.cc', '.cpp', '.cxx', '.C', '.hh', '.hpp', '.h', '.mac', '.cmake', '.png', '.txt', '.sh', '.py', '.md'}):
                    paths.add(f.relative_to(REPO))
        for f in (REPO / 'analysis').rglob('*.py'):
            paths.add(f.relative_to(REPO))
        paths.update(f.relative_to(REPO) for f in (REPO / 'analysis/reference/pmt_r2076').glob('*.csv'))
        for rel in sorted(paths):
            dest = source / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(REPO / rel, dest)
        m['source_sha256'] = {str(x): sha(source / x) for x in sorted(paths)}
        m['environment'] = {k: v for k, v in os.environ.items() if k in ENV_KEYS or (k.startswith('G4') and k.endswith('DATA'))}
        m['build_commands'] = [['cmake', '-S', str(source), '-B', str(bundle / 'build'), '-DCMAKE_EXE_LINKER_FLAGS=-Wl,--enable-new-dtags'],
                               ['cmake', '--build', str(bundle / 'build'), '--target', 'CBDsim', '-j2']]
        if build_type:
            m['build_commands'][0].append('-DCMAKE_BUILD_TYPE=' + build_type)
        m['build_type'] = build_type
        save(bundle / 'manifest.json', m)
        for cmd, log in zip(m['build_commands'], ['configure.log', 'build.log']):
            run_logged(cmd, bundle, bundle / log)
        check_files(source, m['source_sha256'])
        binary = bundle / 'build/CBDsim/CBDsim'
        raw, libs = libraries(binary, os.environ)
        (bundle / 'ldd.txt').write_text(raw)
        artifacts = [binary, bundle / 'build/rootIO/librootIO.so']
        artifacts += list((bundle / 'build/rootIO').glob('*.pcm'))
        m['artifact_sha256'] = {str(x.relative_to(bundle)): sha(x) for x in artifacts}
        m['linked_libraries_sha256'] = libs
        m['state'] = 'ready'
        m['finished_utc'] = now()
        save(bundle / 'manifest.json', m)
    except Exception as e:
        m.update(state='failed', error=str(e), finished_utc=now())
        save(bundle / 'manifest.json', m)
        raise


def execute(bundle, run_id, mode, seed, events, position=(0.,0.,0.), direction=(0.,0.,1.), scan_s=None, primary_entry_audit=False):
    if not re.fullmatch(r'[A-Za-z0-9_-]+', run_id) or events < 1 or seed < 1:
        raise ValueError('Use a simple run ID, positive seed and positive event count')
    if not all(math.isfinite(v) for v in (*position,*direction)) or math.hypot(*direction)==0:
        raise ValueError('Beam coordinates must be finite and direction nonzero')
    if scan_s is not None and (not math.isfinite(scan_s) or not 0 <= scan_s <= 29.5):
        raise ValueError('Assembly scan s must be in [0,29.5] mm')
    direction=tuple(v/math.hypot(*direction) for v in direction)
    m = json.loads((bundle / 'manifest.json').read_text())
    if m['state'] != 'ready':
        raise RuntimeError('Bundle is not ready')
    check_files(bundle / 'source', m['source_sha256'])
    check_files(bundle, m['artifact_sha256'])
    for name, expected in m['linked_libraries_sha256'].items():
        if sha(name) != expected:
            raise RuntimeError('Linked library changed: ' + name)
    run = bundle / 'runs' / run_id
    run.mkdir(parents=True, exist_ok=False)
    r = {'state': 'preparing', 'start_utc': now(), 'mode': mode, 'seed': seed, 'requested_events': events,
         'primary_request': [60000., *position, *direction],
         'assembly_scan_s_mm': scan_s, 'primary_entry_audit': primary_entry_audit,
         'threads': 1, 'bundle_manifest_sha256': sha(bundle / 'manifest.json'), 'analysis_cli': None,
         'purpose': 'reproducibility smoke; final physics baseline waits for B01/M01/M02/W01/W02/A01'}
    save(run / 'manifest.json', r)
    try:
        (run / 'lib').mkdir()
        shutil.copy2(bundle / 'build/CBDsim/CBDsim', run / 'CBDsim')
        for f in (bundle / 'build/rootIO').iterdir():
            if f.name == 'librootIO.so' or f.suffix == '.pcm':
                shutil.copy2(f, run / 'lib' / f.name)
        (run / QE).parent.mkdir(parents=True)
        shutil.copy2(bundle / 'source' / QE, run / QE)
        macro = '/vis/disable\n/run/numberOfThreads 1\n/run/initialize\n/run/verbose 0\n/gun/particle e-\n/gun/energy 60 GeV\n/run/beamOn ' + str(events) + '\n'
        macro=macro.replace('/run/beamOn ', '/gun/position '+ ' '.join(format(v,'.17g') for v in position)+' mm\n/gun/direction '+ ' '.join(format(v,'.17g') for v in direction)+'\n/run/beamOn ')
        (run / 'run.mac').write_text(macro)
        env = os.environ.copy()
        for k in list(env):
            if k.startswith('CBDsim_') or k.startswith('CBDSIM_') or k in {'G4RUN_MANAGER_TYPE', 'G4FORCE_RUN_MANAGER_TYPE', 'G4FORCENUMBEROFTHREADS', 'LD_PRELOAD'}:
                env.pop(k)
        env.update(CBDsim_PROTO_NO_LG='1' if mode == 'nolg' else '0', CBDsim_OPTICAL_DIAG='1',
                   CBDsim_OPTICAL_DIAG_HIST='1', CBDsim_OPTICAL_DIAG_OUT=str(run / 'diagnostics.txt'),
                   CBDsim_OPTICAL_DIAG_HIST_OUT=str(run / 'path.hist'), G4RUN_MANAGER_TYPE='MT')
        if scan_s is not None:
            env['CBDsim_PROTO_SCAN_S_MM'] = format(scan_s, '.17g')
        if primary_entry_audit:
            env['CBDsim_PRIMARY_ENTRY_AUDIT'] = '1'
        env['LD_LIBRARY_PATH'] = str(run / 'lib') + ':' + env.get('LD_LIBRARY_PATH', '')
        raw, libs = libraries(run / 'CBDsim', env)
        if str((run / 'lib/librootIO.so').resolve()) not in libs:
            raise RuntimeError('Run does not resolve its private rootIO library')
        (run / 'ldd.txt').write_text(raw)
        r['environment'] = {k: v for k, v in env.items() if k in ENV_KEYS or k.startswith('CBDsim_') or k == 'G4RUN_MANAGER_TYPE' or (k.startswith('G4') and k.endswith('DATA'))}
        expected_libs = set(m['linked_libraries_sha256'].values())
        if any(value not in expected_libs for value in libs.values()):
            raise RuntimeError('Runtime library set differs from frozen build environment')
        expected_data = {k: v for k, v in m['environment'].items() if k.startswith('G4') and k.endswith('DATA')}
        if any(env.get(k) != v for k, v in expected_data.items()):
            raise RuntimeError('Geant4 data environment differs from frozen bundle')
        r['libraries_sha256'] = libs
        inputs = [run / 'CBDsim', run / QE, run / 'run.mac'] + list((run / 'lib').iterdir())
        r['input_sha256'] = {str(f.relative_to(run)): sha(f) for f in inputs}
        r['command'] = [str(run / 'CBDsim'), 'run.mac', str(seed), 'events']
        r['validation_cli'] = [sys.executable, str(bundle / 'source/scripts/validate_frozen_run.py'), str(run), str(events), str(seed)]
        r['state'] = 'running'
        save(run / 'manifest.json', r)
        run_logged(r['command'], run, run / 'simulation.log', env)
        log = (run / 'simulation.log').read_text()
        if scan_s is not None:
            actual = re.findall(r'\[Proto scan\] s_mm=([0-9.eE+-]+)', log)
            if not actual or any(abs(float(v)-scan_s)>1e-6 for v in actual):
                raise RuntimeError('Actual assembly translation not confirmed')
        expected = '=> no-LG' if mode == 'nolg' else '=> LG'
        if expected not in log or '=== [end optical diagnostics] ===' not in log:
            raise RuntimeError('Actual geometry or run completion not confirmed in log')
        if any(x in log for x in ['GeomNav', 'GeomVol', 'GeomSolids', 'COMMAND NOT FOUND', '***** Illegal']):
            raise RuntimeError('Geometry/navigation/macro error in simulation log')
        run_logged(r['validation_cli'], run, run / 'validation.log', env)
        check_files(run, r['input_sha256'])
        for name, expected_sha in libs.items():
            if sha(name) != expected_sha:
                raise RuntimeError('Runtime library changed during run: ' + name)
        check_files(bundle / 'source', m['source_sha256'])
        r['validation'] = json.loads((run / 'validation.json').read_text())
        r['completed_events'] = r['validation']['events']
        r['output_sha256'] = {f.name: sha(f) for f in run.iterdir() if f.is_file() and f.name not in {'manifest.json', 'CBDsim', 'run.mac'}}
        r.update(state='complete', exit_code=0, finish_utc=now())
        save(run / 'manifest.json', r)
        print(json.dumps({'run': str(run), 'state': r['state'], **r['validation']}))
    except Exception as e:
        r.update(state='failed', error=str(e), finish_utc=now())
        save(run / 'manifest.json', r)
        raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='action', required=True)
    prep = sub.add_parser('prepare'); prep.add_argument('bundle', type=Path)
    prep.add_argument('--build-type', choices=['', 'Release'], default='')
    run = sub.add_parser('run'); run.add_argument('bundle', type=Path); run.add_argument('run_id')
    run.add_argument('--mode', choices=['lg', 'nolg'], required=True)
    run.add_argument('--seed', type=int, default=42); run.add_argument('--events', type=int, default=2)
    run.add_argument('--position-mm', nargs=3, type=float, default=(0.,0.,0.))
    run.add_argument('--direction', nargs=3, type=float, default=(0.,0.,1.))
    run.add_argument('--scan-s-mm', type=float)
    run.add_argument('--primary-entry-audit', action='store_true')
    a = p.parse_args(); bundle = a.bundle.resolve()
    if a.action == 'prepare':
        prepare(bundle, a.build_type)
    else:
        execute(bundle, a.run_id, a.mode, a.seed, a.events, a.position_mm, a.direction, a.scan_s_mm, a.primary_entry_audit)


if __name__ == '__main__':
    main()
