"""Run identical frozen inputs with two CLI builds and the common exchange probe.

Each executable runs once, sequentially, in its own directory. Timings include
output and are observations, not a performance benchmark. Requires Python 3.10+.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous', type=Path, required=True)
    parser.add_argument('--current', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    provenance = {'previous_commit': 'b9dcaf3', 'current_commit': '2bc4d3c',
                  'mpi_ranks': 1, 'runs': []}
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    for revision in ('previous', 'current'):
        binary = getattr(args, revision).resolve()
        for case in ('exchange', 'ferron', 'optical', 'air', 'timing'):
            run = work / case / revision
            run.mkdir(parents=True, exist_ok=True)
            if case == 'exchange':
                executable = binary.with_name('phase_exchange_probe' + binary.suffix)
                command = [str(executable), str(run / 'inventory.csv')]
                config_hash = digest(HERE / 'phase_exchange_probe.cpp')
            else:
                executable = binary
                config = (HERE / f'{case}.txt').read_text(encoding='utf-8')
                config_hash = hashlib.sha256(config.encode('utf-8')).hexdigest()
                config = re.sub(r'^output_dir\s*=.*$',
                                'output_dir = ' + (run / 'results').as_posix(),
                                config, flags=re.MULTILINE)
                (run / 'config.txt').write_text(config, encoding='utf-8')
                command = [str(executable), '-c', str(run / 'config.txt')]
            print(f'Running {case}/{revision}', flush=True)
            started = time.perf_counter()
            with (run / 'run.log').open('w', encoding='utf-8') as log:
                result = subprocess.run(command, cwd=run, env=env, stdout=log,
                                        stderr=subprocess.STDOUT, timeout=1800)
            record = dict(case=case, revision=revision, return_code=result.returncode,
                          elapsed_s=time.perf_counter() - started,
                          executable_sha256=digest(executable), input_sha256=config_hash)
            diagnostics = run / 'results' / 'diagnostics.csv'
            if diagnostics.exists():
                rows = list(csv.DictReader(diagnostics.open(encoding='utf-8')))
                record['diagnostic_rows'] = len(rows)
                record['all_reynolds_converged'] = all(row['converged'] == '1' for row in rows)
                checked = [row['coupling_converged'] for row in rows
                           if row.get('coupling_converged', '')]
                record['all_coupling_converged'] = all(value == '1' for value in checked) if checked else None
            provenance['runs'].append(record)
            (work / 'run_manifest.json').write_text(
                json.dumps(provenance, indent=2) + '\n', encoding='utf-8')
            print(record, flush=True)
            if result.returncode:
                raise RuntimeError(f'Failed run: {run}; see run.log')


if __name__ == '__main__':
    main()
