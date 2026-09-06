"""Check the repaired gas transient at 5 and 2.5 us, saving every actual step."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    args = parser.parse_args()
    binary = args.solver.resolve()
    root = args.work.resolve() / 'air_refinement'
    config = (HERE / 'air.txt').read_text(encoding='utf-8')
    results = []
    for name, dt in (('dt5us', 5e-6), ('dt2p5us', 2.5e-6)):
        run = root / name
        run.mkdir(parents=True, exist_ok=True)
        cfg = re.sub(r'^dt\s*=.*$', f'dt = {dt}', config, flags=re.MULTILINE)
        cfg = re.sub(r'^write_interval\s*=.*$', f'write_interval = {dt}', cfg, flags=re.MULTILINE)
        cfg = re.sub(r'^output_dir\s*=.*$', 'output_dir = ' + (run / 'results').as_posix(),
                     cfg, flags=re.MULTILINE)
        (run / 'config.txt').write_text(cfg, encoding='utf-8')
        start = time.perf_counter()
        with (run / 'run.log').open('w', encoding='utf-8') as log:
            result = subprocess.run([str(binary), '-c', str(run / 'config.txt')], cwd=run,
                env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'),
                stdout=log, stderr=subprocess.STDOUT, timeout=180)
        results.append(dict(name=name, dt=dt, return_code=result.returncode,
            elapsed_s=time.perf_counter() - start,
            executable_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
            input_sha256=hashlib.sha256((run / 'config.txt').read_bytes()).hexdigest()))
        (root / 'manifest.json').write_text(json.dumps(results, indent=2) + '\n', encoding='utf-8')
        if result.returncode:
            raise RuntimeError(f'Failed refinement: {run}')


if __name__ == '__main__':
    main()
