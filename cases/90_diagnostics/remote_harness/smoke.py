import sys, pathlib, subprocess, os
d = pathlib.Path('/home/user/pancake/cases/01_ferron_wallbc_fields')
sys.path.insert(0, str(d))
import force_balance as F
tpl = (d/'config_template.txt').read_text()
for g, spec in F.GROUPS.items():
    w = pathlib.Path('/home/user/smoke')/g
    w.mkdir(parents=True, exist_ok=True)
    (w/'config.txt').write_text(tpl.format(
        e='9.0e-05', attitude='-37.0', n_theta=180, n_z=24,
        thermal_block=spec['thermal'], output_fields=spec['fields'],
        q_supply='1e-4', groove_mixing_model='PINNED',
        coupling_max_iters=15, coupling_relaxation='1.0'))
    r = subprocess.run(['mpirun','-n','8','-x','OMP_NUM_THREADS=1','--oversubscribe',
                        '/home/user/pancake/build-linux/pancake','--config','config.txt'],
                       cwd=str(w), capture_output=True, text=True)
    log = (r.stdout or '') + (r.stderr or '')
    bad = [l for l in log.splitlines() if 'ERROR' in l or 'FATAL' in l or 'error:' in l.lower()]
    print("%-22s rc=%d %s" % (g, r.returncode, ('FAIL: '+bad[0][:110]) if bad else 'ok'))
