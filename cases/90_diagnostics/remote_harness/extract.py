import sys, json, math, pathlib
CASE = pathlib.Path('/home/user/pancake/cases/01_ferron_wallbc_fields')
sys.path.insert(0, str(CASE))
import solver_io as S

R, L = 0.05, 0.08
GROOVE = 90.0
eq = json.loads((CASE/'equilibrium.json').read_text())
out = {}

# favg has no cross-film wall fields; its fluxes are post-derived from the
# Newton-cooling law q = h (T_wall - T_bar), positive INTO the fluid.
FAVG_H = {'journal': 10000.0, 'bearing': 7200.0}
T_WALL = 313.15

for g, rec in eq.items():
    d = CASE / rec['flow']['dir']
    a = S._read_pvts(S.final_pvts(d))
    theta = a['theta_rad'][0]; z = a['z_m'][0]
    n_theta = len(set(round(t, 9) for t in theta))
    n_z = len(set(round(v, 9) for v in z))
    dA = (2*math.pi*R/n_theta) * (L/n_z)
    T = a['T'][0]
    have_xf = 'T_bearing' in a
    if have_xf:
        qb = a['q_bearing'][0]; qj = a['q_journal'][0]
        tb = a['T_bearing'][0]
    else:
        qb = [FAVG_H['bearing']*(T_WALL - t) for t in T]
        qj = [FAVG_H['journal']*(T_WALL - t) for t in T]
        tb = None
    # midplane profiles
    zs = sorted(set(z)); zmid = zs[len(zs)//2]
    rows = [(math.degrees(t) % 360.0, i) for i, (t, zv) in enumerate(zip(theta, z)) if zv == zmid]
    rows.sort()
    xs = [(x - GROOVE) % 360.0 for x, _ in rows]
    idx = [i for _, i in rows]
    order = sorted(range(len(xs)), key=lambda k: xs[k])
    prof = {
        'x_deg':   [xs[k] for k in order],
        'T_bar_C': [T[idx[k]] - 273.15 for k in order],
        'p_MPa':   [(a['p'][0][idx[k]] - 101325.0)/1e6 for k in order],
        'q_bearing': [qb[idx[k]] for k in order],
        'q_journal': [qj[idx[k]] for k in order],
    }
    if have_xf:
        prof['T_bearing_C'] = [tb[idx[k]] - 273.15 for k in order]
    out[g] = {
        'has_wall_fields': have_xf,
        'eps': rec['flow']['eps'], 'psi': rec['flow']['psi'],
        'metrics': {k: rec['flow']['metrics'].get(k) for k in
                    ('W','p_max','T_max','h_min','friction_torque','cav_frac',
                     'q_supply','force_angle_deg','coupling')},
        'newton_steps': rec['flow']['newton_steps'],
        # net wall heat, positive INTO the fluid; sign-flip gives film -> wall
        'net_q_bearing_W': sum(v*dA for v in qb),
        'net_q_journal_W': sum(v*dA for v in qj),
        'profile': prof,
    }
(CASE/'extract.json').write_text(json.dumps(out))
print("groups:", list(out))
for g, v in out.items():
    print("%-22s net q_bear=%9.1f W  net q_jour=%9.1f W  wall_fields=%s"
          % (g, v['net_q_bearing_W'], v['net_q_journal_W'], v['has_wall_fields']))
