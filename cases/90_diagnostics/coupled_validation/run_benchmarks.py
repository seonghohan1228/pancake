"""Independent JFO limit and Etsion/Ludwig photographed gas-cavity comparison."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import xml.etree.ElementTree as ET
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'physics_comparison'))
sys.path.insert(0,str(ROOT/'cases/05_etsion_extent'))
sys.path.insert(0,str(ROOT/'cases/00c_jfo_verification'))
from reduce_comparison import read_frame,write_csv
from conditions import render_config
from run_thermal import configure
import exact_1d_jfo as exact
from run_support import add_launch_arguments, machine_provenance, run_case

def extent(void,threshold):
    mask=np.any(void>threshold,axis=0)
    if mask.all(): return [0,0,360]
    if not mask.any(): return [0,0,0]
    starts=np.flatnonzero(mask & ~np.roll(mask,1)); n=len(mask)
    length=[]
    for start in starts:
        count=0
        while count<n and mask[(start+count)%n]: count+=1
        length.append(count)
    k=int(np.argmax(length)); start=int(starts[k]); count=length[k]
    return [start*360/n,(start+count)*360/n%360,count*360/n]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--evidence',type=Path,default=HERE/'benchmark_evidence')
    parser.add_argument('--only',default='',help='Run names containing this substring')
    add_launch_arguments(parser)
    args=parser.parse_args(); args.binary=args.binary.resolve(); args.work=args.work.resolve()
    args.evidence.mkdir(parents=True,exist_ok=True)
    manifest=dict(machine=machine_provenance(args),binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(),
        source_sha256={p.name:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in sorted((ROOT/'src').glob('*')) if p.is_file()},runs=[])
    tr=exact.solve_theta_r(); tf=exact.solve_theta_f(tr)
    theta,pressure=exact.pressure_profile(tr,tf)
    write_csv(args.evidence/'jfo_exact.csv',[dict(theta_deg=float(a*180/math.pi),pressure_Pa=float(p)) for a,p in zip(theta,pressure)])
    manifest['jfo_exact']=dict(rupture_deg=tr*180/math.pi,reformation_deg=tf*180/math.pi%360)
    scenarios=[]
    jfo='''R = .05
c = .000152
e = .0000912
L = .08
attitude_angle_deg = -90
omega = 418.879020479
mu = .0277
rho = 860
p_cav = 100000
bc_z_south_type = NEUMANN
bc_z_north_type = NEUMANN
fluid_property_model = CONSTANT
temperature_model = ISOTHERMAL
temperature_initial = 313.15
temperature_reference = 313.15
inlet_groove = 90 18 171300
mu_gas = 0
mixture_viscosity_model = DUKLER
solution_mode = STEADY_STATE
dt = .001
output_write_3d = false
output_write_flat = true
output_verbosity = DETAILED
output_fields = pressure, film_content, alpha_gas, h, temperature
diagnostics_interval = 1
'''
    for nt in [180,360,720,1440,2880]:
        scenarios.append((f'jfo_{nt}',configure(jfo,dict(n_theta_global=nt,n_z_global=4,cavitation_model='COUPLED'))))
    # The gas content was not measured in Etsion; retain the previously declared
    # atmospheric-charge assumption and compare kinetics without retuning it.
    for nt,nz,retained in [(120,32,False),(240,64,False),(120,32,True)]:
        template=render_config(1840,13.6,'fbns',n_theta=nt,n_z=nz,p_cav=13.33,cavitation_threshold='SCALAR_PCAV')
        values=dict(cavitation_model='COUPLED',coupled_phase_equilibrium='false' if retained else 'true',
            gas_absorption_rate=0 if retained else -1,mixture_viscosity_model='DUKLER',
            output_verbosity='DETAILED',output_fields='pressure, film_content, alpha_gas, h, temperature')
        scenarios.append((f'etsion_{nt}_'+('retained' if retained else 'equilibrium'),configure(template,values)))
    for name,template in scenarios:
        if args.only and args.only not in name: continue
        run=args.work/name; run.mkdir(parents=True,exist_ok=True)
        config=configure(template,dict(output_dir=(run/'results').as_posix()))
        (run/'config.txt').write_text(config)
        start=time.perf_counter()
        return_code=run_case(args,run/'config.txt',run/'run.log')
        rec=dict(name=name,return_code=return_code,elapsed_s=time.perf_counter()-start)
        for filename,source in [(f'{name}.txt',run/'config.txt'),(f'{name}.log',run/'run.log')]:
            (args.evidence/filename).write_bytes(source.read_bytes())
        if return_code==0:
            pvd=run/'results/flat/results.pvd'; item=ET.parse(pvd).findall('.//DataSet')[-1]
            frame=read_frame(pvd.parent/item.attrib['file']); nz,nt=frame['p'].shape
            angle=(np.arange(nt)+.5)*360/nt; j=nz//2
            write_csv(args.evidence/f'{name}_profile.csv',[dict(theta_deg=angle[i],pressure_Pa=frame['p'][j,i],
                liquid_fraction=frame['film_content'][j,i],max_void=float((1-frame['film_content'][:,i]).max())) for i in range(nt)])
            rec.update(p_min_Pa=float(frame['p'].min()),p_max_Pa=float(frame['p'].max()),
                cavity_angles={str(v):extent(1-frame['film_content'],v) for v in [1e-6,.001,.01,.05]})
            (args.evidence/f'{name}_diagnostics.csv').write_bytes((run/'results/diagnostics.csv').read_bytes())
            if name.startswith('jfo'):
                unwrapped=np.where(angle<99,angle+360,angle)
                land=(unwrapped>=99)&(unwrapped<=441)
                p_exact=np.interp(unwrapped[land]*math.pi/180,theta,pressure)
                rec['pressure_relative_l2']=float(np.linalg.norm(frame['p'][j,land]-p_exact)/np.linalg.norm(p_exact-exact.p_cav))
        manifest['runs'].append(rec)
        (args.evidence/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        print(json.dumps(rec),flush=True)
    return 1 if any(r['return_code'] for r in manifest['runs']) else 0

if __name__=='__main__': raise SystemExit(main())
