"""Thermal bearing verification matrix with explicit synthetic material inputs.

These are sensitivity/implementation cases, not an FW68D property fit.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
import time
import xml.etree.ElementTree as ET
import csv

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'physics_comparison'))
from reduce_comparison import read_frame, measures, write_csv
from run_support import add_launch_arguments, machine_provenance, run_case

def configure(template,values):
    for key,value in values.items():
        line=f'{key} = {value}'
        if re.search(rf'^{key}\s*=',template,re.M):
            template=re.sub(rf'^{key}\s*=.*$',line,template,flags=re.M)
        else: template+='\n'+line+'\n'
    return template

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--evidence',type=Path,default=HERE/'thermal_evidence')
    parser.add_argument('--smoke',action='store_true')
    parser.add_argument('--only',default='',help='Run names containing this substring')
    parser.add_argument('--r410a-equilibrium',action='store_true')
    parser.add_argument('--end-time',type=float,help='Override the physical duration (s)')
    add_launch_arguments(parser)
    args=parser.parse_args()
    args.binary=args.binary.resolve(); args.work=args.work.resolve(); args.evidence.mkdir(parents=True,exist_ok=True)
    template=(HERE.parent/'physics_comparison/air.txt').read_text()
    common=dict(cavitation_model='COUPLED',p_cav=100,gas_alpha_max=1,mixture_viscosity_model='DUKLER',
        output_verbosity='DETAILED',omega=600,mu=.03,liquid_viscosity_model='EMPIRICAL',
        viscosity_temperature_coeff=2500,solution_viscosity_gas_coeff=-3,
        property_reference_pressure=1e5,property_reference_temperature=313.15,
        gas_heat_capacity_cv=717.5,dissolved_heat_capacity=2000,desorption_internal_energy=1e5,
        gas_thermal_conductivity=.02,oil_vapor_temperature_coeff=3000,
        dissolved_gas_henry_temp_coeff=1000,gas_mass_transfer_rate=20,gas_absorption_rate=1,
        bc_z_south_thermal='CONSTANT',bc_z_north_thermal='CONSTANT',
        end_t=.002 if args.smoke else .02,write_interval=.0002,dt=.0001,
        output_fields='pressure, film_content, temperature, dissolved_gas, alpha_gas, rho_liquid_solution, h, mu, mass_residual, energy_residual, internal_energy, gas_mass_transfer, fluid_force_x, fluid_force_y, friction_torque')
    if args.end_time is not None:
        if args.end_time <= 0: raise ValueError('end-time must be positive')
        common['end_t']=args.end_time
    scenarios=[]
    for species in ['AIR','R410A']:
        for mode in ['ISOTHERMAL','ENERGY_EQUATION']:
            scenarios.append((f'{species.lower()}_{mode.lower()}',species,mode,313.15,1e5,48,8,.0001))
    if not args.smoke:
        scenarios += [(f'air_thermal_dt{dt:g}','AIR','ENERGY_EQUATION',313.15,1e5,48,8,dt) for dt in [.00005,.000025]]
        scenarios += [(f'air_thermal_grid{nt}','AIR','ENERGY_EQUATION',313.15,1e5,nt,nz,.0001) for nt,nz in [(24,4),(96,16),(192,32),(384,64)]]
        scenarios += [(f'r410a_{t-273.15:.0f}C_20bar','R410A','ENERGY_EQUATION',t,2e6,48,8,.0001) for t in [293.15,313.15,333.15,353.15]]
    manifest=dict(machine=machine_provenance(args),material='Synthetic oil properties; FW68D prediction is not validated',
        binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(),
        source_sha256={p.name:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in sorted((ROOT/'src').glob('*')) if p.is_file()},runs=[])
    for name,species,mode,t,pout,nt,nz,dt in scenarios:
        if args.only and args.only not in name: continue
        run=args.work/name; run.mkdir(parents=True,exist_ok=True)
        values=dict(common,dissolved_gas_species=species,temperature_model=mode,temperature_initial=t,
            coupled_phase_equilibrium='true' if species=='R410A' and args.r410a_equilibrium else 'false',
            temperature_reference=t,journal_wall_temperature=t,bearing_wall_temperature=t,
            bc_z_south_temperature=t,bc_z_north_temperature=t,bc_z_south_val=pout,bc_z_north_val=pout,
            dissolved_gas_initial=.002 if species=='AIR' else .1,
            dissolved_gas_henry_coeff=1e-8 if species=='AIR' else 1e-7,
            inlet_groove=f'90 30 {2e5 if species=="AIR" else 1e6} {t}',
            output_dir=(run/'results').as_posix(),n_theta_global=nt,n_z_global=nz,dt=dt)
        config=configure(template,values); (run/'config.txt').write_text(config)
        start=time.perf_counter()
        return_code=run_case(args,run/'config.txt',run/'run.log')
        record=dict(name=name,species=species,mode=mode,temperature_K=t,outlet_Pa=pout,
            n_theta=nt,n_z=nz,dt_s=dt,phase_equilibrium=species=='R410A' and args.r410a_equilibrium,
            elapsed_s=time.perf_counter()-start,return_code=return_code)
        (args.evidence/f'{name}.txt').write_text(config)
        (args.evidence/f'{name}.log').write_bytes((run/'run.log').read_bytes())
        if return_code==0:
            pvd=run/'results/flat/results.pvd'
            history=[dict(time_s=float(x.attrib['timestep']),**measures(read_frame(pvd.parent/x.attrib['file']))) for x in ET.parse(pvd).findall('.//DataSet')]
            write_csv(args.evidence/f'{name}_history.csv',history); record['final']=history[-1]
            (args.evidence/f'{name}_diagnostics.csv').write_bytes((run/'results/diagnostics.csv').read_bytes())
            rows=list(csv.DictReader((run/'results/diagnostics.csv').open()))
            for key in ['liquid_residual','gas_residual','gas_volume_excess','energy_residual_W']:
                vals=[abs(float(r[key])) for r in rows if r[key]]
                record['max_abs_'+key]=max(vals) if vals else None
        manifest['runs'].append(record)
        (args.evidence/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        print(json.dumps(record),flush=True)
        if return_code: print((run/'run.log').read_text()[-1800:],flush=True)
    return 1 if any(r['return_code'] for r in manifest['runs']) else 0

if __name__=='__main__': raise SystemExit(main())
