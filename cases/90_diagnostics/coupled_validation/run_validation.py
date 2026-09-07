"""Controlled isothermal solver comparison; all oil/air constants are synthetic.

Preserves inputs, binary/source hashes, return codes and reduced VTK histories.
Requires numpy and the existing comparison reader; no CoolProp runtime dependency.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import sys
import time
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / 'physics_comparison'))
from reduce_comparison import read_frame, measures, write_csv
from run_support import add_launch_arguments, machine_provenance, run_case


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--evidence', type=Path, default=HERE / 'evidence')
    add_launch_arguments(parser)
    args = parser.parse_args()
    args.work = args.work.resolve()
    args.evidence.mkdir(parents=True, exist_ok=True)
    template = (HERE.parent / 'physics_comparison/air.txt').read_text()
    common = {'temperature_model':'ISOTHERMAL', 'p_cav':'100',
              'gas_alpha_max':'1', 'mixture_viscosity_model':'DUKLER',
              'write_interval':'1e-5', 'output_verbosity':'DETAILED'}
    manifest = {'machine':machine_provenance(args), 'material':'Synthetic oil-air verification; not FW68D calibration',
                'binary_sha256':hashlib.sha256(args.binary.read_bytes()).hexdigest(),
                'source_sha256':{}, 'runs':[]}
    for path in sorted((ROOT / 'src').glob('*')):
        if path.is_file():
            manifest['source_sha256'][path.name] = hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    for model in ['FBNS','COUPLED']:
        for dt in [1e-5,5e-6,2.5e-6]:
            name=f'{model.lower()}_{dt*1e6:g}us'
            run=args.work/name
            run.mkdir(parents=True, exist_ok=True)
            changes=dict(common, cavitation_model=model, dt=str(dt), output_dir=(run/'results').as_posix())
            config=template
            for key,value in changes.items():
                if re.search(rf'^{key}\s*=',config,re.M):
                    config=re.sub(rf'^{key}\s*=.*$',f'{key} = {value}',config,flags=re.M)
                else: config+=f'\n{key} = {value}\n'
            (run/'config.txt').write_text(config)
            start=time.perf_counter()
            return_code=run_case(args,run/'config.txt',run/'run.log')
            record={'name':name,'model':model,'dt_s':dt,'return_code':return_code,
                    'elapsed_s':time.perf_counter()-start,
                    'input_sha256':hashlib.sha256(config.encode()).hexdigest()}
            if return_code:
                manifest['runs'].append(record)
                (args.evidence/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
                raise RuntimeError(f'{name} failed; see {run / "run.log"}')
            pvd=run/'results/flat/results.pvd'
            history=[dict(time_s=float(x.attrib['timestep']),**measures(read_frame(pvd.parent/x.attrib['file'])))
                     for x in ET.parse(pvd).findall('.//DataSet')]
            write_csv(args.evidence/f'{name}_history.csv',history)
            (args.evidence/f'{name}_diagnostics.csv').write_bytes((run/'results/diagnostics.csv').read_bytes())
            (args.evidence/f'{name}.log').write_bytes((run/'run.log').read_bytes())
            (args.evidence/f'{name}.txt').write_text(config)
            record['final']=history[-1]
            tail=[row['load_N'] for row in history if row['time_s']>=.0008]
            record['last_0p2ms_load_range_N']=[min(tail),max(tail)]
            diagnostics=list(csv.DictReader((run/'results/diagnostics.csv').open()))
            record['all_converged']=all(row['converged']=='1' for row in diagnostics)
            for key in ['liquid_residual','gas_residual','gas_volume_excess','gas_eos_defect']:
                record['max_abs_'+key]=max(abs(float(row[key])) for row in diagnostics)
            manifest['runs'].append(record)
            (args.evidence/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
            print(json.dumps(record),flush=True)


if __name__=='__main__': main()
