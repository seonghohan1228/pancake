"""Screen steady liquid temperature under explicit provisional cooling inputs.

This solves energy with fixed mu, rho and force-balanced bearing positions.
It is not a viscosity-coupled THD prediction. Raw configurations retain all
thermal assumptions and source files for review.
"""
import concurrent.futures
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import numpy as np
import run_scan as flow

ROOT=Path(__file__).resolve().parent
THERMAL={"inlet_temperature_K":293.15,"wall_temperature_K":293.15,
         "cp_J_kg_K":2000.,"conductivity_W_m_K":.13,
         "wall_htc_W_m2_K":{"adiabatic":0.,"cooled":10000.},
         "viscosity_feedback":False,"cavitated_film_model":"STRIATED"}


def evaluate_temperature(point,mode,use_fine=False):
    """Run one fixed-property energy solve and verify thermal and force closure."""
    load,rpm=point
    prefix="fine" if use_fine else "scan"
    row=json.loads((ROOT/f"{prefix}_n{rpm}_w{load}.json").read_text())
    config=(ROOT/row["dir"]/"config.txt").read_text()
    htc=THERMAL["wall_htc_W_m2_K"][mode]
    updates={"temperature_model":"ENERGY_EQUATION",
             "rho_cp":flow.CFG["rho"]*THERMAL["cp_J_kg_K"],
             "thermal_conductivity":THERMAL["conductivity_W_m_K"],
             "temperature_initial":THERMAL["inlet_temperature_K"],
             "temperature_reference":THERMAL["inlet_temperature_K"],
             "journal_wall_temperature":THERMAL["wall_temperature_K"],
             "bearing_wall_temperature":THERMAL["wall_temperature_K"],
             "journal_heat_transfer":htc,"bearing_heat_transfer":htc,
             "cavitated_film_model":THERMAL["cavitated_film_model"],
             "coupling_max_iters":200}
    for key,value in updates.items():
        pattern=rf"(?m)^{re.escape(key)}\s*=.*$"
        line=f"{key} = {value}"
        config=re.sub(pattern,line,config) if re.search(pattern,config) else config+"\n"+line
    config=re.sub(r"(?m)^inlet_groove\s*=.*$",
        f'inlet_groove = {flow.CFG["groove_angle"]} {flow.CFG["groove_width"]} {flow.CFG["p_end"]} {THERMAL["inlet_temperature_K"]}',config)
    identity=hashlib.sha256(config.encode()).hexdigest()[:14]
    name=f"thermal_{mode}{'_fine' if use_fine else ''}_n{rpm}_w{load}"
    run=ROOT/"runs"/name/identity
    result_path=ROOT/f"{name}.json"
    if result_path.exists():
        prior=json.loads(result_path.read_text())
        if prior.get("config_hash")==identity:
            return prior
    run.mkdir(parents=True,exist_ok=True)
    (run/"config.txt").write_text(config)
    start=time.monotonic()
    env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",
             OMPI_MCA_hwloc_base_binding_policy="none")
    with (run/"run.log").open("w") as output:
        p=subprocess.run(["mpiexec","-n",str(flow.CFG["ranks"]),flow.CFG["executable"],
                          "--config","config.txt"],cwd=run,stdout=output,
                          stderr=subprocess.STDOUT,env=env,timeout=240)
    log=(run/"run.log").read_text(errors="replace")
    residuals=re.findall(r"Coupling:.*?R_p=([\deE+.\-]+).*?R_T=([\deE+.\-]+)",log)
    if p.returncode or not residuals or max(map(float,residuals[-1]))>1e-6:
        raise RuntimeError(f"Thermal solve failed {run}: rc={p.returncode}, residual={residuals[-1:]}")
    fields=flow.read_fields(run)
    temperatures=fields["T"]
    if not np.isfinite(temperatures).all() or np.min(temperatures)<THERMAL["inlet_temperature_K"]-.01:
        raise RuntimeError(f"Invalid temperatures {run}")
    xs,zs=np.unique(fields["theta_rad"]),np.unique(fields["z_m"])
    planes={}
    for key in ("T","p","film_content","h","mu"):
        a=np.full((len(zs),len(xs)),np.nan)
        a[np.searchsorted(zs,fields["z_m"]),np.searchsorted(xs,fields["theta_rad"])]=fields[key]
        if not np.isfinite(a).all():
            raise RuntimeError(f"Incomplete {key} plane")
        planes[key]=a
    np.savez_compressed(run/"fields.npz",theta_deg=xs*180/math.pi,
                       z_over_L=zs/flow.CFG["L"]-.5,**planes)
    force=fields["Fp"][:2]
    target=np.array([0.,-load])
    force_error=float(np.linalg.norm(force-target)/load)
    if force_error>.002 or np.max(np.abs(fields["mu"]-flow.CFG["mu"]))>1e-6:
        raise RuntimeError(f"Fixed-property force field changed unexpectedly in {run}")
    record=dict(load_target=load,rpm=rpm,mode=mode,thermal_inputs=THERMAL,use_fine=use_fine,
                dir=run.relative_to(ROOT).as_posix(),config_hash=identity,
                T_min_C=float(temperatures.min()-273.15),T_max_C=float(temperatures.max()-273.15),
                T_mean_C=float(temperatures.mean()-273.15),
                rise_max_K=float(temperatures.max()-THERMAL["inlet_temperature_K"]),
                R_p=float(residuals[-1][0]),R_T=float(residuals[-1][1]),
                force_error_rel=force_error,power_W=float(fields["friction_torque"][0]*rpm*math.pi/30),
                wall_s=time.monotonic()-start)
    (run/"metrics.json").write_text(json.dumps(record,indent=2))
    result_path.write_text(json.dumps(record,indent=2))
    print(f"DONE {name}: Tmax={record['T_max_C']:.2f} C rise={record['rise_max_K']:.2f} K wall={record['wall_s']:.1f}s",flush=True)
    return record


def main():
    """Evaluate the displayed load/speed matrix for both cooling assumptions."""
    if len(sys.argv)>1 and sys.argv[1]=="refine":
        checks=[]
        for mode in ("cooled","adiabatic"):
            fine=evaluate_temperature((3000,2000),mode,use_fine=True)
            coarse=json.loads((ROOT/f"thermal_{mode}_n2000_w3000.json").read_text())
            checks.append(dict(mode=mode,T_max_coarse_C=coarse["T_max_C"],
                               T_max_fine_C=fine["T_max_C"],delta_K=fine["T_max_C"]-coarse["T_max_C"]))
        (ROOT/"temperature_grid_check.json").write_text(json.dumps(checks,indent=2))
        return
    points=[(w,n) for w in (1000,3000,5000) for n in (500,1000,2000)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(evaluate_temperature,point,mode)
                 for mode in ("cooled","adiabatic") for point in points]
        results=[future.result() for future in futures]
    (ROOT/"temperature_summary.json").write_text(json.dumps(results,indent=2))


if __name__=="__main__":
    main()
