"""Gas effect at the compressor condition: the three B0 steady cases with the VC-only model.

VC only = the B0 (GC+VC) configuration with cavitation_threshold = SCALAR_PCAV and
p_cav = 500 Pa (as case P4v of lge-fbns-tester-20261117): "the refrigerant stays dissolved".
The dissolved-gas keys stay; under SCALAR_PCAV with SINGLE_PHASE and constant properties they
have no effect, which each run checks (dissolved fraction uniform, no free gas).

Cases:
  R-mean     GC+VC at the mean load: reproduction of B0 with this binary and machine
  V-min/mean/max   VC only at the three B0 loads (results)
  V-mean-p1  VC only, p_cav = 1 Pa (vapour-pressure sensitivity, nothing else changed)
  V-mean-p0  VC only, p_cav = 0 Pa (the solver accepts 0 only without dissolved gas, so the
             dissolved fraction is 0 in this run)
Newton on (eps, psi) over STATIC solves as B0 (|F|/W - 1 < 1e-4, angle < 0.01 deg), started from
the B0 solution of the same load.

Usage: python vc_steady.py --work NEW_CAMPAIGN --ref B0_CAMPAIGN --solver EXE [--ranks 24 --mpi mpirun]
Writes summary.json, gc_effect.csv and runs/<case>/... in NEW_CAMPAIGN. Reads REF only.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

import common
import part_b

VC = dict(cavitation_threshold="SCALAR_PCAV", p_cav=500.0)
FIELDS = ("pressure,film_content,h,mu,inlet_indicator,dissolved_gas,alpha_gas,fluid_force_x,fluid_force_y,"
          "fluid_force_z,friction_torque,bearing_x,bearing_y,bearing_z")
REPRO = dict(eps=0.943750, h_min_m=1.68766e-6, torque_Nm=0.077335)   # B0 mean values (prompt)


def details(run_dir, accepted, params, load):
    """Quantities of the gc_effect table for one accepted STATIC solve."""
    run_dir = Path(run_dir)
    fields = common.read_fields(run_dir)
    diag = common.steady_driver.last_diagnostic(run_dir)
    cav = part_b.cavity_report(fields)
    torque = accepted["torque_Nm"]
    out = dict(
        load_N=load, eps=accepted["eps"], attitude_deg=accepted["psi_deg"],
        h_min_um=accepted["h_min_m"] * 1e6, h_min_angle_deg=cav["theta_hmin_deg"],
        h_min_angle_from_hole_deg=common.wrap180(cav["theta_hmin_deg"] - part_b.HOLE_THETA),
        torque_Nm=torque, friction_power_W=torque * part_b.OMEGA,
        p_max_bar=accepted["p_max_Pa"] / 1e5, p_min_bar=accepted["p_min_Pa"] / 1e5, p_min_Pa=accepted["p_min_Pa"],
        cavitated_area_pct=accepted["cavitated_fraction"] * 100,
        cavity_start_from_hole_deg=cav.get("bubble_start_from_hole_deg"),
        cavity_end_from_hole_deg=cav.get("bubble_end_from_hole_deg"),
        film_content_hole=1 - cav["void_hole_mean"], film_content_ring_mean=1 - cav["void_ring_mean"],
        film_content_ring_min=1 - cav["void_ring_max"], bubble_touches_hole_ring=cav["bubble_touches_hole_ring"],
        bubble_gap_to_hole_edge_mm=cav.get("bubble_gap_to_hole_edge_mm"),
        hole_mass_flow_g_s=float(diag["inlet_mass_source"]) * 1e3,
        liquid_residual_rel=float(diag["liquid_residual_rel"]), valid=accepted.get("valid", True),
        force_mag_err=accepted.get("mag_err"), force_angle_err_deg=accepted.get("ang_err_deg"),
        run_dir=str(run_dir))
    if "dissolved_gas" in fields:
        d = fields["dissolved_gas"]
        out.update(dissolved_min=float(d.min()), dissolved_max=float(d.max()))
    if "alpha_gas" in fields:
        out["alpha_gas_max"] = float(np.max(fields["alpha_gas"]))
    out["vapour_cavity"] = bool(out["p_min_Pa"] <= accepted.get("p_cav_Pa", 0.0) + 1.0)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--ref", required=True, type=Path)
    parser.add_argument("--solver", required=True)
    parser.add_argument("--ranks", type=int, default=24)
    parser.add_argument("--mpi", default="mpirun")
    parser.add_argument("--rankfile")
    args = parser.parse_args()
    args.work = args.work.resolve()
    args.work.mkdir(parents=True, exist_ok=True)
    params = json.loads((args.ref / "parameters.json").read_text())
    common.save_json(args.work / "parameters.json", dict(params, source=str(args.ref / "parameters.json")))
    b0 = json.loads((args.ref / "b0" / "summary.json").read_text())
    fields_out = dict(output_fields=FIELDS)
    no_gas = dict(dissolved_gas_initial=0.0, bc_z_south_dissolved_gas=0.0, bc_z_north_dissolved_gas=0.0,
                  inlet_species="0 0 0")
    cases = [
        ("R-mean", "mean", dict(fields_out), None),
        ("V-min", "min", dict(fields_out, **VC), 500.0),
        ("V-mean", "mean", dict(fields_out, **VC), 500.0),
        ("V-max", "max", dict(fields_out, **VC), 500.0),
        ("V-mean-p1", "mean", dict(fields_out, cavitation_threshold="SCALAR_PCAV", p_cav=1.0), 1.0),
        ("V-mean-p0", "mean", dict(fields_out, cavitation_threshold="SCALAR_PCAV", p_cav=0.0, **no_gas), 0.0),
    ]
    summary = dict(reference=str(args.ref), solver=str(args.solver), solver_sha256=common.steady_driver.sha256(args.solver),
                   cases={})
    for name, load_key, overrides, p_cav in cases:
        ref = b0[load_key]
        load = ref["load_N"]
        seed = (ref["accepted"]["eps"], ref["accepted"]["psi_deg"])
        accepted, history = part_b.balance(args, params, name, load, *seed, overrides=overrides,
                                           root=args.work / "runs")
        accepted["p_cav_Pa"] = p_cav if p_cav is not None else params["p_sat_Pa"]
        row = details(accepted["run_dir"], accepted, params, load)
        row.update(case=name, model="GC+VC" if p_cav is None else f"VC p_cav={p_cav:g} Pa", solves=len(history))
        summary["cases"][name] = dict(result=row, history=history)
        common.save_json(args.work / "summary.json", summary)
        print(f"[{name}] eps={row['eps']:.6f} h_min={row['h_min_um']:.5f} um T={row['torque_Nm']:.6f} "
              f"p_min={row['p_min_Pa']:.1f} Pa p_max={row['p_max_bar']:.1f} bar cav={row['cavitated_area_pct']:.1f} % "
              f"flow={row['hole_mass_flow_g_s']:.4f} g/s solves={row['solves']}", flush=True)
        if name == "R-mean":
            b0_mean = b0["mean"]["accepted"]
            target = dict(eps=b0_mean["eps"], h_min_m=b0_mean["h_min_m"], torque_Nm=b0_mean["torque_Nm"])
            repro = {k: row_value / target[k] - 1 for k, row_value in
                     (("eps", row["eps"]), ("h_min_m", row["h_min_um"] * 1e-6), ("torque_Nm", row["torque_Nm"]))}
            summary["reproduction"] = dict(target_b0=target, prompt_rounded=REPRO, rel_diff=repro,
                                           passed=all(abs(v) < 1e-4 for v in repro.values()))
            common.save_json(args.work / "summary.json", summary)
            if not summary["reproduction"]["passed"]:
                raise SystemExit(f"Reproduction check failed: {repro}")
    # GC+VC rows from the B0 accepted solves (reference campaign, read only).
    table = []
    for load_key, vc_name in (("min", "V-min"), ("mean", "V-mean"), ("max", "V-max")):
        ref = b0[load_key]
        gc = details(ref["accepted"]["run_dir"], dict(ref["accepted"], p_cav_Pa=params["p_sat_Pa"]), params, ref["load_N"])
        gc.update(case=f"B0-{load_key}", model="GC+VC (B0)")
        vc = summary["cases"][vc_name]["result"]
        for row in (gc, vc):
            table.append(dict(row, load_case=load_key))
        diff = {f"{q}_gc_minus_vc_pct": (gc[q] / vc[q] - 1) * 100 for q in ("h_min_um", "torque_Nm", "p_max_bar",
                                                                              "hole_mass_flow_g_s")}
        table.append(dict(load_case=load_key, model="(GC+VC - VC)/VC [%]", load_N=ref["load_N"], **diff))
        summary.setdefault("gc_effect", {})[load_key] = diff
    p1 = summary["cases"]["V-mean-p1"]["result"]["h_min_um"]
    p0 = summary["cases"]["V-mean-p0"]["result"]["h_min_um"]
    vm = summary["cases"]["V-mean"]["result"]["h_min_um"]
    summary["pcav_sensitivity"] = dict(h_min_p500_um=vm, h_min_p1_um=p1, h_min_p0_um=p0,
                                       rel_p1=p1 / vm - 1, rel_p0=p0 / vm - 1,
                                       passed=abs(p1 / vm - 1) < 1e-3 and abs(p0 / vm - 1) < 1e-3)
    common.save_json(args.work / "summary.json", summary)
    keys = ["load_case", "model", "load_N", "eps", "attitude_deg", "h_min_um", "h_min_angle_deg", "h_min_angle_from_hole_deg",
            "torque_Nm", "friction_power_W", "p_max_bar", "p_min_bar", "cavitated_area_pct",
            "cavity_start_from_hole_deg", "cavity_end_from_hole_deg", "film_content_hole", "film_content_ring_mean",
            "film_content_ring_min", "bubble_touches_hole_ring", "bubble_gap_to_hole_edge_mm", "hole_mass_flow_g_s",
            "h_min_um_gc_minus_vc_pct", "torque_Nm_gc_minus_vc_pct", "p_max_bar_gc_minus_vc_pct",
            "hole_mass_flow_g_s_gc_minus_vc_pct", "liquid_residual_rel", "valid", "dissolved_min", "dissolved_max",
            "alpha_gas_max", "run_dir"]
    with (args.work / "gc_effect.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(table)
    print(json.dumps(dict(gc_effect=summary["gc_effect"], pcav=summary["pcav_sensitivity"],
                          reproduction=summary["reproduction"]), indent=1))


if __name__ == "__main__":
    main()
