"""Verify and summarize completed A-F free-expansion runs without fitting inputs."""

import argparse
import csv
import json
from pathlib import Path
import re
import sys
from types import SimpleNamespace

import numpy as np
import run_steady as driver
from thermal_expansion import solid_observables, parse_config, geometry_change, expanded_geometry, starting_geometry


def verify_solve(root, case, solve):
    """Recheck raw native fields, provenance, stricter gates and cooling closure."""
    result = solve["result"]
    run = Path(result["run_dir"])
    setup_path = run.parents[2] / "setup.json"
    command = result["command"]
    rankfile = command[command.index("--rankfile") + 1]
    args = SimpleNamespace(solver=str(root / "build-linux/pancake"), work=str(run.parents[1]),
        materials=str(setup_path), ranks=4, rankfile=rankfile, trial=None,
        timeout=7200, max_thermal_iterations=30, coupling_max_iters=int(solve["config"]["coupling_max_iters"]),
        coupling_relaxation=1., initial_state_root=None)
    campaign = driver.Campaign(args)
    checked = campaign.recheck(run, result, campaign.cases[case])
    if checked is None or result["status"] != "load_and_thermal_closed":
        raise ValueError(f"unverified equilibrium {case}/{solve['label']}")
    config = parse_config((run / "config.txt").read_text())
    if config != solve["config"] or driver.sha256(run / "results/steady_state.dat") != result["output_state_sha256"]:
        raise ValueError("saved configuration or checkpoint changed")
    raw = driver.read_fields(run)
    baseline = json.loads((root / "baseline" / f"{case}_config.json").read_text())
    allowed = {"R", "c", "L", "bush_outer_radius", "shaft_ambient_area_ratio", "e", "attitude_angle_deg",
               "bush_outer_heat_transfer", "shaft_ambient_heat_transfer", "output_dir", "solution_table_file",
               "steady_initial_state_file", "steady_force_tolerance", "coupling_tolerance", "fbns_tol"}
    differences = {key for key in set(config) | set(baseline) if config.get(key) != baseline.get(key)}
    if differences - allowed:
        raise ValueError(f"unrelated baseline changes: {differences - allowed}")
    if not np.all(raw["friction_torque"] == result["torque_Nm"]):
        raise ValueError("integrated torque differs across output cells")
    radius = solve["observables"]["sensor_radius_m"]
    actual = solid_observables(config, raw, radius)
    for key in ("Ts_C", "Tb_C", "TC_C"):
        if not np.allclose(actual[key], solve["observables"][key], rtol=0, atol=1e-11):
            raise ValueError("solid postprocessing changed")
    proposed = expanded_geometry(actual["Ts_C"], actual["Tb_C"])
    if proposed != solve["proposed_geometry"]:
        raise ValueError("geometry feedback was not computed from reference dimensions")
    if solve["geometry"] and geometry_change(solve["geometry"], proposed) != solve["feedback"]:
        raise ValueError("geometry feedback stopping check changed")
    material = driver.case_materials(campaign.cases[case], campaign.materials)
    ambient = campaign.cases[case]["tc0_boundary_temperature_k"]
    hb, _ = driver.housing_coefficient(result["T_outer_mean_K"], ambient, campaign.points, material)
    hs, _ = driver.shaft_coefficient(result["T_shaft_mean_K"], ambient, float(config["omega"]), campaign.points, material)
    mismatch = max(abs(hb/result["h_bush_W_m2_K"]-1), abs(hs/result["h_shaft_W_m2_K"]-1))
    if mismatch >= campaign.h_tolerance or result["T_change_K"] >= campaign.wall_tolerance:
        raise ValueError("exterior cooling is not closed")
    log = (run / "run.log").read_text()
    fbns = re.findall(r"fbns iter=\d+ \|G\|=" + driver.NUMBER + r" \|FB\|=" + driver.NUMBER, log)
    fbns_final = max(map(float, fbns[-1]))
    thermal = float(config["coupling_tolerance"])
    if (fbns_final >= float(config["fbns_tol"]) or
        checked["force_error_rel"] >= float(config["steady_force_tolerance"]) or
        float(checked["diagnostics"]["coupling_R_T"]) >= thermal or
        checked["cross_film_profile_residual"] >= thermal or checked["energy_balance_residual_rel"] >= thermal):
        raise ValueError("requested numerical tolerance failed")
    return {"fbns_final": fbns_final, "rechecked_cooling_mismatch": mismatch,
            "fields_sha256": driver.sha256(run.parents[2] / "fields.npz"),
            "config_sha256": driver.sha256(run / "config.txt"),
            "checkpoint_sha256": driver.sha256(run / "results/steady_state.dat")}


def compare(case, experiment, baseline, history, root):
    if history["status"] not in ("complete", "feedback_unresolved_at_cap"):
        raise ValueError(f"case {case} is not completed: {history['status']}")
    solves = history["solves"]
    checks = {s["label"]: verify_solve(root, case, s) for s in solves}
    nominal = solves[0]
    normal = [s for s in solves if s["geometry"] and not s["tight"]]
    tight = [s for s in solves if s["geometry"] and s["tight"]]
    first, settled, final = normal[0], normal[-1], tight[-1]
    starts = json.loads((root / "inputs/LGE-2026-10-17-thermal-expansion-cases.json").read_text())["cases"]
    if first["geometry"] != starting_geometry(next(row for row in starts if row["case"] == case)):
        raise ValueError("first-pass dimensions differ from the supplied full-precision inputs")
    m_exp = float(experiment["new_exp_torque_nm"])
    tc_exp = np.array([float(experiment[f"tc{i}_exp_c"]) for i in (1, 2, 3)])
    m0 = baseline["torque_Nm"]
    tc0 = np.array([baseline[f"TC{i}_CFD_C"] for i in (1, 2, 3)])
    variants = {"archived_baseline": {"torque_Nm": m0, "TC_C": tc0.tolist()}}
    for name, solve in (("first_pass", first), ("normal_settled", settled),
                         ("tight_nominal", nominal), ("final_expanded", final)):
        variants[name] = {"torque_Nm": solve["result"]["torque_Nm"],
                          "TC_C": solve["observables"]["TC_C"]}
    for value in variants.values():
        value["torque_error_pct"] = 100*(value["torque_Nm"]/m_exp-1)
        value["TC_error_K"] = (np.array(value["TC_C"])-tc_exp).tolist()
    delta_normal = variants["normal_settled"]["torque_Nm"]-m0
    delta_tight = variants["final_expanded"]["torque_Nm"]-variants["tight_nominal"]["torque_Nm"]
    discrepancy = abs(delta_normal-delta_tight)
    dt_normal = np.array(variants["normal_settled"]["TC_C"])-tc0
    dt_tight = np.array(variants["final_expanded"]["TC_C"])-np.array(variants["tight_nominal"]["TC_C"])
    tc_discrepancy = np.abs(dt_normal-dt_tight)
    before = variants["archived_baseline"]
    after = variants["final_expanded"]
    residual = m0-m_exp
    result = final["result"]
    geometry = final["geometry"]
    return {"case": case, "experimental_torque_Nm": m_exp, "experimental_TC_C": tc_exp.tolist(),
        "S_reference_TC0": float(experiment["S_TC0"]), "S_reference_TC2": float(experiment["S_TC2"]),
        "variants": variants, "final_observables": final["observables"], "final_geometry": geometry,
        "dimensions": {"Ds_m": 2*geometry["journal_radius_m"],
            "Db_m": 2*(geometry["journal_radius_m"]+geometry["radial_clearance_m"]),
            "Do_m": 2*geometry["bush_outer_radius_m"], "L_m": geometry["bearing_length_m"],
            "c_m": geometry["radial_clearance_m"],
            "shaft_end_area_ratio": geometry["journal_radius_m"]/geometry["bearing_length_m"],
            "shaft_exposed_end_area_m2": 2*np.pi*geometry["journal_radius_m"]**2,
            "physical_bore_radius_m": geometry["journal_radius_m"]+geometry["radial_clearance_m"],
            "physical_wall_thickness_m": geometry["bush_outer_radius_m"]-geometry["journal_radius_m"]-geometry["radial_clearance_m"],
            "groove_reference_arc_scaled_m": .005*geometry["bearing_length_m"]/.04,
            "groove_angle_deg": float(final["config"]["inlet_groove"].split()[1]),
            "shaft_surface_speed_m_s": float(final["config"]["omega"])*geometry["journal_radius_m"]},
        "torque_absolute_error_improved": abs(after["torque_error_pct"]) < abs(before["torque_error_pct"]),
        "torque_original_absolute_residual_removed_pct": 100*(1-abs(after["torque_Nm"]-m_exp)/abs(residual)),
        "torque_gap_removed_by_paired_geometry_pct": -100*delta_tight/residual,
        "TC_absolute_error_improved": (np.abs(after["TC_error_K"]) < np.abs(before["TC_error_K"])).tolist(),
        "paired_torque_absolute_error_improved": abs(after["torque_error_pct"]) < abs(variants["tight_nominal"]["torque_error_pct"]),
        "paired_TC_absolute_error_improved": (np.abs(after["TC_error_K"]) < np.abs(variants["tight_nominal"]["TC_error_K"])).tolist(),
        "numerical_resolution": {"normal_delta_torque_Nm": delta_normal, "tight_delta_torque_Nm": delta_tight,
            "paired_delta_torque_pct": 100*delta_tight/variants["tight_nominal"]["torque_Nm"],
            "paired_delta_difference_Nm": discrepancy,
            "torque_sign_resolved": bool(delta_normal*delta_tight > 0 and abs(delta_tight)>discrepancy),
            "nominal_torque_repeat_difference_Nm": variants["tight_nominal"]["torque_Nm"]-m0,
            "paired_TC_change_K": dt_tight.tolist(), "paired_TC_change_difference_K": tc_discrepancy.tolist(),
            "TC_sign_resolved": ((dt_normal*dt_tight>0) & (np.abs(dt_tight)>tc_discrepancy)).tolist(),
            "note": "Tolerance sensitivity is an empirical stopping-error check, not a statistical confidence interval or mesh-independence claim."},
        "final_diagnostics": {"force_relative": result["force_error_rel"],
            "mass_relative": float(result["diagnostics"]["liquid_residual_rel"]),
            "energy_relative": result["energy_balance_residual_rel"],
            "temperature_relative": float(result["diagnostics"]["coupling_R_T"]),
            "layer_relative": result["cross_film_profile_residual"],
            "h_min_m": result["h_min_m"], "eccentricity_ratio": result["eps"],
            "attitude_deg": result["psi_deg"], "p_max_Pa": result["p_max_Pa"],
            "cavity_area_fraction": result["cavitated_area_fraction"],
            "cooling_temperature_change_K": result["T_change_K"], **checks[final["label"]]},
        "feedback_status": history["status"], "geometry_updates": history["geometry_updates"],
        "final_proposed_geometry": final["proposed_geometry"], "final_feedback": final["feedback"],
        "history": [{"label": s["label"], "tight": s["tight"], "geometry": s["geometry"],
            "torque_Nm": s["result"]["torque_Nm"], "observables": s["observables"],
            "feedback": s.get("feedback"), "elapsed_seconds": s["elapsed_seconds"],
            "run_dir": s["result"]["run_dir"]} for s in solves], "verification": checks}


def plot(rows, out, style_dir):
    sys.path.insert(0, str(style_dir))
    from plot_style import apply_style, save_figure, MARKERS, blues, smart_legend
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.5), constrained_layout=True)
    x = np.arange(6)
    colors = blues(3)[::-1]
    for n, (name, label) in enumerate((("archived_baseline", "Baseline"), ("first_pass", "First expansion"),
                                       ("final_expanded", "Final expansion"))):
        axes[0, 0].plot(x+(n-1)*.13, [r["variants"][name]["torque_error_pct"] for r in rows],
                        ls="none", marker=MARKERS[n], color=colors[n], label=label)
    axes[0, 0].set_ylabel(r"$100(M/M_{\mathrm{exp}}-1)$ (%)")
    smart_legend(axes[0, 0], fontsize=12)
    for n, label in enumerate(("First vs baseline", "Settled vs baseline", "Tight matched pair")):
        values=[]
        for r in rows:
            v=r["variants"]
            if n<2:
                values.append(100*(v[("first_pass", "normal_settled")[n]]["torque_Nm"]/v["archived_baseline"]["torque_Nm"]-1))
            else:
                values.append(r["numerical_resolution"]["paired_delta_torque_pct"])
        axes[1, 0].plot(x+(n-1)*.13, values, ls="none", marker=MARKERS[n], color=colors[n], label=label)
    axes[1, 0].set_ylabel(r"$\Delta M/M_{\mathrm{nominal}}$ (%)")
    smart_legend(axes[1, 0], fontsize=12)
    sensor_colors=["#0072B2", "#D55E00", "#009E73"]
    for i in range(3):
        for n, name in enumerate(("archived_baseline", "final_expanded")):
            axes[0, 1].plot(x+(i-1)*.23+(n-.5)*.075,
                [r["variants"][name]["TC_error_K"][i] for r in rows], ls="none",
                marker=MARKERS[i], color=sensor_colors[i], ms=8 if n else 5,
                label=f"TC{i+1}, {'final' if n else 'baseline'}")
        axes[1, 1].plot(x+(i-1)*.16,
            [r["numerical_resolution"]["paired_TC_change_K"][i]*1000 for r in rows],
            ls="none", marker=MARKERS[i], color=sensor_colors[i], label=f"TC{i+1}")
    axes[0, 1].set_ylabel(r"$T_{\mathrm{CFD}}-T_{\mathrm{exp}}$ (K)")
    axes[1, 1].set_ylabel(r"$\Delta T$, tight matched pair (mK)")
    smart_legend(axes[0, 1], fontsize=10, ncol=2)
    smart_legend(axes[1, 1], fontsize=12)
    for ax in axes.flat:
        ax.set_xticks(x, list("ABCDEF"))
        ax.set_xlim(-.5, 5.5)
        ax.set_xlabel("Case")
    save_figure(fig, str(out / "comparison"))
    plt.close(fig)


def report(rows, out):
    changes=[r["numerical_resolution"]["paired_delta_torque_pct"] for r in rows]
    torque_improved=", ".join(r["case"] for r in rows if r["paired_torque_absolute_error_improved"])
    temperature_improved=sum(sum(r["paired_TC_absolute_error_improved"]) for r in rows)
    removal=max(0., max(r["torque_gap_removed_by_paired_geometry_pct"] for r in rows))
    lines=["# A-F homogeneous thermal-expansion comparison", "",
        f"Tight matched-pair torque changes range from {min(changes):+.5f}% to {max(changes):+.5f}%. Torque agreement improves in {torque_improved or 'none of the cases'}; the largest removal of an original torque residual by expansion is {removal:.4f}%. Expansion reduces {temperature_improved} of the {3*len(rows)} absolute sensor errors. Numerical resolution is assessed separately below.", "",
        "Torque uses the corrected RPM-specific `new_exp_torque_nm`. Temperature errors are CFD minus unchanged experimental TC1-3. The final expansion uses tighter closure; the separate tighter nominal repeat isolates the geometry response from stopping error.", "",
        "| Case | Baseline torque (N m) | First pass (N m) | Final (N m) | Baseline error (%) | Final error (%) | Original absolute torque residual removed (%) |",
        "|---|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        v=r["variants"]
        lines.append(f"| {r['case']} | {v['archived_baseline']['torque_Nm']:.9f} | {v['first_pass']['torque_Nm']:.9f} | {v['final_expanded']['torque_Nm']:.9f} | {v['archived_baseline']['torque_error_pct']:+.5f} | {v['final_expanded']['torque_error_pct']:+.5f} | {r['torque_original_absolute_residual_removed_pct']:+.4f} |")
    lines += ["", "Negative removal means a larger discrepancy. These archived-to-final changes include the tighter stopping correction; the matched-pair table below isolates expansion.", "",
        "| Case | TC1 error, baseline -> final (K) | TC2 error, baseline -> final (K) | TC3 error, baseline -> final (K) | Torque / TC1 / TC2 / TC3 improved |", "|---|---:|---:|---:|---|"]
    for r in rows:
        a,b=(r["variants"][k]["TC_error_K"] for k in ("archived_baseline","final_expanded"))
        cells=" | ".join(f"{x:+.6f} -> {y:+.6f}" for x,y in zip(a,b))
        outcomes=" / ".join("yes" if x else "no" for x in [r["torque_absolute_error_improved"],*r["TC_absolute_error_improved"]])
        lines.append(f"| {r['case']} | {cells} | {outcomes} |")
    lines += ["", "| Case | Final Ts / Tb (C) | Ds / Db (mm) | c (micrometres) | L / Do (mm) | Updates |", "|---|---:|---:|---:|---:|---:|"]
    for r in rows:
        o,d=r["final_observables"],r["dimensions"]
        lines.append(f"| {r['case']} | {o['Ts_C']:.7f} / {o['Tb_C']:.7f} | {d['Ds_m']*1e3:.8f} / {d['Db_m']*1e3:.8f} | {d['c_m']*1e6:.7f} | {d['L_m']*1e3:.8f} / {d['Do_m']*1e3:.8f} | {r['geometry_updates']} |")
    lines += ["", "Reported dimensions are those actually solved. The proposed next dimensions and residuals are retained in results.json. Geometry closure requires delta c <=0.001 micrometre and delta Ds/Db/L/Do <=0.01 micrometre. Displayed precision is for reproducibility, not dimensional measurement accuracy.", "",
        "| Case | Matched expansion torque change (%) | Nominal repeat shift (mN m) | Normal/tight response difference (mN m) | Torque sign resolved | Matched TC1 / TC2 / TC3 change (mK) |", "|---|---:|---:|---:|---|---|"]
    for r in rows:
        n=r["numerical_resolution"]
        tc=" / ".join(f"{v*1e3:+.4f}" for v in n["paired_TC_change_K"])
        lines.append(f"| {r['case']} | {n['paired_delta_torque_pct']:+.6f} | {n['nominal_torque_repeat_difference_Nm']*1e3:+.6f} | {n['paired_delta_difference_Nm']*1e3:.6f} | {n['torque_sign_resolved']} | {tc} |")
    lines += ["", "The response is called sign-resolved only when normal and tighter paired changes have the same sign and their difference is smaller than the tighter response. This is a convergence sensitivity check, not an uncertainty interval. Individual sensor sign-resolution flags are in results.json.", "",
        "| Case | Original torque gap removed by matched expansion (%) | Matched torque / TC1 / TC2 / TC3 improved | TC1 / TC2 / TC3 sign resolved | TC response difference, normal vs tight (mK) |",
        "|---|---:|---|---|---|"]
    for r in rows:
        n=r["numerical_resolution"]
        improved=" / ".join("yes" if value else "no" for value in [r["paired_torque_absolute_error_improved"], *r["paired_TC_absolute_error_improved"]])
        resolved=" / ".join("yes" if value else "no" for value in n["TC_sign_resolved"])
        differences=" / ".join(f"{value*1e3:.4f}" for value in n["paired_TC_change_difference_K"])
        lines.append(f"| {r['case']} | {r['torque_gap_removed_by_paired_geometry_pct']:+.5f} | {improved} | {resolved} | {differences} |")
    lines += ["", "A nominal improvement whose sign is unresolved is not evidence of a physical improvement. Responses comparable to the normal/tight difference have unresolved magnitude at that scale, even when their signs agree.", "",
        "| Case | Force | Mass | Energy | Temperature / layer | FBNS | Cooling h / temperature (K) |", "|---|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        d=r["final_diagnostics"]
        lines.append(f"| {r['case']} | {d['force_relative']:.3e} | {d['mass_relative']:.3e} | {d['energy_relative']:.3e} | {d['temperature_relative']:.3e} / {d['layer_relative']:.3e} | {d['fbns_final']:.3e} | {d['rechecked_cooling_mismatch']:.3e} / {d['cooling_temperature_change_K']:.3e} |")
    lines += ["", "| Case | Minimum gap (micrometres) | Eccentricity ratio | Attitude (degrees) | Maximum pressure (MPa) | Cavity area fraction | Feedback status |", "|---|---:|---:|---:|---:|---:|---|"]
    for r in rows:
        d=r["final_diagnostics"]
        lines.append(f"| {r['case']} | {d['h_min_m']*1e6:.6f} | {d['eccentricity_ratio']:.8f} | {d['attitude_deg']:.6f} | {d['p_max_Pa']/1e6:.6f} | {d['cavity_area_fraction']:.8f} | {r['feedback_status']} |")
    lines += ["", "## Geometry history", "", "| Case | Solve | Clearance (micrometres) | Ts / Tb (C) | Torque (N m) | Proposed clearance update (micrometres) |", "|---|---|---:|---:|---:|---:|"]
    for r in rows:
        for s in r["history"]:
            c=s["geometry"]["radial_clearance_m"]*1e6 if s["geometry"] else 42.5
            delta=f"{s['feedback']['delta_m']['c']*1e6:+.7f}" if s["feedback"] else "not applied"
            o=s["observables"]
            lines.append(f"| {r['case']} | {s['label']} | {c:.7f} | {o['Ts_C']:.7f} / {o['Tb_C']:.7f} | {s['torque_Nm']:.10f} | {delta} |")
    lines += ["", "## Scope and evidence", "",
        "The model uses alpha=16e-6 K^-1 for both 304 analogues, with all dimensions referenced to 20 C. The shaft temperature and volume-mean monolithic-solid temperature set free expansion. It excludes nonuniform thermoelastic distortion, restraint, elastic load deformation, contact and different stainless grades. No measured input or torque subtraction was tuned.", "",
        "The original 720 x 192 film plane, 16 film layers and 16 x 16 solid radial/axial cells are unchanged. This campaign tests closure precision on that grid. It does not establish new mesh independence or experimental uncertainty. Reference Sommerfeld coordinates are retained in results.json; no hot-geometry friction conversion is used.", "",
        "[Full results and history](results.json), [comparison figure](comparison.png), and the adjacent raw production, input, source and verification archives preserve the evidence. The first MPI launch failed on CRLF rankfiles; corrected LF rankfiles used new output directories. Original Windows runtime-path and missing remote test-fixture failures are retained along with successful retries.", "",
        "Figure: top left shows torque errors; bottom left separates first-pass, settled and tighter matched torque responses. Top right compares baseline and final sensor errors; bottom right shows the tighter matched sensor response in mK. Small categorical offsets separate markers. No experimental uncertainty was supplied."]
    (out / "RESULTS.md").write_text("\n".join(lines)+"\n")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,required=True)
    parser.add_argument("--production",default="production_attempt_02")
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--style-dir",type=Path,required=True)
    args=parser.parse_args()
    root=args.root.resolve()
    with (root/"inputs/LGE-2026-10-02-stribeck-rpm-drag-cases.csv").open(newline="") as f:
        experiments={r["case_id"]:r for r in csv.DictReader(f)}
    baseline={r["case"]:r for r in json.loads((root/"baseline/case_results.json").read_text())}
    rows=[]
    for case in "ABCDEF":
        history=json.loads((root/args.production/case/"history.json").read_text())
        rows.append(compare(case,experiments[case],baseline[case],history,root))
    args.out.mkdir(parents=True,exist_ok=False)
    driver.save_json(args.out/"results.json",{"cases":rows,"reporter_sha256":driver.sha256(__file__)})
    report(rows,args.out)
    plot(rows,args.out,args.style_dir)
    print(args.out/"RESULTS.md")


if __name__ == "__main__":
    main()
