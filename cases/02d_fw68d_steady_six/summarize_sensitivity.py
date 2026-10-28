"""Recheck completed thermal scenarios and report changes against the saved baseline.

Run on the workstation with the exact executable used for the campaign. Missing,
failed, capped, or stale runs cannot appear as accepted sensitivity predictions.
"""

import argparse
import csv
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import run_steady as driver
from prepare_sensitivity import parse_config



def execution_rankfile(execution, scenario, case):
    """Use the recorded per-case binding without weakening provenance checks."""
    jobs = [job for job in execution["jobs"] if (job["scenario"], job["case"]) == (scenario, case)]
    if len(jobs) != 1 or jobs[0]["status"] != "completed" or jobs[0]["returncode"] != 0:
        raise ValueError(f"{scenario}/{case} has no unique successful execution")
    job = jobs[0]
    rankfile = Path(job["rankfile"])
    if driver.sha256(rankfile) != job["rankfile_sha256"]:
        raise ValueError(f"{scenario}/{case} rankfile changed after execution")
    return rankfile


def field_metrics(config, raw, points, materials):
    """Mirror the gas-free FBNS fluxes and shaft-shear sum on the output grid."""
    for key, expected in {"cavitation_model": "FBNS", "fluid_property_model": "SINGLE_PHASE",
                          "cavitated_film_model": "STRIATED", "solution_mode": "STEADY_STATE",
                          "bush_thermal_model": "CONDUCTION_3D"}.items():
        if config[key] != expected:
            raise ValueError(f"reporter only supports {key}={expected}")
    nt, nz = int(config["n_theta_global"]), int(config["n_z_global"])
    order = np.lexsort((raw["theta_rad"], raw["z_m"]))
    fields = {key: value[order].reshape(nz, nt) for key, value in raw.items() if len(value) == nt * nz}
    radius, length, omega = (float(config[k]) for k in ("R", "L", "omega"))
    dtheta, dz = 2 * np.pi / nt, length / nz
    area = radius * dtheta * dz
    p, theta, rho = (fields[k] for k in ("p", "film_content", "rho_liquid_solution"))
    f0, f1, f2, hc = (fields[k] for k in ("cross_film_f0", "cross_film_f1", "cross_film_f2", "cross_film_couette_h"))
    full = theta >= 1.
    dp = (np.roll(p, -1, axis=1) - np.roll(p, 1, axis=1)) / (2 * dtheta)
    couette = float(np.sum(omega * radius / f0 * np.where(full, 1., np.maximum(theta, 0.))) * radius * area)
    pressure = float(np.sum(full * f1 / f0 / radius * dp) * radius * area)
    mobility = rho * f2
    def harmonic(a, b):
        return 2 * a * b / (a + b)
    east_rho, east_hc = np.roll(rho, -1, axis=1), np.roll(hc, -1, axis=1)
    mc = .5 * (rho + east_rho) * omega * radius * .5 * (hc + east_hc) * .5 * dz
    east = harmonic(mobility, np.roll(mobility, -1, axis=1)) * dz / (radius * dtheta) * (p - np.roll(p, -1, axis=1))
    east += mc * np.where(mc >= 0, theta, np.roll(theta, -1, axis=1))
    zflux = np.zeros((nz + 1, nt))
    zflux[0] = 2 * radius * dtheta / dz * mobility[0] * (float(config["bc_z_south_val"]) - p[0])
    zflux[-1] = 2 * radius * dtheta / dz * mobility[-1] * (p[-1] - float(config["bc_z_north_val"]))
    zflux[1:-1] = radius * dtheta / dz * harmonic(mobility[:-1], mobility[1:]) * (p[:-1] - p[1:])
    center, width, _, _ = map(float, config["inlet_groove"].split())
    angle = np.degrees(fields["theta_rad"][0])
    groove = np.abs((angle - center + 180) % 360 - 180) <= width / 2
    # Groove rows have no external axial link, even where they overlap the ends.
    zflux[0, groove] = 0
    zflux[-1, groove] = 0
    feed = (east - np.roll(east, 1, axis=1) + zflux[1:] - zflux[:-1])[:, groove]
    fluxes = {}
    for name, incoming in (("south", zflux[0]), ("north", -zflux[-1]), ("groove", feed)):
        fluxes[name + "_in_kg_s"] = float(np.maximum(incoming, 0).sum())
        fluxes[name + "_out_kg_s"] = float(np.maximum(-incoming, 0).sum())
    nr, ns = int(config["bush_radial_cells"]), int(config["bush_axial_cells"])
    ro, kb = float(config["bush_outer_radius"]), float(config["bush_thermal_conductivity"])
    dr, solid_dz = (ro - radius) / nr, length / ns
    rr = radius + (np.arange(nr) + .5) * dr
    solid = np.array([fields[f"cross_film_solid_{k}"] for k in range(nr)])
    mid = .5 * (solid[:, nz // 2 - 1] + solid[:, nz // 2])
    curve = np.array([np.interp(np.log(.026), np.log(rr), mid[:, i]) for i in range(nt)])
    sensors = np.interp([240, 270, 300], angle, curve, period=360)
    hb, hs = float(config["bush_outer_heat_transfer"]), float(config["shaft_ambient_heat_transfer"])
    tb, ts = float(config["bush_outer_temperature"]), float(config["shaft_ambient_temperature"])
    q_bush = float(-fields["q_bearing"].sum() * area)
    q_shaft = float(-fields["q_journal"].sum() * area)
    q_outer = float(hb * (fields["temperature_housing_outer"] - tb).sum() * ro * dtheta * dz)
    q_disk = float(hs * float(config["shaft_ambient_area_ratio"]) * 2 * np.pi * radius * length * (fields["T_journal"].mean() - ts))
    # Same local end-face correlation and half-cell resistance as the native solver.
    nodes = solid[:, [0, -1], :]
    surface = nodes.copy()
    heights = 2 * ro
    tab_t = np.array([v["temperature_c"] for v in points])
    tab_nu = np.log([v["nu_mm2_s"] for v in points])
    for _ in range(20):
        film_c = (surface + tb) / 2 - 273.15
        if film_c.min() < tab_t[0] or film_c.max() > tab_t[-1]:
            raise ValueError("end-cooling properties outside viscosity table")
        nu = np.exp(np.interp(film_c, tab_t, tab_nu)) * 1e-6
        density = driver.density(film_c)
        alpha = float(config["thermal_conductivity"]) / (density * float(config["liquid_specific_heat"]))
        pr = nu / alpha
        ra = 9.80665 * .7284 / density * np.abs(surface - tb) * heights**3 / (nu * alpha)
        he = float(config["thermal_conductivity"]) / heights * (.825 + .387 * ra**(1/6) / (1 + (.492/pr)**(9/16))**(8/27))**2
        surface = tb + (nodes - tb) / (1 + he * solid_dz / (2 * kb))
    q_ends = float(np.sum(rr[:, None, None] * dr * dtheta / (solid_dz/(2*kb) + 1/he) * (nodes - tb)))
    return {"couette_Nm": couette, "pressure_shear_Nm": pressure,
        "torque_reconstruction_error_Nm": couette + pressure - float(fields["friction_torque"][0, 0]),
        "TC_C": (sensors - 273.15).tolist(), "film_mean_K": float(fields["T"].mean()),
        "viscosity_mean_Pa_s": float(fields["mu"].mean()), "viscosity_max_Pa_s": float(fields["mu"].max()),
        "midplane_cavity_arc_deg": float(np.count_nonzero(theta[nz//2] < .999) * 360/nt),
        **fluxes, "bush_heat_in_W": q_bush, "shaft_heat_in_W": q_shaft,
        "outer_cylinder_heat_out_W": q_outer, "annular_ends_heat_out_W": q_ends,
        "shaft_disk_heat_out_W": q_disk,
        "bush_heat_closure_W": q_bush - q_outer - q_ends,
        "shaft_heat_closure_W": q_shaft - q_disk}


def plot_results(rows, baseline, destination):
    import matplotlib
    matplotlib.use("Agg")
    sys.path.insert(0, str(driver.HERE.parents[1] / "validation/journal_bearing"))
    from plot_style import apply_style, save_figure, COLORS
    import matplotlib.pyplot as plt
    apply_style()
    # Match the current lab preset when rerunning with the archived style module.
    matplotlib.rcParams.update({"lines.markeredgecolor":"black", "lines.markeredgewidth":.6,
        "legend.handletextpad":.4, "legend.handlelength":1., "legend.labelspacing":.3,
        "legend.borderaxespad":.4, "legend.borderpad":.3, "legend.numpoints":1})
    names = {"supply_plus_1K":"Supply +1 K", "cooling_plus_1K":"Baths +1 K", "supply_TC2":"TC2 supply"}
    labels = [r["case_id"] + "\n" + names[r["scenario"]] for r in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), constrained_layout=True)
    axes[0].plot(x-.12, [baseline[r["case_id"]]["torque_Nm"] for r in rows], "s", label="Baseline")
    axes[0].plot(x, [r["torque_Nm"] for r in rows], "o", label="Perturbed")
    axes[0].plot(x+.12, [r["new_exp_torque_nm"] for r in rows], "^", label="Experiment")
    axes[0].set_ylabel(r"$M$ (N m)")
    axes[0].legend(loc="best", fontsize=12)
    for i, (old_marker, new_marker) in enumerate(zip(("s", "o", "^"), ("D", "v", ">"))):
        offset = .18 * (i-1)
        axes[1].plot(x+offset-.045, [baseline[r["case_id"]][f"TC{i+1}_error_K"] for r in rows],
                     old_marker, color=COLORS[i], label=f"TC{i+1}, baseline")
        axes[1].plot(x+offset+.045, [r["TC_error_K"][i] for r in rows], new_marker,
                     color=COLORS[i], label=f"TC{i+1}, perturbed")
    axes[1].set_ylabel(r"$T_{\mathrm{CFD}}-T_{\mathrm{exp}}$ (K)")
    axes[1].legend(loc="best", fontsize=11, ncol=2)
    for ax in axes:
        ax.set_xticks(x, labels, fontsize=12)
        ax.set_xlim(-.5, len(rows)-.5)
    save_figure(fig, str(destination / "baseline_vs_perturbed"))
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    axes = axes.ravel()
    axes[0].plot(x, [1e3*r["torque_sensitivity_Nm_K"] for r in rows], "s")
    axes[0].set_ylabel(r"$\Delta M/\Delta T$ (mN m/K)")
    axes[1].plot(x, [1e3*r["torque_absolute_error_change_Nm"] for r in rows], "s")
    axes[1].set_ylabel(r"$\Delta |M_{\mathrm{CFD}}-M_{\mathrm{exp}}|$ (mN m)")
    for i, marker in enumerate(("s", "o", "^")):
        offset = .14 * (i-1)
        axes[2].plot(x+offset, [r["TC_sensitivity_K_K"][i] for r in rows], marker, label=f"TC{i+1}")
        axes[3].plot(x+offset, [abs(r["TC_error_K"][i])-abs(baseline[r["case_id"]][f"TC{i+1}_error_K"]) for r in rows], marker, label=f"TC{i+1}")
    axes[2].set_ylabel(r"$\Delta T_{\mathrm{sensor}}/\Delta T$")
    axes[3].set_ylabel(r"$\Delta |T_{\mathrm{CFD}}-T_{\mathrm{exp}}|$ (K)")
    axes[2].legend(loc="best", fontsize=12)
    axes[3].legend(loc="best", fontsize=12)
    for ax in axes:
        ax.set_xticks(x, labels, rotation=25, ha="right", fontsize=11)
        ax.set_xlim(-.5, len(rows)-.5)
    save_figure(fig, str(destination / "sensitivities_and_error_changes"))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--execution-manifest", type=Path, help="parallel/status.json with successful per-case bindings")
    args = parser.parse_args()
    execution = json.loads(args.execution_manifest.read_text()) if args.execution_manifest else None
    if execution and execution["executable_sha256"] != driver.sha256(args.solver):
        raise ValueError("parallel executable differs from the reporting executable")
    manifest = json.loads((args.prepared / "manifest.json").read_text())
    for name, key in (("experimental_comparison.csv", "experiment_csv_sha256"),
                      ("baseline_comparison.json", "baseline_comparison_sha256")):
        if driver.sha256(args.prepared / name) != manifest[key]:
            raise ValueError(f"comparison evidence changed: {name}")
    baseline = {r["case"]: r for r in json.loads((args.prepared / "baseline_comparison.json").read_text())}
    with (args.prepared / "experimental_comparison.csv").open(newline="", encoding="utf-8-sig") as stream:
        experiment = {r["case_id"]: r for r in csv.DictReader(stream)}
    rows = []
    for scenario in manifest["scenarios"]:
        name = scenario["name"]
        setup = args.prepared / (name + ".json")
        if driver.sha256(setup) != scenario["setup_sha256"]:
            raise ValueError("scenario inputs changed since preparation")
        for key in scenario["cases"]:
            rankfile = execution_rankfile(execution, name, key) if execution else args.prepared / "four_cores.rankfile"
            campaign = driver.Campaign(SimpleNamespace(materials=setup, solver=args.solver,
                work=args.work_root / name, trial=None, rankfile=rankfile,
                ranks=4, coupling_max_iters=2000, coupling_relaxation=1.))
            candidates = sorted((campaign.work / key).glob("result_*.json"), key=lambda p: p.stat().st_mtime)
            if not candidates:
                raise ValueError(f"{name}/{key} has no accepted equilibrium; no sensitivity plot produced")
            record = json.loads(candidates[-1].read_text())
            run = Path(record["run_dir"])
            case = campaign.cases[key]
            checked = campaign.recheck(run, record, case)
            if checked is None or record["status"] != "load_and_thermal_closed":
                raise ValueError(f"{name}/{key} failed independent recheck/provenance")
            config = parse_config((run / "config.txt").read_text())
            expected = campaign.config_for(case, record["seed_eps"], record["seed_psi_deg"],
                record["h_bush_W_m2_K"], record["h_shaft_W_m2_K"], run / "results")
            ignore = {"output_dir", "steady_initial_state_file"}
            if driver.config_identity(expected, ignore) != driver.config_identity((run / "config.txt").read_text(), ignore):
                raise ValueError("accepted record does not match the scenario configuration")
            temp = driver.boundary_temperatures(case)
            hb = driver.housing_coefficient(checked["T_outer_mean_K"], temp["bush_outer_temperature"], campaign.points, campaign.materials)[0]
            hs = driver.shaft_coefficient(checked["T_shaft_mean_K"], temp["shaft_ambient_temperature"], case["angular_speed_magnitude_rad_s"], campaign.points, campaign.materials)[0]
            if max(abs(hb/record["h_bush_W_m2_K"]-1), abs(hs/record["h_shaft_W_m2_K"]-1)) >= 1e-3 or record["T_change_K"] >= .01:
                raise ValueError("exterior coefficient closure failed")
            extra = field_metrics(config, driver.read_fields(run), campaign.points, campaign.materials)
            if abs(extra["torque_reconstruction_error_Nm"]) > 1e-8:
                raise ValueError("torque components do not reproduce native torque")
            diag = checked["diagnostics"]
            axial = sum(extra[f"{end}_in_kg_s"]-extra[f"{end}_out_kg_s"] for end in ("north", "south"))
            groove = extra["groove_in_kg_s"]-extra["groove_out_kg_s"]
            scale = float(diag["liquid_balance_scale_kg_s"])
            if max(abs(axial-float(diag["liquid_boundary_flux"])), abs(groove-float(diag["inlet_mass_source"]))) > 1e-8 * scale:
                raise ValueError("split mass fluxes do not reproduce native diagnostics")
            exp = experiment[key]
            torque_exp = float(exp["new_exp_torque_nm"])
            delta_t = max(temp.values()) - case["tc0_boundary_temperature_k"]
            rows.append({**checked, **extra, "scenario": name, "case_id": key,
                "new_exp_torque_nm": torque_exp,
                "torque_error_percent": 100 * (checked["torque_Nm"]/torque_exp-1),
                "TC_error_K": [extra["TC_C"][i-1]-float(exp[f"tc{i}_exp_c"]) for i in (1,2,3)],
                "TC2_independent_validation": name != "supply_TC2",
                "torque_sensitivity_Nm_K": (checked["torque_Nm"]-baseline[key]["torque_Nm"])/delta_t,
                "TC_sensitivity_K_K": [(extra["TC_C"][i-1]-float(exp[f"tc{i}_cfd_c"]))/delta_t for i in (1,2,3)],
                "torque_absolute_error_change_Nm": abs(checked["torque_Nm"]-torque_exp)-abs(baseline[key]["torque_Nm"]-torque_exp)})
    args.out.mkdir(parents=True, exist_ok=False)
    driver.save_json(args.out / "results.json", rows)
    plot_results(rows, baseline, args.out)
    (args.out / "CAPTIONS.txt").write_text(
        "Baseline and perturbed predictions with the selected RPM-specific experimental subtraction. "
        "All perturbation amplitudes are hypotheses, not measured uncertainty. No uncertainty bars "
        "are available in the selected comparison table. TC2-supply cases use TC2 as an input: their "
        "TC2 residual is not independent validation. Per-kelvin TC2 responses are finite-interval "
        "secants. Negative absolute-error changes mean improvement of that output.\n", encoding="utf-8")
    print(f"Rechecked and reported {len(rows)} equilibria in {args.out}")


if __name__ == "__main__":
    main()
