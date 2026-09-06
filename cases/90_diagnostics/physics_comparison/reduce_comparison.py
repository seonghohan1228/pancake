"""Reduce ASCII VTK snapshots into portable CSVs and checked comparison metrics."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET

import numpy as np

HERE = Path(__file__).resolve().parent
REVISIONS = ('previous', 'current')


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def datasets(run):
    pvd = run / 'results' / 'flat' / 'results.pvd'
    return [(float(item.attrib['timestep']), pvd.parent / item.attrib['file'])
            for item in ET.parse(pvd).findall('.//DataSet')]


def read_frame(pvts):
    """Assemble theta-decomposed physical cells; flat grids have zero third extent."""
    strips = []
    for reference in ET.parse(pvts).findall('.//Piece'):
        piece = ET.parse(pvts.parent / reference.attrib['Source']).find('.//Piece')
        extent = list(map(int, piece.attrib['Extent'].split()))
        nx, nz = extent[1] - extent[0], extent[3] - extent[2]
        arrays = {}
        for data in piece.findall('./CellData/DataArray'):
            components = int(data.attrib.get('NumberOfComponents', '1'))
            shape = (nz, nx) if components == 1 else (nz, nx, components)
            values = np.fromstring(data.text or '', sep=' ').reshape(shape)
            assert np.all(np.isfinite(values)), f'Nonfinite {pvts}: {data.attrib}'
            arrays[data.attrib['Name']] = values
        strips.append((extent[0], arrays))
    strips.sort(key=lambda item: item[0])
    return {name: np.concatenate([strip[1][name] for strip in strips], axis=1)
            for name in strips[0][1]}


def measures(frame):
    result = {'temperature_max_K': float(frame['T'].max())}
    if 'p' in frame:
        result['pressure_max_Pa'] = float(frame['p'].max())
        result['pressure_min_Pa'] = float(frame['p'].min())
    if 'F' in frame:
        force = frame['F'][0, 0]
        result.update(load_N=float(np.linalg.norm(force)),
                      force_x_N=float(force[0]), force_y_N=float(force[1]))
    if 'Fp' in frame:
        force = frame['Fp'][0, 0]
        result.update(pressure_load_N=float(np.linalg.norm(force)),
                      pressure_force_x_N=float(force[0]), pressure_force_y_N=float(force[1]))
    if 'film_content' in frame:
        j = frame['film_content'].shape[0] // 2
        result['liquid_deficit_midplane_arc_deg_theta_below_0p999'] = float(
            np.mean(frame['film_content'][j] < 0.999) * 360)
        result['liquid_deficit_any_axial_arc_deg_theta_below_0p999'] = float(
            np.mean(np.any(frame['film_content'] < 0.999, axis=0)) * 360)
        result['film_content_min'] = float(frame['film_content'].min())
    for field in ('alpha_gas', 'dissolved_gas', 'T_bearing', 'mass_residual'):
        if field in frame:
            result[field + '_max'] = float(frame[field].max())
    if 'xB' in frame:
        result['bearing_z_m'] = float(frame['xB'][0, 0, 2])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=HERE / 'evidence')
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.work / 'run_manifest.json').read_text(encoding='utf-8'))
    assert len(manifest['runs']) == 10, 'Run set incomplete'
    assert all(r['return_code'] == 0 and r.get('all_reynolds_converged', True)
               and r.get('all_coupling_converged') is not False for r in manifest['runs'])
    for case in ('exchange', 'ferron', 'optical', 'air', 'timing'):
        hashes = {r['input_sha256'] for r in manifest['runs'] if r['case'] == case}
        assert len(hashes) == 1, f'Inputs differ: {case}'
    shutil.copy2(args.work / 'run_manifest.json', output / 'run_manifest.json')
    summary = {'revisions': {r: manifest[r + '_commit'] for r in REVISIONS},
               'claims': {'experimental_oil_air_validated': False,
                          'experimental_oil_r410a_validated': False,
                          'gas_transient_load_timestep_accepted': False,
                          'direct_closed_exchange_conservation_verified': True}}
    for case in ('ferron', 'optical', 'air', 'timing'):
        histories = {}
        final_frames = {}
        for revision in REVISIONS:
            run = args.work / case / revision
            snapshots = datasets(run)
            history = []
            for reported_time, path in snapshots:
                # Static air case only: the old driver advances once at label t=0.
                physical_time = reported_time + (1e-5 if case == 'air' and revision == 'previous' else 0)
                frame = read_frame(path)
                history.append(dict(reported_time_s=reported_time,
                                    comparison_time_s=physical_time, **measures(frame)))
                if case != 'air' or physical_time <= 0.001 + 1e-12:
                    final_frames[revision] = frame
            histories[revision] = history
            write_csv(output / f'{case}_{revision}_history.csv', history)
            shutil.copy2(run / 'run.log', output / f'{case}_{revision}.log')
            shutil.copy2(run / 'results' / 'diagnostics.csv',
                         output / f'{case}_{revision}_diagnostics.csv')
            if case in ('ferron', 'optical'):
                j = frame['p'].shape[0] // 2
                # Exact uniform-grid centers avoid VTK's rounded coordinate field.
                angle = (np.arange(frame['p'].shape[1]) + 0.5) * 360 / frame['p'].shape[1]
                fields = ['p', 'T', 'film_content', 'mu']
                if 'T_bearing' in frame:
                    fields.append('T_bearing')
                write_csv(output / f'{case}_{revision}_profile.csv',
                          [dict(angle_deg=float(angle[i]),
                                **{key: float(frame[key][j, i]) for key in fields})
                           for i in range(len(angle))])
        summary[case] = {revision: measures(final_frames[revision]) for revision in REVISIONS}
        summary[case]['snapshot_count'] = {r: len(histories[r]) for r in REVISIONS}
        if case in ('ferron', 'optical', 'air'):
            previous, current = [final_frames[r] for r in REVISIONS]
            summary[case]['max_abs_field_change'] = {
                name: float(np.max(np.abs(current[name] - previous[name])))
                for name in current if name in previous}
        if case == 'air':
            old = {round(r['comparison_time_s'], 10): r for r in histories['previous']
                   if r['comparison_time_s'] <= .001 + 1e-12}
            new = {round(r['comparison_time_s'], 10): r for r in histories['current']
                   if r['comparison_time_s'] > 0}
            assert old.keys() == new.keys() and len(old) == 100
            summary[case]['time_alignment'] = 'previous physical time = label + 10 us; current = label'
            summary[case]['compared_final_time_s'] = .001
            summary[case]['unshifted_previous_at_1ms_label'] = histories['previous'][-1]
            summary[case]['max_abs_aligned_history_change'] = {
                key: max(abs(new[t][key] - old[t][key]) for t in new)
                for key in ('pressure_max_Pa', 'load_N', 'temperature_max_K', 'alpha_gas_max')}
            for revision in REVISIONS:
                rows = list(csv.DictReader((output / f'air_{revision}_diagnostics.csv').open(encoding='utf-8')))
                rows = [r for r in rows if float(r['time']) + (1e-5 if revision == 'previous' else 0) <= .001 + 1e-12]
                summary[case][revision]['diagnostic_max_abs_gas_residual'] = max(abs(float(r['gas_residual'])) for r in rows)
                if revision == 'current':
                    for key in ('gas_volume_excess', 'gas_eos_defect'):
                        summary[case][revision]['diagnostic_max_' + key] = max(float(r[key]) for r in rows)
    summary['exchange'] = {}
    for revision in REVISIONS:
        source = args.work / 'exchange' / revision / 'inventory.csv'
        shutil.copy2(source, output / f'exchange_{revision}.csv')
        data = np.genfromtxt(source, names=True, delimiter=',')
        record = {'final_alpha_gas': float(data['alpha_gas'][-1])}
        for quantity in ('total_gas_area', 'oil_mass_area', 'total_mass_area'):
            error = data[quantity] / data[quantity][0] - 1
            record[quantity + '_final_relative_error'] = float(error[-1])
            record[quantity + '_max_relative_error'] = float(np.max(np.abs(error)))
        if revision == 'current':
            assert record['total_gas_area_max_relative_error'] < 1e-12
            assert record['oil_mass_area_max_relative_error'] < 1e-12
            # Independent frozen p,T analytic law: dMd/dt=-k[(1-ceq)Md-ceq Moil].
            ceq, rate = .001, 20
            oil = data['oil_mass_area'][0]
            equilibrium = ceq * oil / (1 - ceq)
            exact = equilibrium + (data['dissolved_mass_area'][0] - equilibrium) * np.exp(
                -rate * (1 - ceq) * data['time'])
            error = float(np.max(np.abs(data['dissolved_mass_area'] - exact)))
            assert error < 1e-15
            record['analytic_dissolved_max_abs_error_kg_m2'] = error
        summary['exchange'][revision] = record
    timing = summary['timing']['current']
    assert abs(timing['temperature_max_K'] - (300 + .05 * (100 * .02 / .0001)**2 * .015 / 1.9e6)) < .0005
    assert abs(timing['bearing_z_m'] - .015**2) < 1e-12
    refinement = args.work / 'air_refinement'
    if refinement.exists():
        shutil.copy2(refinement / 'manifest.json', output / 'refinement_manifest.json')
        summary['air_refinement'] = {}
        for name in ('dt5us', 'dt2p5us'):
            run = refinement / name
            history = [dict(time_s=t, **measures(read_frame(path))) for t, path in datasets(run)]
            write_csv(output / f'air_{name}_history.csv', history)
            for source, target in (('run.log', '.log'), ('config.txt', '_config.txt'),
                                   ('results/diagnostics.csv', '_diagnostics.csv')):
                shutil.copy2(run / source, output / f'air_{name}{target}')
            rows = list(csv.DictReader((run / 'results/diagnostics.csv').open(encoding='utf-8')))
            assert all(r['converged'] == '1' for r in rows)
            tail = [h for h in history if h['time_s'] >= .0008]
            summary['air_refinement'][name] = {
                'final': history[-1],
                'load_range_last_0p2ms_N': [min(h['load_N'] for h in tail), max(h['load_N'] for h in tail)],
                'max_gas_volume_excess': max(float(r['gas_volume_excess']) for r in rows),
                'max_abs_gas_residual': max(abs(float(r['gas_residual'])) for r in rows)}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
