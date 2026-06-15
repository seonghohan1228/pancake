#!/usr/bin/env bash
# Grid-convergence study (AUDIT §C1 / WP-3 deliverable): runs the steady
# journal case on 60x20, 120x40, and 240x80 and tabulates cavitation extent
# and the conservation residual per grid.
#
# Usage (from the build directory): ../validation/grid_convergence/run.sh [ranks]

set -euo pipefail
ranks="${1:-2}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
base_config="$script_dir/config_base.txt"

summary="grid_convergence_summary.csv"
echo "grid,cavitated_fraction,liquid_residual,outer_iters" > "$summary"

for grid in "60 20" "120 40" "240 80"; do
    set -- $grid
    ntheta=$1; nz=$2; tag="${ntheta}x${nz}"
    out_dir="results/grid_convergence_$tag"
    config="$(mktemp)"
    sed -e "s/^n_theta_global *=.*/n_theta_global = $ntheta/" \
        -e "s/^n_z_global *=.*/n_z_global = $nz/" \
        -e "s#^output_dir *=.*#output_dir = $out_dir#" \
        "$base_config" > "$config"

    echo "Running $tag ..."
    mpirun -n "$ranks" ./pancake -c "$config" > /dev/null

    last="$(tail -n 1 "$out_dir/diagnostics.csv")"
    fraction="$(echo "$last" | cut -d, -f16)"
    residual="$(echo "$last" | cut -d, -f11)"
    outers="$(echo "$last" | cut -d, -f3)"
    echo "$tag,$fraction,$residual,$outers" | tee -a "$summary"
done

echo "Summary written to $summary"
