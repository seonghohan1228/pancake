# Grid-convergence study (AUDIT §C1 / WP-3 deliverable): runs the steady
# journal case on 60x20, 120x40, and 240x80 and tabulates cavitation extent,
# fluid force, friction torque, and the conservation residual per grid.
#
# Usage (from a shell with pancake.exe and mpiexec on PATH, e.g. the build
# directory): powershell -File validation\grid_convergence\run.ps1 [-Ranks 2]

param(
    [int]$Ranks = 2,
    [string]$Pancake = "pancake.exe"
)

$ErrorActionPreference = "Stop"
$script_dir = Split-Path -Parent $MyInvocation.MyCommand.Path
$base_config = Join-Path $script_dir "config_base.txt"
$grids = @(@(60, 20), @(120, 40), @(240, 80))
$rows = @()

foreach ($grid in $grids) {
    $ntheta = $grid[0]
    $nz = $grid[1]
    $tag = "${ntheta}x${nz}"
    $out_dir = "results/grid_convergence_$tag"
    $config = Join-Path $env:TEMP "grid_convergence_$tag.txt"

    (Get-Content $base_config) |
        ForEach-Object {
            $_ -replace '^n_theta_global\s*=.*', "n_theta_global = $ntheta" `
               -replace '^n_z_global\s*=.*', "n_z_global = $nz" `
               -replace '^output_dir\s*=.*', "output_dir = $out_dir"
        } | Set-Content -Encoding ascii $config

    Write-Host "Running $tag ..."
    $stdout = & mpiexec -n $Ranks $Pancake -c $config 2>&1 | Out-String

    $force = "n/a"; $torque = "n/a"
    if ($stdout -match 'F=\(([^)]*)\) N \| M=([-0-9.eE+]+) Nm') {
        $force = $Matches[1]; $torque = $Matches[2]
    }
    $csv = Join-Path $out_dir "diagnostics.csv"
    $last = Import-Csv $csv | Select-Object -Last 1
    $rows += [pscustomobject]@{
        grid = $tag
        cavitated_fraction = $last.cavitated_fraction
        liquid_residual = $last.liquid_residual
        outer_iters = $last.outer_iters
        fluid_force_N = $force
        torque_Nm = $torque
    }
}

$rows | Format-Table -AutoSize
$rows | Export-Csv -NoTypeInformation -Encoding ascii "grid_convergence_summary.csv"
Write-Host "Summary written to grid_convergence_summary.csv"
