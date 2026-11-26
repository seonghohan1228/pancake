# Collects the files an engineer needs to run pancake and its GUI into one folder.
#
#   powershell -File tools/make_release.ps1 [-Build build-windows-mingw] [-Out release]
#
# The folder runs on a Windows machine without MSYS2, CMake or PETSc:
#   pancake_gui.exe   GUI (static; Windows DLLs only)
#   pancake.exe       solver
#   *.dll             MinGW, PETSc, BLAS, METIS runtime deployed by the build
#   fluids/           oil-refrigerant pairs (the synthetic example only)
#   examples/         example cases and the example load table
#   README.txt        what is needed besides this folder
# Not copied: data/ (contains LGE test data) and tests/data/lge.
# Microsoft MPI is a separate install (msmpi.dll; mpiexec for more than 1 core).

param(
    [string]$Build = "build-windows-mingw",
    [string]$Out = "release"
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
$buildDir = Join-Path $root $Build
$outDir = if ([System.IO.Path]::IsPathRooted($Out)) { $Out } else { Join-Path $root $Out }

foreach ($exe in "pancake.exe", "pancake_gui.exe") {
    if (-not (Test-Path (Join-Path $buildDir $exe))) { throw "$exe not found in $buildDir; build first." }
}
New-Item -ItemType Directory -Force $outDir | Out-Null
Copy-Item (Join-Path $buildDir "pancake.exe"), (Join-Path $buildDir "pancake_gui.exe") $outDir -Force
Get-ChildItem (Join-Path $buildDir "*.dll") | Copy-Item -Destination $outDir -Force
foreach ($folder in "fluids", "examples") {
    $target = Join-Path $outDir $folder
    if (Test-Path $target) { Remove-Item $target -Recurse -Force }
    Copy-Item (Join-Path $root $folder) $target -Recurse
}
# Example cases refer to ../fluids from examples/; that layout is kept.

@"
pancake - journal bearings with refrigerant-oil mixtures

Start: pancake_gui.exe. Open examples\example_steady.txt to try a case.
Every case is a text file; the same result from a terminal:
    pancake.exe -c <case.txt>
    mpiexec -n <cores> pancake.exe -c <case.txt>

Needed besides this folder:
- Microsoft MPI runtime (msmpisetup.exe from Microsoft). It provides msmpi.dll,
  which pancake.exe needs even on 1 core, and mpiexec for more than 1 core.
- Windows 10 or 11 with Direct3D 11 (any GPU or the WARP software renderer).

fluids\example_vg68_rx is SYNTHETIC, example only, not for design.
Add your own pairs under fluids\ (format: fluids\README.md).
"@ | Set-Content -Encoding ascii (Join-Path $outDir "README.txt")

Get-ChildItem $outDir | Select-Object Name, Length | Format-Table -AutoSize
