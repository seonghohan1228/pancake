#!/bin/bash
# usage: run.sh <exe> <dir> <ranks> [petsc opts...]
EXE=$1; D=$2; N=$3; shift 3
cd "$D"
rm -rf results
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
t0=$(date +%s)
mpirun -n $N --bind-to core "$EXE" -c config.txt "$@" > run.log 2>&1
t1=$(date +%s)
if grep -q "Simulation finished successfully" run.log; then st=OK; else st=FAIL; fi
echo "$(basename $D) n=$N $st wall=$((t1-t0))s solves=$(grep -c CONVERGED run.log) recover=$(grep -c "recovered on fallback" run.log) iters=$(grep -o "iterations [0-9]*" run.log | awk "{s+=\$2} END {print s+0}")"
