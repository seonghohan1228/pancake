#!/bin/bash
# Strong-scaling probe, MPI-only (OMP_NUM_THREADS=1).
# PETSc here links OpenBLAS+libgomp; unpinned, every rank spawns ~30 threads and
# the node thrashes. All timings below are with threading pinned to 1/rank.
BASE=~/pancake/cases/00_grid_dt_independence
BIN=~/pancake/build-linux/pancake
OUT=$BASE/scaling_results2.txt
PIN="-x OMP_NUM_THREADS=1 -x OPENBLAS_NUM_THREADS=1"
: > $OUT

for n in 4 8 16 24 32 48 56; do
  if [ $n -gt 28 ]; then LAUNCH="--use-hwthread-cpus"; else LAUNCH=""; fi
  d=$BASE/scaling_probe_n$n
  rm -rf $d && mkdir -p $d
  sed -e "s/^end_t = .*/end_t = 2.0e-4/" -e "s/^write_interval = .*/write_interval = 1.0/" \
      $BASE/dt_T2_fbns/config.txt > $d/config10.txt
  sed -e "s/^end_t = .*/end_t = 1.0e-3/" $BASE/dt_T2_fbns/config.txt > $d/config.txt
  cd $d

  s=$(date +%s.%N); mpirun $LAUNCH $PIN -n $n $BIN -c config10.txt > p10.log 2>&1; rc10=$?
  e=$(date +%s.%N); t10=$(echo "$e-$s" | bc); rm -rf results

  s=$(date +%s.%N); mpirun $LAUNCH $PIN -n $n $BIN -c config.txt > p50.log 2>&1; rc50=$?
  e=$(date +%s.%N); t50=$(echo "$e-$s" | bc)

  steps=$(( $(wc -l < results/diagnostics.csv) - 1 ))
  W=$(python3 ~/gate.py . 2>/dev/null | awk '$1=="W"{print $2}')
  echo "n=$n t10=$t10 t50=$t50 steps=$steps W=$W rc10=$rc10 rc50=$rc50 hwt=${LAUNCH:-none}" >> $OUT
  cd $BASE
done

# threading-effect control at the production rank count
d=$BASE/scaling_probe_n24; cd $d
s=$(date +%s.%N); mpirun -n 24 $BIN -c config10.txt > p10_unpinned.log 2>&1
e=$(date +%s.%N)
echo "n=24 UNPINNED t10=$(echo "$e-$s" | bc)" >> $OUT
echo "SCALING2 DONE" >> $OUT
