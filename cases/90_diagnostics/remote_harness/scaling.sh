#!/bin/bash
# Strong-scaling probe on the real dt_T2_fbns config (G6 720x192, dt=2e-5).
# 50-step window (end_t = 1e-3) plus a 10-step companion so fixed MPI/PETSc
# startup can be separated from the marginal per-step cost.
# The box has 28 PHYSICAL cores; n > 28 needs --use-hwthread-cpus and is NOT
# on the same footing (2 ranks share a physical core).
BASE=~/pancake/cases/00_grid_dt_independence
BIN=~/pancake/build-linux/pancake
OUT=$BASE/scaling_results.txt
: > $OUT

for n in 4 8 16 24 32 48 56; do
  if [ $n -gt 28 ]; then LAUNCH="--use-hwthread-cpus"; else LAUNCH=""; fi
  d=$BASE/scaling_probe_n$n
  rm -rf $d && mkdir -p $d

  # 10-step companion
  sed -e "s/^end_t = .*/end_t = 2.0e-4/" -e "s/^write_interval = .*/write_interval = 1.0/" \
      $BASE/dt_T2_fbns/config.txt > $d/config10.txt
  cd $d
  s=$(date +%s.%N); mpirun $LAUNCH -n $n $BIN -c config10.txt > p10.log 2>&1; rc10=$?
  e=$(date +%s.%N); t10=$(echo "$e-$s" | bc)
  rm -rf results

  # 50-step timed window, keeps write_interval so a final frame lands at t=1e-3
  sed -e "s/^end_t = .*/end_t = 1.0e-3/" $BASE/dt_T2_fbns/config.txt > $d/config.txt
  s=$(date +%s.%N); mpirun $LAUNCH -n $n $BIN -c config.txt > p50.log 2>&1; rc50=$?
  e=$(date +%s.%N); t50=$(echo "$e-$s" | bc)

  steps=$(( $(wc -l < results/diagnostics.csv) - 1 ))
  echo "n=$n t10=$t10 t50=$t50 steps=$steps rc10=$rc10 rc50=$rc50 hwt=${LAUNCH:-none}" >> $OUT
  echo "n=$n done t50=$t50" 
  cd $BASE
done
echo "SCALING DONE" >> $OUT
