#!/bin/bash
# Case 00b dt-sweep tail T3 then T4, sequential, n=24, threads pinned.
BASE=~/pancake/cases/00_grid_dt_independence
BIN=~/pancake/build-linux/pancake
PIN="-x OMP_NUM_THREADS=1 -x OPENBLAS_NUM_THREADS=1"
STATUS=$BASE/sweep_t3t4_status.log
echo "SWEEP T3/T4 START $(date -Is) host=$(hostname) np=24 pinned" > $STATUS
for d in dt_T3_fbns dt_T4_fbns; do
  cd $BASE/$d
  rm -rf results
  echo "$d STARTED $(date -Is)" >> $STATUS
  s=$(date +%s)
  mpirun $PIN -n 24 $BIN -c config.txt > run.log 2>&1
  rc=$?
  e=$(date +%s)
  steps=$(( $(wc -l < results/diagnostics.csv) - 1 ))
  echo "$d FINISHED $(date -Is) rc=$rc wall_s=$((e-s)) steps=$steps" >> $STATUS
  tar czf ~/${d}_results.tar.gz results run.log config.txt
  echo "$d TARRED $(date -Is) $(du -h ~/${d}_results.tar.gz | cut -f1)" >> $STATUS
done
echo "SWEEP T3/T4 DONE $(date -Is)" >> $STATUS
