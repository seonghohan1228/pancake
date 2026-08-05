#!/bin/bash
# Sequential Case 00b dt-sweep tail on afml-r6. One case at a time, n=24.
BASE=~/pancake/cases/00_grid_dt_independence
BIN=~/pancake/build-linux/pancake
NP=24
STATUS=$BASE/sweep_status.log

echo "SWEEP START $(date -Is) host=$(hostname) np=$NP" > $STATUS
for d in dt_T2_fbns dt_T3_fbns dt_T4_fbns; do
  cd $BASE/$d
  rm -rf results
  echo "$d STARTED $(date -Is)" >> $STATUS
  s=$(date +%s)
  mpirun -n $NP $BIN -c config.txt > run.log 2>&1
  rc=$?
  e=$(date +%s)
  echo "$d FINISHED $(date -Is) rc=$rc wall_s=$((e-s))" >> $STATUS
done
echo "SWEEP DONE $(date -Is)" >> $STATUS
