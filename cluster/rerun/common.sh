# Sourced by the job scripts of the rerun. Watchdog local to this job: it sums the resident
# memory of the processes in the session of the job script only (other jobs on the same node,
# also of the same user, are not counted), and cancels the job above 80% of the request.
hostname; lscpu | grep 'Model name'
SID=$(ps -o sid= -p $$ | tr -d ' ')
LIMIT_KB=$(( ${SLURM_MEM_PER_NODE:-$(( SLURM_CPUS_PER_TASK * ${SLURM_MEM_PER_CPU:-1024} ))} * 1024 * 8 / 10 ))
( while sleep 30; do
    rss=$(ps -o rss= -s "$SID" | awk '{s+=$1} END{print s+0}')
    if [ "$rss" -gt "$LIMIT_KB" ]; then echo "WATCHDOG: RSS ${rss} KB > ${LIMIT_KB} KB, cancelling"; scancel $SLURM_JOB_ID; fi
  done ) &
CODE=$HOME/dft/rerun/code
PY=$HOME/conda/envs/singular/bin/python
export SINGULAR=$HOME/conda/envs/singular/bin/Singular PROGRESS_EVERY=1800
echo "code: $(cat $CODE/REVISION)"
