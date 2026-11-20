#!/bin/bash --login
# Rerun of the searches over C with the strict backend protocol (one release revision).
# usage (from ~/dft/rerun):
#   G=12 P=2048 D=6 OPTS="--sym2 --galois" sbatch -c 32 --output=$HOME/dft/rerun/g12/slurm_%j.out \
#       run_exact.sh PART [PART ...]
# Parts whose log already has a RESULT line are skipped (resubmission after a time-out).
# Memory: 1 GB per core; 2 GB of virtual memory per process; job-local watchdog (common.sh).
#SBATCH --job-name=f5rerun
#SBATCH --partition=cpu
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
: "${G:?}" "${P:?}" "${D:?}"
cd $HOME/dft/rerun/g$G
. $HOME/dft/rerun/common.sh
ulimit -v 2097152
export G P D OPTS PY CODE
printf '%s\n' "$@" | xargs -P $SLURM_CPUS_PER_TASK -I{} sh -c '
  if [ -f part{}.log ] && grep -q "^RESULT:" part{}.log; then exit 0; fi
  /usr/bin/time -f "part {} maxrss=%MKB elapsed=%es user=%Us sys=%Ss" $PY $CODE/groebner/exact_small.py \
    dft5 $G $OPTS --split {} $P --split-depth $D --show 0 > part{}.log 2>&1'
echo 'job done'
