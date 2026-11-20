#!/bin/bash --login
# Rerun of the local improvement of the F5 circuits of Lev and Winograd over GF(q), windows of
# at most 10 gates and 5 independent entering signals, node limit 2e9 for every window problem.
# usage (from ~/dft/rerun): Q=31 sbatch -c 16 --output=$HOME/dft/rerun/win$Q/slurm_%j.out run_window.sh
#SBATCH --job-name=winrerun
#SBATCH --partition=cpu
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=512M
#SBATCH --time=3-00:00:00
: "${Q:?}"
cd $HOME/dft/rerun/win$Q
. $HOME/dft/rerun/common.sh
ulimit -v 2097152
/usr/bin/time -f "window q=$Q maxrss=%MKB elapsed=%es user=%Us sys=%Ss" $PY $CODE/src/window_search.py \
  $Q 5 5 10 2e9 $SLURM_CPUS_PER_TASK > win_q$Q.log 2>&1
echo "exit $?"
echo 'job done'
