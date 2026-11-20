#!/bin/bash --login
# Rerun of the exact search over GF(11) for F5 with 12 gates (|W| = 7), 256 parts.
# usage (from ~/dft/rerun): sbatch -c 32 --output=$HOME/dft/rerun/fs11/slurm_%j.out run_fs11.sh
#SBATCH --job-name=fs11rerun
#SBATCH --partition=cpu
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=256M
#SBATCH --time=1-00:00:00
cd $HOME/dft/rerun/fs11
. $HOME/dft/rerun/common.sh
ulimit -v 1048576
export CODE
seq 0 255 | xargs -P $SLURM_CPUS_PER_TASK -I{} sh -c '
  if [ -f part{}.log ] && grep -q "nodes=" part{}.log; then exit 0; fi
  /usr/bin/time -f "part {} maxrss=%MKB elapsed=%es user=%Us sys=%Ss" $CODE/src/fs_search 5 11 7 7 dft split {} 256 all 1000 > part{}.log 2>&1'
echo 'job done'
