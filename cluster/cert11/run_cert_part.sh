#!/bin/bash --login
# Certificates for one part of the F5 search over Q(omega) (letter's split), then the independent
# check of that part. Code: release revision of the letter's runs (exact_small.py and
# check_certificates.py, unmodified) plus the wrappers cert_split.py and check_certificates_split.py,
# in $HOME/dft/cert11/code. Job array over the parts.
# usage (from ~/dft/cert11): G=11 P=48 D=4 DIR=G11 sbatch --array=0-47 --output=$HOME/dft/cert11/G11/slurm_%A_%a.out run_cert_part.sh
#SBATCH --job-name=cert11
#SBATCH --partition=cpu
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
: "${G:?}" "${P:?}" "${D:?}" "${DIR:?}"
K=${SLURM_ARRAY_TASK_ID:?}
mkdir -p $HOME/dft/cert11/$DIR && cd $HOME/dft/cert11/$DIR
. $HOME/dft/rerun/common.sh
CODE=$HOME/dft/cert11/code
echo "code: $(cat $CODE/REVISION) + wrappers"
ulimit -v $(( ${SLURM_MEM_PER_NODE:-8192} * 1024 * 5 / 4 ))
export RESTART_EVERY=50000 PROGRESS_EVERY=3600
/usr/bin/time -f "generation part $K maxrss=%MKB elapsed=%es" $PY $CODE/src/explore/cert_split.py dft5 $G \
  cert$K.jsonl.gz --split $K $P --split-depth $D --sym2 > gen$K.log 2>&1 || { echo "generation failed"; exit 1; }
/usr/bin/time -f "check part $K maxrss=%MKB elapsed=%es" $PY $CODE/src/explore/check_certificates_split.py dft5 $G \
  cert$K.jsonl.gz --split $K $P --split-depth $D --sym2 --stream > check$K.log 2>&1
ST=$?; echo "check exit $ST"; (exit $ST)
ST=$?
echo "job done (exit $ST)"
exit $ST  # exit status (code review of 17 November 2026, finding 15)
