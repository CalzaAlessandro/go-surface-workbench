#!/bin/bash -l
#SBATCH --job-name=GO-DFT-bonded
#SBATCH --partition=ciamician_cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --array=0-7%2
#SBATCH --output=dft-%A_%a.out
#SBATCH --error=dft-%A_%a.err
set -eo pipefail
module load gaussian/16
set -u
export OMP_NUM_THREADS=8
model_id=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" tasks.txt)
[[ "$model_id" == pyrene_* ]] || exit 2
cd "$model_id"
export GAUSS_SCRDIR
GAUSS_SCRDIR=$(mktemp -d "${SLURM_TMPDIR:-/tmp}/go-dft-${SLURM_JOB_ID}.XXXXXX")
printf '%s\n' "$GAUSS_SCRDIR" > scratch_path.txt
set +e
g16 < reference.gjf > gaussian.log 2> gaussian.stderr
status=$?
set -e
printf '%s\n' "$status" > gaussian_exit_code.txt
[[ "$status" -eq 0 ]] || exit "$status"
[[ $(grep -c 'Normal termination of Gaussian' gaussian.log) -ge 2 ]] || exit 3
if grep -q 'Error termination' gaussian.log; then exit 4; fi
formchk reference.chk harmonic.fchk > formchk_harmonic.out 2> formchk_harmonic.err
formchk stability.chk stability.fchk > formchk_stability.out 2> formchk_stability.err
grep -q 'Cartesian Force Constants' harmonic.fchk
printf 'DFT_FINISHED_NEEDS_GEOMETRY_FREQUENCY_SPIN_AND_HESSIAN_QA\n'
