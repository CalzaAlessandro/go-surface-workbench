#!/bin/bash
#SBATCH --job-name=GO-Espaloma-control
#SBATCH --partition=ciamician_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=5-00:00:00
#SBATCH --output=benchmark-%j.out
#SBATCH --error=benchmark-%j.err
set -euo pipefail
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
export XDG_CACHE_HOME="$PWD/cache"
export TORCH_HOME="$PWD/cache/torch"
export TMPDIR="$PWD/tmp"
export DGLBACKEND=pytorch
export DGLDEFAULTDIR="$PWD/cache/dgl"
.venv-release/bin/python -u run.py --backend dgl
