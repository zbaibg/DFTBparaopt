#!/usr/bin/env bash
#SBATCH --job-name=asm_test
#SBATCH --output=asm_test.out
#SBATCH --error=asm_test.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --partition=batch
#SBATCH --exclude=compute-0-[0-40,44]
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./run.sh
