#!/usr/bin/env bash
#SBATCH --job-name=xpid-resbins
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=7-00:00:00
#SBATCH --output=xpid_slurm_%j.out
#SBATCH --error=xpid_slurm_%j.err

set -euo pipefail

WORK_DIR="/y/people/bql506/xhpi_paper/paper_section_04_resolution_dependence"
cd "$WORK_DIR"

# A second Slurm submission must not write the same XPID outputs concurrently.
exec 9>"$WORK_DIR/.xpid_resolution.lock"
if ! flock -n 9; then
    printf '[ERROR] Another Slurm XPID resolution job holds %s/.xpid_resolution.lock\n' \
        "$WORK_DIR" >&2
    exit 75
fi

if ! command -v conda >/dev/null 2>&1; then
    printf '[ERROR] conda is unavailable in the Slurm login environment.\n' >&2
    exit 127
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate xhpi_env

export XPID_JOBS="$SLURM_CPUS_PER_TASK"
export PDB_MIRROR="/vault/xhpi_pdb_July2026/pdb/mmCIF"
export REDO_MIRROR="/vault/xhpi_pdb_July2026/pdb-redo"

printf '[SLURM] Job ID    : %s\n' "$SLURM_JOB_ID"
printf '[SLURM] Node      : %s\n' "$SLURMD_NODENAME"
printf '[SLURM] CPU cores : %s\n' "$SLURM_CPUS_PER_TASK"
printf '[SLURM] Started   : %s\n' "$(date --iso-8601=seconds)"

srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" \
    bash "$WORK_DIR/run_all_resolution_bins.sh"

printf '[SLURM] Finished  : %s\n' "$(date --iso-8601=seconds)"
