#!/usr/bin/env bash
# Submit this file with sbatch to run one QM/MM calculation. Set the job resources below.
#SBATCH --job-name=qmmm
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out

set -euo pipefail

usage() {
    printf '%s\n' 'sbatch run.sh [system.json]       Run a calculation (default: water pair)' \
                  'bash run.sh check [system.json]   Check inputs without Slurm' \
                  'bash run.sh status [run_dir]      Show calculation progress'
}

python_bin=${QMMM_PYTHON:-${PYTHON:-python3}}
case ${1:-} in
    check|status)
        exec "$python_bin" -B "$(dirname -- "$0")/qmmm.py" "$@" ;;
    -h|--help) usage; exit 0 ;;
esac
if (( $# > 1 )) || [[ -z ${SLURM_JOB_ID:-} || -z ${SLURM_SUBMIT_DIR:-} ]]; then
    usage >&2
    exit 2
fi

# Slurm executes a spooled copy of this script. Relative input paths refer to
# the submission directory; the installed CLI supplies its package directory.
cd -- "$SLURM_SUBMIT_DIR"
project_dir=${QMMM_PROJECT_DIR:-$PWD}
if [[ ! -f "$project_dir/qmmm.py" ]]; then
    echo 'Submit run.sh from the package directory, or use qmmm submit.' >&2
    exit 2
fi
QMMM_PYTHON=$("$python_bin" -c 'import sys; print(sys.executable)')
export QMMM_PYTHON
run=$("$QMMM_PYTHON" -B "$project_dir/qmmm.py" prepare "$@")
printf 'Run: %s\n' "$run"
printf 'calculation\t%s\n' "$SLURM_JOB_ID" > "$run/jobs.tsv"
exec bash "$run/tools/worker.sh" "$run"
