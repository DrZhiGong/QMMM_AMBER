#!/usr/bin/env bash
# Once Slurm starts the job, this script checks the software and runs the calculation.
set -euo pipefail
run=$1
source "$run/tools/runtime.sh"

# The allocation is already running, so report configuration failures in its directory.
fail_setup() {
    printf '%s\n' "$1" | tee "$run/logs/submission_error.txt" >&2
    exit 1
}
ranks=$("$QMMM_PYTHON" -c 'import json, sys; print(json.load(open(sys.argv[1]))["orca_ranks"])' "$run/config.json")
allocated=${SLURM_NTASKS:-1}
if [[ ! $allocated =~ ^[1-9][0-9]*$ ]] || (( ranks > allocated )); then
    fail_setup "ORCA needs $ranks ranks but Slurm allocated $allocated. Adjust --ntasks or orca_ranks."
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if ! { (load_orca) && load_amber; } > "$run/logs/modules.log" 2>&1; then
    cat "$run/logs/modules.log" >&2
    fail_setup "Software setup failed. See $run/logs/modules.log"
fi
exec "$QMMM_PYTHON" -B -E -s "$run/tools/custom_run.py" "$run"
