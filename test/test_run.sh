#!/usr/bin/env bash
# Try the Slurm job script with stand-ins for AMBER and ORCA so the tests can run locally.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf -- "$tmp"' EXIT
export SLURM_JOB_ID=123 SLURM_SUBMIT_DIR="$tmp/submit" QMMM_PROJECT_DIR="$root"
export QMMM_PYTHON="$tmp/bin/python" TEST_BIN="$tmp/bin"
export PATH="$TEST_BIN:$PATH" AMBER_MODULES="" ORCA_MODULES=""
unset ORCA_BIN AMBER_ORCA_SHIM_DIR AMBER_ORCA_RUN
mkdir -p "$TEST_BIN" "$SLURM_SUBMIT_DIR"
for name in sander cpptraj orca; do
    printf '#!/bin/sh\nexit 0\n' > "$TEST_BIN/$name"
done

cat > "$QMMM_PYTHON" <<'SH'
#!/usr/bin/env bash
set -eu
if [[ $1 == -c && $2 == *sys.executable* ]]; then
    echo "$QMMM_PYTHON"
elif [[ $1 == -c && $2 == *orca_ranks* ]]; then
    echo 2
elif [[ $3 == prepare ]]; then
    [[ $PWD == "$SLURM_SUBMIT_DIR" ]]
    mkdir -p "$TEST_RUN"/{tools,logs,bin}
    cp "$QMMM_PROJECT_DIR/src/runtime/"{runtime,worker}.sh "$TEST_RUN/tools/"
    cp "$QMMM_PROJECT_DIR/src/runtime/orca_wrapper.sh" "$TEST_RUN/bin/orca"
    printf '%s\n' "${4-default}" > "$TEST_RUN/system"
    echo "$TEST_RUN"
elif [[ $4 == */custom_run.py ]]; then
    [[ $AMBER_SANDER == "$TEST_BIN/sander" ]]
    export AMBER_ORCA_RUN="$TEST_RUN" AMBER_ORCA_CONTEXT="$TEST_RUN/context.json"
    export AMBER_ORCA_SHIM_DIR="$TEST_RUN/shim"
    mkdir "$AMBER_ORCA_SHIM_DIR"
    cp "$TEST_RUN/bin/orca" "$AMBER_ORCA_SHIM_DIR/orca"
    export PATH="$AMBER_ORCA_SHIM_DIR:$PATH"
    cd "$TEST_RUN"
    bash "$AMBER_ORCA_SHIM_DIR/orca" orca.inp
elif [[ $4 == */orca_bridge.py ]]; then
    [[ $ORCA_BIN == "$TEST_BIN/orca" && $5 == orca.inp ]]
    echo complete > "$TEST_RUN/result"
else
    exit 99
fi
SH
chmod +x "$TEST_BIN/"*

export SLURM_NTASKS=2 TEST_RUN="$tmp/default"
bash "$root/run.sh" > "$tmp/stdout"
grep -Fxq "Run: $TEST_RUN" "$tmp/stdout"
grep -Fxq default "$TEST_RUN/system"
grep -Fxq $'calculation\t123' "$TEST_RUN/jobs.tsv"
grep -Fxq complete "$TEST_RUN/result"

export TEST_RUN="$tmp/explicit"
bash "$root/run.sh" 'my system/system.json' > "$tmp/stdout"
grep -Fxq 'my system/system.json' "$TEST_RUN/system"
grep -Fxq complete "$TEST_RUN/result"

export TEST_RUN="$tmp/too-small" SLURM_NTASKS=1
if bash "$root/run.sh" > "$tmp/stdout" 2>&1; then exit 1; fi
grep -Fq 'Slurm allocated 1' "$TEST_RUN/logs/submission_error.txt"
[[ ! -e "$TEST_RUN/result" ]]

export TEST_RUN="$tmp/missing" SLURM_NTASKS=2 ORCA_BIN="$tmp/missing-orca"
if bash "$root/run.sh" > "$tmp/stdout" 2>&1; then exit 1; fi
grep -Fq "$ORCA_BIN" "$TEST_RUN/logs/modules.log"
grep -Fq 'Software setup failed' "$TEST_RUN/logs/submission_error.txt"
[[ ! -e "$TEST_RUN/result" ]]
echo 'Slurm script tests passed.'
