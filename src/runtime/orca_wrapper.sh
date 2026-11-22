#!/bin/bash
# Sander calls this script at each QM step. It loads ORCA and passes the input to the Python bridge.
set -euo pipefail
if [[ $# == 1 && "$1" == --check-launcher ]]; then
    printf 'QMMM_LAUNCHER_OK\n'
    exit 0
fi
: "${AMBER_ORCA_RUN:?Missing run directory}"
: "${AMBER_ORCA_CONTEXT:?Missing bridge context}"
: "${QMMM_PYTHON:?Set QMMM_PYTHON to the absolute Python 3.11+ interpreter used by the submission script}"
if [[ "$QMMM_PYTHON" != /* || ! -x "$QMMM_PYTHON" ]]; then
    printf 'Expected an executable absolute Python path, got %s\n' "$QMMM_PYTHON" >&2
    exit 2
fi
source "$AMBER_ORCA_RUN/tools/runtime.sh"
load_orca >> orca_runtime.log 2>&1
exec "$QMMM_PYTHON" -B -E -s \
    "$AMBER_ORCA_RUN/tools/orca_bridge.py" "$@"
