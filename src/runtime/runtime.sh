#!/usr/bin/env bash
# This is where you tell the job how to find AMBER and ORCA on your cluster.
# Put module names below, or export AMBER_MODULES and ORCA_MODULES before you submit.
amber_modules=""
orca_modules=""

load_modules() {
    local names=$1
    [[ -n $names ]] || return 0
    if ! type module >/dev/null 2>&1 && [[ -r /etc/profile.d/modules.sh ]]; then
        source /etc/profile.d/modules.sh || return 1
    fi
    if ! type module >/dev/null 2>&1; then
        echo "Cannot load $names: the module command is unavailable in this job." >&2
        return 1
    fi
    local -a modules
    read -r -a modules <<< "$names"
    module purge && module load "${modules[@]}"
}

load_amber() {
    load_modules "${AMBER_MODULES-$amber_modules}" || return 1
    AMBER_SANDER=$(type -P sander) || { echo 'Cannot find sander. Check AMBER_MODULES or PATH.' >&2; return 1; }
    type -P cpptraj >/dev/null || { echo 'Cannot find cpptraj. Check AMBER_MODULES or PATH.' >&2; return 1; }
    export AMBER_SANDER
    printf 'Sander: %s\n' "$AMBER_SANDER"
}

load_orca() {
    # Sander prepends its private ORCA wrapper to PATH. Its child must run the real driver.
    if [[ -n ${AMBER_ORCA_SHIM_DIR:-} ]]; then
        PATH=${PATH#"$AMBER_ORCA_SHIM_DIR":}
        export PATH
    fi
    load_modules "${ORCA_MODULES-$orca_modules}" || return 1
    local driver=${ORCA_BIN:-}
    if [[ -z $driver ]]; then
        driver=$(type -P orca) || { echo 'Cannot find orca. Check ORCA_MODULES, ORCA_BIN or PATH.' >&2; return 1; }
    fi
    if [[ $driver != /* || ! -f $driver || ! -x $driver ||
          ( -n ${AMBER_ORCA_RUN:-} && $driver -ef "$AMBER_ORCA_RUN/bin/orca" ) ||
          ( -n ${AMBER_ORCA_SHIM_DIR:-} && $driver -ef "$AMBER_ORCA_SHIM_DIR/orca" ) ]]; then
        echo "Expected the absolute path of the ORCA executable, got: $driver" >&2
        return 1
    fi
    export ORCA_BIN="$driver"
    printf 'ORCA: %s\n' "$ORCA_BIN"
}
