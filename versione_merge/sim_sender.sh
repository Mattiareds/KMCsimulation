#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------------
# EDIT THIS: name/path of your executable, including any flags it needs.
# It is run from the current directory, so input files must already be here.
# ----------------------------------------------------------------------------
EXECUTABLE="./Nanopyramid"

# Number of simulations to run
N=100

# Prefix for the output subfolders (run_001, run_002, ...)
PREFIX="run"

# Max wall-clock time per simulation, in seconds (1 hour). This is real
# elapsed time, regardless of how much CPU the process actually uses.
MAX_WALL_SECONDS=1200

for i in $(seq -w 1 "$N"); do
    folder="${PREFIX}_${i}"
    echo "=== Simulation ${i}/${N} ==="

    # 1) start the simulation under `timeout`, which kills it (SIGTERM by
    #    default) if it's still running after MAX_WALL_SECONDS of real time.
    #    `timeout` blocks until the process exits or is killed, so this also
    #    covers "wait until it's finished".
    if timeout "${MAX_WALL_SECONDS}" "${EXECUTABLE}"; then
        echo "Run ${i} finished normally."
    else
        status=$?
        if [ "${status}" -eq 124 ]; then
            echo "Run ${i} exceeded ${MAX_WALL_SECONDS}s of wall time and was killed." >&2
        else
            echo "Run ${i} exited with status ${status}." >&2
        fi
    fi

    # 2) create the subfolder for this run
    mkdir -p "${folder}"

    # 3) move all *.out files produced by the run into that subfolder
    if compgen -G "*.out" > /dev/null; then
        mv ./*.out "${folder}/"
    else
        echo "WARNING: no .out files found after run ${i}" >&2
    fi
done

echo "Done: ${N} simulations completed."
