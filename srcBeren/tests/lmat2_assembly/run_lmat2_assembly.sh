#!/usr/bin/env bash
# Run the lmat2_assembly test under several OMP thread counts.
set -u

BIN_PATH="${1:-$PWD/_build_release/bin/lmat2_assembly}"

if [[ ! -x "$BIN_PATH" ]]; then
    echo "error: test binary not found: $BIN_PATH" >&2
    exit 2
fi

overall=0
for threads in 1 2 3 7 16; do
    echo ""
    echo "############################################################"
    echo "# OMP_NUM_THREADS = ${threads}"
    echo "############################################################"
    OMP_NUM_THREADS="${threads}" "$BIN_PATH"
    rc=$?
    echo ">>> lmat2_assembly with OMP_NUM_THREADS=${threads}: exit=${rc}"
    if [[ $rc -ne 0 ]]; then
        overall=1
    fi
done

echo ""
if [[ $overall -eq 0 ]]; then
    echo "ALL THREAD-COUNT RUNS PASSED"
else
    echo "SOME THREAD-COUNT RUNS FAILED" >&2
fi
exit $overall