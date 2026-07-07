#!/usr/bin/env bash
# run.sh — execute a Python script (or module) inside the xsuite_env conda
# environment without needing to activate it first.
#
# Usage
# -----
#   ./run.sh Animations/save_sequence_SPS.py
#   ./run.sh Animations/phaseSpaceAnimation.py
#   ./run.sh -m pytest
#   ./run.sh -m pytest tests/test_crystal_extraction.py
#
# The first argument is passed directly to the Python interpreter, so any
# valid 'python <args>' invocation works here.

set -euo pipefail

ENV_NAME="xsuite_env"

if ! conda env list 2>/dev/null | grep -qE "^${ENV_NAME}[[:space:]]"; then
    echo "ERROR: conda environment '${ENV_NAME}' not found."
    echo "       Run ./install.sh first to create it."
    exit 1
fi

exec conda run --no-capture-output -n "$ENV_NAME" python "$@"
