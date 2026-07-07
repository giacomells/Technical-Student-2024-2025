#!/usr/bin/env bash
# install.sh — set up a Miniforge environment with all project dependencies.
#
# Supported platforms
# -------------------
#   macOS  : Intel (x86_64) and Apple Silicon (arm64)
#   Linux  : x86_64 and aarch64
#   Windows: not natively supported; install WSL and run this script inside
#            the resulting Linux shell (see instructions printed below).
#
# What this script does
# ---------------------
#   1. Detects the OS and CPU architecture.
#   2. Downloads and installs Miniforge if conda / mamba is not already present.
#   3. Creates a conda environment named 'xsuite_env' (skipped if it exists).
#   4. Installs C compilers (required for Xsuite CPU tracking kernels).
#      On Intel Macs, also installs llvm-openmp for OpenMP support.
#   5. Installs all Python dependencies from requirements-core.txt via pip.
#
# Usage
# -----
#   chmod +x install.sh
#   ./install.sh
#
# After installation, activate the environment with:
#   conda activate xsuite_env

set -euo pipefail

# ── OS / architecture detection ──────────────────────────────────────────────

OS="$(uname -s)"
ARCH="$(uname -m)"
ENV_NAME="xsuite_env"
PYTHON_VERSION="3.11"

echo "Detected OS   : ${OS}"
echo "Architecture  : ${ARCH}"

case "$OS" in
    Darwin)
        MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-${ARCH}.sh"
        ;;
    Linux)
        MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-${ARCH}.sh"
        ;;
    MINGW*|CYGWIN*|MSYS*)
        echo ""
        echo "Windows is not natively supported by Xsuite."
        echo "Please install Windows Subsystem for Linux (WSL) and re-run this"
        echo "script from inside the WSL terminal."
        echo ""
        echo "  Step 1 — open PowerShell or CMD as Administrator and run:"
        echo "             wsl --install"
        echo "  Step 2 — restart your machine."
        echo "  Step 3 — open the WSL terminal and navigate to this directory:"
        echo "             cd /mnt/c/path/to/this/repo"
        echo "  Step 4 — run this script:"
        echo "             bash install.sh"
        echo ""
        echo "WSL guide: https://learn.microsoft.com/en-us/windows/wsl/install"
        exit 1
        ;;
    *)
        echo "ERROR: unsupported OS '${OS}'. Aborting."
        exit 1
        ;;
esac

# ── install Miniforge if conda / mamba is not yet available ──────────────────

if command -v conda &>/dev/null; then
    echo "conda found   : $(conda --version)"
elif command -v mamba &>/dev/null; then
    echo "mamba found   : $(mamba --version)"
else
    echo ""
    echo "conda/mamba not found — downloading Miniforge..."
    INSTALLER="$(basename "$MINIFORGE_URL")"
    curl -fsSL "$MINIFORGE_URL" -o "$INSTALLER"
    bash "$INSTALLER" -b -p "$HOME/miniforge3"
    rm -f "$INSTALLER"
    echo "Miniforge installed at ~/miniforge3"
fi

# ── activate conda in this shell session ─────────────────────────────────────

CONDA_BASE="$(conda info --base 2>/dev/null || echo "$HOME/miniforge3")"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# ── create the project environment ───────────────────────────────────────────

if conda env list | grep -qE "^${ENV_NAME}[[:space:]]"; then
    echo ""
    echo "Environment '${ENV_NAME}' already exists — skipping creation."
else
    echo ""
    echo "Creating conda environment '${ENV_NAME}' (Python ${PYTHON_VERSION})..."
    conda create -y -n "$ENV_NAME" "python=${PYTHON_VERSION}"
fi

conda activate "$ENV_NAME"

# ── compilers (needed to build Xsuite CPU kernels) ───────────────────────────

echo ""
echo "Installing C compilers..."
conda install -y compilers

# Intel Macs need a separate OpenMP package.
if [[ "$OS" == "Darwin" && "$ARCH" == "x86_64" ]]; then
    echo "Intel Mac detected — installing llvm-openmp for OpenMP support..."
    conda install -y llvm-openmp
fi

# ── Python dependencies ───────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REQUIREMENTS="${SCRIPT_DIR}/requirements-core.txt"

echo ""
if [[ -f "$REQUIREMENTS" ]]; then
    echo "Installing project dependencies from requirements-core.txt..."
    pip install -r "$REQUIREMENTS"
else
    echo "requirements-core.txt not found — installing core packages directly..."
    pip install xsuite xcoll cpymad numpy scipy matplotlib pandas pytest
fi

# ── VS Code interpreter setting ───────────────────────────────────────────────
# Write .vscode/settings.json so VS Code automatically selects the right
# Python interpreter without the user having to pick it manually.

PYTHON_BIN="$(conda run -n "$ENV_NAME" python -c 'import sys; print(sys.executable)')"
VSCODE_DIR="${SCRIPT_DIR}/.vscode"
VSCODE_SETTINGS="${VSCODE_DIR}/settings.json"

mkdir -p "$VSCODE_DIR"

if [[ -f "$VSCODE_SETTINGS" ]]; then
    echo ""
    echo "Note: ${VSCODE_SETTINGS} already exists — not overwriting."
    echo "      To auto-select the interpreter, add this line manually:"
    echo "        \"python.defaultInterpreterPath\": \"${PYTHON_BIN}\""
else
    cat > "$VSCODE_SETTINGS" <<EOF
{
    "python.defaultInterpreterPath": "${PYTHON_BIN}"
}
EOF
    echo ""
    echo "VS Code interpreter set to: ${PYTHON_BIN}"
fi

# ── summary ───────────────────────────────────────────────────────────────────

echo ""
echo "┌─────────────────────────────────────────────────────┐"
echo "│  Installation complete!                             │"
echo "│                                                     │"
echo "│  Option A — VS Code (no action needed):             │"
echo "│    The interpreter is already configured.           │"
echo "│    Open any .py file and press ▶ Run.               │"
echo "│                                                     │"
echo "│  Option B — terminal (one-time per session):        │"
printf  "│    conda activate %-33s│\n" "${ENV_NAME}"
echo "│                                                     │"
echo "│  Option C — run a script without activating:        │"
echo "│    ./run.sh Animations/save_sequence_SPS.py         │"
echo "│                                                     │"
echo "│  Verify with the test suite:                        │"
echo "│    pytest           (after conda activate)          │"
echo "│    ./run.sh -m pytest                               │"
echo "│                                                     │"
echo "│  Full Xsuite guide:                                 │"
echo "│    https://xsuite.readthedocs.io/en/latest/         │"
echo "│    installation.html                                │"
echo "└─────────────────────────────────────────────────────┘"
