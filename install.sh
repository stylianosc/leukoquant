#!/usr/bin/env bash
# =============================================================================
# install.sh - leukoquant Python package installer
# Creates a dedicated venv and permanently adds it to PATH via shell rc files.
# After this runs, 'leukoquant' works in every new terminal - no activation.
# =============================================================================
set -euo pipefail

CYAN='\033[0;36m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
log()  { echo -e "${CYAN}[install] $*${NC}"; }
ok()   { echo -e "${GREEN}[install] ✔  $*${NC}"; }
warn() { echo -e "${YELLOW}[install] ⚠  $*${NC}"; }
fail() { echo -e "${RED}[install] ✘  $*${NC}"; exit 1; }

VENV_DIR="$HOME/.venvs/leukoquant"
VENV_BIN="$VENV_DIR/bin"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$SCRIPT_DIR/leukoquant"

# ── Checks ────────────────────────────────────────────────────────────────────
# Find a stable Python >= 3.8 - check versioned binaries and current python3
PYTHON_BIN=""
for candidate in python3.13 python3.12 python3.11 python3.10 python3.9 python3.8 python3; do
    if command -v "$candidate" &>/dev/null; then
        _rl=$("$candidate" -c "import sys; print(sys.version_info.releaselevel)" 2>/dev/null)
        _maj=$("$candidate" -c "import sys; print(sys.version_info.major)" 2>/dev/null)
        _min=$("$candidate" -c "import sys; print(sys.version_info.minor)" 2>/dev/null)
        if [[ "$_rl" == "final" && "$_maj" -eq 3 && "$_min" -ge 8 ]]; then
            PYTHON_BIN=$(command -v "$candidate")
            break
        fi
    fi
done

# If no stable Python found, try to install one via conda automatically
if [[ -z "$PYTHON_BIN" ]]; then
    if command -v conda &>/dev/null; then
        warn "No stable Python >= 3.8 found. Installing Python 3.11 via conda..."
        conda install python=3.11 -y -q
        # hash -r may not update PATH in all shells - resolve directly from the conda prefix
        CONDA_PREFIX_BIN=$(conda info --base 2>/dev/null)/bin
        for candidate in "$CONDA_PREFIX_BIN/python3.11" "$CONDA_PREFIX_BIN/python3"; do
            if [[ -x "$candidate" ]]; then
                _rl=$("$candidate" -c "import sys; print(sys.version_info.releaselevel)" 2>/dev/null)
                if [[ "$_rl" == "final" ]]; then
                    PYTHON_BIN="$candidate"
                    break
                fi
            fi
        done
    fi
fi

if [[ -z "$PYTHON_BIN" ]]; then
    fail "No stable Python >= 3.8 found and conda is not available.\nInstall Python manually: https://www.python.org/downloads/"
fi

PY_VER=$("$PYTHON_BIN" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
ok "Python $PY_VER detected ($PYTHON_BIN)"

[[ -d "$PACKAGE_DIR" ]] || fail "Package not found at $PACKAGE_DIR.\nRun this script from the repo root."

if ! command -v apptainer &>/dev/null && ! command -v singularity &>/dev/null; then
    warn "No container runtime found. Run setup_platform.sh first."
    echo ""
fi

# ── Create venv ───────────────────────────────────────────────────────────────
if [[ -d "$VENV_DIR" ]]; then
    # Check both that pip exists and that the venv Python matches the selected stable Python
    VENV_PY_VER=$("$VENV_BIN/python" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')" 2>/dev/null || echo "")
    if [[ -x "$VENV_BIN/pip" && "$VENV_PY_VER" == "$PY_VER" ]]; then
        ok "Existing venv found at $VENV_DIR (Python $VENV_PY_VER) - reusing it."
    else
        warn "venv at $VENV_DIR is outdated or broken (has Python $VENV_PY_VER, need $PY_VER) - recreating..."
        rm -rf "$VENV_DIR"
        log "Creating virtual environment at $VENV_DIR ..."
        "$PYTHON_BIN" -m venv "$VENV_DIR"
    fi
else
    log "Creating virtual environment at $VENV_DIR ..."
    "$PYTHON_BIN" -m venv "$VENV_DIR"
fi

# ── Install package ───────────────────────────────────────────────────────────
log "Upgrading pip..."
"$VENV_BIN/pip" install --upgrade pip -q

# Install greenlet from binary wheel only - the system g++ may be too old to compile it
"$VENV_BIN/pip" install --only-binary=greenlet greenlet -q

log "Installing leukoquant and all dependencies..."
# pip compares installed versions against pyproject.toml and only installs what has changed
"$VENV_BIN/pip" install -e "$SCRIPT_DIR" -q

ok "Installed: $( "$VENV_BIN/leukoquant" --version 2>&1 )"

# ── Add venv bin to PATH permanently via shell rc (Option A) ──────────────────
log "Adding $VENV_BIN to PATH in shell rc files..."
EXPORT_LINE="export PATH=\"$VENV_BIN:\$PATH\"  # leukoquant"
ADDED=false

for RC in "$HOME/.zshrc" "$HOME/.bashrc" "$HOME/.bash_profile"; do
    if [[ -f "$RC" ]]; then
        # Fix permission issues if the file is not writable
        if [[ ! -w "$RC" ]]; then
            chmod u+w "$RC" 2>/dev/null || warn "Could not fix permissions for $RC; skipping"
            continue
        fi
        # Match the actual venv bin path, not the bare word "leukoquant" --
        # an unrelated PATH line from another tool's setup step (e.g. one
        # that happens to end in a "# leukoquant" comment) can otherwise
        # false-positive this check and make the script wrongly skip adding
        # our own PATH entry.
        if ! grep -qF "$VENV_BIN" "$RC"; then
            { echo ""; echo "$EXPORT_LINE"; } >> "$RC" 2>/dev/null || warn "Could not write to $RC; skipping"
            ok "Updated $RC"
            ADDED=true
        else
            ok "$RC already has leukoquant PATH entry"
            ADDED=true
        fi
    fi
done

if [[ "$ADDED" == false ]]; then
    echo "$EXPORT_LINE" >> "$HOME/.bashrc"
    ok "Created ~/.bashrc with PATH entry"
fi

export PATH="$VENV_BIN:$PATH"

# ── Final verification ────────────────────────────────────────────────────────
echo ""
log "Verifying..."
if command -v leukoquant &>/dev/null; then
    ok "leukoquant is available: $(leukoquant --version 2>&1)"
else
    warn "leukoquant not yet in PATH for this session."
    warn "Run: export PATH=\"$VENV_BIN:\$PATH\""
fi

echo ""
ok "Installation complete!"
echo "  ┌──────────────────────────────────────────────────────┐"
echo "  │  Open a NEW terminal and type: leukoquant --help │"
echo "  └──────────────────────────────────────────────────────┘"
echo ""