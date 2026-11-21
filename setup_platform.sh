#!/usr/bin/env bash
# =============================================================================
# setup_platform.sh - leukoquant platform setup
# Installs Apptainer (or detects Singularity) on Linux, macOS, and Windows.
# Run once per machine before install.sh.
# =============================================================================
set -euo pipefail

CYAN='\033[0;36m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
log()  { echo -e "${CYAN}[setup] $*${NC}"; }
ok()   { echo -e "${GREEN}[setup] ✔  $*${NC}"; }
warn() { echo -e "${YELLOW}[setup] ⚠  $*${NC}"; }
fail() { echo -e "${RED}[setup] ✘  $*${NC}"; exit 1; }

_add_to_path() {
    local dir="$1"
    local export_line="export PATH=\"$dir:\$PATH\"  # leukoquant"
    for RC in "$HOME/.zshrc" "$HOME/.bashrc" "$HOME/.bash_profile"; do
        if [[ -f "$RC" ]]; then
            # Fix permission issues if the file is not writable
            if [[ ! -w "$RC" ]]; then
                chmod u+w "$RC" 2>/dev/null || true
            fi
            # Only add if not already present
            if ! grep -qF "leukoquant" "$RC"; then
                echo "" >> "$RC"
                echo "$export_line" >> "$RC"
                ok "Added $dir to PATH in $RC"
            fi
        fi
    done
    export PATH="$dir:$PATH"
}

detect_platform() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        echo "mac"
    elif grep -qi microsoft /proc/version 2>/dev/null; then
        echo "wsl2"
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        echo "linux"
    else
        fail "Unsupported OS. Install Apptainer manually: https://apptainer.org"
    fi
}

# ── Linux ─────────────────────────────────────────────────────────────────────
install_linux() {
    log "Linux detected"

    if command -v apptainer &>/dev/null; then
        ok "apptainer already available: $(apptainer --version)"; return 0
    fi
    if command -v singularity &>/dev/null; then
        ok "singularity already available: $(singularity --version)"; return 0
    fi

    warn "Neither apptainer nor singularity found in PATH."

    if ! sudo -n true 2>/dev/null; then
        # No sudo - try conda-forge
        CONDA=$(command -v mamba 2>/dev/null || command -v conda 2>/dev/null || true)
        if [[ -n "$CONDA" ]]; then
            log "No sudo - installing Apptainer via conda-forge..."
            $CONDA install -y -c conda-forge apptainer
            ok "Apptainer $(apptainer --version) installed via conda-forge"
            return 0
        fi
        fail "Cannot install Apptainer automatically (no sudo, no conda).\nPlease ask your sysadmin to install Apptainer, or install conda/mamba first:\n  https://conda-forge.org/download/\nThen re-run: bash setup_platform.sh"
    else
        # Has sudo - apt-get via CIQ PPA
        log "Installing Apptainer via CIQ PPA..."
        sudo add-apt-repository -y ppa:apptainer/ppa
        sudo apt-get update -y
        sudo apt-get install -y apptainer
        ok "Apptainer $(apptainer --version) installed"
    fi
}

# ── WSL2 ──────────────────────────────────────────────────────────────────────
install_wsl2() {
    log "WSL2/Ubuntu detected"

    if command -v apptainer &>/dev/null; then
        ok "apptainer already available: $(apptainer --version)"
    elif command -v singularity &>/dev/null; then
        ok "singularity already available: $(singularity --version)"
    else
        warn "Neither apptainer nor singularity found in PATH."
        if ! sudo -n true 2>/dev/null; then
            CONDA=$(command -v mamba 2>/dev/null || command -v conda 2>/dev/null || true)
            if [[ -n "$CONDA" ]]; then
                log "No sudo - installing Apptainer via conda-forge..."
                $CONDA install -y -c conda-forge apptainer
                ok "Apptainer installed via conda-forge"
            else
                fail "Cannot install Apptainer (no sudo, no conda).\nAsk your sysadmin or install conda/mamba first: https://conda-forge.org/download/"
            fi
        else
            log "Installing Apptainer via CIQ PPA..."
            sudo add-apt-repository -y ppa:apptainer/ppa
            sudo apt-get update -y
            sudo apt-get install -y apptainer
            ok "Apptainer $(apptainer --version) installed"
        fi
    fi

    # ── Drop Windows-side .bat shim ──────────────────────────────────────────
    WIN_HOME=$(wslpath "$(powershell.exe -NoProfile -Command 'echo $env:USERPROFILE' \
        2>/dev/null | tr -d '\r')" 2>/dev/null || true)

    if [[ -n "$WIN_HOME" && -d "$WIN_HOME" ]]; then
        WIN_BIN="$WIN_HOME/bin"
        mkdir -p "$WIN_BIN"

        cat > "$WIN_BIN/leukoquant.bat" << 'BAT'
@echo off
setlocal enabledelayedexpansion
set "ARGS="
for %%A in (%*) do (
    set "ARG=%%~A"
    set "FIRST2=!ARG:~1,1!"
    if "!FIRST2!" == ":" (
        set "DRIVE=!ARG:~0,1!"
        set "REST=!ARG:~2!"
        set "REST=!REST:\=/!"
        set "ARG=/mnt/!DRIVE!/!REST!"
    ) else (
        set "ARG=!ARG:\=/!"
    )
    set "ARGS=!ARGS! "!ARG!""
)
wsl bash -lc "leukoquant !ARGS!"
BAT
        ok "Windows shim created: $WIN_BIN/leukoquant.bat"

        # Add to PowerShell profile (Windows equivalent of .bashrc)
        WIN_PROFILE_WIN=$(powershell.exe -NoProfile -Command 'echo $PROFILE' \
            2>/dev/null | tr -d '\r' || true)
        if [[ -n "$WIN_PROFILE_WIN" ]]; then
            WIN_PROFILE=$(wslpath "$WIN_PROFILE_WIN")
            mkdir -p "$(dirname "$WIN_PROFILE")"
            PROFILE_LINE='$env:PATH += ";$env:USERPROFILE\bin"  # leukoquant'
            if ! grep -qF "leukoquant" "$WIN_PROFILE" 2>/dev/null; then
                echo "$PROFILE_LINE" >> "$WIN_PROFILE"
                ok "Added to PowerShell profile: $WIN_PROFILE_WIN"
            else
                ok "PowerShell profile already configured"
            fi
        fi
        ok "Open a new PowerShell window after running install.sh"
    else
        warn "Could not detect Windows home - skipping .bat shim"
    fi
}

# ── macOS ──────────────────────────────────────────────────────────────────────
install_mac() {
    log "macOS detected"

    if command -v apptainer &>/dev/null; then
        ok "apptainer already available: $(apptainer --version)"; return 0
    fi
    if command -v singularity &>/dev/null; then
        ok "singularity already available: $(singularity --version)"; return 0
    fi

    if ! command -v brew &>/dev/null; then
        log "Installing Homebrew..."
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        [[ -f /opt/homebrew/bin/brew ]] && eval "$(/opt/homebrew/bin/brew shellenv)"
    fi

    log "Installing Lima..."
    brew install lima

    if ! limactl list 2>/dev/null | grep -q "apptainer"; then
        log "Creating Lima VM 'apptainer'..."
        limactl start --name=apptainer template://apptainer 2>/dev/null || \
        limactl start --name=apptainer template://ubuntu --plain 2>/dev/null || \
        fail "Could not start Lima VM. Check Lima docs: https://lima-vm.io"
    else
        limactl start apptainer 2>/dev/null || true
    fi

    if ! limactl shell apptainer -- which apptainer &>/dev/null 2>&1; then
        log "Installing Apptainer inside Lima VM..."
        limactl shell apptainer -- bash -c \
            "sudo add-apt-repository -y ppa:apptainer/ppa && \
             sudo apt-get update -y && sudo apt-get install -y apptainer"
    fi

    WRAPPER_DIR="$HOME/.local/bin"
    mkdir -p "$WRAPPER_DIR"
    cat > "$WRAPPER_DIR/apptainer" << 'WRAPPER'
#!/usr/bin/env bash
exec limactl shell apptainer -- apptainer "$@"
WRAPPER
    chmod +x "$WRAPPER_DIR/apptainer"
    ln -sf "$WRAPPER_DIR/apptainer" "$WRAPPER_DIR/singularity"
    _add_to_path "$WRAPPER_DIR"
    ok "Lima wrapper installed: $WRAPPER_DIR/apptainer"
}

# ── Verify ─────────────────────────────────────────────────────────────────────
verify() {
    echo ""
    log "Verifying container runtime..."
    if command -v apptainer &>/dev/null; then
        ok "apptainer: $(apptainer --version)"
    elif command -v singularity &>/dev/null; then
        ok "singularity: $(singularity --version)"
    else
        fail "No container runtime found after setup."
    fi
}

# ── Main ───────────────────────────────────────────────────────────────────────
PLATFORM=$(detect_platform)
log "Platform: $PLATFORM"
case "$PLATFORM" in
    linux) install_linux ;;
    wsl2)  install_wsl2  ;;
    mac)   install_mac   ;;
esac
verify
echo ""
ok "Platform setup complete. Next step: bash install.sh"