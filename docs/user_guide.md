# leukoquant - Installation & User Guide

---

## Before You Start

| Requirement        | Check                 | Notes                                                                                                                                       |
| ------------------ | --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| Python ≥ 3.8      | `python3 --version` | Required                                                                                                                                    |
| Git                | `git --version`     | Required                                                                                                                                    |
| 16 GB RAM          | -                    | Recommended                                                                                                                                 |
| 50 GB free disk    | -                    | Containers + outputs                                                                                                                        |
| FreeSurfer licence | -                    | Required for `process-all` / `process-recon`. Free from [surfer.nmr.mgh.harvard.edu](https://surfer.nmr.mgh.harvard.edu/registration.html) |

**Windows only - install WSL2 first (one reboot required):**

```powershell
# Run in PowerShell as Administrator, then reboot
wsl --install
```

After reboot, open the **Ubuntu** app from the Start Menu for all steps below.

---

## Step 1 - Get the Code

```bash
git clone https://github.com/stylianosc/leukoquant.git
cd leukoquant
```

**What this does:** Downloads the LeukoQuant repository to your machine and enters the directory. All subsequent scripts must be run from this directory.

---

## Step 2 - Platform Setup (one-time per machine)

```bash
bash setup_platform.sh
```

### Linux

1. Checks if `apptainer` or `singularity` are already in PATH - if yes, exits immediately
2. **No sudo + conda/mamba available:** installs Apptainer via `conda-forge` (no root needed)
3. **Has sudo:** installs Apptainer from the official CIQ PPA via `apt-get`
4. **No sudo, no conda:** prints a clear message - ask your sysadmin or install conda first

> No sudo and no conda? Install Miniforge (minimal conda) first:
>
> ```bash
> curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh | bash
> # open a new terminal, then re-run setup_platform.sh
> ```

### macOS

1. Checks if `apptainer` or `singularity` are already in PATH - if yes, exits immediately
2. Installs **Homebrew** if not present (the macOS package manager)
3. Installs **Lima** - a lightweight Linux VM manager for macOS
4. Creates a Lima VM named `apptainer` using an Ubuntu template
5. Installs Apptainer inside the VM
6. Creates a transparent wrapper at `~/.local/bin/apptainer` that silently proxies all commands into the VM - you never interact with Lima directly
7. Symlinks `~/.local/bin/singularity` → `apptainer`
8. Adds `~/.local/bin` to PATH in `~/.zshrc` / `~/.bashrc`

> Your home directory is automatically mounted inside the Lima VM - all your data is accessible with no extra configuration.

#### Lima VM Setup Prompt

During setup, you may see a Lima configuration prompt:

```
[setup] Creating Lima VM 'apptainer'...
? Creating an instance "appainer" [Use arrows to move, type to filter]
> Proceed with the current configuration
  Open an editor to review or modify the current configuration
  Choose another template (docker, podman, archlinux, fedora, ...)
  Exit
```

**Action:** Select "Proceed with the current configuration" and press Enter. This will:
- Create a lightweight Ubuntu Linux VM inside Lima
- Install Apptainer automatically
- Set up transparent command proxying so you don't interact with Lima directly

The setup runs once and takes a few minutes. After completion, `apptainer` commands work transparently on your Mac.

### Windows (run inside the Ubuntu WSL2 app)

1. Checks if `apptainer` or `singularity` are already in PATH - if yes, exits immediately
2. Installs Apptainer inside WSL2 via `apt-get` (or conda-forge if no sudo)
3. Creates `%USERPROFILE%\bin\leukoquant.bat` - a Windows batch file that routes all `leukoquant` commands to WSL2, automatically converting Windows paths (`C:\data\file.nii.gz`) to WSL2 paths (`/mnt/c/data/file.nii.gz`)
4. Writes one line to your **PowerShell profile** (`$PROFILE`) - the Windows equivalent of `.bashrc` - so `leukoquant` is available in every new PowerShell window automatically

---

## Step 3 - Install the Python Package (one-time per machine)

```bash
bash install.sh
```

**What this installs:**

1. Creates a Python virtual environment at `~/.venvs/leukoquant/`
2. Upgrades pip
3. Installs all Python dependencies into the venv:
   - `snakemake 9.20.0` - workflow engine
   - `nibabel 5.3.2` - NIfTI file I/O
   - `numpy`, `pandas`, `scipy`, `scikit-image` - scientific computing
   - `dipy 1.10.0` - diffusion MRI
   - `click 8.3.0` - CLI framework
   - `PyYAML 6.0.2` - config file support
   - `dcm2niix` - DICOM conversion
4. Adds `~/.venvs/leukoquant/bin` to PATH permanently in `~/.zshrc` / `~/.bashrc` - **no activation needed in any future session**

**Open a new terminal after this step**, then verify:

```bash
leukoquant --version
```

---

## Step 4 - FreeSurfer Licence

Place your FreeSurfer `license.txt` inside the repo's `leukoquant/licenses/` directory (this folder exists in the repo and is never committed to git):

```
leukoquant/
└── leukoquant/
    └── licenses/
        └── license.txt   ← put your FreeSurfer license here
```

A FreeSurfer licence is free from [surfer.nmr.mgh.harvard.edu](https://surfer.nmr.mgh.harvard.edu/registration.html). No environment variables need to be set - LeukoQuant picks up the licence automatically.

---

## Step 5 - Run

All commands are **identical on macOS, Linux, and Windows** (PowerShell or terminal).

### Full pipeline

```bash
leukoquant process-all \
  --subject    sub-01 \
  --t1         /data/sub-01_T1.nii.gz \
  --flair      /data/sub-01_FLAIR.nii.gz \
  --dwi        /data/sub-01_DWI.nii.gz \
  --bvecs      /data/sub-01.bvec \
  --bvals      /data/sub-01.bval \
  --output-dir /data/output \
  --cores      4
```

**What runs internally:**

1. FreeSurfer `recon-all` (cortical reconstruction)
2. GIF parcellation (brain atlas)
3. BaMoS WMH segmentation
4. TRACULA tractography
5. DTI fitting (FA, MD, AD, RD)
6. NODDI fitting (NDI, ODI, ISOVF)
7. WMH-informed tract metrics

### Batch processing

```bash
leukoquant process-all \
  --subject    /data/subjects.txt \
  --t1         "/data/{subject}_T1.nii.gz" \
  --flair      "/data/{subject}_FLAIR.nii.gz" \
  --dwi        "/data/{subject}_DWI.nii.gz" \
  --output-dir /data/output \
  --cores      8
```

### Individual steps

```bash
leukoquant process-gif   --subject sub-01 --t1 /data/sub-01_T1.nii.gz --output-dir /data/output
leukoquant process-bamos --subject sub-01 --flair /data/FLAIR.nii.gz --t1 /data/T1.nii.gz --gif-results-dir /data/gif --output-dir /data/output
leukoquant process-dti   --subject sub-01 --dwi /data/DWI.nii.gz --bvecs /data/sub.bvec --bvals /data/sub.bval --mask /data/mask.nii.gz --output-dir /data/output
leukoquant process-noddi --subject sub-01 --dwi /data/DWI.nii.gz --bvecs /data/sub.bvec --bvals /data/sub.bval --mask /data/mask.nii.gz --output-dir /data/output
leukoquant process-metrics --subject sub-01 --tractography-path "/data/tracula:/dpath/*/path.pd.nii.gz:dwi" --output-dir /data/metrics
leukoquant process-zscore  --healthy-list /data/healthy.txt --target-list /data/patients.txt --output-dir /data/zscores
```

### HPC (SGE)

```bash
leukoquant process-all ... --scheduler sge --cores 16
```

### Using a YAML config file

Every command accepts `--config-yaml <path>` instead of (or alongside) its
CLI flags - useful for long invocations you want to save and re-run.
**CLI flags always override the equivalent YAML value** when both are given.

```yaml
# process_all_config.yaml
subject: /data/subjects.txt
t1: "/data/{subject}_T1.nii.gz"
flair: "/data/{subject}_FLAIR.nii.gz"
dwi: "/data/{subject}_DWI.nii.gz"
output_dir: /data/output
scheduler: sge
cores: 8
```

```bash
leukoquant process-all --config-yaml process_all_config.yaml

# override just one value for this run; everything else still comes from the YAML
leukoquant process-all --config-yaml process_all_config.yaml --cores 16
```

See each command's page in the [API Reference](api/index.md) for its exact
recognized YAML keys.

---

## All Commands

| Command                                | Description                             |
| -------------------------------------- | --------------------------------------- |
| `process-all`                        | Full pipeline                           |
| `process-recon`                      | FreeSurfer recon-all                    |
| `process-gif`                        | GIF brain parcellation                  |
| `process-bamos`                      | BaMoS WMH segmentation                  |
| `process-dti`                        | DTI fitting (FA, MD, AD, RD)            |
| `process-noddi`                      | NODDI fitting (NDI, ODI, ISOVF)         |
| `process-tracula`                    | TRACULA tractography                    |
| `process-metrics`                    | WMH-informed lesion/tract metrics       |
| `process-tract-qc`                   | Tractography QC                         |
| `process-zscore`                     | Z-score normalisation vs healthy cohort |
| `process-atlas-conversion`           | Atlas format conversion                 |

Use `leukoquant <command> --help` for full options.

---

## Troubleshooting

**`leukoquant: command not found`**
Open a new terminal, or run: `export PATH="$HOME/.venvs/leukoquant/bin:$PATH"`

**`No container runtime found in PATH`**

```bash
bash setup_platform.sh  # then open a new terminal
```

**`FreeSurfer licence not found`**

Place `license.txt` in `leukoquant/licenses/` (the directory exists in the repo):

```
leukoquant/licenses/license.txt
```

**macOS: `limactl: command not found`**

```bash
brew install lima && bash setup_platform.sh
```

**Windows: `leukoquant is not recognised` in PowerShell**
Open a new PowerShell window. If still failing:

```powershell
cat $PROFILE  # should contain a leukoquant line
```

**Linux: no sudo, no conda**

```bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh | bash
# open new terminal, then:
bash setup_platform.sh
```
