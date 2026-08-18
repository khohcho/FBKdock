# FBKdock — Fully Automated End-to-End Molecular Docking Pipeline

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![DOI](https://zenodo.org/badge/1337801934.svg)](https://doi.org/10.5281/zenodo.21987527)

FBKdock is a cross-platform (Windows & Linux) pipeline that automates the
full molecular docking workflow: **protein preparation**, **co-crystallized
ligand redocking validation**, and **batch virtual screening**, using
**Vina-GPU 2.1** (OpenCL) or **AutoDock Vina** (CPU).

---

## What it does — and why

Docking a ligand library is only meaningful if the receptor and the docking
setup are correct. FBKdock therefore runs a controlled experiment before the
screening, and only then docks your compounds:

```
Protein PDB
    │
    ├─ STEP 1: PROTEIN PREPARATION
    │    clean residues → PDB2PQR (AMBER + PROPKA @ pH 7.4)
    │    → ADFRsuite prepare_receptor → PDBQT
    │
    ├─ STEP 2: REDOCKING VALIDATION (the built-in sanity check)
    │    dock the co-crystallized ligand back into its own pocket
    │    → symmetry-aware RMSD vs. the crystal pose
    │    → PASS (< 2.0 A) ⇒ the grid + receptor are trustworthy
    │    → FAIL ⇒ pipeline stops, fix the grid and re-run
    │
    └─ STEP 3: VIRTUAL SCREENING
         charge detection → GFN2-xTB minimization → Meeko PDBQT
         → one docking run over all ligands
```

**Why each step exists:**

| Step | What would go wrong without it |
|---|---|
| 1. Protein preparation | Missing hydrogens, wrong protonation states or formal charges produce meaningless scores. Metals, waters and non-protein residues are handled interactively. |
| 2. Redocking validation | A wrong grid box or a broken receptor silently produces garbage poses. Docking the known ligand back and measuring the RMSD proves the setup works before you spend hours screening. |
| 3. Screening | Each ligand is minimized at the GFN2-xTB level (implicit water) so bond lengths, angles and charges are physically sensible before docking. |

### Strengths

- **End-to-end automation** — no manual PDB → PDBQT conversions.
- **Built-in quality gate** — RMSD validation refuses to screen with a bad setup.
- **ADFRsuite auto-detection** on Windows and Linux (see Installation).
- **Dual engine** — Vina-GPU (OpenCL 1.2 / 3.0) with automatic fallback to
  AutoDock Vina (CPU) when the GPU fails.
- **Isolated Conda environment** (`FBKdock_env`) built from `environment.yml`.
- **Re-runs never overwrite results** — new workspaces get `_2`, `_3`, ... suffixes.

### Limitations

- The **grid box must be defined manually** in the config templates for each
  receptor; there is no automatic pocket detection (the redocking step is
  what tells you whether your box is right).
- **GPU docking requires OpenCL drivers**; without them you are limited to
  CPU AutoDock Vina.
- **ADFRsuite is an external dependency** (not installable via Conda).
- GFN2-xTB minimization is accurate but slow — large libraries take hours.
- Paths must be **ASCII-only** (a Vina-GPU OpenCL compiler limitation).

---

## Requirements

| Requirement | Notes |
|---|---|
| Conda (Miniconda) | [Install instructions](https://docs.conda.io/en/latest/miniconda.html) |
| ADFRsuite 1.0 | External; see Installation |
| GPU + OpenCL drivers | Only needed for Vina-GPU |
| ~3 GB free disk space | Conda environment + engine files |

> **Python:** `environment.yml` pins Python 3.10 — always run from the
> `FBKdock_env` Conda environment.
>
> **Paths:** Vina-GPU's OpenCL kernel compiler rejects non-ASCII characters
> in file paths. Keep the project in a path made of English letters, digits,
> hyphens and underscores only.

---

## Installation (no setup scripts — follow these steps)

### 1. Install Miniconda

**Windows**
1. Download the 64-bit installer from
   [docs.conda.io](https://docs.conda.io/en/latest/miniconda.html).
2. Run it, accept defaults.
3. Open **Anaconda Prompt (Miniconda3)** from the Start Menu.
4. Verify: `conda --version`

**Linux**
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh        # answer "yes" to init
source ~/.bashrc
conda --version
```

### 2. Create the Conda environment

```bash
cd /path/to/FBKdock
conda env create -f environment.yml     # first time only (10–30 min)
conda activate FBKdock_env              # every terminal session
```

Verify the activation worked — your prompt must show `(FBKdock_env)` and
`python --version` must print `Python 3.10.x`.

### 3. Install ADFRsuite (receptor PDBQT conversion)

Download the installer for your platform from the official page:
<https://ccsb.scripps.edu/adfr/downloads/>

- **Windows:** run the installer `.exe` (double-click).
- **Linux:** download the installer app, then
  `chmod a+x ADFRsuite_Linux-x86_64_1.0_install && ./ADFRsuite_Linux-x86_64_1.0_install`.

**You do not need to tell FBKdock where ADFR is.** At startup it
auto-detects `prepare_receptor` in the standard locations:

| Platform | Searched locations |
|---|---|
| Windows | `C:\Program Files*\ADFRsuite*\bin\prepare_receptor.bat`, `C:\ADFRsuite*\...`, `%LOCALAPPDATA%\ADFRsuite*\...`, `PATH` |
| Linux | `/opt/ADFRsuite*/bin/prepare_receptor`, `/usr/local/ADFRsuite*/...`, `~/ADFRsuite*/...`, `PATH` |

If ADFRsuite lives somewhere non-standard, set the environment variable
`FBKDOCK_PREPARE_RECEPTOR` to the full path of `prepare_receptor`
(`prepare_receptor.bat` on Windows).

---

## Quick start (using the bundled example)

1. Copy the example receptor and ligands into the runtime folders:
   ```
   copy examples\protein\2Y9X.pdb Protein\        (Windows)
   cp   examples/protein/2Y9X.pdb  Protein/        (Linux)

   copy examples\ligands\*.sdf    Ligand\         (Windows)
   cp   examples/ligands/*.sdf    Ligand/          (Linux)
   ```
   `Protein/` must contain **exactly one** PDB file.
2. The grid templates (`config_vinagpu.txt` / `config_autodockvina.txt`) are
   already tuned to the 2Y9X binding site — nothing to edit.
3. Launch:
   ```bash
   conda activate FBKdock_env
   python main.py
   ```

Suggested answers for the 2Y9X walkthrough:

| Prompt | Answer |
|---|---|
| Select docking engine | 1 (Vina-GPU) or AutoDock Vina if you have no GPU |
| Select the co-crystallized ligand | the entry `Chain A \| 0TR 410` |
| Is the ligand covalently bound? | n |
| Residues to delete | `394,791-794,1192-1195,1591-1594,396-399,795-800,1196-1999,1595-1599,1736,1873,2010-2012,2149,(B),(C),(D),(F),(G),(H)` |
| Charge for CU | Enter for default `+2.0` |

If the redocking RMSD passes (< 2.0 A), screening of `3a–3e.sdf` starts
automatically. Results land in `OpenCL-1.2_docking/output/` (or the
equivalent engine folder) and a full transcript is in `fbkdock.log`.

For your own target: replace the files in `Protein/` and `Ligand/` and
**re-center the grid box** on your binding site (see `examples/README.md`).

### Command-line options

| Flag | Description |
|---|---|
| `--skip-interactive` | Run without prompts using safe defaults (good for HPC/scripts) |
| `--verbose`, `-v` | Stream xTB / pdb2pqr / docking output live to the terminal |
| `--version` | Print the version and exit |

Without `--verbose` the console stays clean: step headers, prompts, `[OK]` /
`[FAIL]` status and the final RMSD table only. Full detail is always written
to `fbkdock.log`.

---

## Directory structure

```
FBKdock/
├── main.py                         # CLI entry point
├── fbkdock/                        # Python package (pipeline modules)
├── environment.yml                 # Conda environment definition
├── config_vinagpu.txt              # Vina-GPU grid template (user-editable)
├── config_autodockvina.txt         # AutoDock Vina grid template (user-editable)
├── Vina-GPU/                       # Pre-compiled Vina-GPU 2.1 (OpenCL 1.2 / 3.0)
├── AutodockVina/AutodockVina/      # Pre-compiled AutoDock Vina (CPU)
├── Protein/                        # Put YOUR receptor PDB here (exactly one file)
├── Ligand/                         # Put YOUR screening ligands here (.sdf)
├── cocligand/                      # Optional: ideal ligand SDF / covalent PDBQT
├── examples/                       # Example receptor + ligands + grid reference
├── LICENSE                         # GNU GPLv3
└── CITATION.cff                    # Citation metadata
```

**Runtime folders** (auto-created; re-runs get `_2`, `_3`, ... suffixes —
nothing is ever overwritten):

| Path | Purpose |
|---|---|
| `{engine}_docking/` | Virtual-screening workspace (`input/`, `output/`, `logs/`) |
| `{engine}_redock/` | Redocking workspace + `RMSD.txt` |
| `fbkdock.log` | Detailed pipeline log (always written) |

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `'conda' is not recognized` | Restart the terminal or run `conda init`; on Windows use the Anaconda Prompt. |
| `ModuleNotFoundError: No module named 'rdkit'` | You are not inside the environment: `conda activate FBKdock_env`. |
| `ADFRsuite prepare_receptor not found` | Install ADFRsuite (Installation step 3) or set `FBKDOCK_PREPARE_RECEPTOR`. |
| OpenCL kernel compilation fails | The project path contains non-ASCII characters. Move FBKdock to an ASCII-only path. |
| Vina-GPU fails during docking | FBKdock offers an automatic fallback to AutoDock Vina (CPU); your prepared files are preserved. |
| Redocking RMSD fails | The grid box is wrong for this receptor. Re-center `center_x/y/z` and `size_*` in the config templates and re-run. |
| `No files found in Protein/` | Put exactly one PDB file into `Protein/`. |

---

## Releases

The docking engine binaries are included in this repository, so a fresh
clone works out of the box. **ADFRsuite is not redistributed with FBKdock** —
download the installer for your platform from the official page:
<https://ccsb.scripps.edu/adfr/downloads/>. The GitHub release notes point to
the same official download page.

> Note: ADFRsuite has its own license (<https://ccsb.scripps.edu/adfr/license/>)
> and is not covered by the GPLv3 that applies to the FBKdock source code.

---

## Citation

If you use FBKdock in your research, please cite this repository and the
underlying tools:

- **Vina-GPU 2.1** — Tang, S. et al. (2022). *J. Chem. Inf. Model.*
- **AutoDock Vina** — Trott, O. & Olson, A. J. (2010). *J. Comput. Chem.*, 31, 455–461.
- **RDKit** — Landrum, G. et al. RDKit: Open-Source Cheminformatics Software.
- **xTB / GFN2-xTB** — Bannwarth, C. et al. (2019). *WIREs Comput. Mol. Sci.*, 9, e1493.
- **Meeko** — Forli, S. et al. (2016). *Nat. Protoc.*, 11, 905–919.
- **PDB2PQR / PROPKA** — Dolinsky, T. J. et al. (2004). *Nucleic Acids Res.*, 32, W665–W667.
- **ADFRsuite** — Ravindranath, P. A. et al. (2015). *PLoS Comput. Biol.*, 11, e1004586.
- **BioPython** — Cock, P. J. A. et al. (2009). *Bioinformatics*, 25, 1422–1423.

---

## Author & License

**Author:** Fahrettin Buğra Kılıç — [@khohcho](https://github.com/khohcho)

**Contact:** eczfbkilic@gmail.com

FBKdock is free software licensed under the
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0)
(see [LICENSE](LICENSE)). Bundled third-party binaries are documented in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

```
Copyright (C) 2026 Fahrettin Buğra Kılıç
```
