# Copyright (C) 2026 Fahrettin Buğra Kılıç
# SPDX-License-Identifier: GPL-3.0-or-later
"""
FBKdock Virtual Screening / Docking Module.

Batch-processes all ligands found in the Ligand/ directory:
  - RDKit formal charge
  - GFN2-xTB energy minimisation
  - Meeko PDBQT conversion
  - Docking engine launched via run.bat / run.sh (LOOP=3 internally).
  - On GPU failure, offers template-recovery fallback to AutoDock Vina CPU.
"""

import os
import shutil
import subprocess
import logging
from typing import List

from .utils import (
    ensure_dir, find_sdf_files, run_runner_script,
    recover_workspace_to_cpu, prompt_yes_no,
)
from .ligand import prepare_screening_ligand

logger = logging.getLogger("FBKdock.Docking")


def prepare_all_screening_ligands(
    ligand_dir: str,
    docking_dir: str,
) -> List[str]:
    sdf_files = find_sdf_files(ligand_dir)
    if not sdf_files:
        logger.warning("No ligand files found in %s", ligand_dir)
        return []

    input_dir = ensure_dir(os.path.join(docking_dir, "input"))
    prepared: List[str] = []

    for ligand_path in sdf_files:
        base_name = os.path.splitext(os.path.basename(ligand_path))[0]

        work_dir = os.path.join(docking_dir, ".work_" + base_name)
        ligand_pdbqt = os.path.join(work_dir, f"{base_name}.pdbqt")

        try:
            prepare_screening_ligand(
                ligand_path=ligand_path,
                output_pdbqt=ligand_pdbqt,
                work_dir=work_dir,
            )
            dest = os.path.join(input_dir, os.path.basename(ligand_pdbqt))
            shutil.copy2(ligand_pdbqt, dest)
            prepared.append(base_name)
            print(f"  [OK] {base_name}")
        except Exception as exc:
            logger.error("Ligand %s preparation failed: %s", base_name, exc)
            print(f"  [FAIL] {base_name}: {exc}")

    return prepared


def run_virtual_screening(
    ligand_dir: str,
    docking_dir: str,
    engine: str = "vinagpu",
    skip_interactive: bool = False,
    autodockvina_source: str = "",
    project_root: str = "",
) -> None:
    sdf_files = find_sdf_files(ligand_dir)
    if not sdf_files:
        logger.warning("No ligand files found in %s", ligand_dir)
        print(f"\nWARNING: No files found in {ligand_dir}. Screening aborted.")
        return

    print(f"\n{'='*60}")
    print(f"  VIRTUAL SCREENING - {len(sdf_files)} ligand(s) detected")
    print(f"{'='*60}")
    for i, p in enumerate(sdf_files, 1):
        print(f"  [{i}] {os.path.basename(p)}")
    print(f"{'='*60}")

    if not skip_interactive:
        if not prompt_yes_no("\nProceed with virtual screening?"):
            print("Screening cancelled.")
            return

    print("\n--- Phase A: Preparing all ligands ---")
    prepared = prepare_all_screening_ligands(
        ligand_dir=ligand_dir,
        docking_dir=docking_dir,
    )

    if not prepared:
        print("ERROR: No ligands could be prepared. Aborting.")
        return

    print(f"\n  {len(prepared)}/{len(sdf_files)} ligand(s) prepared successfully.")
    print("  All PDBQTs are in {0}/input/".format(docking_dir))

    print("\n--- Phase B: Launching docking engine ---")
    try:
        run_runner_script(cwd=docking_dir, timeout=14400)
    except (subprocess.CalledProcessError, OSError, FileNotFoundError) as exc:
        logger.warning("Docking engine failed: %s", exc)
        if engine == "vinagpu" and autodockvina_source and project_root:
            print(f"\n  WARNING: Vina-GPU failed.")
            if not skip_interactive:
                if prompt_yes_no("  Fall back to CPU-based AutoDock Vina?"):
                    logger.info("Recovering workspace to AutoDock Vina CPU")
                    recover_workspace_to_cpu(
                        crashed_workspace=docking_dir,
                        autodockvina_source=autodockvina_source,
                        project_root=project_root,
                    )
                    run_runner_script(cwd=docking_dir, timeout=14400)
                else:
                    raise
            else:
                raise
        else:
            raise

    print(f"\n{'='*50}")
    print(f"  SCREENING COMPLETE")
    print(f"  Output:   {docking_dir}/output/")
    print(f"  Logs:     {docking_dir}/logs/")
    print(f"{'='*50}")
