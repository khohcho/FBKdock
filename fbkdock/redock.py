# Copyright (C) 2026 Fahrettin Buğra Kılıç
# SPDX-License-Identifier: GPL-3.0-or-later
"""
FBKdock Redocking Module.

Orchestrates the co-crystallized ligand docking experiment:
  1. Places the prepared ligand PDBQT into the redock input/ directory.
  2. Runs the docking engine via run.bat / run.sh (LOOP=3 internally).
  3. Collects all poses from the output/ directory for RMSD validation.
  4. On GPU failure, offers template-recovery fallback to AutoDock Vina CPU.
"""

import os
import shutil
import subprocess
import logging
from typing import List

from .utils import ensure_dir, run_runner_script, recover_workspace_to_cpu, prompt_yes_no
from .validation import split_pdbqt_models

logger = logging.getLogger("FBKdock.Redock")


def run_redocking(
    redock_dir: str,
    coclig_pdbqt_path: str,
    engine: str = "vinagpu",
    autodockvina_source: str = "",
    project_root: str = "",
) -> List[str]:
    """
    Execute the redocking experiment.

    Parameters
    ----------
    redock_dir : str
        Path to the redock workspace directory.
    coclig_pdbqt_path : str
        Path to the prepared coclig PDBQT file.
    engine : str
        'vinagpu' or 'autodockvina'.
    autodockvina_source : str
        Path to AutoDock Vina source directory (for GPU->CPU fallback).
    project_root : str
        Project root directory (for template_recovery).

    Returns a flat list of MODEL ... ENDMDL blocks.
    """
    input_dir = ensure_dir(os.path.join(redock_dir, "input"))
    output_dir = ensure_dir(os.path.join(redock_dir, "output"))

    input_dest = os.path.join(input_dir, os.path.basename(coclig_pdbqt_path))
    if os.path.abspath(coclig_pdbqt_path) != os.path.abspath(input_dest):
        shutil.copy2(coclig_pdbqt_path, input_dest)
    logger.info("Co-crystallized ligand placed in %s", input_dest)

    try:
        logger.info("Launching redocking via run.bat/run.sh ...")
        run_runner_script(cwd=redock_dir, timeout=10800)
    except (subprocess.CalledProcessError, OSError, FileNotFoundError) as exc:
        logger.warning("Redocking failed: %s", exc)
        if engine == "vinagpu" and autodockvina_source and project_root:
            print(f"\n  WARNING: Vina-GPU kernel compilation or execution failed.")
            if prompt_yes_no("  Fall back to CPU-based AutoDock Vina?"):
                logger.info("Recovering workspace to AutoDock Vina CPU")
                recover_workspace_to_cpu(
                    crashed_workspace=redock_dir,
                    autodockvina_source=autodockvina_source,
                    project_root=project_root,
                )
                run_runner_script(cwd=redock_dir, timeout=10800)
            else:
                raise
        else:
            raise
    logger.info("Redocking complete.")

    pose_blocks = _collect_pose_blocks_from_dir(output_dir)
    logger.info("Redocking - %d total pose blocks collected", len(pose_blocks))
    return pose_blocks


def _collect_pose_blocks_from_dir(directory: str) -> List[str]:
    blocks: List[str] = []
    for root, _dirs, files in os.walk(directory):
        for fname in sorted(files):
            if fname.lower().endswith(".pdbqt"):
                full = os.path.join(root, fname)
                blocks.extend(split_pdbqt_models(full))
    return blocks
