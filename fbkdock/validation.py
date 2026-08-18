# Copyright (C) 2026 Fahrettin Buğra Kılıç
# SPDX-License-Identifier: GPL-3.0-or-later
"""
FBKdock Validation Module.

Symmetry-aware RMSD calculation for redocking validation and the
decision gate that determines whether to proceed with virtual screening.
Exports per-pose RMSD values to RMSD.txt.
"""

import os
import logging
import tempfile
from typing import List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

logger = logging.getLogger("FBKdock.Validation")


def _strip_hydrogens(mol: Chem.Mol) -> Chem.Mol:
    """Return a copy of *mol* with all hydrogen atoms removed."""
    return Chem.RemoveHs(mol, sanitize=True)


def _pdbqt_to_pdb_lines(block: str) -> str:
    """
    Convert ATOM/HETATM records of a PDBQT block to plain PDB lines.

    Only columns 1-54 are kept so RDKit's PDB parser never sees the
    AutoDock atom-type column (e.g. 'A', 'NA', 'OA'), which it would
    misread as chemical elements.
    """
    out: List[str] = []
    for line in block.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            out.append(line[:54])
    out.append("END")
    return "\n".join(out) + "\n"


def _mol_from_pdbqt_block(pdbqt_block: str, sanitize: bool = True) -> Chem.Mol:
    """
    Convert a PDBQT text block to an RDKit Mol.
    """
    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".pdb", prefix="fbkdock_pose_")
    try:
        with os.fdopen(tmp_fd, "w") as fh:
            fh.write(_pdbqt_to_pdb_lines(pdbqt_block))
        mol = Chem.MolFromPDBFile(tmp_path, removeHs=False, sanitize=sanitize)
        if mol is None:
            raise ValueError("RDKit could not parse PDBQT model block")
        return mol
    finally:
        if os.path.isfile(tmp_path):
            os.unlink(tmp_path)


def _pdb_from_pdbqt(pdbqt_path: str) -> str:
    """
    Extract a ligand structure from a PDBQT reference file and save it as
    a temporary PDB for use as an RMSD reference.  Returns the path to
    the temporary PDB.

    Handles both docking output files (MODEL ... ENDMDL blocks) and
    input-style PDBQTs (ROOT/BRANCH/ATOM records without MODEL blocks).

    Used for covalent ligands where the reference is a manually
    prepared PDBQT, not an extracted PDB.
    """
    if not os.path.isfile(pdbqt_path):
        raise FileNotFoundError(f"PDBQT reference not found: {pdbqt_path}")

    models = split_pdbqt_models(pdbqt_path)
    if models:
        block = models[0]
    else:
        # Input-style PDBQT: keep only the atom records.
        with open(pdbqt_path, "r", encoding="utf-8", errors="replace") as fh:
            atom_lines = [
                line.rstrip("\n") for line in fh
                if line.startswith(("ATOM", "HETATM"))
            ]
        if not atom_lines:
            raise ValueError(f"No atom records found in PDBQT: {pdbqt_path}")
        block = "\n".join(atom_lines)

    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".pdb", prefix="fbkdock_ref_")
    with os.fdopen(tmp_fd, "w") as fh:
        fh.write(_pdbqt_to_pdb_lines(block))
    logger.info("PDB reference extracted from PDBQT: %s", tmp_path)
    return tmp_path


def split_pdbqt_models(pdbqt_path: str) -> List[str]:
    """
    Parse a multi-model PDBQT file and return a list of individual
    MODEL ... ENDMDL blocks (as strings).
    """
    if not os.path.isfile(pdbqt_path):
        logger.warning("PDBQT output file not found: %s", pdbqt_path)
        return []

    with open(pdbqt_path, "r", encoding="utf-8", errors="replace") as fh:
        content = fh.read()

    models: List[str] = []
    for part in content.split("MODEL"):
        part = part.strip()
        if not part:
            continue
        if "ENDMDL" not in part:
            continue
        pose_block = "MODEL " + part.split("ENDMDL")[0].strip() + "\nENDMDL\n"
        models.append(pose_block)
    return models


def calculate_rmsd_heavy(
    ref_pdb_path: str,
    pose_blocks: List[str],
) -> List[float]:
    """
    Compute symmetry-aware heavy-atom RMSD between *ref_pdb_path*
    and each docked pose in *pose_blocks*.
    """
    ref_mol = Chem.MolFromPDBFile(ref_pdb_path, removeHs=True, sanitize=True)
    if ref_mol is None:
        raise ValueError(f"Could not read reference ligand from {ref_pdb_path}")
    ref_heavy = _strip_hydrogens(ref_mol)

    rmsd_values: List[float] = []
    for i, block in enumerate(pose_blocks):
        try:
            pose_mol = _mol_from_pdbqt_block(block)
            pose_heavy = _strip_hydrogens(pose_mol)

            if ref_heavy.GetNumAtoms() != pose_heavy.GetNumAtoms():
                logger.debug(
                    "Pose %d: heavy-atom count mismatch (ref=%d, pose=%d) - skipping",
                    i + 1, ref_heavy.GetNumAtoms(), pose_heavy.GetNumAtoms(),
                )
                rmsd_values.append(float("inf"))
                continue

            rmsd = AllChem.GetBestRMS(pose_heavy, ref_heavy)
            rmsd_values.append(rmsd)
            logger.debug("Pose %d: RMSD = %.3f A", i + 1, rmsd)
        except Exception as exc:
            logger.warning("Pose %d: RMSD calculation failed: %s", i + 1, exc)
            rmsd_values.append(float("inf"))

    return rmsd_values


def validate_redocking(
    ref_pdb_path: str,
    pose_blocks: List[str],
    redock_dir: str,
    rmsd_threshold: float = 2.0,
    ref_pdbqt_path: Optional[str] = None,
) -> Tuple[bool, float, int]:
    """
    Decision gate for the redocking experiment.

    If *ref_pdbqt_path* is provided (covalent ligand), converts it to
    a temporary PDB for RMSD reference.  Otherwise uses *ref_pdb_path*.

    Writes per-pose RMSD values to ``{redock_dir}/RMSD.txt`` in a
    single atomic write regardless of pass/fail status.

    Returns (passed, best_rmsd, best_pose_index).
    """
    if not pose_blocks:
        logger.error("No poses available for validation.")
        return False, float("inf"), -1

    effective_ref = ref_pdb_path
    if ref_pdbqt_path and os.path.isfile(ref_pdbqt_path):
        effective_ref = _pdb_from_pdbqt(ref_pdbqt_path)

    logger.info(
        "Validating %d docking poses against reference (threshold = %.1f A)...",
        len(pose_blocks), rmsd_threshold,
    )

    rmsd_list = calculate_rmsd_heavy(effective_ref, pose_blocks)

    finite_rmsds = [(i, v) for i, v in enumerate(rmsd_list) if v != float("inf")]

    if not finite_rmsds:
        best_rmsd = float("inf")
        best_idx = -1
        passed = False
    else:
        best_idx, best_rmsd = min(finite_rmsds, key=lambda x: x[1])
        passed = best_rmsd < rmsd_threshold

    rmsd_path = os.path.join(redock_dir, "RMSD.txt")
    lines: List[str] = []
    lines.append("# Pose  RMSD (A)\n")
    for i, v in enumerate(rmsd_list):
        status = "%.3f" % v if v != float("inf") else "N/A"
        lines.append("%-6d %s\n" % (i + 1, status))
    if not finite_rmsds:
        lines.append("\nNo valid RMSD values could be computed.\n")
    else:
        lines.append(
            "\nBest RMSD: %.3f A  (Pose %d)  -  %s\n"
            % (best_rmsd, best_idx + 1, "PASSED" if passed else "FAILED")
        )
    with open(rmsd_path, "w") as fh:
        fh.writelines(lines)
    logger.info("RMSD values written to %s", rmsd_path)

    logger.info(
        "Best RMSD = %.3f A  (pose %d/%d)  -  %s",
        best_rmsd, best_idx + 1, len(pose_blocks),
        "PASSED" if passed else "FAILED",
    )

    print(f"\n{'='*50}")
    print(f"  REDOCKING VALIDATION RESULTS")
    print(f"{'='*50}")
    print(f"  Total poses evaluated : {len(pose_blocks)}")
    print(f"  Valid RMSD values     : {len(finite_rmsds)}")
    print(f"  Best RMSD             : {best_rmsd:.3f} A")
    print(f"  Threshold             : {rmsd_threshold:.1f} A")
    print(f"  Status                : {'PASSED' if passed else 'FAILED'}")
    print(f"  RMSD file             : {rmsd_path}")
    print(f"{'='*50}")

    return passed, best_rmsd, best_idx + 1
