# Copyright (C) 2026 Fahrettin Buğra Kılıç
# SPDX-License-Identifier: GPL-3.0-or-later
"""
FBKdock Ligand Preparation Module.

Handles co-crystallized ligand extraction, RDKit-based alignment,
and batch ligand preparation for virtual screening (charge + xTB + Meeko).
"""

import os
import sys
import time
import shutil
import logging
import tempfile
import subprocess
import sysconfig
from typing import List, Optional, Tuple

from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
from rdkit import Chem
from rdkit.Chem import AllChem

from .utils import (
    OS_TYPE, run_subprocess, ensure_dir,
)

logger = logging.getLogger("FBKdock.Ligand")

WATER_RESNAMES = {"HOH", "WAT", "H2O", "TIP", "TIP3", "TIP4", "TIP5", "SPC", "SOL"}
ION_RESNAMES = {
    "NA", "K", "CL", "CA", "MG", "ZN", "FE", "MN", "CU", "CO",
    "CD", "HG", "NI", "SO4", "PO4", "NO3", "ACT", "EDO", "GOL",
}


def list_ligand_residues(pdb_path: str) -> List[Tuple[str, str, int, str]]:
    """
    Scan *pdb_path* for HETATM residues that are NOT waters or common ions.
    Returns a list of (chain_id, resname, resid, insertion_code).
    """
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    candidates: List[Tuple[str, str, int, str]] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                hetfield, resid, ins_code = residue.get_id()
                if hetfield == " ":
                    continue
                resname = residue.get_resname().strip().upper()
                if resname in WATER_RESNAMES or resname in ION_RESNAMES:
                    continue
                candidates.append((chain.get_id(), resname, resid, ins_code.strip()))
    return candidates


def extract_ligand_residue(
    pdb_path: str,
    chain_id: str,
    resname: str,
    resid: int,
    insertion: str,
    output_path: str,
) -> str:
    """
    Extract a single HETATM residue from *pdb_path* and write it to *output_path* as PDB.
    Matching is performed on (chain_id, resname, resid, insertion).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    class LigandSelector(Select):
        def accept_residue(self, residue):
            try:
                cid = residue.get_parent().get_id()
            except AttributeError:
                cid = ""
            _, r_id, ins = residue.get_id()
            r_name = residue.get_resname().strip()
            return (cid == chain_id and r_name.upper() == resname.upper()
                    and r_id == resid and ins.strip() == insertion.strip())

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path, LigandSelector())
    if not os.path.isfile(output_path) or os.path.getsize(output_path) == 0:
        raise RuntimeError(
            f"Failed to extract ligand residue "
            f"(chain={chain_id}, {resname} {resid}:{insertion}) to {output_path}"
        )
    logger.info("Extracted co-crystallized ligand to %s", output_path)
    return output_path


def add_hydrogens_rdkit(input_path: str, output_path: str) -> str:
    """
    Read a molecule from *input_path* (PDB, SDF), add hydrogens while
    preserving heavy-atom 3D coordinates, write as SDF.
    """
    fmt = os.path.splitext(input_path)[1].lower()
    if fmt == ".pdb":
        mol = Chem.MolFromPDBFile(input_path, removeHs=False, sanitize=True)
    elif fmt in (".mol", ".sdf"):
        suppl = Chem.SDMolSupplier(input_path, removeHs=False, sanitize=True)
        mol = next(suppl) if suppl else None
    else:
        raise ValueError(f"Unsupported input format for RDKit: {fmt}")

    if mol is None:
        raise ValueError(f"RDKit could not read molecule from {input_path}")

    mol = Chem.AddHs(mol, addCoords=True)

    out_fmt = os.path.splitext(output_path)[1].lower()
    if out_fmt == ".pdb":
        w = Chem.PDBWriter(output_path)
        w.write(mol)
        w.close()
    elif out_fmt in (".mol", ".sdf"):
        w = Chem.SDWriter(output_path)
        w.write(mol)
        w.close()
    else:
        raise ValueError(f"Unsupported output format: {out_fmt}")

    logger.info("Hydrogens added (%d atoms total) -> %s", mol.GetNumAtoms(), output_path)
    return output_path


def align_coclig_sdf(
    sdf_path: str,
    ref_pdb_path: str,
    output_path: str,
) -> str:
    """
    Produce a ligand SDF whose bond orders come from the ideal SDF template
    and whose 3D coordinates come from the reference (co-crystallized) PDB.

    Hydrogens are stripped from both molecules before matching: the ideal
    SDF often carries explicit hydrogens while the PDB-extracted residue
    does not, which makes ``AssignBondOrdersFromTemplate`` report
    "No matching found".  After the bond orders are assigned, hydrogens
    are added back explicitly (Meeko requires them) with 3D coordinates.
    """
    template_mol = Chem.SDMolSupplier(sdf_path, removeHs=True, sanitize=True)
    template_mol = next(template_mol) if template_mol else None
    if template_mol is None:
        raise ValueError(f"Could not read SDF template: {sdf_path}")

    ref_mol = Chem.MolFromPDBFile(ref_pdb_path, removeHs=True, sanitize=True)
    if ref_mol is None:
        raise ValueError(f"Could not read reference PDB: {ref_pdb_path}")

    logger.info("Assigning bond orders from SDF template (%d atoms) -> PDB (%d atoms) ...",
                template_mol.GetNumAtoms(), ref_mol.GetNumAtoms())

    try:
        matched_mol = AllChem.AssignBondOrdersFromTemplate(template_mol, ref_mol)
    except Exception as exc:
        logger.error("AssignBondOrdersFromTemplate failed: %s", exc)
        raise RuntimeError(
            f"Could not assign bond orders from template to reference: {exc}"
        )

    conf = ref_mol.GetConformer()
    matched_mol.RemoveAllConformers()
    matched_mol.AddConformer(conf, assignId=True)

    # Meeko requires explicit hydrogens; add them on the now-correct bond
    # orders (heavy-atom coordinates come from the reference conformer).
    matched_mol = Chem.AddHs(matched_mol, addCoords=True)

    logger.info("Template bonds + PDB coords applied; RMSD = 0.0 A (exact).")

    w = Chem.SDWriter(output_path)
    w.write(matched_mol)
    w.close()
    logger.info("Aligned coclig SDF written to %s", output_path)
    return output_path


_xtb_version_logged = False


def _log_xtb_version_once(xtb_exe: str) -> None:
    global _xtb_version_logged
    if _xtb_version_logged:
        return
    _xtb_version_logged = True
    try:
        proc = subprocess.run(
            [xtb_exe, "--version"], capture_output=True, text=True, timeout=15
        )
        raw = (proc.stdout or proc.stderr or "").strip().splitlines()
        logger.info("xTB binary: %s (version: %s)", xtb_exe, raw[0] if raw else "?")
    except Exception:
        logger.info("xTB binary: %s", xtb_exe)


def run_xtb_optimization(
    input_path: str,
    output_dir: str,
    gfn: int = 2,
    net_charge: int = 0,
    alpb: bool = True,
) -> str:
    """
    Perform GFN2-xTB geometry optimisation on *input_path*.

    Parameters
    ----------
    input_path : str
        Input file (SDF, PDB, XYZ - xTB auto-detects the format).
    output_dir : str
        Working directory for xTB output.
    gfn : int
        GFN level (1 or 2). Default 2.
    net_charge : int
        Net molecular charge passed to xTB via ``--chrg``.
    alpb : bool
        If True, use the ALPB implicit solvation model (water).

    Returns the path to the optimised output file.
    """
    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"xTB input not found: {input_path}")
    ensure_dir(output_dir)

    xtb_exe = shutil.which("xtb")
    if xtb_exe is None:
        raise FileNotFoundError(
            "xTB executable not found on PATH. Install it with: "
            "conda install -n FBKdock_env -c conda-forge xtb"
        )

    cmd = ["xtb", input_path, "--gfn", str(gfn), "--opt"]
    if alpb:
        cmd.append("--alpb")
        cmd.append("water")
    if net_charge != 0:
        cmd.extend(["--chrg", str(net_charge)])

    _log_xtb_version_once(xtb_exe)
    threads = max(1, os.cpu_count() or 1)
    logger.info(
        "Running GFN%d-xTB optimisation (charge=%d, alpb=%s, OMP_NUM_THREADS=%d) on %s ...",
        gfn, net_charge, alpb, threads, input_path,
    )
    xtb_env = {"OMP_NUM_THREADS": str(threads)}
    if OS_TYPE == "linux":
        xtb_env["OPENBLAS_NUM_THREADS"] = "1"
    t0 = time.monotonic()
    run_subprocess(
        cmd, cwd=output_dir, timeout=3600,
        env=xtb_env,
    )
    dt = time.monotonic() - t0
    logger.info("xTB optimisation finished in %.1f s (%s)", dt, xtb_exe)

    ext = os.path.splitext(input_path)[1].lower()
    if ext == ".sdf":
        optimized = os.path.join(output_dir, "xtbopt.sdf")
    else:
        optimized = os.path.join(output_dir, "xtbopt.xyz")
    if not os.path.isfile(optimized):
        raise RuntimeError(f"xTB did not produce output in {output_dir}")
    logger.info("xTB optimization complete: %s", optimized)
    return optimized


def run_meeko(input_path: str, output_pdbqt: str) -> str:
    """
    Convert a ligand file (MOL2, SDF, PDB) to PDBQT using Meeko.
    """
    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"Meeko input not found: {input_path}")
    ensure_dir(os.path.dirname(output_pdbqt) or ".")

    try:
        ligand_script = os.path.join(
            sysconfig.get_paths()["scripts"], "mk_prepare_ligand.py"
        )
        if not os.path.isfile(ligand_script):
            ligand_script = os.path.join(sys.prefix, "Scripts", "mk_prepare_ligand.py")
        run_subprocess(
            [sys.executable, ligand_script, "-i", input_path, "-o", output_pdbqt],
            timeout=300,
        )
    except Exception:
        logger.warning("mk_prepare_ligand.py failed; trying Python API fallback ...")
        try:
            from meeko import MoleculePreparation, PDBQTWriterLegacy
            ext = os.path.splitext(input_path)[1].lower()
            if ext in (".mol2",):
                mol = Chem.MolFromMol2File(input_path, removeHs=False, sanitize=True)
            elif ext in (".sdf", ".mol"):
                suppl = Chem.SDMolSupplier(input_path, removeHs=False, sanitize=True)
                mol = next(suppl) if suppl else None
            elif ext == ".pdb":
                mol = Chem.MolFromPDBFile(input_path, removeHs=False, sanitize=True)
            else:
                mol = Chem.MolFromMol2File(input_path, removeHs=False, sanitize=True)
            if mol is None:
                raise RuntimeError("Meeko fallback: cannot read input file")
            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(mol)
            pdbqt_string = ""
            for setup in mol_setups:
                pdbqt_str, is_ok, err_msg = PDBQTWriterLegacy.write_string(setup)
                if not is_ok or not pdbqt_str.strip():
                    raise RuntimeError(
                        f"Meeko fallback: PDBQT writing failed: {err_msg}"
                    )
                pdbqt_string += pdbqt_str + "\n"
            ensure_dir(os.path.dirname(output_pdbqt) or ".")
            with open(output_pdbqt, "w") as fh:
                fh.write(pdbqt_string)
        except Exception as fallback_err:
            logger.error("Meeko fallback also failed: %s", fallback_err)
            raise RuntimeError(
                f"Meeko failed for both CLI and Python API approaches: {fallback_err}"
            )

    if not os.path.isfile(output_pdbqt) or os.path.getsize(output_pdbqt) == 0:
        raise RuntimeError(f"Meeko did not produce output: {output_pdbqt}")
    logger.info("Meeko PDBQT written to %s", output_pdbqt)
    return output_pdbqt


def _mol_is_2d(mol: Chem.Mol) -> bool:
    """Return True when the molecule is flat (all z equal) or has no conformer.

    PubChem and drawing-tool exports are 2D (z = 0 for every atom); xTB
    cannot build a sane 3D structure from such a geometry and its SCF
    diverges.
    """
    try:
        conf = mol.GetConformer()
    except ValueError:
        return True
    zs = [conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())]
    return abs(max(zs) - min(zs)) < 1e-3


def _embed_3d(mol: Chem.Mol) -> Chem.Mol:
    """Generate a 3D conformer for a flat/2D molecule.

    Hydrogens are added first so the embedding produces a complete,
    protonated 3D structure.  ETKDG embedding (deterministic random seeds,
    random-coordinate fallback) is followed by MMFF94 minimisation with a
    UFF fallback, which covers elements MMFF94 lacks (P, Si, metals, ...).
    """
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xA710
    cid = -1
    for attempt in range(5):
        params.randomSeed = 0xA710 + attempt
        cid = AllChem.EmbedMolecule(mol, params)
        if cid == 0:
            break
    if cid != 0:
        params.useRandomCoords = True
        cid = AllChem.EmbedMolecule(mol, params)
    if cid != 0:
        raise ValueError("RDKit could not generate a 3D conformer")
    if AllChem.MMFFOptimizeMolecule(mol) != 0:
        AllChem.UFFOptimizeMolecule(mol)
    return mol


def prepare_screening_ligand(
    ligand_path: str,
    output_pdbqt: str,
    work_dir: Optional[str] = None,
) -> str:
    """
    Full ligand-preparation pipeline for a single virtual-screening ligand:

    1. Read SDF with RDKit, compute net formal charge.
    2. If the input is 2D (flat, e.g. PubChem export), generate a 3D
       conformer (ETKDG embedding + MMFF/UFF minimisation).
    3. GFN2-xTB energy minimisation (--opt --alpb water --chrg <charge>).
    4. Complete any missing explicit hydrogens on the optimized frame.
    5. Meeko: xtbopt_H.sdf -> PDBQT (direct, no RDKit re-read).

    Parameters
    ----------
    ligand_path : str
        Path to input ligand SDF file.
    output_pdbqt : str
        Desired output PDBQT path.
    work_dir : str or None
        Temporary working directory (auto-created if None).
    """
    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix="fbkdock_ligand_")
    ensure_dir(work_dir)

    base = os.path.splitext(os.path.basename(ligand_path))[0]
    logger.info("=== Preparing screening ligand: %s ===", base)

    suppl = Chem.SDMolSupplier(ligand_path, removeHs=False, sanitize=True)
    mol = next(suppl) if suppl else None
    if mol is None:
        raise ValueError(f"RDKit could not read SDF: {ligand_path}")

    charge = Chem.GetFormalCharge(mol)
    logger.info("Net formal charge (RDKit): %d", charge)

    xtb_input = ligand_path
    if _mol_is_2d(mol):
        logger.info("Input SDF is 2D (flat); generating a 3D conformer ...")
        mol3d = _embed_3d(mol)
        xtb_input = os.path.join(work_dir, f"{base}_3d.sdf")
        with Chem.SDWriter(xtb_input) as writer:
            writer.write(mol3d)

    xtb_opt_path = run_xtb_optimization(
        xtb_input, work_dir, gfn=2, net_charge=charge, alpb=True,
    )

    # xTB writes back exactly the atom set of the input SDF, and screening
    # SDFs are often only partially protonated.  Meeko requires fully
    # explicit hydrogens, so complete any missing ones on the optimized
    # frame.  Idempotent: molecules that already carry all their hydrogens
    # are passed through unchanged.
    suppl = Chem.SDMolSupplier(xtb_opt_path, removeHs=False, sanitize=True)
    opt_mol = next(suppl) if suppl else None
    if opt_mol is None:
        raise ValueError(f"RDKit could not read xTB output: {xtb_opt_path}")
    opt_h = Chem.AddHs(opt_mol, addCoords=True)
    xtb_h_path = os.path.join(work_dir, "xtbopt_H.sdf")
    with Chem.SDWriter(xtb_h_path) as writer:
        writer.write(opt_h)

    run_meeko(xtb_h_path, output_pdbqt)

    logger.info("=== Ligand %s prepared -> %s ===", base, output_pdbqt)
    return output_pdbqt
