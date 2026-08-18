#!/usr/bin/env python3
# Copyright (C) 2026 Fahrettin Buğra Kılıç
# SPDX-License-Identifier: GPL-3.0-or-later
"""
FBKdock - Fully Automated End-to-End Molecular Docking Pipeline
================================================================

Orchestrates protein preparation, co-crystallized ligand redocking
validation, and batch virtual screening using Vina-GPU 2.1 or
AutoDock Vina (CPU).
"""

import os
import sys
import shutil
import glob as _g
import logging
import argparse
from typing import Optional, Set

from fbkdock import __version__
from fbkdock import utils
from fbkdock.utils import (
    OS_TYPE, ensure_dir, copy_directory_contents,
    ensure_executable,
    parse_config, prompt_choice, prompt_yes_no,
    find_fbkdock_files, find_pdbqt_files,
)
from fbkdock.protein import prepare_protein, analyze_pdb_structure
from fbkdock.ligand import (
    list_ligand_residues, extract_ligand_residue,
    add_hydrogens_rdkit, run_meeko, align_coclig_sdf,
)
from fbkdock.redock import run_redocking
from fbkdock.validation import validate_redocking
from fbkdock.docking import run_virtual_screening

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
VINA_GPU_DIR = os.path.join(PROJECT_ROOT, "Vina-GPU")
AUTODOCKVINA_DIR = os.path.join(PROJECT_ROOT, "AutodockVina", "AutodockVina")
PROTEIN_DIR = os.path.join(PROJECT_ROOT, "Protein")
LIGAND_DIR = os.path.join(PROJECT_ROOT, "Ligand")
COCLIGAND_DIR = os.path.join(PROJECT_ROOT, "cocligand")
CONFIG_VINAGPU = os.path.join(PROJECT_ROOT, "config_vinagpu.txt")
CONFIG_AUTODOCKVINA = os.path.join(PROJECT_ROOT, "config_autodockvina.txt")


def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    fmt = logging.Formatter(
        "[%(asctime)s] %(levelname)-8s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh = logging.FileHandler("fbkdock.log", mode="w")
    fh.setLevel(level)
    fh.setFormatter(fmt)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.WARNING if not verbose else logging.DEBUG)
    ch.setFormatter(fmt)
    root_logger = logging.getLogger("FBKdock")
    root_logger.setLevel(level)
    root_logger.addHandler(fh)
    root_logger.addHandler(ch)


def _setup_console_encoding() -> None:
    """Force UTF-8 output with lossy fallback so any locale can print."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    try:
        sys.stdin.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


def select_engine() -> tuple:
    options = []
    mappings = []

    for ver in ("OpenCL-1.2", "OpenCL-3"):
        path = os.path.join(VINA_GPU_DIR, ver)
        if os.path.isdir(path):
            options.append(f"Vina-GPU  ({ver})")
            mappings.append(("vinagpu", path, CONFIG_VINAGPU, ver))

    if os.path.isdir(AUTODOCKVINA_DIR):
        options.append("AutoDock Vina  (CPU)")
        mappings.append(("autodockvina", AUTODOCKVINA_DIR, CONFIG_AUTODOCKVINA, "AutodockVina"))

    if not options:
        print("ERROR: No docking engine found.")
        sys.exit(1)

    if len(options) == 1:
        print(f"Only one engine available: {options[0]}")
        return mappings[0]

    idx = prompt_choice("Select docking engine:", options)
    return mappings[idx - 1]


def _next_free_workspace_pair(prefix: str) -> tuple:
    """Return (docking_dir, redock_dir) names that do not exist yet.

    First run uses the base names; later runs get '_2', '_3', ... so
    earlier results are never deleted or overwritten.
    """
    base_d = os.path.join(PROJECT_ROOT, f"{prefix}_docking")
    base_r = os.path.join(PROJECT_ROOT, f"{prefix}_redock")
    if not (os.path.isdir(base_d) or os.path.isdir(base_r)):
        return base_d, base_r
    k = 2
    while (
        os.path.isdir(os.path.join(PROJECT_ROOT, f"{prefix}_docking_{k}"))
        or os.path.isdir(os.path.join(PROJECT_ROOT, f"{prefix}_redock_{k}"))
    ):
        k += 1
    return (
        os.path.join(PROJECT_ROOT, f"{prefix}_docking_{k}"),
        os.path.join(PROJECT_ROOT, f"{prefix}_redock_{k}"),
    )


def setup_workspace(engine_type: str, source_dir: str, config_template: str, prefix: str) -> tuple:
    docking_dir, redock_dir = _next_free_workspace_pair(prefix)

    for ddir in (docking_dir, redock_dir):
        copy_directory_contents(source_dir, ddir)
        ensure_executable(ddir)
        shutil.copy2(config_template, os.path.join(ddir, "config.txt"))
        print(f"  Workspace created: {ddir}")

    return docking_dir, redock_dir


def run_protein_prep(
    docking_dir: str,
    redock_dir: str,
    pdb_path: str,
    config_template: str,
    skip_interactive: bool = False,
    auto_remove_indices: Optional[Set[int]] = None,
) -> tuple:
    print("\n" + "=" * 60)
    print("  STEP 1: PROTEIN PREPARATION")
    print("=" * 60)

    receptor_pdbqt = prepare_protein(
        pdb_path=pdb_path,
        docking_dir=docking_dir,
        redock_dir=redock_dir,
        config_template_path=config_template,
        skip_interactive=skip_interactive,
        auto_remove_indices=auto_remove_indices,
    )

    cleaned_pdb = os.path.join(
        os.path.dirname(pdb_path),
        os.path.splitext(os.path.basename(pdb_path))[0] + "_cleaned.pdb",
    )
    return receptor_pdbqt, cleaned_pdb


def run_redock_pipeline(
    redock_dir: str,
    coclig_pdb_path: str,
    coclig_sdf_path: str,
    engine: str = "vinagpu",
    ref_pdbqt_path: Optional[str] = None,
    autodockvina_source: str = "",
    project_root: str = "",
) -> tuple:
    print("\n" + "=" * 60)
    print("  STEP 2: CO-CRYSTALLIZED LIGAND REDOCKING")
    print("=" * 60)

    coclig_work_dir = ensure_dir(os.path.join(redock_dir, "cocligand"))
    output_pdbqt = os.path.join(coclig_work_dir, "coclig.pdbqt")

    if ref_pdbqt_path and os.path.isfile(ref_pdbqt_path):
        dest = os.path.join(redock_dir, "input", os.path.basename(ref_pdbqt_path))
        ensure_dir(os.path.dirname(dest))
        shutil.copy2(ref_pdbqt_path, dest)
        print(f"  Covalent ligand PDBQT placed: {dest}")
    elif coclig_sdf_path and os.path.isfile(coclig_sdf_path):
        aligned_sdf = os.path.join(coclig_work_dir, "coclig_aligned.sdf")
        align_coclig_sdf(coclig_sdf_path, coclig_pdb_path, aligned_sdf)
        run_meeko(aligned_sdf, output_pdbqt)
        print(f"  Coclig SDF aligned and prepared: {output_pdbqt}")
    elif coclig_pdb_path and os.path.isfile(coclig_pdb_path):
        coclig_h_sdf = os.path.join(coclig_work_dir, "coclig_H.sdf")
        add_hydrogens_rdkit(coclig_pdb_path, coclig_h_sdf)
        run_meeko(coclig_h_sdf, output_pdbqt)
        print(f"  Coclig prepared from PDB: {output_pdbqt}")
    else:
        print(f"  WARNING: No co-crystallized ligand to prepare.")

    if ref_pdbqt_path:
        print(f"  Reference PDBQT: {ref_pdbqt_path}")
    elif coclig_pdb_path:
        print(f"  Reference PDB: {coclig_pdb_path}")

    effective_ligand = output_pdbqt
    if ref_pdbqt_path and os.path.isfile(ref_pdbqt_path):
        effective_ligand = os.path.join(redock_dir, "input", os.path.basename(ref_pdbqt_path))

    print("\nLaunching redocking (3 runs via run.bat/run.sh) ...")
    pose_blocks = run_redocking(
        redock_dir=redock_dir,
        coclig_pdbqt_path=effective_ligand,
        engine=engine,
        autodockvina_source=autodockvina_source,
        project_root=project_root,
    )
    print(f"  Collected {len(pose_blocks)} pose blocks from redocking runs.")

    print("\nValidating RMSD...")
    passed, best_rmsd, best_idx = validate_redocking(
        ref_pdb_path=coclig_pdb_path if coclig_pdb_path and os.path.isfile(coclig_pdb_path) else "",
        pose_blocks=pose_blocks,
        redock_dir=redock_dir,
        rmsd_threshold=2.0,
        ref_pdbqt_path=ref_pdbqt_path if ref_pdbqt_path else None,
    )

    return passed, best_rmsd, pose_blocks


def main() -> None:
    _setup_console_encoding()
    parser = argparse.ArgumentParser(
        description="FBKdock - Automated molecular docking pipeline",
    )
    parser.add_argument(
        "--skip-interactive", action="store_true",
        help="Run without interactive prompts.",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Stream subprocess output (xTB, pdb2pqr, docking) live to the terminal.",
    )
    parser.add_argument(
        "--version", action="version", version=f"FBKdock {__version__}",
    )
    args = parser.parse_args()

    utils.LIVE_OUTPUT = args.verbose
    setup_logging(verbose=args.verbose)
    logger = logging.getLogger("FBKdock")

    print("\n" + "=" * 60)
    print("  FBKdock v{0}".format(__version__))
    print("  Automated Molecular Docking Pipeline")
    print("  OS: {0}".format(OS_TYPE))
    print("=" * 60)

    # --- Engine selection ---
    engine_type, source_dir, config_template, prefix = select_engine()
    print(f"\nEngine: {engine_type}")

    # --- Config validation ---
    if not os.path.isfile(config_template):
        print(f"ERROR: Config template not found: {config_template}")
        sys.exit(1)
    cfg = parse_config(config_template)
    for key in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z"):
        if key in cfg:
            print(f"  {key} = {cfg[key]}")

    # --- Workspace setup ---
    print("\nSetting up docking workspaces...")
    docking_dir, redock_dir = setup_workspace(engine_type, source_dir, config_template, prefix)

    # --- Locate protein PDB (any extension) ---
    pdb_files = find_fbkdock_files(PROTEIN_DIR)
    if pdb_files:
        pdb_path = pdb_files[0]
        print(f"Using protein file: {pdb_path}")
    else:
        print(f"ERROR: No files found in {PROTEIN_DIR}/")
        sys.exit(1)

    _, descriptors = analyze_pdb_structure(pdb_path)

    # --- Co-crystallized ligand selection ---
    coclig_info: Optional[tuple] = None
    auto_remove_indices: Set[int] = set()
    coclig_sdf_path: Optional[str] = None
    is_covalent: bool = False
    ref_pdbqt_path: Optional[str] = None

    candidates = list_ligand_residues(pdb_path)
    if candidates:
        print("\n" + "=" * 60)
        print("  CO-CRYSTALLIZED LIGAND DETECTION")
        print("=" * 60)
        print("\nDetected hetero-ligand candidates in the original protein PDB:")
        for idx, (chain, resn, resid, ins) in enumerate(candidates, 1):
            ins_str = f":{ins}" if ins else ""
            print(f"  [{idx}] Chain {chain} | {resn} {resid}{ins_str}")

        if args.skip_interactive:
            chosen = candidates[0]
            print(f"Auto-selected first candidate: {chosen}")
        else:
            idx = prompt_choice("Select the co-crystallized ligand (or 0 to skip redocking):", [
                f"Chain {c[0]} | {c[1]} {c[2]}{':' + c[3] if c[3] else ''}"
                for c in candidates
            ])
            chosen = candidates[idx - 1]

        chain_id, resname, resid, insertion = chosen

        coclig_desc_index = None
        for di, d in enumerate(descriptors):
            if (d.chain_id == chain_id and d.resname.upper() == resname.upper()
                    and d.resid == resid and d.insertion == insertion):
                coclig_desc_index = di
                break

        if coclig_desc_index is not None:
            auto_remove_indices.add(coclig_desc_index)
            coclig_info = (chain_id, resname, resid, insertion, coclig_desc_index)

            if args.skip_interactive:
                is_covalent = False
            else:
                is_covalent = prompt_yes_no(
                    f"Is the co-crystallized ligand ({resname}) covalently bound?"
                )

            if is_covalent:
                pdbqt_files = find_pdbqt_files(COCLIGAND_DIR)
                if pdbqt_files:
                    ref_pdbqt_path = pdbqt_files[0]
                    print(f"\n  Covalent ligand PDBQT found: {ref_pdbqt_path}")
                else:
                    print(f"\n  WARNING: No PDBQT file found in {COCLIGAND_DIR}/")
                    print(f"  Place your manually prepared PDBQT there and re-run.")
                    if not args.skip_interactive:
                        if not prompt_yes_no("Continue to virtual screening without redocking?"):
                            sys.exit(0)
            else:
                cocligand_dir = ensure_dir(os.path.join(redock_dir, "cocligand"))
                coclig_pdb_path = os.path.join(cocligand_dir, "coclig.pdb")
                extract_ligand_residue(pdb_path, chain_id, resname, resid, insertion, coclig_pdb_path)
                print(f"\n  Co-crystallized ligand extracted to: {coclig_pdb_path}")

                sdf_files_root = find_fbkdock_files(COCLIGAND_DIR)
                for sf in sdf_files_root:
                    ext = os.path.splitext(sf)[1].lower()
                    if ext in (".sdf", ".mol"):
                        coclig_sdf_path = sf
                        print(f"  Ideal SDF found: {coclig_sdf_path} (will be aligned)")
                        break

            print(f"  This ligand will be automatically removed from the receptor.")
        else:
            print("\nWARNING: Could not locate the selected ligand in the residue list.")
    else:
        print("\nNo non-water, non-ion HETATM ligands detected in the protein PDB.")
        print("Redocking validation will be skipped.")

    # --- STEP 1: Protein preparation ---
    receptor_pdbqt_path, cleaned_pdb_path = run_protein_prep(
        docking_dir, redock_dir, pdb_path, config_template,
        skip_interactive=args.skip_interactive,
        auto_remove_indices=auto_remove_indices,
    )

    print("\n" + "=" * 60)
    print("  PROTEIN PDBQT READY")
    print(f"  {receptor_pdbqt_path}")
    print("=" * 60)

    # --- STEP 2: Redocking validation (UNIFIED - covalent + non-bonded) ---
    has_coclig = coclig_info is not None

    if has_coclig:
        coclig_pdb_path = os.path.join(redock_dir, "cocligand", "coclig.pdb")

        passed, best_rmsd, pose_blocks = run_redock_pipeline(
            redock_dir=redock_dir,
            coclig_pdb_path=coclig_pdb_path if not is_covalent else "",
            coclig_sdf_path=coclig_sdf_path if not is_covalent else "",
            engine=engine_type,
            ref_pdbqt_path=ref_pdbqt_path if is_covalent else None,
            autodockvina_source=AUTODOCKVINA_DIR,
            project_root=PROJECT_ROOT,
        )

        if not passed:
            print("\n" + "!" * 60)
            print("  RMSD VALIDATION FAILED")
            print(f"  Best RMSD: {best_rmsd:.3f} A  (threshold: 2.0 A)")
            print("  Please re-evaluate the grid center/size in the config")
            print(f"  template ({config_template}) and re-run.")
            print("!" * 60)
            sys.exit(1)
        else:
            print(f"\nRMSD validation PASSED (best = {best_rmsd:.3f} A).")
            print("Proceeding to virtual screening...")
    else:
        print("\nWARNING: Skipping redocking validation.")
        if not args.skip_interactive:
            if not prompt_yes_no("Continue directly to virtual screening?"):
                sys.exit(0)

    # --- STEP 3: Virtual screening ---
    print("\n" + "=" * 60)
    print("  STEP 3: VIRTUAL SCREENING")
    print("=" * 60)

    run_virtual_screening(
        ligand_dir=LIGAND_DIR,
        docking_dir=docking_dir,
        engine=engine_type,
        skip_interactive=args.skip_interactive,
        autodockvina_source=AUTODOCKVINA_DIR,
        project_root=PROJECT_ROOT,
    )

    print("\n" + "=" * 60)
    print("  FBKdock PIPELINE COMPLETE")
    print(f"  Docking output: {docking_dir}/output/")
    print(f"  Docking logs:   {docking_dir}/logs/")
    print(f"  Pipeline log:   fbkdock.log")
    print("  NOTE: Re-running keeps earlier workspaces; new runs get _2, _3, ...")
    print("=" * 60)


if __name__ == "__main__":
    main()
