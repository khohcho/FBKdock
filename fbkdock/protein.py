# Copyright (C) 2026 Fahrettin Buğra Kılıç
# SPDX-License-Identifier: GPL-3.0-or-later
"""
FBKdock Protein Preparation Module.

Receptor PDBQT preparation via ADFRsuite prepare_receptor4.py:
  1. Interactive residue cleaning (BioPython explorer).
  2. Split cleaned PDB into hetatm.pdb (non-water HETATM, chain preserved)
     and prot_wat.pdb / prot.pdb (ATOM + water lines, chain preserved).
  3. PDB2PQR (AMBER ff + PROPKA @ pH 7.4, --keep-chain) on the protein PDB.
  4. ADFRsuite prepare_receptor.bat converts PDBQ -> PDBQT
     (-A checkhydrogens -C -U nphs_lps) and hetatm.pdb -> PDBQT
     (-A bonds_hydrogens -U nphs_lps).
  5. Merge both PDBQTs with TER between chains and TER before the
     HETATM block, then interactive metal charge patching.
"""

import os
import sys
import shutil
import glob
import logging
import subprocess
from typing import Dict, List, Set, Tuple, Optional

from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
from Bio.PDB.Structure import Structure

from .utils import (
    run_subprocess, read_grid_params, update_config_file,
    ensure_dir, prompt_yes_no, detect_protein_key,
)

logger = logging.getLogger("FBKdock.Protein")

WATER_RESNAMES = {"HOH", "WAT", "H2O", "TIP", "TIP3", "TIP4", "TIP5", "SPC", "SOL"}
ION_RESNAMES = {
    "NA", "K", "CL", "CA", "MG", "ZN", "FE", "MN", "CU", "CO",
    "CD", "HG", "NI", "SO4", "PO4", "NO3", "ACT", "EDO", "GOL",
}

METAL_ELEMENTS = {"ZN", "MG", "CA", "FE", "MN", "CU", "CO", "NI", "K", "NA"}


class ResidueDesc:
    """Lightweight descriptor for one residue presented to the user."""

    def __init__(self, chain_id: str, resname: str, resid: int, insertion: str,
                 is_hetero: bool, is_water: bool, num_atoms: int):
        self.chain_id = chain_id
        self.resname = resname
        self.resid = resid
        self.insertion = insertion.strip()
        self.is_hetero = is_hetero
        self.is_water = is_water
        self.num_atoms = num_atoms

    def __str__(self) -> str:
        flag = ""
        if self.is_water:
            flag = " [WATER]"
        elif self.is_hetero:
            flag = " [HETERO]"
        ins = f":{self.insertion}" if self.insertion else ""
        return (f"Chain {self.chain_id:2s} | {self.resname:4s} {self.resid:5d}{ins}"
                f" | {self.num_atoms:4d} atoms{flag}")


def analyze_pdb_structure(pdb_path: str) -> Tuple[Structure, List[ResidueDesc]]:
    """
    Parse *pdb_path* with BioPython and return the Structure object together
    with a human-readable list of all residue descriptors.
    """
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"Protein PDB file not found: {pdb_path}")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    descriptors: List[ResidueDesc] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                hetfield, resid, ins_code = residue.get_id()
                is_hetero = (hetfield != " ")
                resname = residue.get_resname().strip()
                is_water = resname.upper() in WATER_RESNAMES
                atoms = list(residue.get_atoms())
                descriptors.append(ResidueDesc(
                    chain_id=chain.get_id(),
                    resname=resname,
                    resid=resid,
                    insertion=ins_code.strip(),
                    is_hetero=is_hetero,
                    is_water=is_water,
                    num_atoms=len(atoms),
                ))
    return structure, descriptors


def _interactive_explorer(
    descriptors: List[ResidueDesc],
    auto_remove_indices: Optional[Set[int]] = None,
) -> Set[int]:
    """
    Interactive hierarchical PDB explorer.

    Main menu lists each chain separately plus Heteroatoms / Waters / Finish.
    Returns the union of *auto_remove_indices* and user-selected indices.
    """
    if auto_remove_indices is None:
        auto_remove_indices = set()

    chains: Dict[str, List[int]] = {}
    hetero: List[int] = []
    waters: List[int] = []

    for idx, d in enumerate(descriptors):
        if idx in auto_remove_indices:
            continue
        if d.is_water:
            waters.append(idx)
        elif d.is_hetero:
            hetero.append(idx)
        else:
            chains.setdefault(d.chain_id, []).append(idx)

    chain_names = sorted(chains.keys())

    while True:
        print("\n" + "=" * 55)
        print("  PROTEIN EXPLORER - Main Menu")
        print("=" * 55)
        menu_idx = 1
        chain_map: Dict[int, str] = {}
        for cid in chain_names:
            label = f"Chain {cid} ({len(chains[cid])} residues)"
            print(f"  [{menu_idx}] {label}")
            chain_map[menu_idx] = cid
            menu_idx += 1
        print(f"  [{menu_idx}] Heteroatoms ({len(hetero)} residue(s))")
        hetero_menu = menu_idx
        menu_idx += 1
        print(f"  [{menu_idx}] Water Molecules ({len(waters)} residue(s))")
        waters_menu = menu_idx
        menu_idx += 1
        print(f"  [{menu_idx}] Finish & Delete")
        finish_menu = menu_idx
        if auto_remove_indices:
            auto_labels = ", ".join(str(i + 1) for i in sorted(auto_remove_indices))
            print(f"\n  >>> Residue(s) {auto_labels} will be removed automatically (ligand).")
        print("=" * 55)

        choice = input("Select [1-{}]: ".format(finish_menu)).strip()

        try:
            ci = int(choice)
        except ValueError:
            print("Invalid choice.")
            continue

        if ci in chain_map:
            _explore_list(chains[chain_map[ci]], f"Chain {chain_map[ci]}", descriptors)
        elif ci == hetero_menu:
            _explore_list(hetero, "HETEROATOMS", descriptors)
        elif ci == waters_menu:
            _explore_list(waters, "WATER MOLECULES", descriptors)
        elif ci == finish_menu:
            return _prompt_deletion(descriptors, auto_remove_indices)
        else:
            print(f"Invalid choice. Enter 1-{finish_menu}.")


def _explore_list(
    indices: List[int],
    title: str,
    descriptors: List[ResidueDesc],
) -> None:
    """Display a flat list of residues."""
    print(f"\n--- {title} ({len(indices)} residue(s)) ---")
    for j in indices:
        d = descriptors[j]
        print(f"  [{j + 1:3d}] {d}")
    print(f"  [b] Back")
    while True:
        inp = input("Enter 'b' to go back: ").strip().lower()
        if inp == "b":
            return


def _prompt_deletion(
    descriptors: List[ResidueDesc],
    auto_remove_indices: Set[int],
) -> Set[int]:
    """Prompt user for residue indices to delete.

    Input may mix residue numbers, ranges, and whole chains, e.g.
    ``1,2-4,(A),(E),100,187-190`` deletes residues 1, 2-4, 100, 187-190
    plus every residue of chains A and E.
    """
    while True:
        print("\n" + "=" * 55)
        print("  SELECT RESIDUES TO DELETE")
        print("=" * 55)
        print("  Enter indices, ranges, or whole chains, e.g. 1,2-4,(A),(E),100,187-190.")
        print("  A chain token like (A) deletes every residue of that chain.")
        print("  Enter 'none' to keep protein as-is (aside from auto-remove).")
        print("  Enter 'all' to delete everything.")

        raw = input("\nIndices to delete: ").strip()
        if raw.lower() == "none":
            return auto_remove_indices.copy()
        if raw.lower() == "all":
            return set(range(len(descriptors)))

        selected: Set[int] = set()
        selected_chains: Set[str] = set()
        invalid_chains: List[str] = []
        valid_chain_ids = {d.chain_id for d in descriptors}
        try:
            parts = raw.replace(",", " ").split()
            for part in parts:
                if part.startswith("(") and part.endswith(")"):
                    chain = part[1:-1].strip()
                    if not chain or chain not in valid_chain_ids:
                        invalid_chains.append(chain or part)
                        continue
                    selected_chains.add(chain)
                    selected.update(
                        i for i, d in enumerate(descriptors) if d.chain_id == chain
                    )
                elif "-" in part:
                    lo_s, hi_s = part.split("-", 1)
                    lo, hi = int(lo_s.strip()), int(hi_s.strip())
                    for i in range(lo, hi + 1):
                        if 1 <= i <= len(descriptors):
                            selected.add(i - 1)
                else:
                    i = int(part)
                    if 1 <= i <= len(descriptors):
                        selected.add(i - 1)
        except ValueError:
            print("Invalid format. Use numbers, ranges, or chains like (A).")
            continue

        if invalid_chains:
            print(
                "  Warning: unknown chain(s): %s"
                % ", ".join("(%s)" % c for c in invalid_chains)
            )

        final = auto_remove_indices | selected
        if not selected and not selected_chains:
            print("No additional residues selected.")
            if prompt_yes_no("Proceed with only auto-remove residues?"):
                return final
            continue

        print("\nResidues to be deleted (auto-remove + your selection):")
        for chain in sorted(selected_chains):
            count = sum(
                1 for i, d in enumerate(descriptors)
                if d.chain_id == chain and i in final
            )
            print(f"  Chain {chain}  -  {count} residue(s)")
        for i in sorted(final):
            marker = " [AUTO-REMOVE]" if i in auto_remove_indices else ""
            print(f"  [{i + 1:3d}] {descriptors[i]}{marker}")
        if prompt_yes_no("Confirm deletion of these residues?"):
            return final
        print("Redo selection.\n")


def clean_protein_pdb(
    structure: Structure,
    descriptors: List[ResidueDesc],
    remove_indices: Set[int],
    output_path: str,
) -> str:
    """Write a cleaned PDB to *output_path* with the selected residues removed."""
    residues_to_remove: Set[Tuple[int, str, Tuple[str, int, str]]] = set()
    idx = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                if idx in remove_indices:
                    key = (model.get_id(), chain.get_id(), residue.get_id())
                    residues_to_remove.add(key)
                idx += 1

    io = PDBIO()
    io.set_structure(structure)

    class CleanSelect(Select):
        def accept_residue(self, residue):
            try:
                model_id = residue.get_parent().get_parent().get_id()
            except AttributeError:
                model_id = 0
            try:
                chain_id = residue.get_parent().get_id()
            except AttributeError:
                chain_id = ""
            key = (model_id, chain_id, residue.get_id())
            return key not in residues_to_remove

    io.save(output_path, CleanSelect())
    logger.info("Cleaned PDB written to %s  (%d residues removed)", output_path, len(remove_indices))
    return output_path


# ---------------------------------------------------------------------------
# ADFRsuite Receptor Preparation Pipeline (PDB2PQR + prepare_receptor4.py)
# ---------------------------------------------------------------------------

ADFR_PREPARE_RECEPTOR = r"C:\Program Files (x86)\ADFRsuite-1.0\bin\prepare_receptor.bat"


def _find_prepare_receptor() -> str:
    """Locate the ADFRsuite prepare_receptor launcher.

    Resolution order: FBKDOCK_PREPARE_RECEPTOR env var, the known Windows
    ADFRsuite path, version-independent directory globs, then PATH lookup.
    """
    candidates: List[str] = []
    env_path = os.environ.get("FBKDOCK_PREPARE_RECEPTOR")
    if env_path:
        candidates.append(env_path)
    if os.name == "nt":
        candidates.append(ADFR_PREPARE_RECEPTOR)
        candidates.extend(
            glob.glob(r"C:\Program Files*\ADFRsuite*\bin\prepare_receptor.bat")
        )
        candidates.extend(glob.glob(r"C:\ADFRsuite*\bin\prepare_receptor.bat"))
        candidates.extend(
            glob.glob(
                os.path.expandvars(
                    r"%LOCALAPPDATA%\ADFRsuite*\bin\prepare_receptor.bat"
                )
            )
        )
        candidates.append("prepare_receptor")
    else:
        for base in ("/opt", "/usr/local", os.path.expanduser("~")):
            candidates.extend(
                glob.glob(os.path.join(base, "ADFRsuite*", "bin", "prepare_receptor"))
            )
        candidates.append("prepare_receptor")
        candidates.append("prepare_receptor.py")
    for cand in candidates:
        if os.path.isfile(cand):
            return cand
        resolved = shutil.which(cand)
        if resolved:
            return resolved
    raise RuntimeError(
        "ADFRsuite prepare_receptor not found. Install ADFRsuite or set "
        "FBKDOCK_PREPARE_RECEPTOR to the full path of prepare_receptor(.bat)."
    )


def _resolve_direct_python_invocation(
    launcher: str,
) -> Optional[Tuple[List[str], Dict[str, str]]]:
    """Resolve a Windows ADFRsuite .bat launcher to a shell-free invocation.

    The standard .bat calls '<bin>/python.exe' with the bundled
    AutoDockTools prepare_receptor4.py script and sets BABEL_DATADIR/PATH.
    Invoking the interpreter directly avoids cmd.exe argument parsing
    entirely.  Returns ([python_exe, script_path], env) or None.
    """
    if not launcher.lower().endswith(".bat"):
        return None
    bin_dir = os.path.dirname(os.path.abspath(launcher))
    root = os.path.dirname(bin_dir)
    python_exe = os.path.join(bin_dir, "python.exe")
    if not os.path.isfile(python_exe):
        python_exe = os.path.join(root, "python.exe")
    script = os.path.join(
        root, "Lib", "site-packages", "AutoDockTools",
        "Utilities24", "prepare_receptor4.py",
    )
    if not (os.path.isfile(python_exe) and os.path.isfile(script)):
        return None
    env: Dict[str, str] = {}
    babel_data = os.path.join(root, "OpenBabel-2.4.1", "data")
    if os.path.isdir(babel_data):
        env["BABEL_DATADIR"] = babel_data
    ob_dir = os.path.join(root, "OpenBabel-2.4.1")
    env["PATH"] = bin_dir + os.pathsep + ob_dir + os.pathsep + os.environ.get("PATH", "")
    return [python_exe, script], env


def _parse_bat_launcher(bat_path: str) -> Optional[Tuple[List[str], Dict[str, str]]]:
    """Parse an ADFRsuite-style .bat into a shell-free invocation.

    The standard launcher only sets environment variables and invokes
    '"<python.exe>" "<prepare_receptor4.py>" %*'.  Reproducing that
    invocation directly avoids cmd.exe entirely (no argument injection).
    Returns ([python_exe, script_path], env) or None if unparseable.
    """
    try:
        with open(bat_path, "r") as fh:
            content = fh.read()
    except OSError:
        return None

    env: Dict[str, str] = {}
    invocation: Optional[Tuple[str, str]] = None
    for raw_line in content.splitlines():
        line = raw_line.strip()
        if line[:4].lower() == "set ":
            key, _, value = line[4:].partition("=")
            env[key.strip()] = value.strip().strip('"')
        elif '"' in line and "%*" in line:
            parts = line.split('"')
            if len(parts) >= 5:
                exe, script = parts[1], parts[3]
                if os.path.isfile(exe) and os.path.isfile(script):
                    invocation = (exe, script)
    if not invocation:
        return None
    exe, script = invocation
    if env.get("path"):
        env["PATH"] = env.pop("path") + os.pathsep + os.environ.get("PATH", "")
    if not env.get("BABEL_DATADIR"):
        env.pop("BABEL_DATADIR", None)
    return [exe, script], env


def _run_prepare_receptor(
    input_path: str,
    output_path: str,
    repairs: str,
    preserve_charges: bool,
    cleanup: str,
) -> str:
    """Convert a PDB/PQR to PDBQT via ADFRsuite prepare_receptor4.py.

    On Windows the .bat launcher is never executed through cmd.exe; either
    the bundled python.exe + prepare_receptor4.py are resolved from the
    install layout or parsed out of the .bat itself.
    """
    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")
    ensure_dir(os.path.dirname(output_path) or ".")

    launcher = _find_prepare_receptor()
    args = ["-r", input_path, "-o", output_path, "-A", repairs]
    if preserve_charges:
        args.append("-C")
    args.extend(["-U", cleanup])

    env: Optional[Dict[str, str]] = None
    if os.name == "nt" and launcher.lower().endswith(".bat"):
        direct = _resolve_direct_python_invocation(launcher)
        if not direct:
            direct = _parse_bat_launcher(launcher)
        if direct:
            cmd, env = direct
            cmd = cmd + args
        else:
            raise RuntimeError(
                "ADFRsuite prepare_receptor.bat could not be resolved to its "
                "bundled Python interpreter. Set FBKDOCK_PREPARE_RECEPTOR to "
                "the prepare_receptor.bat of a standard ADFRsuite install or "
                "point it directly at the prepare_receptor4.py script path."
            )
    else:
        cmd = [launcher] + args

    print(f"  prepare_receptor ({repairs}, -U {cleanup}): {input_path}")
    logger.info("ADFRsuite prepare_receptor on %s", input_path)
    run_subprocess(cmd, timeout=1800, env=env)
    if not os.path.isfile(output_path):
        raise RuntimeError(f"prepare_receptor did not produce output: {output_path}")
    logger.info("PDBQT written to %s", output_path)
    return output_path


def _split_cleaned_pdb(
    cleaned_pdb: str,
    hetatm_path: str,
    prot_wat_path: str,
    prot_path: str,
) -> Tuple[str, Optional[str]]:
    """Split the cleaned PDB into a HETATM file and a protein+water file.

    hetatm.pdb receives all non-water HETATM lines.  The protein file
    receives all ATOM lines plus water HETATM lines (original order) and is
    named prot_wat.pdb when any water remains, otherwise prot.pdb.
    Chain information is preserved verbatim in both files.
    """
    hetatm_lines: List[str] = []
    prot_lines: List[str] = []
    has_water = False
    with open(cleaned_pdb, "r") as fh:
        for line in fh:
            if line.startswith("HETATM"):
                resname = line[17:20].strip().upper()
                if resname in WATER_RESNAMES:
                    has_water = True
                    prot_lines.append(line)
                else:
                    hetatm_lines.append(line)
            elif line.startswith("ATOM"):
                prot_lines.append(line)

    prot_output = prot_wat_path if has_water else prot_path
    with open(prot_output, "w") as fh:
        fh.writelines(prot_lines)
        fh.write("END\n")
    logger.info("Protein/water PDB written to %s (%d atom lines)", prot_output, len(prot_lines))

    hetatm_output: Optional[str] = None
    if hetatm_lines:
        with open(hetatm_path, "w") as fh:
            fh.writelines(hetatm_lines)
            fh.write("END\n")
        hetatm_output = hetatm_path
        logger.info("HETATM PDB written to %s (%d lines)", hetatm_path, len(hetatm_lines))
    else:
        logger.info("No non-water HETATM lines; hetatm.pdb skipped.")
    return prot_output, hetatm_output


# ---------------------------------------------------------------------------
# Standalone metal/ion handling (bypass ADFRsuite's isolated-atom crash)
# ---------------------------------------------------------------------------

def _list_hetatm_residues(
    hetatm_pdb: str,
) -> List[Tuple[str, str, str, str, List[str]]]:
    """Group the HETATM lines of *hetatm_pdb* into ordered residues.

    Returns (chain, resname, resSeq, iCode, [lines]) tuples in file order.
    """
    residues: List[Tuple[str, str, str, str, List[str]]] = []
    index: Dict[Tuple[str, str, str, str], int] = {}
    with open(hetatm_pdb, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.startswith(("ATOM", "HETATM")):
                continue
            key = (
                line[21:22],
                line[17:20].strip(),
                line[22:26].strip(),
                line[26:27].strip(),
            )
            if key not in index:
                index[key] = len(residues)
                residues.append((key[0], key[1], key[2], key[3], []))
            residues[index[key]][4].append(line)
    return residues


def _hetatm_residue_element(lines: List[str]) -> str:
    """Element symbol of a hetatm residue (first line, PDB columns 77-78)."""
    for line in lines:
        if len(line) >= 78:
            element = line[76:78].strip().upper()
            if element:
                return element
    return ""


def _select_isolated_metals(
    residues: List[Tuple[str, str, str, str, List[str]]],
    skip_interactive: bool,
) -> Set[int]:
    """Let the user pick which hetatm residues are standalone atoms.

    ADFRsuite crashes on isolated atoms whose element symbol starts with
    C/N/O/H/S/P (Cu, Na, Ni, Cl, Co, ...).  Selected residues bypass ADFR
    and are written directly as PDBQT; their charges are prompted at merge
    time.  In skip-interactive mode every single-atom residue is selected.
    Returns a set of 0-based residue indices.
    """
    single_atom_ids = {
        idx for idx, res in enumerate(residues) if len(res[4]) == 1
    }
    if skip_interactive:
        for idx in sorted(single_atom_ids):
            chain, resname, resseq, icode, _ = residues[idx]
            logger.info(
                "Auto-selecting standalone atom %s %s%s chain %s for direct PDBQT writing.",
                resname, resseq, icode, chain or "?",
            )
        return single_atom_ids

    print("\n" + "=" * 60)
    print("  SELECT STANDALONE METALS / IONS IN hetatm.pdb")
    print("=" * 60)
    print("  ADFRsuite crashes on isolated atoms whose element symbol starts")
    print("  with C/N/O/H/S/P (Cu, Na, Ni, Cl, Co, Cd, ...).")
    print("  Select those residues below; they are written directly as PDBQT")
    print("  (coordinates preserved, ADFR metal atom types) and their charge")
    print("  is asked during the merge step.")
    for idx, (chain, resname, resseq, icode, lines) in enumerate(residues):
        element = _hetatm_residue_element(lines)
        tag = " [TEK ATOM]" if len(lines) == 1 else f" ({len(lines)} atoms)"
        print(f"  [{idx + 1:2d}] Chain {chain or ' '} | {resname} {resseq}{icode} | {element}{tag}")
    print("  Enter indices (e.g. 1,3 or 1-3), 'all', or 'none' (default: none).")

    while True:
        raw = input("\nStandalone metal indices: ").strip()
        if raw.lower() in ("", "none"):
            return set()
        if raw.lower() == "all":
            selected = set(range(len(residues)))
        else:
            try:
                selected = set()
                for part in raw.replace(",", " ").split():
                    if "-" in part:
                        lo_s, hi_s = part.split("-", 1)
                        lo, hi = int(lo_s.strip()), int(hi_s.strip())
                        for i in range(lo, hi + 1):
                            if 1 <= i <= len(residues):
                                selected.add(i - 1)
                    else:
                        i = int(part)
                        if 1 <= i <= len(residues):
                            selected.add(i - 1)
            except ValueError:
                print("Invalid format. Use indices, ranges, 'all', or 'none'.")
                continue
        if not selected:
            continue

        multi = [i for i in sorted(selected) if len(residues[i][4]) > 1]
        for i in multi:
            chain, resname, resseq, icode, lines = residues[i]
            print(
                f"  WARNING: {resname} {resseq}{icode} (chain {chain}) has "
                f"{len(lines)} atoms - it will skip ADFR preparation "
                "(no hydrogens added, charges stay 0.000)."
            )
        print("\nSelected for direct PDBQT writing:")
        for i in sorted(selected):
            chain, resname, resseq, icode, lines = residues[i]
            element = _hetatm_residue_element(lines)
            print(f"  [{i + 1:2d}] Chain {chain or ' '} | {resname} {resseq}{icode} | {element}")
        if prompt_yes_no("Confirm this selection?"):
            return selected
        print("Redo selection.\n")


def _write_direct_hetatm_pdbqt(lines: List[str], output_path: str) -> str:
    """Write HETATM lines directly as PDBQT, bypassing ADFRsuite.

    Coordinates, chain and residue info are preserved verbatim.  Charges
    are placeholders (+0.000) that are overwritten at merge time / by
    _patch_metal_charges.  The AutoDock type is the title-cased element
    symbol (matches ADFRsuite's metal types: Cu, Zn, Na, Cl, ...) written
    in ADFR's column layout (charge 71-76, type 78-79).  TER is inserted
    on chain changes and the file ends with TER (no END - Vina-GPU 2.x).
    """
    out: List[str] = []
    serial = 0
    last_chain: Optional[str] = None
    for raw_line in lines:
        line = raw_line.rstrip("\n")
        if len(line) < 54:
            continue
        serial += 1
        chain = line[21:22]
        if last_chain is not None and chain != last_chain:
            out.append("TER")
        last_chain = chain
        element = line[76:78].strip() if len(line) >= 78 else ""
        ion_type = element[:1].upper() + element[1:].lower() if element else "X"
        occupancy = line[54:66] if len(line) >= 66 else " " * 12
        out.append(
            f"{line[:6]}{serial:>5}{line[11:54]}{occupancy}"
            f"{' ' * 4}{0.0:>+6.3f} {ion_type:>2}"
        )
    if out and not out[-1].startswith("TER"):
        out.append("TER")
    with open(output_path, "w") as fh:
        fh.write("\n".join(out) + "\n")
    logger.info("Direct PDBQT written for %d hetatm atom(s): %s", serial, output_path)
    return output_path


def _next_free_run_dir(parent: str, base: str) -> str:
    """First free directory of '<base>', '<base>_2', '<base>_3', ..."""
    if not os.path.isdir(os.path.join(parent, base)):
        return os.path.join(parent, base)
    k = 2
    while os.path.isdir(os.path.join(parent, f"{base}_{k}")):
        k += 1
    return os.path.join(parent, f"{base}_{k}")


def _renumber_large_resseq(pdb_path: str) -> int:
    """Renumber residues with resSeq >= 1000 to values <= 999.

    ADFRsuite's MolKit PQRParser loses chain/resSeq information when a PQR
    line has chain and resSeq concatenated (10 whitespace tokens, which
    happens for resSeq >= 1000 such as RCSB waters).  Renumbering these
    residues in the protein PDB before PDB2PQR keeps chain IDs intact.
    Only columns 23-26 (resSeq) are rewritten; chain IDs are untouched.
    A mapping (chain, resSeq, iCode) -> new_resSeq guarantees that every
    atom of a multi-atom residue receives the same new number.
    """
    used_nums: Dict[str, Set[int]] = {}
    renumber_map: Dict[Tuple[str, int, str], int] = {}
    renumbered = 0
    with open(pdb_path, "r") as fh:
        lines = fh.readlines()

    for line in lines:
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 27:
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        if resseq < 1000:
            used_nums.setdefault(line[21:22], set()).add(resseq)

    for i, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 27:
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        if resseq < 1000:
            continue
        chain = line[21:22]
        icode = line[26:27].strip()
        key = (chain, resseq, icode)
        if key not in renumber_map:
            nums = used_nums.setdefault(chain, set())
            new_seq = max(nums) + 1 if nums else 901
            while new_seq in nums:
                new_seq += 1
            if new_seq >= 1000:
                logger.warning(
                    "Could not renumber resSeq %s in chain %s below 1000; leaving as-is.",
                    resseq, chain or "?",
                )
                continue
            nums.add(new_seq)
            renumber_map[key] = new_seq
        lines[i] = line[:22] + f"{renumber_map[key]:>4}" + line[26:]
        renumbered += 1

    with open(pdb_path, "w") as fh:
        fh.writelines(lines)
    if renumbered:
        logger.info("Renumbered %d residue(s) with resSeq >= 1000 in %s", renumbered, pdb_path)
    return renumbered


def _detect_metals(pdb_path: str) -> Set[str]:
    """Scan HETATM lines for metal elements (PDB columns 76-78)."""
    metals: Set[str] = set()
    with open(pdb_path, "r") as fh:
        for line in fh:
            if not line.startswith("HETATM"):
                continue
            if len(line) < 78:
                continue
            element = line[76:78].strip().upper()
            if element in METAL_ELEMENTS:
                metals.add(element)
    if metals:
        print(f"  Metals detected: {', '.join(sorted(metals))}")
        logger.info("Metals detected: %s", sorted(metals))
    else:
        print("  No metals detected.")
        logger.info("No metals detected.")
    return metals


def _pdb2pqr_prepare(input_pdb: str, output_pqr: str) -> str:
    """Run PDB2PQR with AMBER ff + PROPKA at pH 7.4."""
    if not os.path.isfile(input_pdb):
        raise FileNotFoundError(f"Input PDB not found: {input_pdb}")
    ensure_dir(os.path.dirname(output_pqr) or ".")
    print("  PDB2PQR (AMBER ff + PROPKA @ pH 7.4) ...")
    logger.info("pdb2pqr on %s", input_pdb)
    run_subprocess(
        [sys.executable, "-m", "pdb2pqr",
         "--ff", "AMBER", "--titration-state-method", "propka",
         "--with-ph", "7.4", "--keep-chain",
         input_pdb, output_pqr],
        timeout=1800,
    )
    if not os.path.isfile(output_pqr):
        raise RuntimeError(f"PDB2PQR did not produce output: {output_pqr}")
    logger.info("PQR written to %s", output_pqr)
    return output_pqr


def _pqr_to_pdbq_for_adfr(pqr_path: str, pdbq_path: str) -> str:
    """Convert pdb2pqr PQR to fixed-column PDBQ for ADFRsuite.

    ADFRsuite's MolKit PQRParser decides between token-based and fixed-column
    parsing by counting whitespace tokens, which corrupts lines whose
    coordinate fields abut (negative / full-width values merge into one
    token, dropping chain IDs).  The PDBQ parser (PdbqParser) is purely
    fixed-column: it reads the charge from columns 71-76 and derives the
    element from the atom name.  Lines are therefore re-emitted in standard
    PDB columns with pdb2pqr-style name alignment (' CA ' not 'CA  ') so
    element derivation yields single-letter elements for protein atoms.
    Atom serials are wrapped into 1-9999 (only used for CONECT records,
    absent here).
    """
    with open(pqr_path, "r") as fh:
        lines = [line.rstrip("\n") for line in fh]

    out: List[str] = []
    serial = 0
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            out.append(line)
            continue
        if len(line) < 69:
            logger.warning("Skipping short PQR line: %s", line)
            out.append(line)
            continue
        record = line[:6].strip() or "ATOM"
        name = line[12:16].strip()
        resname = line[16:20].strip()
        chain = line[21:22]
        try:
            resseq_int = int(line[22:26])
        except ValueError:
            resseq_int = 0
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            charge = float(line[54:62])
        except ValueError:
            logger.warning("Skipping unparsable PQR line: %s", line)
            out.append(line)
            continue
        serial += 1
        wrapped_serial = ((serial - 1) % 9999) + 1
        if resseq_int >= 1000:
            resseq_int = 900 + (wrapped_serial % 100)
            logger.warning(
                "Renumbering resSeq >= 1000 to %d in PDBQ (chain %s).",
                resseq_int, chain,
            )
        if len(name) == 4:
            name_field = name.ljust(4)
        else:
            name_field = " " + name.ljust(3)
        resname_field = " " + resname.ljust(3)
        out.append(
            f"{record:<6}{wrapped_serial:>5} {name_field}{resname_field} "
            f"{chain}{resseq_int:>4}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"{' ' * 16}{charge:>+6.3f}"
        )
    with open(pdbq_path, "w") as fh:
        fh.write("\n".join(out) + "\n")
    logger.info("PDBQ written for ADFRsuite: %s", pdbq_path)
    return pdbq_path


def _merge_pdbqt(
    prot_pdbqt: str,
    hetatm_pdbqt: Optional[str],
    metals_pdbqt: Optional[str],
    final_path: str,
    selected_elements: Set[str],
    skip_interactive: bool = False,
) -> Tuple[str, Set[str]]:
    """Merge protein, ADFR-hetatm and directly-written metal PDBQTs.

    Layout: protein block (chains with TER between them), a TER line before
    the HETATM lines begin, the ADFR hetatm block (internal TERs kept),
    then the standalone-metal block.  Chain information in every line is
    preserved verbatim.  Charges for the selected standalone metals are
    prompted interactively here (defaults: +2.0; K/NA +1.0).  No trailing
    END record is written (Vina-GPU 2.x rejects it).

    Returns (final_path, charged_ion_types) so the later metal patching
    step skips these already-charged types.
    """
    if not os.path.isfile(prot_pdbqt):
        raise FileNotFoundError(f"Protein PDBQT not found: {prot_pdbqt}")

    # ---- resolve charges for the selected standalone metals ----
    metal_charges: Dict[str, float] = {}
    for element in sorted(selected_elements):
        ion_type = element.capitalize()
        default = "2.0"
        if element in ("K", "NA"):
            default = "1.0"
        if skip_interactive:
            metal_charges[ion_type] = float(default)
        else:
            raw = input(
                f"  Charge for {element} (default +{default}): "
            ).strip()
            metal_charges[ion_type] = float(raw) if raw else float(default)

    def _patch_block(lines: List[str]) -> List[str]:
        patched: List[str] = []
        for line in lines:
            if line.startswith("HETATM") and len(line) >= 78:
                fields = line.split()
                if fields and fields[-1] in metal_charges:
                    charge_str = f"{metal_charges[fields[-1]]:>+6.3f}"
                    line = line[:70] + charge_str + line[76:]
            patched.append(line)
        return patched

    with open(prot_pdbqt, "r") as fh:
        prot_lines = [line.rstrip("\n") for line in fh]
    while prot_lines and not prot_lines[-1].strip():
        prot_lines.pop()

    het_lines: List[str] = []
    if hetatm_pdbqt and os.path.isfile(hetatm_pdbqt):
        with open(hetatm_pdbqt, "r") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(("ATOM", "HETATM", "TER")):
                    het_lines.append(line)
    if hetatm_pdbqt and not het_lines:
        logger.warning("hetatm.pdbqt contains no atom records; skipping HETATM merge.")

    metals_lines: List[str] = []
    if metals_pdbqt and os.path.isfile(metals_pdbqt):
        with open(metals_pdbqt, "r") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(("ATOM", "HETATM", "TER")):
                    metals_lines.append(line)
    metals_lines = _patch_block(metals_lines)

    merged: List[str] = []
    if het_lines or metals_lines:
        while prot_lines and prot_lines[-1].startswith("TER"):
            prot_lines.pop()
        merged.extend(prot_lines)
        merged.append("TER")
        merged.extend(het_lines)
        if metals_lines:
            merged.extend(metals_lines)
    else:
        merged.extend(prot_lines)

    with open(final_path, "w") as fh:
        fh.write("\n".join(merged) + "\n")
    logger.info("Merged PDBQT written to %s (%d lines)", final_path, len(merged))
    return final_path, set(metal_charges.keys())


def _patch_metal_charges(
    pdbqt_path: str,
    metals: Set[str],
    skip_interactive: bool = False,
    skip_types: Optional[Set[str]] = None,
) -> str:
    """Interactively patch partial charges for metal ions in the merged PDBQT.

    ADFRsuite assigns metals 0.000 Gasteiger charge (no parameters), so the
    user is prompted for the oxidation-state charge per detected element.
    Only HETATM lines whose AutoDock type token matches the title-cased
    ion type (e.g. 'Na', 'Zn') are patched, so acceptor nitrogens typed
    'NA' in the protein block are never touched.

    *skip_types* contains ion types whose charges were already assigned at
    merge time (standalone metals); they are neither prompted nor patched.
    """
    skip_types = skip_types or set()
    metals = {m for m in metals if m.capitalize() not in skip_types}
    if not metals:
        logger.info("No metals to patch (all handled at merge time).")
        return pdbqt_path

    metal_charges: Dict[str, float] = {}
    for element in sorted(metals):
        ion_type = element.capitalize()
        default = "2.0"
        if element in ("K", "NA"):
            default = "1.0"
        if skip_interactive:
            metal_charges[ion_type] = float(default)
        else:
            raw = input(
                f"  Charge for {element} (default +{default}): "
            ).strip()
            metal_charges[ion_type] = float(raw) if raw else float(default)

    repaired: List[str] = []
    patched = 0
    with open(pdbqt_path, "r") as fh:
        for line in fh:
            if line.startswith("HETATM") and len(line) >= 78:
                fields = line.split()
                if fields and fields[-1] in metal_charges:
                    charge_val = metal_charges[fields[-1]]
                    charge_str = f"{charge_val:>+6.3f}"
                    line = line[:70] + charge_str + line[76:]
                    patched += 1
            repaired.append(line)

    with open(pdbqt_path, "w") as fh:
        fh.writelines(repaired)

    print(f"  Patched {patched} metal atom(s).")
    logger.info("Metal charge patching complete (%d atoms)", patched)
    return pdbqt_path


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def prepare_protein(
    pdb_path: str,
    docking_dir: str,
    redock_dir: str,
    config_template_path: str,
    skip_interactive: bool = False,
    remove_indices: Optional[Set[int]] = None,
    auto_remove_indices: Optional[Set[int]] = None,
) -> str:
    """
    Full protein-preparation pipeline:

    1. Parse PDB, interactive explorer for residue selection.
    2. Clean PDB (remove selected residues).
    3. Split (into Protein/prepared[_n]/) hetatm.pdb (non-water HETATM) and
       prot_wat.pdb / prot.pdb (ATOM + water lines); chains preserved.
    4. Standalone metal/ion selection: chosen hetatm residues are written
       directly as PDBQT (coordinates preserved, ADFR metal atom types);
       the rest goes through ADFRsuite.
    5. Renumber resSeq >= 1000 in the protein PDB (ADFR parser limits).
    6. PDB2PQR (AMBER ff + PROPKA @ pH 7.4, --keep-chain), then convert the
       PQR to fixed-column PDBQ for ADFRsuite's parser.
    7. ADFRsuite prepare_receptor:
       protein PDBQ -> prot_wat.pdbqt (-A checkhydrogens -C -U nphs_lps);
       remaining hetatm -> hetatm.pdbqt (-A bonds_hydrogens -U nphs_lps),
       with a direct-write fallback if ADFR crashes.
    8. Merge PDBQTs: TER between chains, TER before the HETATM block; the
       standalone metals' charges are prompted at merge time.
       (No trailing END - Vina-GPU 2.x rejects it.)
    9. Interactive metal charge patching for ADFR-processed metals.
    10. Copy protein PDBQT into *docking*/protein/ and *redock*/protein/.
    11. Update config.txt in both directories.

    Returns the path to the final protein PDBQT file.
    """
    protein_dir = os.path.dirname(pdb_path)
    base_name = os.path.splitext(os.path.basename(pdb_path))[0]

    structure, descriptors = analyze_pdb_structure(pdb_path)

    if remove_indices is None:
        if skip_interactive:
            remove_indices = (auto_remove_indices or set()).copy()
        else:
            remove_indices = _interactive_explorer(descriptors, auto_remove_indices)
    else:
        if auto_remove_indices:
            remove_indices = remove_indices | auto_remove_indices

    cleaned_pdb = os.path.join(protein_dir, f"{base_name}_cleaned.pdb")
    clean_protein_pdb(structure, descriptors, remove_indices, cleaned_pdb)

    # ---- split cleaned PDB into hetatm and protein(+water) streams ----
    # Intermediates live in a prepared[_n]/ subdirectory so they never
    # overwrite user inputs or results of earlier runs.
    work_dir = ensure_dir(_next_free_run_dir(protein_dir, "prepared"))
    hetatm_pdb = os.path.join(work_dir, "hetatm.pdb")
    prot_wat_pdb = os.path.join(work_dir, "prot_wat.pdb")
    prot_pdb = os.path.join(work_dir, "prot.pdb")
    prot_input, hetatm_input = _split_cleaned_pdb(
        cleaned_pdb, hetatm_pdb, prot_wat_pdb, prot_pdb,
    )

    # ---- standalone metal/ion selection + direct PDBQT writing ----
    metals_pdbqt: Optional[str] = None
    selected_elements: Set[str] = set()
    hetatm_pdbqt: Optional[str] = None
    if hetatm_input and os.path.isfile(hetatm_input):
        residues = _list_hetatm_residues(hetatm_input)
        chosen = _select_isolated_metals(residues, skip_interactive)
        selected_lines: List[str] = []
        rest_lines: List[str] = []
        for idx, res in enumerate(residues):
            target = selected_lines if idx in chosen else rest_lines
            target.extend(res[4])
        if selected_lines:
            metals_pdbqt = os.path.join(work_dir, "metals.pdbqt")
            _write_direct_hetatm_pdbqt(selected_lines, metals_pdbqt)
            selected_elements = {
                _hetatm_residue_element(res[4])
                for idx, res in enumerate(residues) if idx in chosen
            }
            selected_elements.discard("")
        if rest_lines:
            hetatm_rest_pdb = os.path.join(work_dir, "hetatm_rest.pdb")
            with open(hetatm_rest_pdb, "w") as fh:
                for line in rest_lines:
                    fh.write(line + "\n")
                fh.write("END\n")
            hetatm_pdbqt = os.path.join(work_dir, "hetatm.pdbqt")
            try:
                _run_prepare_receptor(
                    hetatm_rest_pdb, hetatm_pdbqt,
                    repairs="bonds_hydrogens", preserve_charges=False,
                    cleanup="nphs_lps",
                )
            except (subprocess.CalledProcessError, RuntimeError) as exc:
                logger.warning(
                    "ADFRsuite failed on the hetatm stream (%s); writing the "
                    "remaining residues directly with 0.000 charges.",
                    exc,
                )
                _write_direct_hetatm_pdbqt(rest_lines, hetatm_pdbqt)

    # ---- ADFR PQR parser cannot handle resSeq >= 1000 -> renumber ----
    _renumber_large_resseq(prot_input)

    prot_stem = os.path.splitext(os.path.basename(prot_input))[0]
    protein_pqr = os.path.join(work_dir, f"{prot_stem}.pqr")
    _pdb2pqr_prepare(prot_input, protein_pqr)
    protein_pdbq = os.path.join(work_dir, f"{prot_stem}.pdbq")
    _pqr_to_pdbq_for_adfr(protein_pqr, protein_pdbq)

    prot_pdbqt = os.path.join(work_dir, f"{prot_stem}.pdbqt")
    _run_prepare_receptor(
        protein_pdbq, prot_pdbqt,
        repairs="checkhydrogens", preserve_charges=True, cleanup="nphs_lps",
    )

    # ---- merge with TER records; standalone-metal charges asked here ----
    pdbqt_path = os.path.join(work_dir, f"{base_name}.pdbqt")
    _, charged_types = _merge_pdbqt(
        prot_pdbqt, hetatm_pdbqt, metals_pdbqt, pdbqt_path,
        selected_elements, skip_interactive=skip_interactive,
    )

    metals = _detect_metals(cleaned_pdb)
    _patch_metal_charges(
        pdbqt_path, metals,
        skip_interactive=skip_interactive, skip_types=charged_types,
    )

    # ---- Copy to workspace & update config ----
    protein_basename = os.path.basename(pdbqt_path)
    protein_key = detect_protein_key(config_template_path)

    for target_dir in (docking_dir, redock_dir):
        protein_target = os.path.join(target_dir, "protein", protein_basename)
        ensure_dir(os.path.dirname(protein_target))
        shutil.copy2(pdbqt_path, protein_target)
        logger.info("Protein PDBQT copied to %s", protein_target)

    grid_params = read_grid_params(config_template_path)
    config_updates = {
        protein_key: f"protein/{protein_basename}",
        "center_x": str(grid_params["center_x"]),
        "center_y": str(grid_params["center_y"]),
        "center_z": str(grid_params["center_z"]),
        "size_x": str(grid_params["size_x"]),
        "size_y": str(grid_params["size_y"]),
        "size_z": str(grid_params["size_z"]),
    }
    for target_dir in (docking_dir, redock_dir):
        cfg_path = os.path.join(target_dir, "config.txt")
        update_config_file(cfg_path, config_updates)
        logger.info("Config updated: %s", cfg_path)

    return pdbqt_path
