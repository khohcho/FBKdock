# Copyright (C) 2026 Fahrettin Buğra Kılıç
# SPDX-License-Identifier: GPL-3.0-or-later
"""
FBKdock Utilities - OS detection, config handling, subprocess execution,
and file-system operations shared across all pipeline modules.
"""

import os
import sys
import shutil
import glob
import subprocess
import logging
import threading
from pathlib import Path
from typing import Dict, Optional, List, Tuple

logger = logging.getLogger("FBKdock")

# When True, subprocess output is streamed live to the terminal.
# Set from main.py based on the --verbose flag.
LIVE_OUTPUT: bool = False


def detect_os() -> str:
    """Return 'windows' or 'linux' based on the host platform."""
    if sys.platform.startswith("win"):
        return "windows"
    return "linux"


OS_TYPE: str = detect_os()


def get_runner_script() -> str:
    """Return the appropriate launcher script name for the current OS."""
    if OS_TYPE == "windows":
        return "run.bat"
    return "run.sh"


def get_vina_exe_name(engine: str = "vinagpu") -> str:
    """Return the docking executable filename for the current OS and engine."""
    if engine == "autodockvina":
        return "vina.exe" if OS_TYPE == "windows" else "vina"
    else:
        return "Vina-GPU.exe" if OS_TYPE == "windows" else "Vina-GPU"


def run_subprocess_live(
    cmd: List[str],
    cwd: Optional[str] = None,
    timeout: Optional[int] = None,
    stdin: Optional[int] = None,
    env: Optional[Dict[str, str]] = None,
) -> subprocess.CompletedProcess:
    """
    Execute a command via subprocess with REAL-TIME stdout display.

    Lines are printed to terminal as they arrive AND captured for logging.
    Raises subprocess.CalledProcessError on non-zero exit.
    """
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)

    logger.info("Running: %s  (cwd=%s)", " ".join(cmd), cwd or ".")
    process = subprocess.Popen(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        encoding="utf-8",
        errors="replace",
        stdin=stdin,
        env=merged_env,
        bufsize=1,
        universal_newlines=True,
    )

    output_lines: List[str] = []

    def _read_stdout() -> None:
        if process.stdout is None:
            return
        for line in iter(process.stdout.readline, ""):
            if LIVE_OUTPUT:
                print(line, end="")
            output_lines.append(line)
        process.stdout.close()

    reader = threading.Thread(target=_read_stdout, daemon=True)
    reader.start()

    try:
        process.wait(timeout=timeout)
        reader.join(timeout=5)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()
        raise subprocess.TimeoutExpired(cmd, timeout)

    if process.returncode != 0:
        logger.error("Command failed (rc=%d): %s", process.returncode, " ".join(cmd))
        logger.error("STDOUT/STDERR: %s", "".join(output_lines[-20:]).strip())
        raise subprocess.CalledProcessError(
            process.returncode, cmd,
            output="".join(output_lines),
            stderr="",
        )

    return subprocess.CompletedProcess(
        args=cmd,
        returncode=0,
        stdout="".join(output_lines),
        stderr="",
    )


def run_subprocess(
    cmd: List[str],
    cwd: Optional[str] = None,
    timeout: Optional[int] = None,
    stdin: Optional[int] = None,
    env: Optional[Dict[str, str]] = None,
) -> subprocess.CompletedProcess:
    """Alias for run_subprocess_live - all subprocess calls use live logging."""
    return run_subprocess_live(cmd, cwd=cwd, timeout=timeout, stdin=stdin, env=env)


def run_runner_script(cwd: str, timeout: Optional[int] = None) -> subprocess.CompletedProcess:
    """
    Execute the launcher script (run.bat / run.sh) in *cwd*.
    On Windows, stdin is piped from DEVNULL so that ``pause`` inside
    run.bat does not block.
    """
    script = get_runner_script()
    if OS_TYPE == "windows":
        return run_subprocess(
            ["cmd.exe", "/c", script],
            cwd=cwd,
            stdin=subprocess.DEVNULL,
            timeout=timeout,
        )
    else:
        return run_subprocess(
            ["bash", script],
            cwd=cwd,
            timeout=timeout,
        )


def run_vina_directly(
    cwd: str,
    engine: str = "vinagpu",
    config_path: str = "config.txt",
    split_log: bool = True,
    log_dir: Optional[str] = None,
    timeout: Optional[int] = None,
) -> subprocess.CompletedProcess:
    """
    Call the docking binary directly (bypassing run.bat/run.sh).
    """
    exe = os.path.join(cwd, get_vina_exe_name(engine))
    if not os.path.isfile(exe):
        raise FileNotFoundError(f"Docking binary not found: {exe}")
    cmd = [exe, "--config", config_path]
    if split_log:
        cmd.append("--split_log")
    if log_dir:
        cmd.extend(["--log_dir", log_dir])
    return run_subprocess(cmd, cwd=cwd, timeout=timeout)


def parse_config(config_path: str) -> Dict[str, str]:
    """Parse a key = value configuration file into a dict."""
    config: Dict[str, str] = {}
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if "=" not in stripped:
                continue
            key, _, value = stripped.partition("=")
            config[key.strip()] = value.strip()
    return config


def write_config(config_path: str, config: Dict[str, str]) -> None:
    """Write a dict to a key = value configuration file."""
    os.makedirs(os.path.dirname(config_path) or ".", exist_ok=True)
    with open(config_path, "w", encoding="utf-8") as fh:
        for key, value in config.items():
            fh.write(f"{key} = {value}\n")


def update_config_file(
    config_path: str,
    updates: Dict[str, str],
) -> None:
    """Read an existing config file, apply *updates*, and write it back."""
    cfg = parse_config(config_path)
    cfg.update(updates)
    write_config(config_path, cfg)


def detect_protein_key(config_path: str) -> str:
    """
    Detect which key the config template uses for the protein path.
    Returns 'protein' or 'receptor' based on what's found in the file.
    """
    cfg = parse_config(config_path)
    if "protein" in cfg:
        return "protein"
    if "receptor" in cfg:
        return "receptor"
    return "protein"  # default


def read_grid_params(config_path: str) -> Dict[str, float]:
    """
    Extract the grid centre and size parameters from a config file.
    Returns keys: center_x, center_y, center_z, size_x, size_y, size_z.
    """
    cfg = parse_config(config_path)
    required = ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]
    params: Dict[str, float] = {}
    for key in required:
        if key not in cfg:
            raise KeyError(f"Missing required grid parameter '{key}' in {config_path}")
        params[key] = float(cfg[key])
    return params


def ensure_dir(path: str) -> str:
    """Create a directory (including parents) if it does not exist."""
    os.makedirs(path, exist_ok=True)
    return path


def clean_directory(path: str) -> None:
    """Remove all files and subdirectories inside *path* (does not remove *path* itself)."""
    if not os.path.isdir(path):
        return
    for entry in os.listdir(path):
        full = os.path.join(path, entry)
        if os.path.isfile(full) or os.path.islink(full):
            os.unlink(full)
        elif os.path.isdir(full):
            shutil.rmtree(full)


def ensure_executable(workspace_dir: str) -> None:
    """Set the executable bit on engine binaries (Linux only).

    Files copied from a Windows filesystem (/mnt/c) lose the +x bit on
    Linux; without it run.sh cannot execute ./Vina-GPU or ./vina.
    """
    if OS_TYPE != "linux":
        return
    for name in ("Vina-GPU", "vina", "vina_split"):
        path = os.path.join(workspace_dir, name)
        if os.path.isfile(path):
            os.chmod(path, 0o755)


def copy_directory_contents(src: str, dst: str) -> None:
    """Recursively copy the *contents* of *src* into *dst*."""
    ensure_dir(dst)
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            if os.path.exists(d):
                shutil.rmtree(d)
            shutil.copytree(s, d)
        else:
            shutil.copy2(s, d)


def replace_engine_binary(workspace_dir: str, autodockvina_source: str) -> None:
    """
    Replace Vina-GPU binaries in *workspace_dir* with AutoDock Vina CPU.

    Copies vina.exe, vina_split.exe, boost DLLs, run.bat, and run.sh
    from *autodockvina_source* into *workspace_dir*, overwriting the
    GPU binaries.  Also converts config keys from Vina-GPU format
    (ligand_directory, output_directory) to AutoDock Vina format
    (batch, dir).
    """
    config_path = os.path.join(workspace_dir, "config.txt")
    if os.path.isfile(config_path):
        cfg = parse_config(config_path)
        if "ligand_directory" in cfg:
            cfg["batch"] = cfg.pop("ligand_directory")
        if "output_directory" in cfg:
            cfg["dir"] = cfg.pop("output_directory")
        for key in ("opencl_binary_path", "thread", "search_depth"):
            cfg.pop(key, None)
        if "exhaustiveness" not in cfg:
            cfg["exhaustiveness"] = "8"
        if "scoring" not in cfg:
            cfg["scoring"] = "vina"
        if "verbosity" not in cfg:
            cfg["verbosity"] = "1"
        if "cpu" not in cfg:
            cfg["cpu"] = "0"
        write_config(config_path, cfg)

    for fname in os.listdir(autodockvina_source):
        src_path = os.path.join(autodockvina_source, fname)
        dst_path = os.path.join(workspace_dir, fname)
        if os.path.isfile(src_path):
            shutil.copy2(src_path, dst_path)
        elif os.path.isdir(src_path) and not os.path.exists(dst_path):
            shutil.copytree(src_path, dst_path)
    ensure_executable(workspace_dir)


def recover_workspace_to_cpu(
    crashed_workspace: str,
    autodockvina_source: str,
    project_root: str,
) -> str:
    """
    Save prepared files from *crashed_workspace*, wipe it, recreate from
    *autodockvina_source*, and restore the saved files.  Returns the
    path to the recovered workspace (same as *crashed_workspace*).

    Used when Vina-GPU crashes and the user wants to fall back to
    CPU-based AutoDock Vina without losing prepared ligand PDBQTs,
    protein, config, or logs.
    """
    recovery_dir = os.path.join(project_root, "template_recovery")
    ensure_dir(recovery_dir)

    for sub in ("input", "protein", "output", "logs"):
        src_sub = os.path.join(crashed_workspace, sub)
        dst_sub = os.path.join(recovery_dir, sub)
        if os.path.isdir(src_sub):
            shutil.copytree(src_sub, dst_sub)

    config_src = os.path.join(crashed_workspace, "config.txt")
    if os.path.isfile(config_src):
        shutil.copy2(config_src, os.path.join(recovery_dir, "config.txt"))

    logger.info("Workspace %s saved to template_recovery", crashed_workspace)

    shutil.rmtree(crashed_workspace)
    copy_directory_contents(autodockvina_source, crashed_workspace)

    for sub in ("input", "protein", "output", "logs"):
        src_sub = os.path.join(recovery_dir, sub)
        dst_sub = os.path.join(crashed_workspace, sub)
        if os.path.isdir(src_sub):
            if os.path.exists(dst_sub):
                shutil.rmtree(dst_sub)
            shutil.copytree(src_sub, dst_sub)

    recovered_config = os.path.join(recovery_dir, "config.txt")
    if os.path.isfile(recovered_config):
        shutil.copy2(recovered_config, os.path.join(crashed_workspace, "config.txt"))

    shutil.rmtree(recovery_dir)
    replace_engine_binary(crashed_workspace, autodockvina_source)

    logger.info("Workspace recovered as AutoDock Vina: %s", crashed_workspace)
    return crashed_workspace


def find_files(directory: str, pattern: str) -> List[str]:
    """Return sorted list of file paths matching *pattern* inside *directory*."""
    return sorted(glob.glob(os.path.join(directory, pattern)))


def find_fbkdock_files(directory: str) -> List[str]:
    """Return sorted list of all readable files in *directory* (any extension)."""
    if not os.path.isdir(directory):
        return []
    files = []
    for entry in os.listdir(directory):
        full = os.path.join(directory, entry)
        if os.path.isfile(full):
            files.append(full)
    return sorted(files)


def find_pdb_files(directory: str) -> List[str]:
    """Return sorted list of all files in *directory* (any extension, case-insensitive)."""
    return find_fbkdock_files(directory)


def find_sdf_files(directory: str) -> List[str]:
    """Return sorted list of all files in *directory* (any extension, case-insensitive)."""
    return find_fbkdock_files(directory)


def find_pdbqt_files(directory: str) -> List[str]:
    """Return sorted list of all files in *directory* (any extension, case-insensitive)."""
    return find_fbkdock_files(directory)


def prompt_choice(prompt_text: str, options: List[str]) -> int:
    """Display a numbered list and return the user's 1-based index choice."""
    print(f"\n{prompt_text}")
    for idx, option in enumerate(options, start=1):
        print(f"  [{idx}] {option}")
    while True:
        try:
            choice = input("Enter number (or 0 to abort): ").strip()
            if choice == "0":
                sys.exit(0)
            idx = int(choice)
            if 1 <= idx <= len(options):
                return idx
            print(f"Please enter a number between 1 and {len(options)}.")
        except (ValueError, EOFError):
            print("Invalid input. Enter a number.")


def prompt_yes_no(question: str) -> bool:
    """Ask a yes/no question; return True for 'y'/'yes'."""
    while True:
        ans = input(f"{question} [y/n]: ").strip().lower()
        if ans in ("y", "yes"):
            return True
        if ans in ("n", "no"):
            return False
        print("Please answer 'y' or 'n'.")


def prompt_float(prompt_text: str, default: Optional[float] = None) -> float:
    """Prompt the user for a floating-point number with an optional default."""
    default_str = f" [{default}]" if default is not None else ""
    while True:
        try:
            val = input(f"{prompt_text}{default_str}: ").strip()
            if not val and default is not None:
                return default
            return float(val)
        except ValueError:
            print("Invalid number. Please try again.")
