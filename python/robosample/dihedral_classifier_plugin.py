import os
import subprocess
import sys
import tempfile

from pymol import cmd

# Path to the worker script sitting next to this file
WORKER = "/home/victor/Robosample/python/robosample/dihedral_classifier_worker.py"


def classify_dihedrals(selection="all"):
    tmp_pdb = None

    try:
        fd, tmp_pdb = tempfile.mkstemp(suffix=".pdb")
        os.close(fd)

        print(f"[plugin] saving selection '{selection}' -> {tmp_pdb}")

        cmd.save(tmp_pdb, selection)

        result = subprocess.run(
            [sys.executable, WORKER, tmp_pdb],
            capture_output=True,
            text=True,
            check=False,
            cwd="/home/victor/Robosample/python/robosample",
        )

        print("\n========== WORKER STDOUT ==========")
        print(result.stdout)

        if result.stderr:
            print("\n========== WORKER STDERR ==========")
            print(result.stderr)

        print(f"\n[plugin] return code = {result.returncode}")

    finally:
        if tmp_pdb and os.path.exists(tmp_pdb):
            os.remove(tmp_pdb)


cmd.extend("classify_dihedrals", classify_dihedrals)
