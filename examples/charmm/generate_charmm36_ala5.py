#!/usr/bin/env python3
"""
generate_charmm36_ala5.py

Regenerates ``charmm36_ala5.prmtop`` / ``charmm36_ala5.rst7`` in this
directory: a clean CHARMM36 CHAMBER-format AMBER prmtop/rst7 for a small
capped peptide (ACE-(ALA)5-NME). Exercises the CHAMBER-specific sections of
the prmtop format (CMAP, Urey-Bradley, CHARMM impropers, CHARMM 1-4
nonbonded) end to end, with NO NBFIX and NO box, so it can be validated both
against native OpenMM (``openmm.app.AmberPrmtopFile``/``AmberInpcrdFile``) and
against Robosample's own loader (see
``tests/test_openmm_potential_energy.py::charmm36_ala5``).

Pipeline
--------
1. ``tleap`` (AmberTools, ff14SB) builds an extended-conformation seed PDB for
   the capped pentapeptide. This step only supplies a reasonable, clash-free,
   non-degenerate starting *geometry* -- AMBER's ff14SB and CHARMM36 standard
   backbone/CB bond lengths and angles are close enough that the coordinates
   are a fine seed once the atom names are remapped (step 2).
2. AMBER atom names/residue numbering are remapped to CHARMM36 naming.
   CHARMM applies the ACE/NME caps as *patches* fused onto the flanking ALA
   residues (not standalone residues, unlike AMBER's ACE/NME), so the tleap
   ACE/NME residues are folded into psfgen residues 1 and 5 under the
   PRES ACE / PRES CT3 atom names (see ``_remap_amber_to_charmm``).
3. VMD's ``psfgen`` Tcl plugin (topology ``top_all36_prot.rtf``) builds the
   PSF + coordinates for the CHARMM-typed structure from the remapped seed
   PDB via ``coordpdb`` (every atom is matched by name -- ``guesscoord`` has
   nothing left to guess).
4. A small deterministic Gaussian jitter is added to the coordinates. The
   raw tleap-built chain is exactly planar along the backbone (phi = psi =
   180 deg), which would make the harmonic-improper (CHARMM_IMPROPERS) and
   CMAP correction-map energies come out at *exactly* zero for every term --
   numerically "correct" but a self-defeating regression test (a broken
   improper/CMAP implementation that always returns 0 would pass too). The
   jitter (sigma = 0.1 Angstrom, fixed seed) breaks the exact planarity while
   staying far below any bond length, so every bonded force group ends up
   genuinely non-zero without producing clashes or a non-physical structure.
5. ParmEd (``CharmmPsfFile`` + ``CharmmParameterSet(RTF, PRM)``) loads the
   typed structure and its CHARMM36 parameters, then writes it out as a
   CHAMBER-format AMBER prmtop + rst7 (``Structure.save(..., format='amber')``
   / ``format='rst7'``).

Requires (present in the ``robo_cuda13.0`` conda env / system):
  - ``tleap`` (AmberTools, on PATH)
  - ``vmd`` with the ``psfgen`` Tcl plugin (on PATH; VMD 1.9's psfgen1.4 was
    used to develop this script)
  - ``parmed`` (Python, importable)
  - CHARMM36 protein parameters: ``$CONDA_PREFIX/dat/chamber/{top,par}_all36_prot.*``

Run: ``python3 examples/charmm/generate_charmm36_ala5.py``
"""

from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
N_ALA = 5
JITTER_SIGMA_ANGSTROM = 0.1
JITTER_SEED = 42

# AMBER (tleap ff14SB) -> CHARMM36 (top_all36_prot.rtf) atom-name remap. CHARMM
# applies ACE/NME as *patches* fused onto the flanking ALA residues (not
# separate residues), so the AMBER ACE (resid 1) and NME (resid N_ALA + 2)
# atoms are folded into CHARMM psfgen resid 1 and N_ALA resp.
_ACE_MAP = {"CH3": "CAY", "H1": "HY1", "H2": "HY2", "H3": "HY3", "C": "CY", "O": "OY"}
_NME_MAP = {"N": "NT", "H": "HNT", "C": "CAT", "H1": "HT1", "H2": "HT2", "H3": "HT3"}
_ALA_MAP = {"H": "HN"}  # everything else (N,CA,HA,CB,HB1,HB2,HB3,C,O) is identical


def _require(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        raise RuntimeError(
            f"required tool '{tool}' not found on PATH -- cannot regenerate "
            "the CHARMM36 CHAMBER test asset. See this script's module "
            "docstring for the full tool list."
        )
    return path


def _conda_prefix() -> pathlib.Path:
    prefix = os.environ.get("CONDA_PREFIX")
    if not prefix:
        raise RuntimeError(
            "CONDA_PREFIX is not set -- activate the robo_cuda13.0 (or "
            "robo_cpu) conda env before running this script."
        )
    return pathlib.Path(prefix)


def _build_tleap_seed(tleap: str, workdir: pathlib.Path) -> pathlib.Path:
    """Extended-conformation ACE-(ALA)n-NME seed PDB, AMBER ff14SB naming."""
    seed_pdb = workdir / "seed_amber.pdb"
    tleap_in = workdir / "tleap_seed.in"
    sequence = " ".join(["ACE"] + ["ALA"] * N_ALA + ["NME"])
    tleap_in.write_text(
        "source leaprc.protein.ff14SB\n"
        f"pep = sequence {{ {sequence} }}\n"
        f"savepdb pep {seed_pdb}\n"
        "quit\n"
    )
    subprocess.run(
        [tleap, "-f", str(tleap_in)], cwd=workdir, check=True, capture_output=True
    )
    if not seed_pdb.exists():
        raise RuntimeError("tleap did not produce the expected seed PDB")
    return seed_pdb


def _parse_pdb_atoms(path: pathlib.Path):
    atoms = []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        resname = line[17:20].strip()
        resid = int(line[22:26])
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        atoms.append((name, resname, resid, x, y, z))
    return atoms


def _remap_amber_to_charmm(atoms):
    """AMBER (name, resname, resid) -> CHARMM (name, 'ALA', psfgen resid)."""
    out = []
    for name, resname, resid, x, y, z in atoms:
        if resname == "ACE":
            cname, cresid = _ACE_MAP[name], 1
        elif resname == "NME":
            cname, cresid = _NME_MAP[name], N_ALA
        elif resname == "ALA":
            cname, cresid = _ALA_MAP.get(name, name), resid - 1
        else:
            raise ValueError(f"unexpected residue {resname!r} in tleap seed PDB")
        out.append((cname, "ALA", cresid, x, y, z))
    return out


def _write_pdb(records, path: pathlib.Path) -> None:
    lines = []
    for i, (name, resname, resid, x, y, z) in enumerate(records, start=1):
        atom_field = f" {name:<3s}" if len(name) < 4 else name
        lines.append(
            f"ATOM  {i:5d} {atom_field:<4s} {resname:>3s} P{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def _build_psf(
    vmd: str, rtf: pathlib.Path, charmm_seed_pdb: pathlib.Path, workdir: pathlib.Path
) -> tuple[pathlib.Path, pathlib.Path]:
    """Run VMD's psfgen to build the CHARMM-typed PSF + coordinate PDB."""
    psf = workdir / "ala5.psf"
    pdb = workdir / "ala5.pdb"
    residues = "\n".join(f"    residue {i} ALA" for i in range(1, N_ALA + 1))
    tcl = workdir / "build_psf.tcl"
    tcl.write_text(
        "package require psfgen\n"
        "resetpsf\n"
        f"topology {rtf}\n"
        "segment PEP {\n"
        "    first ACE\n"
        "    last CT3\n"
        f"{residues}\n"
        "}\n"
        f"coordpdb {charmm_seed_pdb} PEP\n"
        "guesscoord\n"
        f"writepsf {psf}\n"
        f"writepdb {pdb}\n"
        "exit\n"
    )
    result = subprocess.run(
        [vmd, "-dispdev", "text", "-e", str(tcl)],
        cwd=workdir,
        check=True,
        capture_output=True,
        text=True,
    )
    if "guessing coordinates for 0 atoms" not in result.stdout:
        raise RuntimeError(
            "psfgen had to guess >=1 atom coordinate (expected 0 -- every "
            "atom should have been matched by coordpdb); the AMBER->CHARMM "
            "name remap is likely incomplete. psfgen output:\n" + result.stdout
        )
    if not (psf.exists() and pdb.exists()):
        raise RuntimeError("psfgen did not produce the expected PSF/PDB")
    return psf, pdb


def main() -> None:
    tleap = _require("tleap")
    vmd = _require("vmd")
    conda_prefix = _conda_prefix()
    rtf = conda_prefix / "dat" / "chamber" / "top_all36_prot.rtf"
    prm = conda_prefix / "dat" / "chamber" / "par_all36_prot.prm"
    if not rtf.exists() or not prm.exists():
        raise RuntimeError(f"CHARMM36 params not found: {rtf} / {prm}")

    import parmed as pmd
    from parmed.charmm import CharmmParameterSet, CharmmPsfFile

    with tempfile.TemporaryDirectory(prefix="charmm36_ala5_") as tmp:
        workdir = pathlib.Path(tmp)

        seed_pdb = _build_tleap_seed(tleap, workdir)
        amber_atoms = _parse_pdb_atoms(seed_pdb)
        charmm_records = _remap_amber_to_charmm(amber_atoms)
        charmm_seed_pdb = workdir / "seed_charmm.pdb"
        _write_pdb(charmm_records, charmm_seed_pdb)

        psf_path, pdb_path = _build_psf(vmd, rtf, charmm_seed_pdb, workdir)

        params = CharmmParameterSet(str(rtf), str(prm))
        psf = CharmmPsfFile(str(psf_path))
        psf.load_parameters(params)

        coord_struct = pmd.load_file(str(pdb_path))
        rng = np.random.default_rng(JITTER_SEED)
        jitter = rng.normal(scale=JITTER_SIGMA_ANGSTROM, size=coord_struct.coordinates.shape)
        psf.coordinates = coord_struct.coordinates + jitter
        psf.box = None  # no periodicity -- genuinely absent, not degenerate

        out_prmtop = HERE / "charmm36_ala5.prmtop"
        out_rst7 = HERE / "charmm36_ala5.rst7"
        psf.save(str(out_prmtop), format="amber", overwrite=True)
        psf.save(str(out_rst7), format="rst7", overwrite=True)

        print(f"wrote {out_prmtop}")
        print(f"wrote {out_rst7}")
        print(f"  {len(psf.atoms)} atoms, {len(psf.cmaps)} CMAP terms, "
              f"{len(psf.urey_bradleys)} Urey-Bradley terms, "
              f"{len(psf.impropers)} CHARMM impropers")


if __name__ == "__main__":
    main()
